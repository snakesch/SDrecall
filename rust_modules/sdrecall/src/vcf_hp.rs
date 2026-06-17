//! HP-support annotation for clean per-island VCFs.
//!
//! Ports `realign_recall/annotate_HP_tag_to_vars.py::annotate_vcf`: for each
//! variant allele, collect HP tags from reads supporting that allele in the
//! HP-tagged clean BAM and write them as FORMAT/HPSUP.

use std::collections::{BTreeSet, HashMap};
use std::path::{Path, PathBuf};

use rust_htslib::bam::pileup::Indel;
use rust_htslib::bam::record::Aux;
use rust_htslib::{bam, bam::Read as BamRead, bcf, bcf::Read as BcfRead};
use sdrecall_utils::{Result, SdError};

const MIN_MAPQ: u8 = 10;
const MIN_BASEQ: u8 = 15;

pub fn annotate_vcf_hp(
    input_vcf: &Path,
    output_vcf: &Path,
    hp_bam: &Path,
    threads: usize,
) -> Result<()> {
    let mut reader = bcf::Reader::from_path(input_vcf)
        .map_err(|e| SdError::Vcf(format!("open {}: {e}", input_vcf.display())))?;
    let input_header = reader.header().clone();

    let mut header = bcf::Header::from_template(reader.header());
    header.push_record(
        b"##FORMAT=<ID=HPSUP,Number=A,Type=String,Description=\"Supporting HP values for each variant allele\">",
    );

    let tmp_vcf = temp_vcf_path(output_vcf);
    let mut writer = bcf::Writer::from_path(&tmp_vcf, &header, true, bcf::Format::Vcf)
        .map_err(|e| SdError::Vcf(format!("create {}: {e}", tmp_vcf.display())))?;

    let mut bam_reader = bam::IndexedReader::from_path(hp_bam)
        .map_err(|e| SdError::Htslib(format!("open {}: {e}", hp_bam.display())))?;
    if threads > 1 {
        bam_reader.set_threads(threads - 1).map_err(hts_err)?;
    }
    let bam_header = bam_reader.header().to_owned();

    for result in reader.records() {
        let mut rec = result.map_err(|e| SdError::Vcf(format!("read record: {e}")))?;
        let Some(rid) = rec.rid() else {
            writer.translate(&mut rec);
            writer.write(&rec).map_err(vcf_write_err)?;
            continue;
        };
        let chrom = input_header
            .rid2name(rid)
            .map(|n| String::from_utf8_lossy(n).to_string())
            .map_err(|e| SdError::Vcf(format!("rid2name({rid}): {e}")))?;
        let pos0 = rec.pos();
        let alleles: Vec<String> = rec
            .alleles()
            .iter()
            .map(|a| String::from_utf8_lossy(a).to_string())
            .collect();
        if alleles.len() < 2 {
            writer.translate(&mut rec);
            writer.write(&rec).map_err(vcf_write_err)?;
            continue;
        }
        let ref_allele = alleles[0].clone();
        let alts: Vec<String> = alleles[1..].to_vec();

        let supporting = supporting_hp_tags(
            &mut bam_reader,
            &bam_header,
            &chrom,
            pos0,
            &ref_allele,
            &alts,
        )?;
        let hpsup = alts
            .iter()
            .map(|alt| {
                supporting
                    .get(alt)
                    .filter(|tags| !tags.is_empty())
                    .map(|tags| tags.iter().cloned().collect::<Vec<_>>().join(";"))
                    .unwrap_or_else(|| ".".to_string())
            })
            .collect::<Vec<_>>()
            .join(",");

        writer.translate(&mut rec);
        rec.push_format_string(b"HPSUP", &[hpsup.as_bytes()])
            .map_err(|e| SdError::Vcf(format!("set HPSUP: {e}")))?;
        writer.write(&rec).map_err(vcf_write_err)?;
    }
    drop(writer);

    crate::tools::bcftools_sort_index(&tmp_vcf, output_vcf, threads)?;
    let _ = std::fs::remove_file(&tmp_vcf);
    Ok(())
}

fn supporting_hp_tags(
    reader: &mut bam::IndexedReader,
    header: &bam::HeaderView,
    chrom: &str,
    pos0: i64,
    ref_allele: &str,
    alts: &[String],
) -> Result<HashMap<String, BTreeSet<String>>> {
    let Some(tid) = header.tid(chrom.as_bytes()) else {
        return Ok(HashMap::new());
    };
    reader
        .fetch((tid as i32, pos0, pos0 + 1))
        .map_err(|e| SdError::Htslib(format!("fetch {chrom}:{pos0}-{}: {e}", pos0 + 1)))?;

    let alt_set: BTreeSet<&str> = alts.iter().map(String::as_str).collect();
    let mut supporting: HashMap<String, BTreeSet<String>> = HashMap::new();

    for pileup_result in reader.pileup() {
        let pileup = pileup_result.map_err(hts_err)?;
        if pileup.tid() != tid || pileup.pos() as i64 != pos0 {
            continue;
        }
        for aln in pileup.alignments() {
            let record = aln.record();
            if record.mapq() < MIN_MAPQ {
                continue;
            }
            let Some(hp_tag) = hp_tag(&record) else {
                continue;
            };

            if aln.is_del() {
                if ref_allele.len() > 1 {
                    let del_alt = &ref_allele[1..];
                    if alt_set.contains(del_alt) {
                        supporting
                            .entry(del_alt.to_string())
                            .or_default()
                            .insert(hp_tag);
                    }
                }
                continue;
            }

            if let Indel::Ins(len) = aln.indel() {
                let Some(qpos) = aln.qpos() else {
                    continue;
                };
                if !baseq_passes(&record, qpos) {
                    continue;
                }
                let seq = record.seq().as_bytes();
                let end = qpos.saturating_add(len as usize).min(seq.len());
                let ins_seq = String::from_utf8_lossy(&seq[qpos..end]);
                for alt in alts {
                    if alt.starts_with(&ref_allele[..1]) && alt.contains(ins_seq.as_ref()) {
                        supporting.entry(alt.clone()).or_default().insert(hp_tag);
                        break;
                    }
                }
                continue;
            }

            if aln.is_refskip() {
                continue;
            }
            let Some(qpos) = aln.qpos() else {
                continue;
            };
            if !baseq_passes(&record, qpos) {
                continue;
            }
            let seq = record.seq().as_bytes();
            let Some(&base) = seq.get(qpos) else {
                continue;
            };
            let read_base = (base as char).to_string();
            if alt_set.contains(read_base.as_str()) {
                supporting.entry(read_base).or_default().insert(hp_tag);
            }
        }
    }

    Ok(supporting)
}

fn hp_tag(record: &bam::Record) -> Option<String> {
    match record.aux(b"HP").ok()? {
        Aux::String(s) => Some(s.to_string()),
        _ => None,
    }
}

fn baseq_passes(record: &bam::Record, qpos: usize) -> bool {
    record.qual().get(qpos).copied().unwrap_or(0) >= MIN_BASEQ
}

fn temp_vcf_path(output_vcf: &Path) -> PathBuf {
    let name = output_vcf
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("output.vcf.gz");
    output_vcf.with_file_name(format!("{}.{}.hp.tmp.vcf", name, std::process::id()))
}

fn hts_err(e: impl std::fmt::Display) -> SdError {
    SdError::Htslib(e.to_string())
}

fn vcf_write_err(e: impl std::fmt::Display) -> SdError {
    SdError::Vcf(format!("write record: {e}"))
}
