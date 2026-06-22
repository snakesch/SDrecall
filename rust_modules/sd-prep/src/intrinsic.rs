//! Intrinsic alignment — counterpart reference sequences mapped against the RG's
//! masked genome, producing the per-RG intrinsic BAM used by the BILC model.
//!
//! Ports `intrinsic_alignment.py::getIntrinsicBam` (l.14-66) + `getRawseq`
//! (seq.py l.5-22) + `filter_intrinsic_alignments` (l.116-228). minimap2 is the
//! only external algorithm (via the [`crate::minimap`] FFI / minimap2-rs htslib
//! bridge); the total-intrinsic-BAM merge routes through `samtools merge/sort/index`
//! (the leaf-subprocess survivor the T8 design sanctions — BAM merge is NOT
//! re-implemented here).
//!
//! ## Scope note (pass-criterion boundary)
//!
//! The T8 pass criteria are the **RG node groupings**, the **paralog-pair set** and
//! the **masked-genome FASTA md5** — the intrinsic BAM is not one of them. This
//! module produces a correct per-RG intrinsic BAM (query seqs → masked genome,
//! self-location + enclosed/duplicate filtering, secondary→primary promotion) so
//! the driver can assemble the total intrinsic BAM, but its byte-level parity with
//! the Python BAM is not asserted (alignment is the same minimap2 algorithm; record
//! ordering/aux differs).

use minimap2::{ffi::MM_F_EQX, Aligner};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{self, header::HeaderRecord, Read as _};
use sdrecall_utils::{GenomicInterval, Result, SdError};
use std::path::Path;

/// A query sequence for intrinsic alignment: the reference subsequence of a
/// counterpart region, named `{chrom}:{start}-{end}:{strand}` (the qname the
/// Python `getRawseq` produces from the BED, used by `filter_intrinsic_alignments`
/// to detect self-location and enclosure).
struct IntrinsicQuery {
    name: String,
    seq: Vec<u8>,
}

/// Extract counterpart reference sequences for intrinsic alignment — the in-process
/// analog of `getRawseq` (seq.py l.5-22). For each interval in `all_homo_regions`,
/// the reference subsequence `[start, end)` is fetched (uppercased) and named
/// `{chrom}:{start}-{end}:{strand}`. The Python adds `padding` slop + 1000bp
/// windows; here we extract the per-interval sequence directly (the windowing is a
/// performance tiling for minimap2 and does not change the alignment targets).
fn get_raw_seqs(
    all_homo_regions: &[GenomicInterval],
    ref_fa: &Path,
) -> Result<Vec<IntrinsicQuery>> {
    let mut reader = bio::io::fasta::IndexedReader::from_file(&ref_fa).map_err(|e| SdError::Io {
        path: ref_fa.display().to_string(),
        source: std::io::Error::other(e.to_string()),
    })?;
    let mut out = Vec::with_capacity(all_homo_regions.len());
    let mut buf = Vec::new();
    for iv in all_homo_regions {
        if iv.end <= iv.start {
            continue;
        }
        reader
            .fetch(&iv.chrom, iv.start.max(0) as u64, iv.end.max(0) as u64)
            .map_err(|e| SdError::Io {
                path: format!("{}:{}-{}", iv.chrom, iv.start, iv.end),
                source: e,
            })?;
        buf.clear();
        reader.read(&mut buf).map_err(|e| SdError::Io {
            path: format!("{}:{}-{}", iv.chrom, iv.start, iv.end),
            source: e,
        })?;
        buf.make_ascii_uppercase();
        let strand = match iv.strand {
            sdrecall_utils::Strand::Reverse => '-',
            _ => '+',
        };
        out.push(IntrinsicQuery {
            name: format!("{}:{}-{}:{}", iv.chrom, iv.start, iv.end, strand),
            seq: buf.clone(),
        });
    }
    Ok(out)
}

/// Build the per-RG intrinsic BAM — port of `getIntrinsicBam` (l.14-66).
///
/// Extracts counterpart sequences ([`get_raw_seqs`]), maps them against the masked
/// genome with minimap2 `asm20` (via the minimap2-rs htslib bridge → BAM Records),
/// then **remaps the masked-genome (`{chrom}:{start}`, local) alignments back to
/// original-genome coordinates** ([`sdrecall_io::remap_masked_bam_to_genomic`]) and
/// finally applies [`filter_intrinsic_alignments`]. This mirrors Python's
/// `getIntrinsicBam`, which runs the remap *inside* `independent_minimap2_masked`
/// (shell_utils.sh l.472-479) **before** `filter_intrinsic_alignments` — so the
/// self-location check (`start == reference_start + 1`) runs against genomic
/// coordinates (in local coordinates it could never fire). The remap also drops
/// `FLAG ≥ 256` (Python's `$2 < 256`), so the filter — like Python — sees a
/// primary-only BAM. Returns the path to the filtered, indexed intrinsic BAM.
///
/// `all_homo_regions` are the counterpart + query regions (the `*_related_homo`
/// BED); `masked_genome` is the RG's masked FASTA (its `.fai` must exist);
/// `ref_fa` is the reference for sequence extraction + the remap @SQ dictionary;
/// `rg_tag` is the RG label appended to each remapped QNAME (e.g. `RG0`).
pub fn intrinsic_bam(
    all_homo_regions: &[GenomicInterval],
    masked_genome: &Path,
    ref_fa: &Path,
    rg_tag: &str,
    intrinsic_bam: &Path,
) -> Result<()> {
    let queries = get_raw_seqs(all_homo_regions, ref_fa)?;

    // Header SQ dictionary from the masked-genome `.fai`, in `.fai` order.
    let sq = read_masked_sq(masked_genome)?;
    let mut header = bam::Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
    let mut tid_of: ahash::AHashMap<String, i32> = ahash::AHashMap::new();
    for (i, (name, len)) in sq.iter().enumerate() {
        header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", name).push_tag(b"LN", len));
        tid_of.insert(name.clone(), i as i32);
    }

    // Build the minimap2 asm20 aligner indexed on the masked genome.
    let unsorted = format!("{}.unsorted.bam", intrinsic_bam.display());
    if queries.is_empty() {
        // No counterparts → header-only BAM so downstream merge has a valid input.
        {
            let _w = bam::Writer::from_path(&unsorted, &header, bam::Format::Bam)
                .map_err(|e| SdError::Htslib(format!("open {unsorted}: {e}")))?;
        }
    } else {
        let aligner = Aligner::builder()
            .asm20()
            .with_cigar()
            .with_index_threads(1)
            .with_index(masked_genome, None)
            .map_err(|e| {
                SdError::Compute(format!("minimap2 index {}: {e}", masked_genome.display()))
            })?;
        let mut writer = bam::Writer::from_path(&unsorted, &header, bam::Format::Bam)
            .map_err(|e| SdError::Htslib(format!("open {unsorted}: {e}")))?;
        for q in &queries {
            // map(seq, cs=false, md=false, extra_flags=MM_F_EQX, ...) — Python's
            // independent_minimap2_masked always passes `--eqx`, and fp-control's
            // hap-vector encoding requires explicit `=`/`X` CIGAR ops.
            let eqx_flags = [MM_F_EQX as u64];
            let mappings = aligner
                .map(
                    &q.seq,
                    false,
                    false,
                    None,
                    Some(&eqx_flags),
                    Some(q.name.as_bytes()),
                )
                .map_err(|e| SdError::Compute(format!("minimap2 map {}: {e}", q.name)))?;
            for m in &mappings {
                if let Some(rec) = mapping_to_record(q, m, &tid_of)? {
                    writer
                        .write(&rec)
                        .map_err(|e| SdError::Htslib(format!("write record: {e}")))?;
                }
            }
        }
    }

    // Remap masked(local `{chrom}:{start}`) → original-genome coordinates, drop
    // FLAG ≥ 256, append the `:rg_tag` QNAME suffix, then coordinate-sort + index
    // (Python's `independent_minimap2_masked`). This must precede the filter so its
    // self-location check runs in genomic coordinates (matching `getIntrinsicBam`).
    sdrecall_io::remap_masked_bam_to_genomic(
        Path::new(&unsorted),
        ref_fa,
        rg_tag,
        intrinsic_bam,
    )?;
    let _ = std::fs::remove_file(&unsorted);

    // Filter self-location + enclosed/duplicate intervals (genomic, primary-only).
    filter_intrinsic_alignments(intrinsic_bam)?;
    samtools_index(intrinsic_bam)?;
    Ok(())
}

/// Read the masked-genome `.fai` into `(name, length)` SQ entries (in `.fai`
/// order), the header dictionary for the intrinsic BAM.
fn read_masked_sq(masked_genome: &Path) -> Result<Vec<(String, i64)>> {
    let fai = format!("{}.fai", masked_genome.display());
    let text = std::fs::read_to_string(&fai).map_err(|e| SdError::Io {
        path: fai.clone(),
        source: e,
    })?;
    let mut out = Vec::new();
    for line in text.lines() {
        let mut f = line.split('\t');
        if let (Some(name), Some(len)) = (f.next(), f.next()) {
            if let Ok(len) = len.parse::<i64>() {
                out.push((name.to_string(), len));
            }
        }
    }
    Ok(out)
}

/// Convert a minimap2 [`minimap2::Mapping`] into a `rust_htslib::bam::Record`,
/// reproducing the SAM fields the intrinsic BAM needs (qname, flag, tid, pos,
/// mapq, CIGAR with soft-clips, the query sequence on the mapping strand). Returns
/// `None` for mappings whose target contig is absent from the header (should not
/// happen — the index IS the masked genome) or with no CIGAR.
fn mapping_to_record(
    q: &IntrinsicQuery,
    m: &minimap2::Mapping,
    tid_of: &ahash::AHashMap<String, i32>,
) -> Result<Option<bam::Record>> {
    let target = match &m.target_name {
        Some(t) => t.as_str(),
        None => return Ok(None),
    };
    let tid = match tid_of.get(target) {
        Some(&t) => t,
        None => return Ok(None),
    };
    let aln = match &m.alignment {
        Some(a) => a,
        None => return Ok(None),
    };
    let cigar_ops = match &aln.cigar {
        Some(c) => c,
        None => return Ok(None),
    };

    // The aligned sequence is the query, reverse-complemented if on the minus
    // strand (so it matches the forward reference, SAM convention).
    let reverse = matches!(m.strand, minimap2::Strand::Reverse);
    let seq: Vec<u8> = if reverse {
        bio::alphabets::dna::revcomp(&q.seq)
    } else {
        q.seq.clone()
    };
    let qlen = seq.len() as i32;
    // query_start/query_end are on the original query orientation; convert to the
    // read orientation for soft-clip lengths.
    let (clip_front, clip_back) = if reverse {
        (qlen - m.query_end, m.query_start)
    } else {
        (m.query_start, qlen - m.query_end)
    };

    // Build the CIGAR: leading/trailing soft-clips + the core ops. With
    // MM_F_EQX, minimap2 emits 7/8 for =/X instead of 0/M.
    let mut cigar: Vec<Cigar> = Vec::with_capacity(cigar_ops.len() + 2);
    if clip_front > 0 {
        cigar.push(Cigar::SoftClip(clip_front as u32));
    }
    for &(len, op) in cigar_ops {
        cigar.push(minimap2_cigar_op(len, op)?);
    }
    if clip_back > 0 {
        cigar.push(Cigar::SoftClip(clip_back as u32));
    }
    let cigar = CigarString(cigar);

    let mut rec = bam::Record::new();
    let qual = vec![255u8; seq.len()];
    rec.set(q.name.as_bytes(), Some(&cigar), &seq, &qual);
    rec.set_tid(tid);
    rec.set_pos(m.target_start as i64);
    rec.set_mapq(m.mapq.min(255) as u8);
    rec.set_mtid(-1);
    rec.set_mpos(-1);
    rec.set_insert_size(0);
    // flags: reverse strand (0x10), secondary (0x100) for non-primary,
    // supplementary (0x800).
    let mut flag: u16 = 0;
    if reverse {
        flag |= 0x10;
    }
    if !m.is_primary && !m.is_supplementary {
        flag |= 0x100;
    }
    if m.is_supplementary {
        flag |= 0x800;
    }
    rec.set_flags(flag);
    Ok(Some(rec))
}

fn minimap2_cigar_op(len: u32, op: u8) -> Result<Cigar> {
    match op {
        0 => Ok(Cigar::Match(len)),
        1 => Ok(Cigar::Ins(len)),
        2 => Ok(Cigar::Del(len)),
        3 => Ok(Cigar::RefSkip(len)),
        4 => Ok(Cigar::SoftClip(len)),
        5 => Ok(Cigar::HardClip(len)),
        6 => Ok(Cigar::Pad(len)),
        7 => Ok(Cigar::Equal(len)),
        8 => Ok(Cigar::Diff(len)),
        _ => Err(SdError::Compute(format!(
            "unexpected minimap2 CIGAR op {op}"
        ))),
    }
}

/// Per-distinct-interval enclosure status — port of `compute_interval_status`
/// (intrinsic_alignment.py l.70-112). An interval `(chrom,start,end)` parsed from a
/// qname is `false` (disallowed) when a larger interval on the same chrom encloses
/// it; for identical/contained intervals only the first (by start asc, end desc) is
/// allowed.
fn compute_interval_status(
    records: &[(String, i64, i64)],
) -> ahash::AHashMap<(String, i64, i64), bool> {
    use std::collections::BTreeMap;
    // intervals per chrom (dedup via set).
    let mut by_chr: BTreeMap<String, std::collections::BTreeSet<(i64, i64)>> = BTreeMap::new();
    for (c, s, e) in records {
        by_chr.entry(c.clone()).or_default().insert((*s, *e));
    }
    let mut allowed = ahash::AHashMap::new();
    for (chrom, set) in by_chr {
        // sort by start asc, end desc.
        let mut sorted: Vec<(i64, i64)> = set.into_iter().collect();
        sorted.sort_by(|a, b| a.0.cmp(&b.0).then(b.1.cmp(&a.1)));
        let mut max_end = -1i64;
        for (i, (s, e)) in sorted.iter().enumerate() {
            if i == 0 {
                allowed.insert((chrom.clone(), *s, *e), true);
                max_end = *e;
            } else if *e <= max_end {
                allowed.insert((chrom.clone(), *s, *e), false);
            } else {
                allowed.insert((chrom.clone(), *s, *e), true);
                max_end = *e;
            }
        }
    }
    allowed
}

/// Parse a qname of shape `chrom:start-end(:label)?` into `(chrom, start, end)`.
fn parse_qname_interval(qname: &str) -> Option<(String, i64, i64)> {
    // Match the LAST ":" group as optional label: take up to the coord block.
    // Format: <chrom>:<start>-<end>[:<label>]. chrom may itself contain no ':'
    // in this pipeline (contig names like chr1), so split on ':' then '-'.
    let mut parts = qname.splitn(3, ':');
    let chrom = parts.next()?;
    let coords = parts.next()?;
    let (s, e) = coords.split_once('-')?;
    let start = s.parse::<i64>().ok()?;
    let end = e.parse::<i64>().ok()?;
    Some((chrom.to_string(), start, end))
}

/// Filter intrinsic alignments — port of `filter_intrinsic_alignments`
/// (intrinsic_alignment.py l.116-228). Drops reads mapping to their own genomic
/// location (`start == reference_start + 1`), enclosed/duplicate intervals, and
/// promotes the relevant secondary alignments to primary. Rewrites `bam` in place.
pub fn filter_intrinsic_alignments(bam: &Path) -> Result<()> {
    // ---- pass 1: collect intervals for enclosure status ----
    let mut intervals: Vec<(String, i64, i64)> = Vec::new();
    {
        let mut reader =
            bam::Reader::from_path(bam).map_err(|e| SdError::Htslib(format!("open {}: {e}", bam.display())))?;
        let mut rec = bam::Record::new();
        while let Some(r) = reader.read(&mut rec) {
            r.map_err(|e| SdError::Htslib(format!("read: {e}")))?;
            if rec.is_unmapped() || rec.seq_len() == 0 {
                continue;
            }
            if let Ok(qn) = std::str::from_utf8(rec.qname()) {
                if let Some(iv) = parse_qname_interval(qn) {
                    intervals.push(iv);
                }
            }
        }
    }
    let allowed = compute_interval_status(&intervals);

    // ---- pass 2: filter + secondary→primary promotion ----
    let tmp = format!("{}.filtered.bam", bam.display());
    {
        let reader =
            bam::Reader::from_path(bam).map_err(|e| SdError::Htslib(format!("open {}: {e}", bam.display())))?;
        let header = bam::Header::from_template(reader.header());
        let mut writer = bam::Writer::from_path(&tmp, &header, bam::Format::Bam)
            .map_err(|e| SdError::Htslib(format!("open {tmp}: {e}")))?;
        let mut reader =
            bam::Reader::from_path(bam).map_err(|e| SdError::Htslib(format!("open {}: {e}", bam.display())))?;

        let mut seen_intervals: ahash::AHashSet<(String, i64, i64)> = ahash::AHashSet::new();
        let mut primary_origin: ahash::AHashSet<String> = ahash::AHashSet::new();
        let mut sec_to_pri: ahash::AHashSet<String> = ahash::AHashSet::new();
        let mut written_qnames: ahash::AHashSet<String> = ahash::AHashSet::new();
        let mut buffer_sec: ahash::AHashMap<String, bam::Record> = ahash::AHashMap::new();

        let mut rec = bam::Record::new();
        while let Some(r) = reader.read(&mut rec) {
            r.map_err(|e| SdError::Htslib(format!("read: {e}")))?;
            if rec.is_unmapped() || rec.seq_len() == 0 {
                continue;
            }
            let qname = match std::str::from_utf8(rec.qname()) {
                Ok(q) => q.to_string(),
                Err(_) => continue,
            };
            if written_qnames.contains(&qname) {
                continue;
            }
            if let Some(iv) = parse_qname_interval(&qname) {
                if !allowed.get(&iv).copied().unwrap_or(true) {
                    continue;
                }
                if seen_intervals.contains(&iv) {
                    continue;
                }
                seen_intervals.insert(iv.clone());

                if primary_origin.contains(&qname) {
                    if rec.is_supplementary() {
                        continue;
                    } else if rec.is_secondary() {
                        sec_to_pri.insert(qname.clone());
                        let new_flag = rec.flags() & !0x100;
                        rec.set_flags(new_flag);
                        writer
                            .write(&rec)
                            .map_err(|e| SdError::Htslib(format!("write: {e}")))?;
                        continue;
                    }
                } else if rec.is_supplementary() {
                    continue;
                } else if rec.is_secondary() {
                    buffer_sec.insert(qname.clone(), rec.clone());
                    continue;
                }
                // self-location check: start_val == reference_start + 1.
                if iv.1 != rec.reference_start() + 1 {
                    writer
                        .write(&rec)
                        .map_err(|e| SdError::Htslib(format!("write: {e}")))?;
                    written_qnames.insert(qname.clone());
                } else {
                    primary_origin.insert(qname.clone());
                }
            } else {
                writer
                    .write(&rec)
                    .map_err(|e| SdError::Htslib(format!("write: {e}")))?;
                written_qnames.insert(qname.clone());
            }
        }
        // promote buffered secondaries whose primary mapped to its own location.
        let promote: Vec<String> = primary_origin.difference(&sec_to_pri).cloned().collect();
        for qname in promote {
            if let Some(mut sec) = buffer_sec.remove(&qname) {
                let new_flag = sec.flags() & !0x100;
                sec.set_flags(new_flag);
                writer
                    .write(&sec)
                    .map_err(|e| SdError::Htslib(format!("write: {e}")))?;
            }
        }
    }
    std::fs::rename(&tmp, bam).map_err(|e| SdError::Io {
        path: bam.display().to_string(),
        source: e,
    })?;
    Ok(())
}

/// Coordinate-sort `input` → `output` via `samtools sort` (sanctioned leaf
/// subprocess; BAM sort is NOT re-implemented).
fn samtools_sort(input: &Path, output: &Path) -> Result<()> {
    let out = std::process::Command::new("samtools")
        .args([
            "sort",
            "-O",
            "bam",
            "-o",
            &output.display().to_string(),
            &input.display().to_string(),
        ])
        .output()
        .map_err(|e| SdError::Io {
            path: "samtools sort".to_string(),
            source: e,
        })?;
    if !out.status.success() {
        return Err(SdError::Compute(format!(
            "samtools sort failed: {}",
            String::from_utf8_lossy(&out.stderr)
        )));
    }
    Ok(())
}

/// Index a coordinate-sorted BAM via `samtools index`.
fn samtools_index(bam: &Path) -> Result<()> {
    let out = std::process::Command::new("samtools")
        .args(["index", &bam.display().to_string()])
        .output()
        .map_err(|e| SdError::Io {
            path: "samtools index".to_string(),
            source: e,
        })?;
    if !out.status.success() {
        return Err(SdError::Compute(format!(
            "samtools index failed: {}",
            String::from_utf8_lossy(&out.stderr)
        )));
    }
    Ok(())
}

/// Merge per-RG intrinsic BAMs (+ the raw intrinsic BAM if present) into the total
/// intrinsic BAM, then sort + index — the analog of build_beds l.96-113. Routes
/// through `samtools merge/sort/index` (the leaf-subprocess survivor the T8 design
/// sanctions; BAM merge is NOT re-implemented). The header is taken from the first
/// input (Python rebuilds SQ lines from the reference via `modify_bam_sq_lines`;
/// here we let `samtools merge` reconcile headers).
pub fn merge_total_intrinsic_bam(inputs: &[&Path], total: &Path) -> Result<()> {
    if inputs.is_empty() {
        return Err(SdError::Compute("no intrinsic BAMs to merge".to_string()));
    }
    let merged = format!("{}.merged.bam", total.display());
    let mut args: Vec<String> = vec!["merge".into(), "-f".into(), merged.clone()];
    for p in inputs {
        args.push(p.display().to_string());
    }
    let out = std::process::Command::new("samtools")
        .args(&args)
        .output()
        .map_err(|e| SdError::Io {
            path: "samtools merge".to_string(),
            source: e,
        })?;
    if !out.status.success() {
        return Err(SdError::Compute(format!(
            "samtools merge failed: {}",
            String::from_utf8_lossy(&out.stderr)
        )));
    }
    samtools_sort(Path::new(&merged), total)?;
    let _ = std::fs::remove_file(&merged);
    samtools_index(total)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_qname_interval_basic() {
        assert_eq!(
            parse_qname_interval("chr1:1000-2000:+"),
            Some(("chr1".to_string(), 1000, 2000))
        );
        assert_eq!(
            parse_qname_interval("chr12:63557340-63559211"),
            Some(("chr12".to_string(), 63557340, 63559211))
        );
        assert_eq!(parse_qname_interval("nocoords"), None);
    }

    #[test]
    fn minimap2_cigar_op_preserves_eqx_ops() {
        assert!(matches!(minimap2_cigar_op(5, 7).unwrap(), Cigar::Equal(5)));
        assert!(matches!(minimap2_cigar_op(2, 8).unwrap(), Cigar::Diff(2)));
        assert!(matches!(minimap2_cigar_op(3, 1).unwrap(), Cigar::Ins(3)));
        assert!(minimap2_cigar_op(1, 9).is_err());
    }

    #[test]
    fn minimap2_eqx_flag_emits_equal_or_diff_ops() {
        let reference = b"ACGGTAGAGAGGAAGAAGAAGGAATAGCGGACTTGTGTATTTTATCGTCATTCGTGGTTATCATATAGTTTATTGATTTGAAGACTACGTAAGTAATTTGAGGACTGATTAAAATTTTCTTTTTTAGCTTAGAGTCAATTAAAGAGGGCAAAATTTTCTCAAAAGACCATGGTGCATATGACGATAGCTTTAGTAGTATGGATTGGGCTCTTCTTTCATGGATGTTATTCAGAAGGAGTGATATATCGAGGTGTTTGAAACACCAGCGACACCAGAAGGCTGTGGATGTTAAATCGTAGAACCTATAGACGAGTTCTAAAATATACTTTGGGGTTTTCAGCGATGCAAAA";
        let mut query = reference.to_vec();
        query[37] = if query[37] == b'A' { b'C' } else { b'A' };
        let aligner = Aligner::builder()
            .asm20()
            .with_cigar()
            .with_seq(reference)
            .expect("build in-memory minimap2 index");
        let eqx_flags = [MM_F_EQX as u64];
        let mappings = aligner
            .map(&query, false, false, None, Some(&eqx_flags), Some(b"query"))
            .expect("map query");
        let cigar = mappings
            .iter()
            .find_map(|m| m.alignment.as_ref()?.cigar.as_ref())
            .expect("mapping with CIGAR");
        assert!(cigar.iter().any(|(_, op)| *op == 7 || *op == 8));
        assert!(!cigar.iter().any(|(_, op)| *op == 0));
    }

    #[test]
    fn enclosure_status_marks_enclosed_false() {
        // [100,1000) encloses [200,500); the larger is allowed, the enclosed not.
        let recs = vec![
            ("chr1".to_string(), 100, 1000),
            ("chr1".to_string(), 200, 500),
            ("chr1".to_string(), 1100, 1500), // disjoint → allowed
        ];
        let allowed = compute_interval_status(&recs);
        assert_eq!(allowed.get(&("chr1".to_string(), 100, 1000)), Some(&true));
        assert_eq!(allowed.get(&("chr1".to_string(), 200, 500)), Some(&false));
        assert_eq!(allowed.get(&("chr1".to_string(), 1100, 1500)), Some(&true));
    }

    #[test]
    fn enclosure_identical_intervals_keep_first_only() {
        // Two identical intervals: first allowed, the second (same end) enclosed.
        let recs = vec![
            ("chrX".to_string(), 10, 50),
            ("chrX".to_string(), 10, 50),
        ];
        let allowed = compute_interval_status(&recs);
        // The set dedups identical intervals → only one entry, allowed.
        assert_eq!(allowed.get(&("chrX".to_string(), 10, 50)), Some(&true));
    }
}
