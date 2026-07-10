use rust_htslib::bam::record::Aux;
use rust_htslib::{bam, bam::Read};
use std::collections::{HashMap, HashSet};
use std::fs::{create_dir_all, File};
use std::io::{BufWriter, Write};
use std::path::Path;

/// Read BED file and return regions as Vec<(chr, start, end)>.
fn read_bed_regions(bed_path: &str) -> anyhow::Result<Vec<(String, u64, u64)>> {
    use std::io::{BufRead, BufReader};

    let file = File::open(bed_path)
        .map_err(|e| anyhow::anyhow!("Failed to open BED file {bed_path}: {e}"))?;
    let reader = BufReader::new(file);

    let mut regions = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 3 {
            let chr = parts[0].to_string();
            let start: u64 = parts[1]
                .parse()
                .map_err(|_| anyhow::anyhow!("Invalid start position: {}", parts[1]))?;
            let end: u64 = parts[2]
                .parse()
                .map_err(|_| anyhow::anyhow!("Invalid end position: {}", parts[2]))?;
            regions.push((chr, start, end));
        }
    }
    Ok(regions)
}

/// Sort + merge overlapping/adjacent regions (book-ended intervals fuse, like
/// `bedtools sort | merge`). Collapsing overlapping input rows means a read pair
/// is not fetched (and written) once per overlapping row, and it cuts BAM seeks.
/// Straddling pairs (R1/R2 in two disjoint regions) are deduped separately by the
/// written-qname set in [`bam_to_fastq`].
fn merge_regions(mut regions: Vec<(String, u64, u64)>) -> Vec<(String, u64, u64)> {
    regions.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
    let mut merged: Vec<(String, u64, u64)> = Vec::with_capacity(regions.len());
    for (chrom, start, end) in regions {
        match merged.last_mut() {
            Some((c, _s, e)) if *c == chrom && start <= *e => {
                if end > *e {
                    *e = end;
                }
            }
            _ => merged.push((chrom, start, end)),
        }
    }
    merged
}

/// Maximum exclusive `|AS - XS|` gap below which a read is treated as a
/// near-tied multi-mapper. BWA folds multi-mapping ambiguity into MAPQ, but
/// the raw AS/XS gap is the direct measure: a small gap means the read's best
/// and second-best alignments are essentially interchangeable, so it is a
/// genuine SD multi-mapper that may carry the alt haplotype from a paralog.
const NFC_AS_XS_MAX_GAP: i64 = 10;

/// Read an integer aux tag as `i64`. Returns `None` for a missing tag or a
/// non-integer (string/array) variant, so a caller can treat absence as
/// "not a near-tied multi-mapper" rather than erroring out of extraction.
fn aux_integer(record: &bam::Record, tag: &[u8]) -> Option<i64> {
    match record.aux(tag).ok()? {
        Aux::I8(v) => Some(v as i64),
        Aux::U8(v) => Some(v as i64),
        Aux::I16(v) => Some(v as i64),
        Aux::U16(v) => Some(v as i64),
        Aux::I32(v) => Some(v as i64),
        Aux::U32(v) => Some(v as i64),
        _ => None,
    }
}

/// `|AS - XS| < NFC_AS_XS_MAX_GAP` — the read's best and second-best
/// alignments are near-tied, so it is a genuine SD multi-mapper. Both tags
/// must be present; minimap2 BAMs (which omit `XS`) simply return `false`.
fn near_tied_multimapper(record: &bam::Record) -> bool {
    match (aux_integer(record, b"AS"), aux_integer(record, b"XS")) {
        (Some(as_score), Some(xs_score)) => (as_score - xs_score).abs() < NFC_AS_XS_MAX_GAP,
        _ => false,
    }
}

/// Check if a read carries NFC evidence.
fn should_include_read(record: &bam::Record, multi_aligned: bool) -> bool {
    if !multi_aligned {
        // Python FC extraction:
        // `samtools view -h -P -L {region_bed} -u {input_bam} | bamtofastq`
        // No MAPQ/tag predicate is applied here.
        return true;
    }

    // NFC extraction recruits non-split reads with explicit alternate
    // alignments (`XA`) or near-tied best/second-best alignment scores. MAPQ is
    // deliberately not a recruitment arm here: it is a downstream confidence
    // summary, while AS/XS is the direct evidence of SD alignment ambiguity.
    if record.aux(b"SA").is_ok() {
        return false;
    }
    record.aux(b"XA").is_ok() || near_tied_multimapper(record)
}

/// Read-pair level extraction rule. If either mate carries NFC evidence, both
/// mates are emitted to preserve the paired FASTQ unit for realignment.
fn should_include_pair(r1: &bam::Record, r2: &bam::Record, multi_aligned: bool) -> bool {
    should_include_read(r1, multi_aligned) || should_include_read(r2, multi_aligned)
}

/// Encode BAM Phred qualities as a Sanger-FASTQ quality line.
///
/// Valid Phred scores are `0..=93` (ASCII `!`..=`~`). `0xFF` is htslib's
/// "qualities unavailable" sentinel and any value `> 93` is not representable in
/// Sanger FASTQ; either is a hard error rather than a wrapped/invented character,
/// because this FASTQ feeds realignment and a silent `0xFF + 33` wrap would emit
/// corrupt quality scores (and panics in debug builds).
fn quality_to_string(qual: &[u8], qname: &str) -> anyhow::Result<String> {
    let mut out = String::with_capacity(qual.len());
    for &q in qual {
        if q > 93 {
            return Err(anyhow::anyhow!(
                "read {qname}: base quality {q} is not encodable as Sanger FASTQ \
                 (expected 0..=93; 0xFF means the BAM has no stored qualities)"
            ));
        }
        out.push((q + 33) as char);
    }
    Ok(out)
}

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
//  Pure-Rust API (no PyO3 dependency)
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

/// Extract paired reads from a BAM file that overlap the given BED regions,
/// writing R1/R2 FASTQ files. Returns `(r1_path, r2_path)`.
///
/// When `multi_aligned` is true, applies the multi-alignment filter:
/// `!SA && (XA || |AS - XS| < 10)`. When false, every fetched read pair is eligible.
/// Singleton reads (only one mate found) are discarded; if either mate of a pair
/// passes, both are written.
pub fn bam_to_fastq(
    input_bam: &str,
    region_bed: &str,
    output_freads: &str,
    output_rreads: &str,
    multi_aligned: bool,
    threads: usize,
) -> anyhow::Result<(String, String)> {
    if let Some(parent) = Path::new(output_freads).parent() {
        create_dir_all(parent)?;
    }
    if let Some(parent) = Path::new(output_rreads).parent() {
        create_dir_all(parent)?;
    }

    // Merge overlapping/adjacent input rows so a pair is not extracted once per
    // overlapping row (and to cut BAM seeks).
    let regions = merge_regions(read_bed_regions(region_bed)?);

    let mut bam_reader = bam::IndexedReader::from_path(input_bam)
        .map_err(|e| anyhow::anyhow!("Failed to open BAM {input_bam}: {e}"))?;
    let header = bam_reader.header().clone();
    if threads > 1 {
        bam_reader.set_threads(threads - 1)?;
    }

    let r1_file = File::create(output_freads)?;
    let r2_file = File::create(output_rreads)?;
    let mut r1_writer = BufWriter::new(r1_file);
    let mut r2_writer = BufWriter::new(r2_file);

    // Global dedup across regions: a pair straddling two disjoint regions (R1 in
    // one, R2 in another) is reconstructed in BOTH via the mate-fetch and would
    // otherwise be written twice, inflating realignment depth. Emit each qname once.
    let mut written: HashSet<Vec<u8>> = HashSet::new();

    for (chr, start, end) in &regions {
        let tid = header
            .tid(chr.as_bytes())
            .ok_or_else(|| anyhow::anyhow!("Chromosome {chr} not found in BAM"))?;

        bam_reader.fetch((tid, *start as i64, *end as i64))?;

        let mut read_pairs: HashMap<Vec<u8>, (Option<bam::Record>, Option<bam::Record>)> =
            HashMap::new();

        for result in bam_reader.records() {
            let record = result?;
            let qname = record.qname().to_vec();
            let entry = read_pairs.entry(qname).or_insert((None, None));
            if record.is_first_in_template() {
                entry.0 = Some(record);
            } else {
                entry.1 = Some(record);
            }
        }

        // Fetch mates for singletons.
        for entry in read_pairs.values_mut() {
            if let (Some(read), None) | (None, Some(read)) = entry {
                if read.is_paired() && !read.is_mate_unmapped() {
                    let mtid = read.mtid();
                    let mpos = read.mpos();
                    let window_start = (mpos - 5).max(0);
                    let window_end = mpos + 5 + read.seq_len() as i64;
                    bam_reader.fetch((mtid, window_start, window_end))?;
                    let mut mate = bam::Record::new();
                    while let Some(Ok(())) = bam_reader.read(&mut mate) {
                        if mate.qname() == read.qname()
                            && mate.is_first_in_template() != read.is_first_in_template()
                        {
                            if read.is_first_in_template() {
                                entry.1 = Some(mate.clone());
                            } else {
                                entry.0 = Some(mate.clone());
                            }
                            break;
                        }
                    }
                }
            }
        }

        // Write pairs where at least one mate passes the filter.
        for (r1_opt, r2_opt) in read_pairs.into_values() {
            let passes = match (&r1_opt, &r2_opt) {
                (Some(r1), Some(r2)) => should_include_pair(r1, r2, multi_aligned),
                _ => false,
            };
            if passes {
                if let (Some(r1), Some(r2)) = (r1_opt, r2_opt) {
                    // Skip if this qname was already emitted from an earlier region.
                    if written.insert(r1.qname().to_vec()) {
                        write_fastq_record(&mut r1_writer, &r1)?;
                        write_fastq_record(&mut r2_writer, &r2)?;
                    }
                }
            }
        }
    }

    r1_writer.flush()?;
    r2_writer.flush()?;
    Ok((output_freads.to_string(), output_rreads.to_string()))
}

fn write_fastq_record(w: &mut impl Write, rec: &bam::Record) -> anyhow::Result<()> {
    let name = std::str::from_utf8(rec.qname())
        .map_err(|e| anyhow::anyhow!("Read name is not valid UTF-8: {e}"))?;
    // Encode sequence/quality first so an invalid score errors before any
    // partial record is written to the buffer. BAM stores reverse-strand
    // alignments in mapping orientation; FASTQ must restore read orientation.
    let seq = fastq_sequence(rec);
    let qual = fastq_quality_string(rec, name)?;
    writeln!(w, "@{name}")?;
    writeln!(w, "{}", String::from_utf8_lossy(&seq))?;
    writeln!(w, "+")?;
    writeln!(w, "{qual}")?;
    Ok(())
}

fn fastq_sequence(rec: &bam::Record) -> Vec<u8> {
    let seq = rec.seq().as_bytes();
    if rec.is_reverse() {
        seq.iter().rev().map(|&b| complement_base(b)).collect()
    } else {
        seq
    }
}

fn fastq_quality_string(rec: &bam::Record, qname: &str) -> anyhow::Result<String> {
    if rec.is_reverse() {
        let reversed: Vec<u8> = rec.qual().iter().rev().copied().collect();
        quality_to_string(&reversed, qname)
    } else {
        quality_to_string(rec.qual(), qname)
    }
}

fn complement_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'D' => b'H',
        b'H' => b'D',
        b'V' => b'B',
        b'N' => b'N',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        b'r' => b'y',
        b'y' => b'r',
        b's' => b's',
        b'w' => b'w',
        b'k' => b'm',
        b'm' => b'k',
        b'b' => b'v',
        b'd' => b'h',
        b'h' => b'd',
        b'v' => b'b',
        b'n' => b'n',
        _ => b'N',
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Aux, Cigar, CigarString};

    #[test]
    fn quality_encodes_valid_phred_range() {
        // 0 → '!' (33), 40 → 'I' (73), 93 → '~' (126): the Sanger FASTQ bounds.
        assert_eq!(quality_to_string(&[0, 40, 93], "r1").unwrap(), "!I~");
        assert_eq!(quality_to_string(&[], "r1").unwrap(), "");
    }

    #[test]
    fn quality_rejects_missing_sentinel() {
        // 0xFF is htslib's "qualities absent" marker — must error (not wrap to ' ').
        let err = quality_to_string(&[30, 0xFF, 30], "read42").unwrap_err();
        assert!(
            err.to_string().contains("read42"),
            "error should name the offending read: {err}"
        );
    }

    #[test]
    fn quality_rejects_out_of_range() {
        // 94 is one past the Sanger cap and would still wrap-corrupt downstream.
        assert!(quality_to_string(&[94], "r1").is_err());
    }

    #[test]
    fn reverse_strand_fastq_restores_read_orientation() {
        let mut rec = bam::Record::new();
        let cigar = CigarString(vec![Cigar::Match(5)]);
        rec.set(b"read1", Some(&cigar), b"ACGTN", &[10, 20, 30, 40, 50]);
        rec.set_flags(0x10);

        assert_eq!(fastq_sequence(&rec), b"NACGT");
        assert_eq!(fastq_quality_string(&rec, "read1").unwrap(), "SI?5+");
    }

    fn read_with_tags(
        mapq: u8,
        xa: bool,
        sa: bool,
        as_score: Option<i32>,
        xs_score: Option<i32>,
    ) -> bam::Record {
        let mut rec = bam::Record::new();
        let cigar = CigarString(vec![Cigar::Match(1)]);
        rec.set(b"read1", Some(&cigar), b"A", &[30]);
        rec.set_mapq(mapq);
        if xa {
            rec.push_aux(b"XA", Aux::String("chr1,+100,1M,0;")).unwrap();
        }
        if sa {
            rec.push_aux(b"SA", Aux::String("chr1,100,+,1M,60,0;"))
                .unwrap();
        }
        if let Some(v) = as_score {
            rec.push_aux(b"AS", Aux::I32(v)).unwrap();
        }
        if let Some(v) = xs_score {
            rec.push_aux(b"XS", Aux::I32(v)).unwrap();
        }
        rec
    }

    #[test]
    fn fc_extraction_keeps_high_mapq_reads() {
        let rec = read_with_tags(60, false, false, None, None);
        assert!(should_include_read(&rec, false));
    }

    #[test]
    fn nfc_extraction_uses_xa_or_near_tied_scores_not_mapq() {
        assert!(should_include_read(
            &read_with_tags(60, true, false, None, None),
            true
        ));
        assert!(!should_include_read(
            &read_with_tags(49, false, false, None, None),
            true
        ));
        assert!(!should_include_read(
            &read_with_tags(50, false, false, None, None),
            true
        ));
        assert!(!should_include_read(
            &read_with_tags(10, true, true, None, None),
            true
        ));
    }

    #[test]
    fn nfc_extraction_recruits_near_tied_multimappers() {
        // High MAPQ, no XA, but |AS - XS| = 5 < 10 marks it as a genuine SD
        // multi-mapper to keep.
        assert!(should_include_read(
            &read_with_tags(60, false, false, Some(148), Some(143)),
            true
        ));
        // |AS - XS| = 0 (perfect tie, e.g. BWA AS:i:148 XS:i:148) is kept.
        assert!(should_include_read(
            &read_with_tags(60, false, false, Some(148), Some(148)),
            true
        ));
        // |AS - XS| = 10 (exactly the cap) is dropped; the bound is strict.
        assert!(!should_include_read(
            &read_with_tags(60, false, false, Some(150), Some(140)),
            true
        ));
        // |AS - XS| = 11 (over the cap), high MAPQ, no XA -> dropped.
        assert!(!should_include_read(
            &read_with_tags(60, false, false, Some(150), Some(139)),
            true
        ));
        // Missing XS (e.g. a minimap2 BAM) → no AS/XS arm, falls back to the
        // XA arm; high MAPQ + no XA -> dropped (no false recruit).
        assert!(!should_include_read(
            &read_with_tags(60, false, false, Some(148), None),
            true
        ));
    }

    #[test]
    fn nfc_as_xs_arm_does_not_override_split_alignments() {
        // A near-tied read that also carries SA (split alignment) is still
        // excluded — the `![SA]` gate wins, matching the Python predicate.
        assert!(!should_include_read(
            &read_with_tags(60, false, true, Some(148), Some(148)),
            true
        ));
    }

    #[test]
    fn nfc_extraction_is_pair_level() {
        let plain_mate = read_with_tags(60, false, false, None, None);
        let near_tied_mate = read_with_tags(60, false, false, Some(148), Some(143));
        assert!(should_include_pair(&plain_mate, &near_tied_mate, true));

        let low_mapq_only_mate = read_with_tags(1, false, false, None, None);
        assert!(!should_include_pair(&plain_mate, &low_mapq_only_mate, true));
    }

    #[test]
    fn merge_regions_collapses_overlaps_and_sorts() {
        let regions = vec![
            ("chr1".to_string(), 50, 100),
            ("chr1".to_string(), 10, 60), // overlaps [50,100) → [10,100)
            ("chr1".to_string(), 100, 120), // book-ended → fuse to [10,120)
            ("chr2".to_string(), 5, 9),
        ];
        assert_eq!(
            merge_regions(regions),
            vec![("chr1".to_string(), 10, 120), ("chr2".to_string(), 5, 9)]
        );
    }
}

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
//  PyO3 bindings (only with `python` feature)
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

#[cfg(feature = "python")]
mod python_bindings {
    use pyo3::exceptions::PyRuntimeError;
    use pyo3::prelude::*;

    #[pyfunction]
    #[pyo3(signature = (input_bam, region_bed, output_freads, output_rreads, multi_aligned=false, threads=1, _tmp_dir="/tmp"))]
    fn bam_to_fastq_biobambam(
        input_bam: &str,
        region_bed: &str,
        output_freads: &str,
        output_rreads: &str,
        multi_aligned: bool,
        threads: usize,
        _tmp_dir: &str,
    ) -> PyResult<(String, String)> {
        crate::bam_to_fastq(
            input_bam,
            region_bed,
            output_freads,
            output_rreads,
            multi_aligned,
            threads,
        )
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    #[pymodule]
    fn rust_read_extraction(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(bam_to_fastq_biobambam, m)?)?;
        Ok(())
    }
}
