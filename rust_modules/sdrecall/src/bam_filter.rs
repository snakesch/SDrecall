//! In-process BAM filtering by qname set — the post-fp-control step that
//! keeps only correctly-aligned read pairs in the per-island clean BAM.
//!
//! Replaces the Python filtering loop in
//! `fp_control/realign_filter_per_cov.py:375-437`.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use rust_htslib::{bam, bam::record::Aux, bam::Read};
use sdrecall_utils::{Result, SdError};

const CLEAN_MIN_MAPQ: u8 = 10;

pub(crate) struct HpTagAssignments<'a> {
    pub correct_qnames: &'a HashSet<String>,
    pub mismap_qnames: &'a HashSet<String>,
    pub lowqual_qnames: &'a HashSet<String>,
    pub qname_hap: &'a HashMap<String, i32>,
}

/// Write a new BAM containing only primary mapped reads whose qname is in
/// `keep_qnames`, excluding Python's duplicate/QC-fail/secondary/supplementary
/// cases. Each retained read is annotated with its fp-control HP tag.
///
/// The output BAM is coordinate-sorted and indexed. The returned count is the
/// number of alignments written before sorting.
pub fn filter_bam_by_qnames(
    input_bam: &Path,
    keep_qnames: &HashSet<String>,
    lowqual_qnames: &HashSet<String>,
    qname_hap: &HashMap<String, i32>,
    output_bam: &Path,
    chunk_id: usize,
    threads: usize,
) -> Result<usize> {
    let mut reader = bam::Reader::from_path(input_bam).map_err(hts_err)?;
    if threads > 1 {
        reader.set_threads(threads - 1).map_err(hts_err)?;
    }
    let header = bam::Header::from_template(reader.header());

    let unsorted_bam = temp_bam_path(output_bam, "clean.unsorted");
    let mut writer =
        bam::Writer::from_path(&unsorted_bam, &header, bam::Format::Bam).map_err(hts_err)?;
    if threads > 1 {
        writer.set_threads(threads - 1).map_err(hts_err)?;
    }

    let mut written = 0usize;
    for result in reader.records() {
        let mut record = result.map_err(hts_err)?;
        let qname = std::str::from_utf8(record.qname())
            .unwrap_or("")
            .to_string();
        if keep_qnames.contains(&qname)
            && !lowqual_qnames.contains(qname.as_str())
            && passes_python_clean_policy(&record)
        {
            set_hp_tag(&mut record, &clean_hp_tag(&qname, qname_hap, chunk_id))?;
            writer.write(&record).map_err(hts_err)?;
            written += 1;
        }
    }

    drop(writer);

    crate::tools::samtools_sort_index(&unsorted_bam, output_bam, threads)?;
    let _ = std::fs::remove_file(&unsorted_bam);

    Ok(written)
}

/// Write a new BAM with HP (haplotype) tag annotations for visualization.
///
/// Reads in `correct_qnames` get `HP:Z:chunk{chunk_id}_hap{hap_id}`.
/// Reads in `mismap_qnames` get `HP:Z:chunk{chunk_id}_hap{hap_id}_HIGHVD`.
/// Other reads get `HP:Z:LOWQUAL`.
pub fn annotate_hp_tags(
    input_bam: &Path,
    assignments: &HpTagAssignments<'_>,
    output_bam: &Path,
    chunk_id: usize,
    threads: usize,
) -> Result<usize> {
    let mut reader = bam::Reader::from_path(input_bam).map_err(hts_err)?;
    if threads > 1 {
        reader.set_threads(threads - 1).map_err(hts_err)?;
    }
    let header = bam::Header::from_template(reader.header());

    let unsorted_bam = temp_bam_path(output_bam, "hp.unsorted");
    let mut writer =
        bam::Writer::from_path(&unsorted_bam, &header, bam::Format::Bam).map_err(hts_err)?;
    if threads > 1 {
        writer.set_threads(threads - 1).map_err(hts_err)?;
    }

    let mut written = 0usize;
    for result in reader.records() {
        let mut record = result.map_err(hts_err)?;
        if record.is_secondary() || record.is_supplementary() {
            continue;
        }
        let qname = std::str::from_utf8(record.qname())
            .unwrap_or("")
            .to_string();

        let hp_tag = raw_hp_tag(
            &qname,
            assignments.correct_qnames,
            assignments.mismap_qnames,
            assignments.lowqual_qnames,
            assignments.qname_hap,
            chunk_id,
        );
        set_hp_tag(&mut record, &hp_tag)?;
        writer.write(&record).map_err(hts_err)?;
        written += 1;
    }

    drop(writer);
    crate::tools::samtools_sort_index(&unsorted_bam, output_bam, threads)?;
    let _ = std::fs::remove_file(&unsorted_bam);
    Ok(written)
}

fn hts_err(e: impl std::fmt::Display) -> SdError {
    SdError::Htslib(e.to_string())
}

fn passes_python_clean_policy(record: &bam::Record) -> bool {
    !record.is_secondary()
        && !record.is_supplementary()
        && !record.is_duplicate()
        && !record.is_quality_check_failed()
        && !record.is_unmapped()
        && record.mapq() > CLEAN_MIN_MAPQ
}

fn set_hp_tag(record: &mut bam::Record, hp_tag: &str) -> Result<()> {
    let _ = record.remove_aux(b"HP");
    record
        .push_aux(b"HP", Aux::String(hp_tag))
        .map_err(hts_err)?;
    Ok(())
}

fn clean_hp_tag(qname: &str, qname_hap: &HashMap<String, i32>, chunk_id: usize) -> String {
    match qname_hap.get(qname) {
        Some(hap_id) => format!("chunk{chunk_id}_{hap_id}"),
        None => format!("chunk{chunk_id}_NA"),
    }
}

fn raw_hp_tag(
    qname: &str,
    correct_qnames: &HashSet<String>,
    mismap_qnames: &HashSet<String>,
    lowqual_qnames: &HashSet<String>,
    qname_hap: &HashMap<String, i32>,
    chunk_id: usize,
) -> String {
    if lowqual_qnames.contains(qname) {
        return format!("chunk{chunk_id}_LOWQUAL");
    }
    let base = clean_hp_tag(qname, qname_hap, chunk_id);
    if mismap_qnames.contains(qname) || !correct_qnames.contains(qname) {
        format!("{base}_HIGHVD")
    } else {
        base
    }
}

fn temp_bam_path(target: &Path, suffix: &str) -> std::path::PathBuf {
    let name = target
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("output.bam");
    target.with_file_name(format!("{}.{}.{}.bam", name, std::process::id(), suffix))
}

#[cfg(test)]
mod tests {
    // Integration tests require BAM fixtures; covered by the T9 HG006
    // differential. Unit tests for the filter logic are trivial (HashSet
    // membership) so we test the error paths instead.

    use super::*;

    #[test]
    fn missing_bam_returns_error() {
        let qs: HashSet<String> = HashSet::new();
        let lowqual: HashSet<String> = HashSet::new();
        let qname_hap: HashMap<String, i32> = HashMap::new();
        let r = filter_bam_by_qnames(
            Path::new("/nonexistent.bam"),
            &qs,
            &lowqual,
            &qname_hap,
            Path::new("/tmp/out.bam"),
            1,
            1,
        );
        assert!(r.is_err());
    }
}
