//! Write HP (haplotype) tags to a BAM file based on phasing results.

use std::collections::HashMap;
use std::process::Command;

use rust_htslib::bam::{self, record::Aux, Read};
use sdrecall_utils::{Result, SdError};

/// Write a new BAM with `HP:Z:hap{N}` tags for each phased read.
///
/// Unphased reads (not present in `vertex_hap`) get `HP:Z:unphased`.
/// The output is coordinate-sorted and indexed via `samtools sort` + `samtools index`.
///
/// Returns the total number of records written.
pub fn write_hp_tagged_bam(
    input_bam: &str,
    output_bam: &str,
    vertex_hap: &HashMap<i32, i32>,
    vertex_qname: &[String],
    threads: u8,
) -> Result<usize> {
    let qname_to_hap: HashMap<&str, i32> = vertex_qname
        .iter()
        .enumerate()
        .filter_map(|(idx, qname)| {
            vertex_hap
                .get(&(idx as i32))
                .map(|&hap| (qname.as_str(), hap))
        })
        .collect();

    let mut reader = bam::Reader::from_path(input_bam)
        .map_err(|e| SdError::Htslib(format!("open {input_bam}: {e}")))?;
    if threads > 1 {
        reader
            .set_threads((threads - 1) as usize)
            .map_err(|e| SdError::Htslib(format!("set reader threads: {e}")))?;
    }
    let header = bam::Header::from_template(reader.header());

    let tmp_unsorted = format!("{output_bam}.unsorted.tmp.bam");
    let mut writer = bam::Writer::from_path(&tmp_unsorted, &header, bam::Format::Bam)
        .map_err(|e| SdError::Htslib(format!("create {tmp_unsorted}: {e}")))?;
    if threads > 1 {
        writer
            .set_threads((threads - 1) as usize)
            .map_err(|e| SdError::Htslib(format!("set writer threads: {e}")))?;
    }

    let mut n_written = 0usize;
    for result in reader.records() {
        let mut record = result.map_err(|e| SdError::Htslib(format!("BAM read error: {e}")))?;
        let qname = std::str::from_utf8(record.qname()).unwrap_or("");

        let hp_tag = match qname_to_hap.get(qname) {
            Some(&hap_id) => format!("hap{hap_id}"),
            None => "unphased".to_string(),
        };

        // Replace any pre-existing HP tag: rust-htslib's push_aux returns
        // BamAuxTagAlreadyPresent on a duplicate tag, which would abort the whole
        // BAM on pre-phased/annotated input. remove_aux errs only when the tag is
        // absent, which we ignore.
        let _ = record.remove_aux(b"HP");
        record
            .push_aux(b"HP", Aux::String(&hp_tag))
            .map_err(|e| SdError::Compute(format!("push HP tag: {e}")))?;
        writer
            .write(&record)
            .map_err(|e| SdError::Htslib(format!("write record: {e}")))?;
        n_written += 1;
    }
    drop(writer);

    samtools_sort_and_index(&tmp_unsorted, output_bam, threads)?;
    let _ = std::fs::remove_file(&tmp_unsorted);

    log::info!("[hp_writer] wrote {n_written} records to {output_bam}");
    Ok(n_written)
}

fn samtools_sort_and_index(unsorted: &str, output: &str, threads: u8) -> Result<()> {
    let status = Command::new("samtools")
        .args(["sort", "-@", &threads.to_string(), "-o", output, unsorted])
        .status()
        .map_err(|e| SdError::Htslib(format!("samtools sort: {e}")))?;
    if !status.success() {
        return Err(SdError::Compute("samtools sort failed".into()));
    }

    let status = Command::new("samtools")
        .args(["index", "-@", &threads.to_string(), output])
        .status()
        .map_err(|e| SdError::Htslib(format!("samtools index: {e}")))?;
    if !status.success() {
        return Err(SdError::Compute("samtools index failed".into()));
    }

    Ok(())
}
