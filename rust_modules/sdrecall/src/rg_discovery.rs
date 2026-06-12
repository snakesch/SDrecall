//! RG (realign-group) discovery — port of
//! `realign_recall/stat_realign_group_regions.py::stat_all_RG_region_size`.
//!
//! Reads the `*_related_homo_regions.bed` files written by sd-prep (Phase 1)
//! for each RG, extracts the FC subgroup IDs, and returns the RG list sorted
//! by total region size (largest first) for load-balanced parallelism.

use std::io::{BufRead, BufReader};
use std::path::Path;

use sdrecall_utils::{Result, SdError};

/// One realign group with its subgroup IDs and total region span (for sorting).
#[derive(Clone, Debug)]
pub struct RgInfo {
    pub label: String,
    pub subgroup_ids: Vec<String>,
    pub total_span: i64,
}

/// Discover RGs and their subgroups from the Phase-1 all-homo-regions BED
/// files. Returns the RGs sorted by total span (largest first), mirroring
/// the Python's pybedtools `total_coverage()` sort.
///
/// `rg_labels` is the list of RG labels produced by sd-prep (from
/// `PrepResult.rg_outputs`). For each, reads
/// `paths.all_homo_regions_bed_path(rg)`.
pub fn stat_all_rg_region_size(
    rg_labels: &[String],
    all_homo_beds: &[impl AsRef<Path>],
) -> Result<Vec<RgInfo>> {
    assert_eq!(rg_labels.len(), all_homo_beds.len());

    let mut rgs: Vec<RgInfo> = Vec::with_capacity(rg_labels.len());

    for (label, bed_path) in rg_labels.iter().zip(all_homo_beds.iter()) {
        let bed_path = bed_path.as_ref();
        let (subgroup_ids, total_span) = parse_all_regions_bed(bed_path)?;
        rgs.push(RgInfo {
            label: label.clone(),
            subgroup_ids,
            total_span,
        });
    }

    // Sort largest-first (load balancing).
    rgs.sort_by(|a, b| b.total_span.cmp(&a.total_span));
    Ok(rgs)
}

/// Parse an all-homo-regions BED file. Returns (FC subgroup IDs, total span).
///
/// The BED6 file has lines like:
///   chr1  1000  2000  .  .  FC:0
///   chr1  3000  4000  .  .  NFC:0
///
/// FC rows indicate fully-covered subgroups; we extract the unique FC
/// subgroup indices as sequential IDs, mirroring the Python `fc_bedf` logic.
fn parse_all_regions_bed(path: &Path) -> Result<(Vec<String>, i64)> {
    let file = std::fs::File::open(path).map_err(|e| SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    let reader = BufReader::new(file);

    let mut fc_count: usize = 0;
    let mut total_span: i64 = 0;

    for line in reader.lines() {
        let line = line.map_err(|e| SdError::Io {
            path: path.display().to_string(),
            source: e,
        })?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 {
            continue;
        }

        let start: i64 = cols[1].parse().unwrap_or(0);
        let end: i64 = cols[2].parse().unwrap_or(0);
        total_span += end - start;

        // The last column contains "FC:N" or "NFC:N" (strand field in BED6).
        if let Some(last) = cols.last() {
            if last.starts_with("FC:") {
                fc_count += 1;
            }
        }
    }

    // Python: sub_ids = [i for i in range(fc_bedf.shape[0])] after drop_duplicates.
    // We use a simple sequential count here.
    let subgroup_ids: Vec<String> = (0..fc_count).map(|i| i.to_string()).collect();

    Ok((subgroup_ids, total_span))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn parses_all_regions_bed() {
        let dir = tempfile::tempdir().unwrap();
        let bed = dir.path().join("RG0_related_homo_regions.bed");
        {
            let mut f = std::fs::File::create(&bed).unwrap();
            writeln!(f, "chr1\t1000\t2000\t.\t.\tFC:0").unwrap();
            writeln!(f, "chr1\t3000\t5000\t.\t.\tNFC:0").unwrap();
            writeln!(f, "chr1\t6000\t7000\t.\t.\tFC:1").unwrap();
        }
        let (ids, span) = parse_all_regions_bed(&bed).unwrap();
        assert_eq!(ids, vec!["0", "1"]);
        assert_eq!(span, 4000); // 1000 + 2000 + 1000
    }

    #[test]
    fn stat_sorts_by_total_span() {
        let dir = tempfile::tempdir().unwrap();
        let bed0 = dir.path().join("RG0.bed");
        let bed1 = dir.path().join("RG1.bed");
        {
            let mut f = std::fs::File::create(&bed0).unwrap();
            writeln!(f, "chr1\t0\t100\t.\t.\tFC:0").unwrap();
        }
        {
            let mut f = std::fs::File::create(&bed1).unwrap();
            writeln!(f, "chr1\t0\t5000\t.\t.\tFC:0").unwrap();
            writeln!(f, "chr1\t6000\t8000\t.\t.\tFC:1").unwrap();
        }
        let rgs = stat_all_rg_region_size(
            &["RG0".into(), "RG1".into()],
            &[&bed0, &bed1],
        )
        .unwrap();
        assert_eq!(rgs[0].label, "RG1"); // larger first
        assert_eq!(rgs[1].label, "RG0");
    }
}
