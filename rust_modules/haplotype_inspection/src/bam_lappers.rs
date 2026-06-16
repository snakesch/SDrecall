/// BAM reading and Lapper construction for interval queries.
/// Replaces Python's NCLS-based `migrate_bam_to_ncls` with Lapper interval trees.
use rust_htslib::bam::{self, Read, Record};
use rust_lapper::{Lapper, Interval};
use rustc_hash::FxHashMap;
use std::collections::{HashMap, HashSet};
#[cfg(unix)]
use std::os::unix::io::AsRawFd;
use std::path::Path;
use std::process::{Command, Stdio};
use tempfile::NamedTempFile;
use log::{info, warn};

/// Result structure from BAM → Lapper construction (replaces Python's migrate_bam_to_ncls).
#[derive(Debug)]
pub struct BamLapperResult {
    /// Per-chromosome Lapper interval trees (chrom → Lapper<qname_idx>)
    pub lapper_dict: HashMap<String, Lapper<u32, u32>>,
    /// Read objects indexed by qname_idx
    pub read_dict: FxHashMap<u32, Vec<Record>>,
    /// qname_idx -> qname mapping
    pub qname_dict: FxHashMap<u32, String>,
    /// qname -> qname_idx mapping
    pub qname_idx_dict: HashMap<String, u32>,
    /// Set of noisy qnames that were filtered out
    pub noisy_qnames: HashSet<String>,
    /// Processing statistics
    pub stats: BamProcessingStats,
}

#[derive(Debug)]
pub struct BamProcessingStats {
    pub total_reads_processed: usize,
    pub skipped_alignments: usize,
    pub qnames_retained: usize,
    pub noisy_qnames_filtered: usize,
}

/// Check whether samtools is available on PATH.
fn samtools_available() -> bool {
    Command::new("samtools")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

/// Spawn `samtools collate` with stdout piped so we can stream BAM records
/// directly from the pipe via `/dev/fd/N`, avoiding temp-file disk I/O.
///
/// Returns `(Reader, Child)` on success.  The caller **must** call
/// `child.wait()` after consuming all records to reap the subprocess.
///
/// Falls back to `None` if samtools is unavailable (caller should use
/// the temp-file or two-pass approach instead).
#[cfg(unix)]
fn spawn_collate_pipe(
    bam_file_path: &str,
    threads: u8,
) -> Result<Option<(bam::Reader, std::process::Child)>, Box<dyn std::error::Error>> {
    if !samtools_available() {
        warn!("[spawn_collate_pipe] samtools not found");
        return Ok(None);
    }

    info!(
        "[spawn_collate_pipe] Running: samtools collate -f -@ {threads} {bam_file_path} -o - (piped)"
    );

    let mut child = Command::new("samtools")
        .args([
            "collate",
            "-f",
            "-@",
            &threads.to_string(),
            bam_file_path,
            "-o",
            "-", // write BAM to stdout
        ])
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;

    // Take the child's stdout pipe and open it via /dev/fd/N
    let stdout = child
        .stdout
        .take()
        .ok_or("Failed to capture samtools stdout")?;
    let fd = stdout.as_raw_fd();
    let fd_path = format!("/dev/fd/{fd}");

    match bam::Reader::from_path(&fd_path) {
        Ok(reader) => {
            // htslib's hts_open already dup'd the pipe fd, so our ChildStdout copy
            // is redundant: drop it to close that fd. `mem::forget` here would leak
            // one pipe fd per call and exhaust the fd table under per-island
            // fan-out. The child stays alive in the returned handle and is reaped
            // by the `child.wait()` at the end of `build_lapper_from_bam`.
            drop(stdout);
            info!("[spawn_collate_pipe] Opened BAM reader from pipe fd {fd}");
            Ok(Some((reader, child)))
        }
        Err(e) => {
            warn!(
                "[spawn_collate_pipe] from_path({fd_path}) failed: {e}, will fall back to temp file"
            );
            let _ = child.kill();
            let _ = child.wait();
            Ok(None)
        }
    }
}

/// Fallback: Run `samtools collate` writing to a temp file.
/// Returns the temp file handle (auto-deleted on drop), or None if samtools
/// is unavailable.
fn collate_bam_file(
    bam_file_path: &str,
    threads: u8,
) -> Result<Option<NamedTempFile>, Box<dyn std::error::Error>> {
    if !samtools_available() {
        warn!("[collate_bam_file] samtools not found, falling back to two-pass approach");
        return Ok(None);
    }

    let temp_file = NamedTempFile::with_suffix_in(".bam", Path::new("."))?;
    let temp_path = temp_file
        .path()
        .to_str()
        .ok_or("Failed to convert temp path to string")?;

    info!(
        "[collate_bam_file] Running: samtools collate -f -@ {threads} {bam_file_path} -o {temp_path}"
    );

    let output = Command::new("samtools")
        .args([
            "collate",
            "-f",
            "-@",
            &threads.to_string(),
            bam_file_path,
            "-o",
            temp_path,
        ])
        .output()?;

    if !output.status.success() {
        warn!(
            "[collate_bam_file] samtools collate failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        return Ok(None);
    }

    info!("[collate_bam_file] Collation complete: {temp_path}");
    Ok(Some(temp_file))
}

/// Build per-chromosome Lapper interval trees from a BAM file.
/// Replaces Python's `migrate_bam_to_ncls` with Lapper-based implementation.
///
/// For paired-end data, uses `samtools collate` to group reads by qname so we can
/// stream one qname group at a time without loading the entire BAM into memory.
/// Falls back to a two-pass in-memory approach if samtools is unavailable.
pub fn build_lapper_from_bam(
    bam_path: &str,
    mapq_filter: u8,
    basequal_median_filter: u8,
    paired: bool,
    filter_noisy: bool,
) -> Result<BamLapperResult, Box<dyn std::error::Error>> {
    // For paired-end data, try pipe-based collation first (zero disk I/O),
    // fall back to temp-file collation, then to two-pass in-memory.
    enum BamSource {
        /// Pipe: Reader is already constructed; Child must be waited on.
        Pipe(bam::Reader, std::process::Child),
        /// Temp file: path + handle kept alive.
        TempFile(String, NamedTempFile),
        /// Original BAM path (no collation).
        Original(String),
    }

    let source = if paired {
        // Try 1: pipe-based streaming (no disk I/O)
        #[cfg(unix)]
        {
            match spawn_collate_pipe(bam_path, 4)? {
                Some((reader, child)) => {
                    info!("[build_lapper_from_bam] Using pipe-based collation (zero disk I/O)");
                    BamSource::Pipe(reader, child)
                }
                None => {
                    // Try 2: temp-file collation
                    match collate_bam_file(bam_path, 4)? {
                        Some(tf) => {
                            let p = tf.path().to_str()
                                .ok_or("Failed to convert temp path to string")?
                                .to_string();
                            info!("[build_lapper_from_bam] Using temp-file collation: {p}");
                            BamSource::TempFile(p, tf)
                        }
                        None => {
                            info!("[build_lapper_from_bam] No collation available, two-pass fallback");
                            BamSource::Original(bam_path.to_string())
                        }
                    }
                }
            }
        }
        #[cfg(not(unix))]
        {
            match collate_bam_file(bam_path, 4)? {
                Some(tf) => {
                    let p = tf.path().to_str()
                        .ok_or("Failed to convert temp path to string")?
                        .to_string();
                    BamSource::TempFile(p, tf)
                }
                None => BamSource::Original(bam_path.to_string()),
            }
        }
    } else {
        BamSource::Original(bam_path.to_string())
    };

    // Track whether we got collated data before destructuring `source`
    let collated = matches!(source, BamSource::Pipe(_, _) | BamSource::TempFile(_, _));

    // Open BAM reader (Pipe variant already has one)
    let (mut bam, mut collate_child) = match source {
        BamSource::Pipe(reader, child) => (reader, Some(child)),
        BamSource::TempFile(ref path, ref _tf) => {
            (bam::Reader::from_path(path)?, None)
        }
        BamSource::Original(ref path) => {
            (bam::Reader::from_path(path)?, None)
        }
    };

    let header = bam.header().clone();
    let chroms: Vec<String> = header
        .target_names()
        .iter()
        .map(|name| String::from_utf8_lossy(name).to_string())
        .collect();

    // Initialize data structures
    let mut read_dict: FxHashMap<u32, Vec<Record>> = FxHashMap::default();
    let mut qname_interval_dict: HashMap<String, Vec<(i64, i64, u32)>> =
        chroms.iter().map(|c| (c.clone(), Vec::new())).collect();
    let mut qname_idx_dict: HashMap<String, u32> = HashMap::new();
    let mut qname_dict: FxHashMap<u32, String> = FxHashMap::default();
    let mut noisy_qnames: HashSet<String> = HashSet::new();
    let mut total_qnames: HashSet<String> = HashSet::new();

    let mut qname_idx_counter: u32 = 0;
    let mut total_reads_processed: usize = 0;
    let mut skipped_alignments: usize = 0;

    // Determine whether reads are collated (pipe or temp-file)
    // (computed above before `source` was moved)

    if collated {
        // --- Streaming approach: BAM is collated, so same-qname reads are adjacent ---
        // Buffer one qname group at a time and process immediately.
        let mut current_qname: Option<String> = None;
        let mut current_group: Vec<Record> = Vec::new();

        for result in bam.records() {
            let record = result?;
            total_reads_processed += 1;

            if record.is_secondary() || record.is_supplementary() || record.is_duplicate() {
                skipped_alignments += 1;
                continue;
            }

            let qname = String::from_utf8_lossy(record.qname()).to_string();

            if current_qname.as_deref() == Some(&qname) {
                // Same qname group — accumulate
                current_group.push(record);
            } else {
                // New qname — flush the previous group
                if let Some(ref prev_qname) = current_qname {
                    if !current_group.is_empty() {
                        qname_idx_counter = process_qname_group(
                            prev_qname,
                            &current_group,
                            &header,
                            &mut read_dict,
                            &mut qname_interval_dict,
                            &mut qname_idx_dict,
                            &mut qname_dict,
                            &mut noisy_qnames,
                            &mut total_qnames,
                            qname_idx_counter,
                            paired,
                            mapq_filter,
                            basequal_median_filter,
                            filter_noisy,
                        );
                    }
                }
                current_qname = Some(qname);
                current_group.clear();
                current_group.push(record);
            }
        }

        // Flush the last group (return value unused: this is the final group)
        if let Some(ref prev_qname) = current_qname {
            if !current_group.is_empty() {
                let _ = process_qname_group(
                    prev_qname,
                    &current_group,
                    &header,
                    &mut read_dict,
                    &mut qname_interval_dict,
                    &mut qname_idx_dict,
                    &mut qname_dict,
                    &mut noisy_qnames,
                    &mut total_qnames,
                    qname_idx_counter,
                    paired,
                    mapq_filter,
                    basequal_median_filter,
                    filter_noisy,
                );
            }
        }
    } else {
        // --- Two-pass fallback: no collation (single-end, or samtools unavailable) ---
        // For single-end, reads are unique per qname so streaming works directly.
        // For paired without collation, we must collect all reads first.
        if paired {
            // Paired but no samtools — must load all into memory
            warn!("[build_lapper_from_bam] Two-pass fallback for paired data (high memory usage)");
            let mut reads_by_qname: HashMap<String, Vec<Record>> = HashMap::new();

            for result in bam.records() {
                let record = result?;
                total_reads_processed += 1;

                if record.is_secondary() || record.is_supplementary() || record.is_duplicate() {
                    skipped_alignments += 1;
                    continue;
                }

                let qname = String::from_utf8_lossy(record.qname()).to_string();
                reads_by_qname.entry(qname).or_default().push(record);
            }

            for (qname, reads) in &reads_by_qname {
                if !reads.is_empty() {
                    qname_idx_counter = process_qname_group(
                        qname,
                        reads,
                        &header,
                        &mut read_dict,
                        &mut qname_interval_dict,
                        &mut qname_idx_dict,
                        &mut qname_dict,
                        &mut noisy_qnames,
                        &mut total_qnames,
                        qname_idx_counter,
                        paired,
                        mapq_filter,
                        basequal_median_filter,
                        filter_noisy,
                    );
                }
            }
        } else {
            // Single-end: each record is its own group, stream directly
            for result in bam.records() {
                let record = result?;
                total_reads_processed += 1;

                if record.is_secondary() || record.is_supplementary() || record.is_duplicate() {
                    skipped_alignments += 1;
                    continue;
                }

                let qname = String::from_utf8_lossy(record.qname()).to_string();
                qname_idx_counter = process_qname_group(
                    &qname,
                    &[record],
                    &header,
                    &mut read_dict,
                    &mut qname_interval_dict,
                    &mut qname_idx_dict,
                    &mut qname_dict,
                    &mut noisy_qnames,
                    &mut total_qnames,
                    qname_idx_counter,
                    paired,
                    mapq_filter,
                    basequal_median_filter,
                    filter_noisy,
                );
            }
        }
    }

    // Reap the samtools collate child process (if pipe-based)
    if let Some(mut child) = collate_child.take() {
        let status = child.wait()?;
        if !status.success() {
            warn!("[build_lapper_from_bam] samtools collate exited with: {status}");
        }
    }

    // Build Lapper for each chromosome
    let mut lapper_dict: HashMap<String, Lapper<u32, u32>> = HashMap::new();

    for chrom in &chroms {
        let chrom_intervals = &qname_interval_dict[chrom];
        if !chrom_intervals.is_empty() {
            let mut intervals: Vec<Interval<u32, u32>> =
                Vec::with_capacity(chrom_intervals.len());
            for (start, end, qname_idx) in chrom_intervals {
                // Checked i64 → u32: human coordinates fit comfortably, but a
                // >4 Gb contig (or a stray negative position) would silently
                // truncate and corrupt the interval tree — fail loudly instead.
                let start = u32::try_from(*start)
                    .map_err(|_| format!("interval start {start} on {chrom} exceeds u32 range"))?;
                let stop = u32::try_from(*end)
                    .map_err(|_| format!("interval end {end} on {chrom} exceeds u32 range"))?;
                intervals.push(Interval { start, stop, val: *qname_idx });
            }

            intervals.sort_by_key(|iv| iv.start);
            lapper_dict.insert(chrom.clone(), Lapper::new(intervals));
        }
    }

    // Check if we have enough data for paired mode
    if paired && total_qnames.len() <= 2 {
        return Err("Insufficient paired-end data".into());
    }

    // Compute stats before moving into result struct
    let qnames_retained = qname_idx_dict.len();
    let noisy_qnames_filtered = noisy_qnames.len();

    info!(
        "[build_lapper_from_bam] Done: {total_reads_processed} reads processed, {skipped_alignments} skipped, {qnames_retained} qnames retained, {noisy_qnames_filtered} noisy filtered"
    );

    Ok(BamLapperResult {
        lapper_dict,
        read_dict,
        qname_dict,
        qname_idx_dict,
        noisy_qnames,
        stats: BamProcessingStats {
            total_reads_processed,
            skipped_alignments,
            qnames_retained,
            noisy_qnames_filtered,
        },
    })
}

/// Process all reads for a single qname
#[allow(clippy::too_many_arguments)]
fn process_qname_group(
    qname: &str,
    reads: &[Record],
    header: &rust_htslib::bam::HeaderView,
    read_dict: &mut FxHashMap<u32, Vec<Record>>,
    qname_interval_dict: &mut HashMap<String, Vec<(i64, i64, u32)>>,
    qname_idx_dict: &mut HashMap<String, u32>,
    qname_dict: &mut FxHashMap<u32, String>,
    noisy_qnames: &mut HashSet<String>,
    total_qnames: &mut HashSet<String>,
    qname_idx_counter: u32,
    paired: bool,
    mapq_filter: u8,
    basequal_median_filter: u8,
    filter_noisy: bool,
) -> u32 {
    total_qnames.insert(qname.to_string());

    // Check if any read in the group is noisy
    let is_noisy = reads.iter().any(|read| {
        is_read_noisy(read, header, paired, mapq_filter, basequal_median_filter, filter_noisy)
    });

    if is_noisy {
        noisy_qnames.insert(qname.to_string());
        return qname_idx_counter;
    }

    // For paired-end mode, ensure we have both read1 and read2
    if paired {
        let has_read1 = reads.iter().any(|r| r.is_first_in_template());
        let has_read2 = reads.iter().any(|r| r.is_last_in_template());

        if !has_read1 || !has_read2 {
            return qname_idx_counter;
        }
    }

    // Assign qname_idx and record reads
    let qname_idx = if let Some(&idx) = qname_idx_dict.get(qname) {
        idx
    } else {
        let idx = qname_idx_counter;
        qname_idx_dict.insert(qname.to_string(), idx);
        qname_dict.insert(idx, qname.to_string());
        idx
    };

    // Append reads for this qname (don't overwrite — in coordinate-sorted BAMs,
    // R1 and R2 of the same qname arrive as separate groups)
    read_dict.entry(qname_idx).or_default().extend(reads.iter().cloned());

    // Store one interval per read (not merged across R1/R2)
    for read in reads {
        if let Some(chrom_name) = get_reference_name(read, header) {
            let start = read.pos();
            let end = read.cigar().end_pos();

            qname_interval_dict
                .get_mut(&chrom_name)
                .unwrap()
                .push((start, end, qname_idx));
        }
    }

    if qname_idx == qname_idx_counter {
        qname_idx_counter + 1
    } else {
        qname_idx_counter
    }
}

/// Check if a read is noisy based on quality criteria
fn is_read_noisy(
    read: &Record,
    header: &rust_htslib::bam::HeaderView,
    paired: bool,
    mapq_filter: u8,
    basequal_median_filter: u8,
    filter_noisy: bool,
) -> bool {
    // Basic quality checks
    if read.is_unmapped() {
        return true;
    }

    // Python (bam_ncls.is_read_noisy) drops QC-fail reads.
    if read.is_quality_check_failed() {
        return true;
    }

    if read.mapq() < mapq_filter {
        return true;
    }

    // Paired-end specific checks
    if paired {
        let aln_len = read.cigar().end_pos() - read.pos();
        if aln_len < 75 {
            return true;
        }

        // Ensure mate on same reference. Python compares reference_name !=
        // next_reference_name, which also flags a read whose mate reference is
        // unset (None) — so compare the Options directly rather than requiring
        // both to be Some.
        if get_reference_name(read, header) != get_next_reference_name(read, header) {
            return true;
        }

        // Check proper pair flag (warning only, not filtering)
        // if !read.is_proper_pair() {
        //     // Python logs warning but doesn't filter
        // }
    } else {
        // Single-end: duplicates already filtered earlier
    }

    // Base quality and soft-clip checks
    if filter_noisy {
        let qualities = read.qual();
        if !qualities.is_empty() {
            let median_q = fast_median(qualities);
            if median_q <= basequal_median_filter as f32 {
                return true;
            }

            let num_low = qualities.iter().filter(|&&q| q < basequal_median_filter).count();
            if num_low >= 75 {
                return true;
            }

            // Soft-clip length from CIGAR
            let softclip: u32 = read
                .cigar()
                .iter()
                .filter_map(|c| {
                    if let rust_htslib::bam::record::Cigar::SoftClip(len) = c {
                        Some(*len)
                    } else {
                        None
                    }
                })
                .sum();

            if softclip >= 75 {
                return true;
            }
        }
    }

    false
}

/// Median of quality scores, matching Python's `np.median`.
///
/// Returns a float: for an even-length input the mean of the two middle
/// elements is NOT integer-truncated (e.g. median of [15, 16] is 15.5, not 15).
/// This parity matters at the `median <= cutoff` boundary in `is_read_noisy`.
fn fast_median(qualities: &[u8]) -> f32 {
    if qualities.is_empty() {
        return 0.0;
    }

    let mut sorted = qualities.to_vec();
    sorted.sort_unstable();
    let mid = sorted.len() / 2;

    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] as f32 + sorted[mid] as f32) / 2.0
    } else {
        sorted[mid] as f32
    }
}

/// Get reference name for a read
fn get_reference_name(read: &Record, header: &rust_htslib::bam::HeaderView) -> Option<String> {
    let tid = read.tid();
    if tid >= 0 {
        Some(String::from_utf8_lossy(header.tid2name(tid as u32)).to_string())
    } else {
        None
    }
}

/// Get next reference name for a read (mate)
fn get_next_reference_name(read: &Record, header: &rust_htslib::bam::HeaderView) -> Option<String> {
    let mtid = read.mtid();
    if mtid >= 0 {
        Some(String::from_utf8_lossy(header.tid2name(mtid as u32)).to_string())
    } else {
        None
    }
}

/// Query overlapping reads by qname_idx
/// Replaces Python's overlap_qname_idx_iterator
pub fn query_overlapping_qname_indices(
    lapper: &Lapper<u32, u32>,
    start: u32,
    end: u32,
) -> Vec<u32> {
    lapper.find(start, end).map(|iv| iv.val).collect()
}

/// Query overlapping reads and return actual read objects
/// Replaces Python's overlapping_reads_iterator
pub fn query_overlapping_reads<'a>(
    lapper_dict: &HashMap<String, Lapper<u32, u32>>,
    read_dict: &'a FxHashMap<u32, Vec<Record>>,
    chrom: &str,
    start: u32,
    end: u32,
) -> Vec<&'a Record> {
    let mut result = Vec::new();

    if let Some(lapper) = lapper_dict.get(chrom) {
        // Get unique qname indices
        let mut seen = HashSet::new();
        for qname_idx in lapper.find(start, end).map(|iv| iv.val) {
            if seen.insert(qname_idx) {
                // Get all reads for this qname
                if let Some(reads) = read_dict.get(&qname_idx) {
                    for read in reads {
                        // Filter to ensure actual overlap
                        let read_start = read.pos() as u32;
                        let read_end = read.cigar().end_pos() as u32;
                        if read_start < end && read_end > start {
                            result.push(read);
                        }
                    }
                }
            }
        }
    }

    result
}

/// Query overlapping read pairs (for paired-end data)
/// Returns pairs of (read1, read2) that overlap the interval
pub fn query_overlapping_read_pairs<'a>(
    lapper_dict: &HashMap<String, Lapper<u32, u32>>,
    read_dict: &'a FxHashMap<u32, Vec<Record>>,
    chrom: &str,
    start: u32,
    end: u32,
) -> Vec<(&'a Record, &'a Record)> {
    let mut result = Vec::new();

    if let Some(lapper) = lapper_dict.get(chrom) {
        let mut seen = HashSet::new();
        for qname_idx in lapper.find(start, end).map(|iv| iv.val) {
            if seen.insert(qname_idx) {
                if let Some(reads) = read_dict.get(&qname_idx) {
                    // Find read1 and read2
                    let mut read1 = None;
                    let mut read2 = None;

                    for read in reads {
                        if read.is_first_in_template() {
                            read1 = Some(read);
                        } else if read.is_last_in_template() {
                            read2 = Some(read);
                        }
                    }

                    // Only include if we have both reads and at least one overlaps
                    if let (Some(r1), Some(r2)) = (read1, read2) {
                        let r1_start = r1.pos() as u32;
                        let r1_end = r1.cigar().end_pos() as u32;
                        let r2_start = r2.pos() as u32;
                        let r2_end = r2.cigar().end_pos() as u32;

                        if (r1_start < end && r1_end > start) || (r2_start < end && r2_end > start) {
                            result.push((r1, r2));
                        }
                    }
                }
            }
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fast_median() {
        // Matches np.median: odd length → middle element; even length → mean
        // of the two middle elements as a float (NOT integer-truncated).
        assert_eq!(fast_median(&[1, 2, 3, 4, 5]), 3.0);
        assert_eq!(fast_median(&[1, 2, 3, 4]), 2.5);
        assert_eq!(fast_median(&[5, 1, 3, 2, 4]), 3.0);
        assert_eq!(fast_median(&[]), 0.0);
    }

    /// Integration test: run build_lapper_from_bam in paired mode on real HG002 BAM.
    /// This exercises the pipe-based collation path (or temp-file fallback).
    #[test]
    fn test_collate_pipe_real_bam() {
        let _ = env_logger::try_init();

        let bam_path = "/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002_hg38_exome_SDrecall/recall_results/HG002.pooled.clean.bam";
        if !std::path::Path::new(bam_path).exists() {
            eprintln!("SKIP: test BAM not found at {bam_path}");
            return;
        }

        let result = build_lapper_from_bam(
            bam_path,
            0,     // mapq_filter
            0,     // basequal_median_filter
            true,  // paired
            false, // filter_noisy
        )
        .expect("build_lapper_from_bam failed");

        eprintln!("=== Collate Pipe Test Results ===");
        eprintln!("Total reads processed: {}", result.stats.total_reads_processed);
        eprintln!("Skipped alignments: {}", result.stats.skipped_alignments);
        eprintln!("Qnames retained: {}", result.stats.qnames_retained);
        eprintln!("Noisy qnames filtered: {}", result.stats.noisy_qnames_filtered);
        eprintln!("Chromosomes with intervals: {}", result.lapper_dict.len());
        eprintln!("read_dict entries: {}", result.read_dict.len());

        let total_records: usize = result.read_dict.values().map(|v| v.len()).sum();
        eprintln!("Total BAM records in read_dict: {total_records}");

        // In paired mode, every retained qname should have exactly 2 records
        let mut singles = 0usize;
        let mut pairs = 0usize;
        let mut triples_plus = 0usize;
        for records in result.read_dict.values() {
            match records.len() {
                1 => singles += 1,
                2 => pairs += 1,
                _ => triples_plus += 1,
            }
        }
        eprintln!("Read grouping: 1-rec={singles} 2-rec={pairs} 3+-rec={triples_plus}");

        // Basic sanity checks
        assert!(result.stats.total_reads_processed > 0, "No reads processed");
        assert!(result.stats.qnames_retained > 0, "No qnames retained");
        assert!(!result.read_dict.is_empty(), "read_dict is empty");
        // Paired mode should produce mostly pairs
        assert!(pairs > singles, "Expected more pairs than singles in paired mode");
    }
}
