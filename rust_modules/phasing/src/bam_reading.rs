use crate::structs::{AlleleDepthMap, ReadPair, ReadPairMap, SortedVecIntervals};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{self, Read, Record};

use ahash::AHashMap; // Faster HashMap for string keys
use clap::ValueEnum;
use log::{debug, error, info, warn};
use std::collections::VecDeque;
use std::io::{BufRead, BufReader, Read as IoRead};
#[cfg(unix)]
use std::os::fd::AsRawFd;
use std::path::Path;
use std::process::{Command, Stdio};
use std::sync::{Arc, Mutex};
use std::thread;
use tempfile::{NamedTempFile, TempDir};

/// How coordinate-sorted BAM records are grouped into read pairs.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub enum PairingEngine {
    /// Retained benchmark control: write a collated BAM and reopen it.
    #[value(name = "temp-file")]
    SamtoolsTempFile,
    /// Stream uncompressed collated BAM records through an OS pipe.
    #[default]
    #[value(name = "samtools-pipe")]
    SamtoolsPipe,
    /// Group primary records directly in Rust without invoking samtools.
    #[value(name = "rust-memory")]
    RustMemory,
}

impl std::fmt::Display for PairingEngine {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let value = match self {
            Self::SamtoolsTempFile => "temp-file",
            Self::SamtoolsPipe => "samtools-pipe",
            Self::RustMemory => "rust-memory",
        };
        formatter.write_str(value)
    }
}

fn htslib_additional_threads(total_threads: u8) -> u8 {
    total_threads.saturating_sub(1)
}

struct CollateScratch {
    _directory: TempDir,
    prefix: String,
}

impl CollateScratch {
    fn new() -> Result<Self, Box<dyn std::error::Error>> {
        let directory = tempfile::Builder::new()
            .prefix("sdrecall-collate-")
            .tempdir()?;
        let prefix = directory
            .path()
            .join("spill")
            .to_string_lossy()
            .into_owned();
        Ok(Self {
            _directory: directory,
            prefix,
        })
    }
}

fn should_skip_alignment(read: &Record) -> bool {
    // Skip secondary, supplementary, and duplicate alignments early
    // These are not quality issues but different types of alignments
    read.is_secondary() || read.is_supplementary() || read.is_duplicate()
}

/// Median Phred score with `np.median` semantics: odd length → the middle
/// element; even length → the mean of the two middle elements (NOT
/// integer-truncated, e.g. median of `[15, 16]` is `15.5`). Mirrors
/// `haplotype_inspection::fast_median` and the Python `numba_operators.fast_median`
/// (`np.median`) so both BAM-read paths apply the `median <= cutoff` noise filter
/// identically.
fn median_phred(quals: &[u8]) -> f32 {
    if quals.is_empty() {
        return 0.0;
    }
    let mut sorted: Vec<u8> = quals.to_vec();
    sorted.sort_unstable();
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] as f32 + sorted[mid] as f32) / 2.0
    } else {
        sorted[mid] as f32
    }
}

fn is_read_noisy(
    read: &Record,
    mapq_filter: u8,
    basequal_median_filter: u8,
    filter_noisy: bool,
) -> bool {
    let qname = String::from_utf8_lossy(read.qname());

    // Only evaluate primary alignments for noise - secondary or supplementary alignments are skipped
    if read.is_secondary() || read.is_supplementary() {
        debug!("[is_read_noisy] skip_secondary_supplementary - {} skipped (secondary/supplementary alignment); not considered noisy", qname);
        return false;
    }

    // Common fast checks for both paired and unpaired
    if read.is_unmapped() {
        debug!(
            "[is_read_noisy] unmapped_check - {} flagged noisy: unmapped read",
            qname
        );
        return true;
    }

    if read.is_quality_check_failed() {
        debug!(
            "[is_read_noisy] qc_fail_check - {} flagged noisy: QC fail flag set",
            qname
        );
        return true;
    }

    if read.mapq() < mapq_filter {
        debug!(
            "[is_read_noisy] mapq_check - {} flagged noisy: MAPQ {} < threshold {}",
            qname,
            read.mapq(),
            mapq_filter
        );
        return true;
    }

    if read.seq_len() == 0 {
        debug!(
            "[is_read_noisy] seq_len_check - {} flagged noisy: missing query_sequence",
            qname
        );
        return true;
    }

    // Paired-end specific checks
    let aln_len = read.reference_end() - read.reference_start();
    if aln_len < 75 {
        debug!(
            "[is_read_noisy] aln_len_check - {} flagged noisy: alignment span {} < 75",
            qname, aln_len
        );
        return true;
    }

    // Ensure mate on same reference for proper pairing in this pipeline
    if read.tid() != read.mtid() {
        debug!("[is_read_noisy] mate_tid_check - {} flagged noisy: different reference chromosomes for mate pair", qname);
        return true;
    }

    if !read.is_proper_pair() {
        debug!(
            "[is_read_noisy] proper_pair_check - warning:{} is not a proper pair",
            qname
        );
    }

    // Base quality and soft-clip based checks (controlled by filter_noisy)
    if filter_noisy && !read.qual().is_empty() {
        // np.median parity: average the two middle values for even-length quality
        // arrays (see `median_phred`); the `order_stat::kth(len/2)` this replaced
        // returned the upper-middle element, diverging from Python at the boundary.
        let median_qual = median_phred(read.qual());

        if median_qual <= basequal_median_filter as f32 {
            debug!("[is_read_noisy] median_qual_check - {} flagged noisy: median baseQ {} <= threshold {}", qname, median_qual, basequal_median_filter);
            return true;
        }

        let low_qual_count = read
            .qual()
            .iter()
            .filter(|&&q| q < basequal_median_filter)
            .count();
        if low_qual_count >= 75 {
            debug!("[is_read_noisy] low_qual_count_check - {} flagged noisy: #bases with Q<{} is {} >= 75", qname, basequal_median_filter, low_qual_count);
            return true;
        }

        let soft_clip_bases: u32 = read
            .cigar()
            .iter()
            .filter(|c| c.char() == 'S')
            .map(|c| c.len())
            .sum();

        if soft_clip_bases >= 75 {
            debug!("[is_read_noisy] soft_clip_check - {} flagged noisy: total soft-clip length {} >= 75\n", qname, soft_clip_bases);
            return true;
        }
    }

    // If none of the noisy conditions triggered
    false
}

/// Retained control engine: collate through a temporary BAM on disk.
fn collate_bam_temp_file(
    bam_file_path: &str,
    threads: u8,
) -> Result<NamedTempFile, Box<dyn std::error::Error>> {
    let temp_file = NamedTempFile::with_suffix_in(".bam", Path::new("."))?;
    let temp_path = temp_file
        .path()
        .to_str()
        .ok_or("failed to convert collated BAM path to UTF-8")?;

    let scratch = CollateScratch::new()?;
    let additional_threads = htslib_additional_threads(threads).to_string();
    let output = Command::new("samtools")
        .args([
            "collate",
            "-f",
            "-T",
            &scratch.prefix,
            "-@",
            &additional_threads,
            bam_file_path,
            "-o",
            temp_path,
        ])
        .output()?;

    if !output.status.success() {
        return Err(format!(
            "samtools collate failed: {}",
            String::from_utf8_lossy(&output.stderr)
        )
        .into());
    }

    Ok(temp_file)
}

#[cfg(unix)]
struct CollateProcess {
    child: std::process::Child,
    stderr_reader: Option<thread::JoinHandle<Vec<u8>>>,
    _scratch: CollateScratch,
    finished: bool,
}

#[cfg(unix)]
impl CollateProcess {
    fn finish(mut self) -> Result<(), Box<dyn std::error::Error>> {
        let status = self.child.wait()?;
        let stderr = self
            .stderr_reader
            .take()
            .ok_or("samtools collate stderr reader was already consumed")?
            .join()
            .map_err(|_| "samtools collate stderr reader panicked")?;
        self.finished = true;

        if status.success() {
            Ok(())
        } else {
            Err(format!(
                "samtools collate failed with {status}: {}",
                String::from_utf8_lossy(&stderr)
            )
            .into())
        }
    }
}

#[cfg(unix)]
impl Drop for CollateProcess {
    fn drop(&mut self) {
        if self.finished {
            return;
        }
        let _ = self.child.kill();
        let _ = self.child.wait();
        if let Some(stderr_reader) = self.stderr_reader.take() {
            let _ = stderr_reader.join();
        }
    }
}

enum PairingReaderGuard {
    None,
    TempFile {
        _file: NamedTempFile,
    },
    #[cfg(unix)]
    Pipe(CollateProcess),
}

impl PairingReaderGuard {
    fn finish(self) -> Result<(), Box<dyn std::error::Error>> {
        match self {
            Self::None | Self::TempFile { .. } => Ok(()),
            #[cfg(unix)]
            Self::Pipe(process) => process.finish(),
        }
    }
}

#[cfg(unix)]
fn collate_bam_pipe(
    bam_file_path: &str,
    threads: u8,
) -> Result<(bam::Reader, PairingReaderGuard), Box<dyn std::error::Error>> {
    let scratch = CollateScratch::new()?;
    let additional_threads = htslib_additional_threads(threads).to_string();
    let mut child = Command::new("samtools")
        .args([
            "collate",
            "-f",
            "-u",
            "-O",
            "-T",
            &scratch.prefix,
            "-@",
            &additional_threads,
            bam_file_path,
        ])
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;

    let stdout = child
        .stdout
        .take()
        .ok_or("failed to capture samtools collate stdout")?;
    let stderr = child
        .stderr
        .take()
        .ok_or("failed to capture samtools collate stderr")?;
    let stderr_reader = thread::spawn(move || {
        let mut bytes = Vec::new();
        let _ = BufReader::new(stderr).read_to_end(&mut bytes);
        bytes
    });

    let fd_path = format!("/dev/fd/{}", stdout.as_raw_fd());
    let reader = match bam::Reader::from_path(&fd_path) {
        Ok(reader) => reader,
        Err(error) => {
            let _ = child.kill();
            let _ = child.wait();
            let _ = stderr_reader.join();
            return Err(format!("failed to open collated BAM pipe {fd_path}: {error}").into());
        }
    };
    drop(stdout);

    Ok((
        reader,
        PairingReaderGuard::Pipe(CollateProcess {
            child,
            stderr_reader: Some(stderr_reader),
            _scratch: scratch,
            finished: false,
        }),
    ))
}

fn open_pairing_reader(
    bam_file_path: &str,
    engine: PairingEngine,
    threads: u8,
) -> Result<(bam::Reader, PairingReaderGuard), Box<dyn std::error::Error>> {
    match engine {
        PairingEngine::SamtoolsTempFile => {
            let temp_file = collate_bam_temp_file(bam_file_path, threads)?;
            let reader = bam::Reader::from_path(temp_file.path())?;
            Ok((reader, PairingReaderGuard::TempFile { _file: temp_file }))
        }
        PairingEngine::SamtoolsPipe => {
            #[cfg(unix)]
            {
                collate_bam_pipe(bam_file_path, threads)
            }
            #[cfg(not(unix))]
            {
                let temp_file = collate_bam_temp_file(bam_file_path, threads)?;
                let reader = bam::Reader::from_path(temp_file.path())?;
                Ok((reader, PairingReaderGuard::TempFile { _file: temp_file }))
            }
        }
        PairingEngine::RustMemory => {
            let mut reader = bam::Reader::from_path(bam_file_path)?;
            let additional_threads = htslib_additional_threads(threads);
            if additional_threads > 0 {
                reader.set_threads(usize::from(additional_threads))?;
            }
            Ok((reader, PairingReaderGuard::None))
        }
    }
}

/// Process BAM file by iterating through qname-grouped reads
/// This is more efficient than the original single-read iteration
pub fn migrate_bam_to_sorted_intervals_grouped(
    bam_file_path: &str,
    mapq_filter: u8,
    basequal_median_filter: u8,
    filter_noisy: bool,
    engine: PairingEngine,
    threads: u8,
) -> Result<(ReadPairMap, rust_htslib::bam::HeaderView), Box<dyn std::error::Error>> {
    info!(
        "[migrate_bam_to_sorted_intervals_grouped] file={bam_file_path} \
         mapq_filter={mapq_filter} basequal_median_filter={basequal_median_filter} \
         filter_noisy={filter_noisy} engine={engine} threads={threads}"
    );

    let (mut bam, reader_guard) = open_pairing_reader(bam_file_path, engine, threads)?;
    let header = bam.header().clone();
    info!("[migrate_bam_to_sorted_intervals_grouped] BAM file opened successfully, {} chromosomes found", header.target_count());

    let mut result = ReadPairMap::new();
    let mut qname_idx_counter = 0usize;

    // Pre-allocate chromosome interval trees
    let chrom_count = header.target_count() as usize;
    result.interval_trees.reserve(chrom_count);

    for tid in 0..header.target_count() {
        let ref_name = header.tid2name(tid);
        let chrom = String::from_utf8_lossy(ref_name).to_string();
        result
            .interval_trees
            .insert(chrom, SortedVecIntervals::new());
    }
    info!(
        "[migrate_bam_to_sorted_intervals_grouped] Interval trees initialized for \
         {chrom_count} chromosomes"
    );

    if engine == PairingEngine::RustMemory {
        let mut qname_order = Vec::new();
        let mut groups: AHashMap<String, Vec<Record>> = AHashMap::new();
        let mut total_reads_processed = 0usize;
        let mut skipped_alignments = 0usize;

        for read_result in bam.records() {
            let read = read_result?;
            total_reads_processed += 1;
            if should_skip_alignment(&read) {
                skipped_alignments += 1;
                continue;
            }

            let qname = String::from_utf8_lossy(read.qname()).to_string();
            if let Some(reads) = groups.get_mut(&qname) {
                reads.push(read);
            } else {
                qname_order.push(qname.clone());
                groups.insert(qname, vec![read]);
            }
        }

        for qname in qname_order {
            let mut reads = groups
                .remove(&qname)
                .expect("qname order and in-memory groups stay synchronized");
            process_qname_group(
                &mut result,
                &header,
                qname,
                &mut reads,
                &mut qname_idx_counter,
                mapq_filter,
                basequal_median_filter,
                filter_noisy,
            )?;
        }

        drop(bam);
        reader_guard.finish()?;
        for interval_tree in result.interval_trees.values_mut() {
            interval_tree.finalize()?;
        }
        info!(
            "[migrate_bam_to_sorted_intervals_grouped] BAM processing complete: {} total reads processed, {} alignments skipped, {} read pairs retained, {} noisy qnames filtered",
            total_reads_processed,
            skipped_alignments,
            result.readpair_dict.len(),
            result.noisy_qnames.len()
        );
        return Ok((result, header));
    }

    // Buffer for collecting reads with the same qname
    let mut current_qname: Option<String> = None;
    let mut current_reads: Vec<Record> = Vec::with_capacity(2);

    let mut total_reads_processed = 0usize;
    let mut skipped_alignments = 0usize;
    info!("[migrate_bam_to_sorted_intervals_grouped] Starting to process BAM records");

    // Process reads grouped by qname
    for read_result in bam.records() {
        let read = read_result?;
        total_reads_processed += 1;

        // Skip secondary, supplementary, and duplicate alignments
        if should_skip_alignment(&read) {
            skipped_alignments += 1;
            continue;
        }

        let qname = String::from_utf8_lossy(read.qname()).to_string();

        // Check if we've moved to a new qname
        // as_ref(): converts &Option<String> to Option<&String> for comparison
        if current_qname.as_ref() != Some(&qname) {
            // Process the previous qname group if any
            if let Some(prev_qname) = current_qname.take() {
                process_qname_group(
                    &mut result,
                    &header,
                    prev_qname,
                    &mut current_reads,
                    &mut qname_idx_counter,
                    mapq_filter,
                    basequal_median_filter,
                    filter_noisy,
                )?; // ? operator: propagates any error from process_qname_group
            }

            current_qname = Some(qname.clone());
            current_reads.clear();
        }

        // Skip if this qname is already marked as noisy
        if result.noisy_qnames.contains_key(&qname) {
            current_reads.clear();
            continue;
        }

        current_reads.push(read);
    }

    // Process the last qname group
    if let Some(qname) = current_qname {
        process_qname_group(
            &mut result,
            &header,
            qname,
            &mut current_reads,
            &mut qname_idx_counter,
            mapq_filter,
            basequal_median_filter,
            filter_noisy,
        )?; // ? operator: propagates any error from process_qname_group
    }

    drop(bam);
    reader_guard.finish()?;

    // Finalize all interval trees
    for interval_tree in result.interval_trees.values_mut() {
        // ? operator: propagates any error from finalize
        interval_tree.finalize()?;
    }

    info!("[migrate_bam_to_sorted_intervals_grouped] BAM processing complete: {} total reads processed, {} alignments skipped, {} read pairs retained, {} noisy qnames filtered",
         total_reads_processed, skipped_alignments, result.readpair_dict.len(), result.noisy_qnames.len());

    Ok((result, header))
}

/// Process all reads for a single qname
fn process_qname_group(
    result: &mut ReadPairMap,
    header: &bam::HeaderView,
    qname: String,
    reads: &mut Vec<Record>,
    qname_idx_counter: &mut usize,
    mapq_filter: u8,
    basequal_median_filter: u8,
    filter_noisy: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    // Check if any read in the group is noisy
    let is_noisy = reads
        .iter()
        .any(|read| is_read_noisy(read, mapq_filter, basequal_median_filter, filter_noisy));

    if is_noisy {
        debug!(
            "[process_qname_group] This qname {} is noisy. Skip it.\n",
            qname
        );
        result.noisy_qnames.insert(qname.clone(), ());
        return Ok(());
    }

    // Filter to get exactly 2 primary reads (read1 and read2)
    let mut read1_opt: Option<Record> = None;
    let mut read2_opt: Option<Record> = None;

    for read in reads.iter() {
        if read.is_first_in_template() {
            read1_opt = Some(read.clone());
        } else if read.is_last_in_template() {
            read2_opt = Some(read.clone());
        }
    }

    // We need both reads for a complete pair
    let (read1, read2) = match (read1_opt, read2_opt) {
        (Some(r1), Some(r2)) => (r1, r2),
        _ => return Ok(()), // Skip incomplete pairs
    };

    // Assign qname_idx and create ReadPair
    // * operator: dereferences the mutable reference to get the actual usize value
    let qname_idx = *qname_idx_counter;
    result.qname_to_idx.insert(qname.clone(), qname_idx);
    result.idx_to_qname.insert(qname_idx, qname.clone());
    // * operator: dereferences the mutable reference to modify the actual usize value
    *qname_idx_counter += 1;

    // Add intervals for both reads
    for read in [&read1, &read2] {
        let tid = read.tid();
        let ref_name = header.tid2name(tid as u32);
        let chrom = String::from_utf8_lossy(ref_name).to_string();
        // ? operator: propagates any error from add_read_interval
        add_read_interval(&mut result.interval_trees, &chrom, read, qname_idx)?;
    }

    // Create and store the complete ReadPair
    let readpair = ReadPair::new_complete(read1, read2, qname, qname_idx);
    result.readpair_dict.insert(qname_idx, readpair);

    Ok(())
}

/// Generalized function to add interval for any read record
fn add_read_interval(
    interval_trees: &mut AHashMap<String, SortedVecIntervals>,
    chrom: &str,
    read: &Record,
    qname_idx: usize,
) -> Result<(), String> {
    let start = read.reference_start();
    let end = read.reference_end();

    if let Some(interval_tree) = interval_trees.get_mut(chrom) {
        interval_tree.add_interval(start, end, qname_idx)
    } else {
        Err(format!("Chromosome {} not found in interval trees", chrom))
    }
}

/// Equivalent to Python's stat_ad_to_dict function
/// Builds allele depth information using bcftools mpileup and processes with Polars
/// Returns a more efficient Rust data structure instead of nested HashMaps
pub fn build_allele_depth_map(
    bam_file: &str,
    reference_genome: &str,
    mapq_filter: u8,
    base_qual_filter: u8,
) -> Result<AlleleDepthMap, Box<dyn std::error::Error>> {
    info!(
        "[build_allele_depth_map] Starting with BAM: {}, MAPQ>={}, BaseQ>={}",
        bam_file, mapq_filter, base_qual_filter
    );

    // Stream mpileup stdout directly into query stdin, and write query stdout to the .ad file
    let mut mpileup_child = Command::new("bcftools")
        .args([
            "mpileup",
            "-Ou",
            "--fasta-ref",
            reference_genome,
            "-a",
            "FORMAT/AD",
            "--indels-2.0",
            "-A",
            "-q",
            &mapq_filter.to_string(),
            "-Q",
            &base_qual_filter.to_string(),
            bam_file,
        ])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;

    // drain mpileup stderr concurrently into a bounded buffer (no logging here to avoid GIL contention)
    let mpileup_stderr_buf: Arc<Mutex<VecDeque<String>>> =
        Arc::new(Mutex::new(VecDeque::with_capacity(500)));
    let mpileup_stderr_buf_reader = Arc::clone(&mpileup_stderr_buf);
    let mpileup_stderr_jh = if let Some(stderr) = mpileup_child.stderr.take() {
        Some(thread::spawn(move || {
            let reader = BufReader::new(stderr);
            for line_res in reader.lines() {
                if let Ok(line) = line_res {
                    let mut buf = mpileup_stderr_buf_reader.lock().unwrap();
                    if buf.len() == buf.capacity() {
                        buf.pop_front();
                    }
                    buf.push_back(line);
                }
            }
        }))
    } else {
        None
    };

    let mpileup_stdout = mpileup_child
        .stdout
        .take()
        .ok_or("Failed to capture mpileup stdout")?;

    let mut query_child = Command::new("bcftools")
        .args(["query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%AD]\\n", "-"])
        .stdin(Stdio::from(mpileup_stdout))
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;

    // drain query stderr concurrently into a bounded buffer (no logging here)
    let query_stderr_buf: Arc<Mutex<VecDeque<String>>> =
        Arc::new(Mutex::new(VecDeque::with_capacity(500)));
    let query_stderr_buf_reader = Arc::clone(&query_stderr_buf);
    let query_stderr_jh = if let Some(stderr) = query_child.stderr.take() {
        Some(thread::spawn(move || {
            let reader = BufReader::new(stderr);
            for line_res in reader.lines() {
                if let Ok(line) = line_res {
                    let mut buf = query_stderr_buf_reader.lock().unwrap();
                    if buf.len() == buf.capacity() {
                        buf.pop_front();
                    }
                    buf.push_back(line);
                }
            }
        }))
    } else {
        None
    };

    // Stream-parse bcftools query output directly into allele_depth_map
    let query_stdout = query_child
        .stdout
        .take()
        .ok_or("Failed to capture bcftools query stdout")?;

    let reader = BufReader::new(query_stdout);
    let mut allele_depth_map = AlleleDepthMap::new();
    let mut line_count = 0usize;
    for line_res in reader.lines() {
        let line = line_res?;
        if line.trim().is_empty() {
            continue;
        }
        line_count += 1;

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }

        let chrom = fields[0];
        // pileup POS is 1-based; `checked_sub(1)` converts to 0-based and skips a
        // malformed POS==0 that would otherwise wrap to u32::MAX in release builds.
        let pos: u32 = match fields[1].parse::<u32>().ok().and_then(|p| p.checked_sub(1)) {
            Some(p) => p,
            None => continue,
        };
        let ref_allele = fields[2];
        let alt_alleles = fields[3];
        let ad_str = fields[4];

        // Even if ALT is empty, still record REF depth/DP (positions with no ALT)

        // Parse AD values. u32 (not u16): the per-allele depths and especially
        // their sum (`total_depth`) can exceed 65535 in ultra-high-coverage
        // pileups, which would overflow u16 (panic in debug, wrap in release).
        let ad_values: Vec<u32> = ad_str.split(',').filter_map(|s| s.parse().ok()).collect();

        if ad_values.is_empty() {
            continue;
        }

        let ref_depth = ad_values[0];
        let total_depth: u32 = ad_values.iter().sum();

        if total_depth == 0 {
            continue;
        }

        // Parse ALT alleles  (filter placeholders and empties)
        let mut alt_list: Vec<&str> = alt_alleles.trim_end_matches(",<*>").split(',').collect();

        // Remove placeholders and empty entries
        alt_list.retain(|&alt| !alt.is_empty() && alt != "<*>");

        // Only report positions with at least one ALT allele after filtering
        if alt_list.is_empty() {
            continue;
        }

        // Create position entry - just a simple array
        let mut position_data = AlleleDepthMap::new_position_data(total_depth);

        // Add reference depth
        AlleleDepthMap::set_allele_depth(
            &mut position_data,
            base_to_index(ref_allele.chars().next().unwrap_or('N')),
            ref_depth,
        );

        // Add alt allele depths
        for (i, &alt_allele) in alt_list.iter().enumerate() {
            if i + 1 < ad_values.len() && !alt_allele.is_empty() {
                let alt_depth = ad_values[i + 1];
                if alt_depth > 0 {
                    // Only handle SNVs for now (same as Python version)
                    if alt_allele.len() == 1 && ref_allele.len() == 1 {
                        let alt_char = alt_allele.chars().next().unwrap();
                        AlleleDepthMap::set_allele_depth(
                            &mut position_data,
                            base_to_index(alt_char),
                            alt_depth,
                        );
                    }
                }
            }
        }

        allele_depth_map.insert(chrom, pos, position_data);
        debug!("[build_allele_depth_map] Inserted position data for {} at position {} with ref {} and alt {} where the ADs are {:?}", chrom, pos, ref_allele, alt_alleles, position_data);
    }
    info!(
        "[build_allele_depth_map] Processed {} lines, created map with {} chromosomes",
        line_count,
        allele_depth_map.chromosome_count()
    );

    // Ensure both processes exited successfully
    let query_status = query_child.wait()?;
    let mpileup_status = mpileup_child.wait()?;
    if let Some(jh) = mpileup_stderr_jh {
        let _ = jh.join();
    }
    if let Some(jh) = query_stderr_jh {
        let _ = jh.join();
    }

    // Flush buffered stderr lines to logger now (on main thread)
    if let Ok(buf) = mpileup_stderr_buf.lock() {
        for line in buf.iter() {
            let lower = line.to_lowercase();
            if lower.contains("error") {
                error!("{}", line);
            } else if lower.contains("warn") {
                warn!("{}", line);
            } else {
                info!("{}", line);
            }
        }
    }
    if let Ok(buf) = query_stderr_buf.lock() {
        for line in buf.iter() {
            let lower = line.to_lowercase();
            if lower.contains("error") {
                error!("{}", line);
            } else if lower.contains("warn") {
                warn!("{}", line);
            } else {
                info!("{}", line);
            }
        }
    }

    if !mpileup_status.success() {
        return Err("bcftools mpileup failed".into());
    }
    if !query_status.success() {
        return Err("bcftools query failed".into());
    }

    Ok(allele_depth_map)
}

/// Convert DNA base to array index (more efficient than HashMap lookups)
#[inline]
/// Canonical base → allele-depth-array index used to BOTH populate and query
/// `PositionAlleleDepth` ([A, T, C, G, N, total]). Keeping a single mapping
/// guarantees the lookup in `is_sequencing_error` indexes the same slot that
/// `build_allele_depth_map` filled.
pub(crate) fn base_to_index(base: char) -> usize {
    match base.to_ascii_uppercase() {
        'A' => 0,
        'T' => 1,
        'C' => 2,
        'G' => 3,
        _ => 4, // N or any other base
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::header::HeaderRecord;
    use rust_htslib::bam::record::{Cigar, CigarString};
    use rust_htslib::bam::{Format, Header, Writer};
    use std::path::Path;

    const PAIRED: u16 = 0x1;
    const PROPER_PAIR: u16 = 0x2;
    const READ1: u16 = 0x40;
    const READ2: u16 = 0x80;
    const SECONDARY: u16 = 0x100;
    const DUPLICATE: u16 = 0x400;
    const SUPPLEMENTARY: u16 = 0x800;

    #[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
    struct RecordPayload {
        tid: i32,
        pos: i64,
        mtid: i32,
        mpos: i64,
        insert_size: i64,
        flags: u16,
        mapq: u8,
        cigar: String,
        sequence: Vec<u8>,
        qualities: Vec<u8>,
    }

    #[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
    struct PairPayload {
        qname: String,
        read1: RecordPayload,
        read2: RecordPayload,
    }

    fn make_record(qname: &[u8], pos: i64, mate_pos: i64, flags: u16, mapq: u8) -> Record {
        let sequence = vec![b'A'; 100];
        let qualities = vec![30; 100];
        let mut record = Record::new();
        record.set(
            qname,
            Some(&CigarString(vec![Cigar::Match(100)])),
            &sequence,
            &qualities,
        );
        record.set_tid(0);
        record.set_pos(pos);
        record.set_mtid(0);
        record.set_mpos(mate_pos);
        record.set_insert_size(mate_pos - pos);
        record.set_flags(flags);
        record.set_mapq(mapq);
        record
    }

    fn write_pairing_fixture(path: &Path) {
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "stream1");
        sq.push_tag(b"LN", 10_000);
        let mut header = Header::new();
        header.push_record(&sq);

        let valid_read1 = PAIRED | PROPER_PAIR | READ1;
        let valid_read2 = PAIRED | PROPER_PAIR | READ2;
        let records = vec![
            make_record(b"alpha", 100, 700, valid_read1, 60),
            make_record(b"beta", 150, 600, valid_read1, 60),
            make_record(b"incomplete", 175, 900, valid_read1, 60),
            make_record(b"noisy", 200, 650, valid_read1, 5),
            make_record(b"duplicate", 225, 750, valid_read1 | DUPLICATE, 60),
            make_record(b"alpha", 250, 700, valid_read1 | SECONDARY, 60),
            make_record(b"alpha", 350, 700, valid_read1 | SUPPLEMENTARY, 60),
            make_record(b"beta", 600, 150, valid_read2, 60),
            make_record(b"noisy", 650, 200, valid_read2, 60),
            make_record(b"alpha", 700, 100, valid_read2, 60),
            make_record(b"duplicate", 750, 225, valid_read2 | DUPLICATE, 60),
        ];

        let mut writer = Writer::from_path(path, &header, Format::Bam).expect("create fixture BAM");
        for record in &records {
            writer.write(record).expect("write fixture record");
        }
    }

    fn record_payload(record: &Record) -> RecordPayload {
        RecordPayload {
            tid: record.tid(),
            pos: record.pos(),
            mtid: record.mtid(),
            mpos: record.mpos(),
            insert_size: record.insert_size(),
            flags: record.flags(),
            mapq: record.mapq(),
            cigar: record.cigar().to_string(),
            sequence: record.seq().as_bytes(),
            qualities: record.qual().to_vec(),
        }
    }

    fn canonical_pairs(result: &ReadPairMap) -> Vec<PairPayload> {
        let mut pairs: Vec<_> = result
            .readpair_dict
            .values()
            .map(|pair| PairPayload {
                qname: pair.qname.clone(),
                read1: record_payload(&pair.read1),
                read2: record_payload(pair.read2.as_ref().expect("complete pair")),
            })
            .collect();
        pairs.sort_unstable();
        pairs
    }

    fn noisy_qnames(result: &ReadPairMap) -> Vec<String> {
        let mut qnames: Vec<_> = result.noisy_qnames.keys().cloned().collect();
        qnames.sort_unstable();
        qnames
    }

    #[test]
    fn median_phred_matches_np_median() {
        // Odd length → the middle element.
        assert_eq!(median_phred(&[10, 20, 30]), 20.0);
        // Even length → mean of the two middles (the parity fix vs order_stat::kth).
        assert_eq!(median_phred(&[15, 16]), 15.5);
        assert_eq!(median_phred(&[30, 10, 20, 40]), 25.0); // unsorted even
        assert_eq!(median_phred(&[]), 0.0);
    }

    #[test]
    fn total_thread_budget_converts_to_htslib_additional_threads() {
        assert_eq!(htslib_additional_threads(0), 0);
        assert_eq!(htslib_additional_threads(1), 0);
        assert_eq!(htslib_additional_threads(2), 1);
        assert_eq!(htslib_additional_threads(4), 3);
    }

    #[test]
    fn collate_scratch_uses_managed_system_temp_directory() {
        let first = CollateScratch::new().expect("first scratch directory");
        let second = CollateScratch::new().expect("second scratch directory");
        assert!(first._directory.path().starts_with(std::env::temp_dir()));
        assert!(second._directory.path().starts_with(std::env::temp_dir()));
        assert_ne!(first.prefix, second.prefix);
        assert!(first._directory.path().is_dir());
        assert!(second._directory.path().is_dir());
    }

    #[test]
    fn pairing_engines_match_on_interleaved_coordinate_sorted_records() {
        let fixture = NamedTempFile::with_suffix(".bam").expect("create fixture path");
        write_pairing_fixture(fixture.path());
        let fixture_path = fixture.path().to_string_lossy();

        let run = |engine| {
            migrate_bam_to_sorted_intervals_grouped(&fixture_path, 10, 15, true, engine, 2)
                .expect("pair fixture records")
                .0
        };
        let temp_file = run(PairingEngine::SamtoolsTempFile);
        let pipe = run(PairingEngine::SamtoolsPipe);
        let rust_memory = run(PairingEngine::RustMemory);

        let expected_qnames = vec!["alpha".to_string(), "beta".to_string()];
        let retained_qnames = |result: &ReadPairMap| {
            let mut qnames: Vec<_> = result.qname_to_idx.keys().cloned().collect();
            qnames.sort_unstable();
            qnames
        };
        assert_eq!(retained_qnames(&temp_file), expected_qnames);
        assert_eq!(retained_qnames(&pipe), expected_qnames);
        assert_eq!(retained_qnames(&rust_memory), expected_qnames);

        assert_eq!(canonical_pairs(&pipe), canonical_pairs(&temp_file));
        assert_eq!(canonical_pairs(&rust_memory), canonical_pairs(&temp_file));
        assert_eq!(noisy_qnames(&temp_file), vec!["noisy"]);
        assert_eq!(noisy_qnames(&pipe), noisy_qnames(&temp_file));
        assert_eq!(noisy_qnames(&rust_memory), noisy_qnames(&temp_file));

        assert_eq!(pipe.qname_to_idx, temp_file.qname_to_idx);
        assert_eq!(pipe.idx_to_qname, temp_file.idx_to_qname);

        for result in [&temp_file, &pipe, &rust_memory] {
            let overlap_ids = result.interval_trees["stream1"]
                .find_overlaps(0, 1_000)
                .expect("query finalized intervals");
            let mut overlap_qnames: Vec<_> = overlap_ids
                .into_iter()
                .map(|idx| result.idx_to_qname[&idx].clone())
                .collect();
            overlap_qnames.sort_unstable();
            assert_eq!(overlap_qnames, expected_qnames);
        }
    }
}
