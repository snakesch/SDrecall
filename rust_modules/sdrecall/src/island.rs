//! Coverage island detection + BAM slicing — port of
//! `realign_recall/slice_bam_by_cov.py::split_bam_by_cov`.
//!
//! Given a deduplicated raw BAM and the target recall BED, this module:
//! 1. Computes sparse covered-position depth via `samtools depth` (subprocess)
//! 2. Extracts continuous coverage blocks (min_depth ≥ 3)
//! 3. Processes target regions: large (>10 kbp) split by coverage, small
//!    (<2 kbp) merged with nearby coverage islands, medium padded
//! 4. Slices the BAM into per-island chunks, either with one-pass in-process
//!    routing or the retained `samtools view -L` compatibility engine
//!
//! The result is a `Vec<IslandPaths>` used by the fp-control per-island
//! fan-out in `pipeline.rs`.

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

use rayon::prelude::*;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::{self, ext::BamRecordExtensions, Read};
use rust_htslib::tpool::ThreadPool;
use sdrecall_utils::{GenomicInterval, Result, SdError};

const SLICE_ENGINE_ENV: &str = "SDRECALL_ISLAND_SLICE_ENGINE";

/// Selects the implementation used to create per-island BAMs.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum IslandSliceEngine {
    /// One sequential read of each source BAM, with in-process routing and coverage.
    Fused,
    /// Original per-island `samtools` and `bedtools` subprocess pipelines.
    Legacy,
}

impl IslandSliceEngine {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Fused => "fused",
            Self::Legacy => "legacy",
        }
    }
}

pub(crate) fn configured_slice_engine() -> Result<IslandSliceEngine> {
    match std::env::var(SLICE_ENGINE_ENV) {
        Ok(value) if value.eq_ignore_ascii_case("fused") => Ok(IslandSliceEngine::Fused),
        Ok(value) if value.eq_ignore_ascii_case("legacy") => Ok(IslandSliceEngine::Legacy),
        Ok(value) => Err(SdError::Compute(format!(
            "invalid {SLICE_ENGINE_ENV}={value:?}; expected 'fused' or 'legacy'"
        ))),
        Err(std::env::VarError::NotPresent) => Ok(IslandSliceEngine::Fused),
        Err(error) => Err(SdError::Compute(format!(
            "cannot read {SLICE_ENGINE_ENV}: {error}"
        ))),
    }
}

/// Paths for one coverage island.
#[derive(Clone, Debug)]
pub struct IslandPaths {
    pub id: usize,
    pub raw_bam: PathBuf,
    pub intrinsic_bam: PathBuf,
    pub coverage_bed: PathBuf,
}

/// Split both the raw BAM and intrinsic BAM into per-island chunks.
///
/// `avg_frag_size` is used for the delimiter: Python uses
/// `ceil(avg_frag_size * 1.5)` as the slop, but the actual
/// `split_bam_by_cov` uses a fixed `delimiter_size = 1000`.
pub(crate) fn split_bams_into_islands(
    raw_bam: &Path,
    intrinsic_bam: &Path,
    target_bed: &Path,
    chrom_sizes: &ahash::AHashMap<String, i64>,
    threads: usize,
    tmp_dir: &Path,
    engine: IslandSliceEngine,
) -> Result<Vec<IslandPaths>> {
    let total_start = Instant::now();
    std::fs::create_dir_all(tmp_dir).map_err(|e| SdError::Io {
        path: tmp_dir.display().to_string(),
        source: e,
    })?;

    // Step 1: coverage depth.
    let depth_start = Instant::now();
    let depth_file = tmp_dir.join("raw.depth");
    crate::tools::samtools_depth(raw_bam, &depth_file, threads)?;
    let depth_time = depth_start.elapsed();

    // Step 2: extract coverage blocks (depth ≥ 3).
    let cov_bed_path = tmp_dir.join("raw.cov.bed");
    let cov_blocks = extract_depth_blocks(&depth_file, 3)?;
    write_intervals(&cov_bed_path, &cov_blocks)?;

    // Step 3: process target regions with coverage data.
    let islands = process_target_regions(target_bed, &cov_blocks, chrom_sizes, 1000)?;
    if islands.is_empty() {
        log::warn!("[island] no coverage islands found — skipping fp-control");
        return Ok(Vec::new());
    }
    log::info!("[island] {} coverage islands detected", islands.len());

    // Step 4: pad each island coverage BED by 1000 bp (Python: slop + merge).
    let padded = sdrecall_io::slop(&islands, 1000, chrom_sizes)?;
    let padded = sdrecall_io::sort_merge_bed(&padded, false);

    // Step 5: materialize the paths and slice both source BAMs.
    let slice_start = Instant::now();
    let plans = prepare_island_plans(raw_bam, intrinsic_bam, &padded, tmp_dir)?;
    let (mut result, slice_jobs) = match engine {
        IslandSliceEngine::Fused => (
            slice_bams_fused(raw_bam, intrinsic_bam, &plans, threads)?,
            1,
        ),
        IslandSliceEngine::Legacy => {
            let jobs = plans.len().clamp(1, 2);
            (
                slice_bams_legacy(raw_bam, intrinsic_bam, &plans, threads, jobs)?,
                jobs,
            )
        }
    };
    let slice_time = slice_start.elapsed();

    // Sort islands by raw BAM file size (largest first) for load balancing.
    result.sort_by(|a, b| {
        let sa = std::fs::metadata(&a.raw_bam).map(|m| m.len()).unwrap_or(0);
        let sb = std::fs::metadata(&b.raw_bam).map(|m| m.len()).unwrap_or(0);
        sb.cmp(&sa)
    });

    // Preserve the legacy file-size filter only for the legacy engine. The fused
    // engine knows every output was created and indexed; applying a compressed
    // byte-size threshold would make retention depend on header compression.
    if engine == IslandSliceEngine::Legacy {
        result.retain(|ip| {
            std::fs::metadata(&ip.raw_bam)
                .map(|m| m.len() > 1000)
                .unwrap_or(false)
        });
    }

    log::warn!(
        concat!(
            "[island_slice_metrics] engine={} islands={} jobs={} threads={} ",
            "t_depth_s={:.3} t_slice_s={:.3} t_total_s={:.3}"
        ),
        engine.as_str(),
        result.len(),
        slice_jobs,
        threads,
        depth_time.as_secs_f64(),
        slice_time.as_secs_f64(),
        total_start.elapsed().as_secs_f64()
    );

    Ok(result)
}

#[derive(Clone, Debug)]
struct IslandPlan {
    id: usize,
    interval: GenomicInterval,
    bed: PathBuf,
    raw_bam: PathBuf,
    intrinsic_bam: PathBuf,
    coverage_bed: PathBuf,
}

fn prepare_island_plans(
    raw_bam: &Path,
    intrinsic_bam: &Path,
    intervals: &[GenomicInterval],
    tmp_dir: &Path,
) -> Result<Vec<IslandPlan>> {
    intervals
        .iter()
        .enumerate()
        .map(|(index, interval)| {
            let id = index + 1;
            let island_bed = tmp_dir.join(format!("island_{id}.bed"));
            write_intervals(&island_bed, std::slice::from_ref(interval))?;
            Ok(IslandPlan {
                id,
                interval: interval.clone(),
                bed: island_bed,
                raw_bam: path_with_island(raw_bam, id),
                intrinsic_bam: path_with_island(intrinsic_bam, id),
                coverage_bed: tmp_dir.join(format!("island_{id}.actual_cov.bed")),
            })
        })
        .collect()
}

fn slice_bams_legacy(
    raw_bam: &Path,
    intrinsic_bam: &Path,
    plans: &[IslandPlan],
    threads: usize,
    jobs: usize,
) -> Result<Vec<IslandPaths>> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(jobs)
        .build()
        .map_err(|error| SdError::Compute(error.to_string()))?;
    let slices: Vec<Result<IslandPaths>> = pool.install(|| {
        plans
            .par_iter()
            .map(|plan| {
                crate::tools::samtools_view_region(raw_bam, &plan.bed, &plan.raw_bam, threads)?;
                crate::tools::samtools_view_region(
                    intrinsic_bam,
                    &plan.bed,
                    &plan.intrinsic_bam,
                    threads,
                )?;
                actual_coverage_bed(&plan.raw_bam, &plan.coverage_bed, threads)?;
                Ok(plan.paths())
            })
            .collect()
    });
    slices.into_iter().collect()
}

fn slice_bams_fused(
    raw_bam: &Path,
    intrinsic_bam: &Path,
    plans: &[IslandPlan],
    threads: usize,
) -> Result<Vec<IslandPaths>> {
    let samtools_version = samtools_version()?;
    let raw_start = Instant::now();
    let raw_outputs: Vec<&Path> = plans.iter().map(|plan| plan.raw_bam.as_path()).collect();
    let raw_partition = partition_bam_one_pass(
        raw_bam,
        &raw_outputs,
        plans,
        threads,
        true,
        &samtools_version,
    )?;
    let raw_time = raw_start.elapsed();

    let intrinsic_start = Instant::now();
    let intrinsic_outputs: Vec<&Path> = plans
        .iter()
        .map(|plan| plan.intrinsic_bam.as_path())
        .collect();
    let intrinsic_partition = partition_bam_one_pass(
        intrinsic_bam,
        &intrinsic_outputs,
        plans,
        threads,
        false,
        &samtools_version,
    )?;
    let intrinsic_time = intrinsic_start.elapsed();

    let coverage_start = Instant::now();
    let coverage = raw_partition.coverage.ok_or_else(|| {
        SdError::Compute("fused raw partition did not return coverage intervals".into())
    })?;
    for (plan, intervals) in plans.iter().zip(&coverage) {
        write_intervals(&plan.coverage_bed, intervals)?;
    }
    let coverage_time = coverage_start.elapsed();

    let index_start = Instant::now();
    let all_outputs: Vec<&Path> = raw_outputs
        .iter()
        .chain(&intrinsic_outputs)
        .copied()
        .collect();
    index_bams(&all_outputs, threads)?;
    let index_time = index_start.elapsed();

    let empty_raw = raw_partition
        .record_counts
        .iter()
        .filter(|&&count| count == 0)
        .count();
    log::warn!(
        concat!(
            "[island_fused_metrics] source_raw_records={} raw_assignments={} ",
            "source_intrinsic_records={} intrinsic_assignments={} empty_raw={} ",
            "t_raw_partition_s={:.3} t_intrinsic_partition_s={:.3} ",
            "t_coverage_write_s={:.3} t_index_s={:.3}"
        ),
        raw_partition.source_records,
        raw_partition.record_counts.iter().sum::<u64>(),
        intrinsic_partition.source_records,
        intrinsic_partition.record_counts.iter().sum::<u64>(),
        empty_raw,
        raw_time.as_secs_f64(),
        intrinsic_time.as_secs_f64(),
        coverage_time.as_secs_f64(),
        index_time.as_secs_f64(),
    );

    Ok(plans.iter().map(IslandPlan::paths).collect())
}

impl IslandPlan {
    fn paths(&self) -> IslandPaths {
        IslandPaths {
            id: self.id,
            raw_bam: self.raw_bam.clone(),
            intrinsic_bam: self.intrinsic_bam.clone(),
            coverage_bed: self.coverage_bed.clone(),
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct RoutedInterval {
    start: i64,
    end: i64,
    island_index: usize,
}

struct IntervalRouter {
    by_tid: Vec<Vec<RoutedInterval>>,
}

impl IntervalRouter {
    fn new(plans: &[IslandPlan], header: &bam::HeaderView) -> Self {
        let mut by_tid = vec![Vec::new(); header.target_count() as usize];
        for (island_index, plan) in plans.iter().enumerate() {
            let Some(tid) = header.tid(plan.interval.chrom.as_bytes()) else {
                log::warn!(
                    "[island] contig '{}' is absent from BAM header",
                    plan.interval.chrom
                );
                continue;
            };
            by_tid[tid as usize].push(RoutedInterval {
                start: plan.interval.start,
                end: plan.interval.end,
                island_index,
            });
        }
        for intervals in &mut by_tid {
            intervals.sort_unstable_by_key(|interval| (interval.start, interval.end));
        }
        Self { by_tid }
    }

    fn intervals_for_tid(&self, tid: i32) -> &[RoutedInterval] {
        usize::try_from(tid)
            .ok()
            .and_then(|tid| self.by_tid.get(tid))
            .map(Vec::as_slice)
            .unwrap_or_default()
    }
}

struct PartitionResult {
    source_records: u64,
    record_counts: Vec<u64>,
    coverage: Option<Vec<Vec<GenomicInterval>>>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct SamtoolsPgLink {
    id: String,
    previous_id: Option<String>,
}

fn partition_bam_one_pass(
    input_bam: &Path,
    output_bams: &[&Path],
    plans: &[IslandPlan],
    threads: usize,
    collect_coverage: bool,
    samtools_version: &str,
) -> Result<PartitionResult> {
    if output_bams.len() != plans.len() {
        return Err(SdError::Compute(format!(
            "island partition output count {} differs from plan count {}",
            output_bams.len(),
            plans.len()
        )));
    }

    let mut reader = bam::Reader::from_path(input_bam).map_err(|error| {
        SdError::Htslib(format!("open source BAM {}: {error}", input_bam.display()))
    })?;
    let header_view = reader.header().clone();
    let pg_links = next_samtools_pg_links(&header_view);
    let router = IntervalRouter::new(plans, &header_view);

    // A threaded BGZF writer owns a dispatcher thread even when many writers
    // share one htslib worker pool. Keeping every island writer threaded would
    // therefore create one OS thread per island while all outputs are open.
    // Parallelize source decompression here; output compression remains serial
    // per writer and preserves the same routed record stream.
    let reader_thread_pool = if threads > 1 {
        let pool = ThreadPool::new(u32::from(sdrecall_utils::clamp_threads_u8(threads - 1)))
            .map_err(|error| SdError::Htslib(format!("create BAM thread pool: {error}")))?;
        reader
            .set_thread_pool(&pool)
            .map_err(|error| SdError::Htslib(format!("set BAM reader thread pool: {error}")))?;
        Some(pool)
    } else {
        None
    };

    let mut writers = Vec::with_capacity(output_bams.len());
    for (plan, output_bam) in plans.iter().zip(output_bams) {
        remove_bam_indexes(output_bam)?;
        let header = samtools_view_header(
            &header_view,
            &pg_links,
            input_bam,
            &plan.bed,
            output_bam,
            threads,
            samtools_version,
        );
        let writer =
            bam::Writer::from_path(output_bam, &header, bam::Format::Bam).map_err(|error| {
                SdError::Htslib(format!(
                    "create island BAM {}: {error}",
                    output_bam.display()
                ))
            })?;
        writers.push(writer);
    }

    let mut source_records = 0u64;
    let mut record_counts = vec![0u64; plans.len()];
    let mut coverage = collect_coverage.then(|| vec![Vec::new(); plans.len()]);
    for result in reader.records() {
        let record = result.map_err(|error| {
            SdError::Htslib(format!("read source BAM {}: {error}", input_bam.display()))
        })?;
        source_records += 1;
        if record.tid() < 0 || record.pos() < 0 {
            continue;
        }
        let contributes_coverage = !record.is_unmapped();
        let start = record.pos();
        let end = record.reference_end().max(start + 1);
        let intervals = router.intervals_for_tid(record.tid());
        for interval in overlapping_intervals(intervals, start, end) {
            writers[interval.island_index]
                .write(&record)
                .map_err(|error| {
                    SdError::Htslib(format!(
                        "write island BAM {}: {error}",
                        output_bams[interval.island_index].display()
                    ))
                })?;
            record_counts[interval.island_index] += 1;
            if let Some(all_coverage) = coverage.as_mut().filter(|_| contributes_coverage) {
                merge_coverage_interval(
                    &mut all_coverage[interval.island_index],
                    &plans[interval.island_index].interval.chrom,
                    start,
                    end,
                );
            }
        }
    }
    drop(reader);
    drop(writers);
    drop(reader_thread_pool);

    Ok(PartitionResult {
        source_records,
        record_counts,
        coverage,
    })
}

fn samtools_version() -> Result<String> {
    let output = Command::new("samtools")
        .arg("--version-only")
        .output()
        .map_err(|source| SdError::Io {
            path: "<samtools --version-only>".into(),
            source,
        })?;
    if !output.status.success() {
        return Err(SdError::Compute(format!(
            "samtools --version-only failed with exit {:?}",
            output.status.code()
        )));
    }
    let full_version = String::from_utf8_lossy(&output.stdout);
    let version = full_version
        .trim()
        .split_once('+')
        .map_or(full_version.trim(), |(version, _)| version);
    if version.is_empty() {
        return Err(SdError::Compute(
            "samtools --version-only returned an empty version".into(),
        ));
    }
    Ok(version.to_string())
}

fn next_samtools_pg_links(header: &bam::HeaderView) -> Vec<SamtoolsPgLink> {
    #[derive(Debug)]
    struct PgRecord {
        id: String,
        previous_id: Option<String>,
    }

    let mut records = Vec::new();
    for line in header.as_bytes().split(|byte| *byte == b'\n') {
        if !line.starts_with(b"@PG\t") {
            continue;
        }
        let mut id = None;
        let mut previous_id = None;
        for field in line.split(|byte| *byte == b'\t').skip(1) {
            if let Some(value) = field.strip_prefix(b"ID:") {
                id = Some(String::from_utf8_lossy(value).into_owned());
            } else if let Some(value) = field.strip_prefix(b"PP:") {
                previous_id = Some(String::from_utf8_lossy(value).into_owned());
            }
        }
        if let Some(id) = id {
            records.push(PgRecord { id, previous_id });
        }
    }

    let mut ids: HashSet<String> = records.iter().map(|record| record.id.clone()).collect();
    let index_by_id: HashMap<&str, usize> = records
        .iter()
        .enumerate()
        .map(|(index, record)| (record.id.as_str(), index))
        .collect();
    let mut is_parent = vec![false; records.len()];
    let mut chain_size = vec![0usize; records.len()];
    for (index, record) in records.iter().enumerate() {
        if let Some(parent) = record
            .previous_id
            .as_deref()
            .and_then(|id| index_by_id.get(id).copied())
        {
            is_parent[parent] = true;
            chain_size[index] = chain_size[parent] + 1;
        }
    }

    let mut tails: Vec<String> = records
        .iter()
        .enumerate()
        .filter(|(index, _)| !is_parent[*index] && chain_size[*index] > 0)
        .map(|(_, record)| record.id.clone())
        .collect();
    if tails.is_empty() {
        if let Some(record) = records.last() {
            tails.push(record.id.clone());
        }
    }

    let link_count = tails.len().max(1);
    let mut links = Vec::with_capacity(link_count);
    for index in 0..link_count {
        let id = if ids.insert("samtools".to_string()) {
            "samtools".to_string()
        } else {
            let mut suffix = 1usize;
            loop {
                let candidate = format!("samtools.{suffix}");
                if ids.insert(candidate.clone()) {
                    break candidate;
                }
                suffix += 1;
            }
        };
        links.push(SamtoolsPgLink {
            id,
            previous_id: tails.get(index).cloned(),
        });
    }
    links
}

fn samtools_view_header(
    input_header: &bam::HeaderView,
    pg_links: &[SamtoolsPgLink],
    input_bam: &Path,
    region_bed: &Path,
    output_bam: &Path,
    threads: usize,
    samtools_version: &str,
) -> bam::Header {
    let command_line = format!(
        "samtools view -@ {threads} -b -L {} -o {} {}",
        region_bed.display(),
        output_bam.display(),
        input_bam.display()
    );
    let mut header = bam::Header::from_template(input_header);
    for link in pg_links {
        let mut record = HeaderRecord::new(b"PG");
        record.push_tag(b"ID", &link.id).push_tag(b"PN", "samtools");
        if let Some(previous_id) = &link.previous_id {
            record.push_tag(b"PP", previous_id);
        }
        record
            .push_tag(b"VN", samtools_version)
            .push_tag(b"CL", &command_line);
        header.push_record(&record);
    }
    header
}

fn overlapping_intervals(
    intervals: &[RoutedInterval],
    start: i64,
    end: i64,
) -> impl Iterator<Item = &RoutedInterval> {
    let first = intervals.partition_point(|interval| interval.end <= start);
    intervals[first..]
        .iter()
        .take_while(move |interval| interval.start < end)
}

fn merge_coverage_interval(merged: &mut Vec<GenomicInterval>, chrom: &str, start: i64, end: i64) {
    match merged.last_mut() {
        Some(previous) if previous.chrom == chrom && start <= previous.end.saturating_add(100) => {
            previous.end = previous.end.max(end);
        }
        _ => merged.push(GenomicInterval::new(chrom, start, end)),
    }
}

fn index_bams(paths: &[&Path], threads: usize) -> Result<()> {
    let jobs = paths.len().min(threads.max(1));
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(jobs.max(1))
        .build()
        .map_err(|error| SdError::Compute(error.to_string()))?;
    pool.install(|| {
        paths.par_iter().try_for_each(|path| {
            bam::index::build(*path, None::<&Path>, bam::index::Type::Bai, 1).map_err(|error| {
                SdError::Htslib(format!("index island BAM {}: {error}", path.display()))
            })
        })
    })
}

fn remove_bam_indexes(bam_path: &Path) -> Result<()> {
    for suffix in ["bai", "csi"] {
        let index_path = PathBuf::from(format!("{}.{}", bam_path.display(), suffix));
        match std::fs::remove_file(&index_path) {
            Ok(()) => {}
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
            Err(source) => {
                return Err(SdError::Io {
                    path: index_path.display().to_string(),
                    source,
                });
            }
        }
    }
    Ok(())
}

// ─────────────────────────── depth blocks ────────────────────────────────

/// Extract continuous coverage blocks from a `samtools depth` file.
/// Coordinate gaps in sparse output represent zero-depth positions and break
/// the current block through the same adjacency check used for dense output.
/// Mirrors `slice_bam_by_cov.py::extract_depth_blocks`.
fn extract_depth_blocks(depth_file: &Path, min_depth: u32) -> Result<Vec<GenomicInterval>> {
    let file = std::fs::File::open(depth_file).map_err(|e| SdError::Io {
        path: depth_file.display().to_string(),
        source: e,
    })?;
    let reader = BufReader::new(file);

    let mut blocks = Vec::new();
    let mut current: Option<(String, i64, i64)> = None; // (chrom, start_0based, last_pos_0based)

    for line in reader.lines() {
        let line = line.map_err(|e| SdError::Io {
            path: depth_file.display().to_string(),
            source: e,
        })?;
        let mut cols = line.split('\t');
        let chrom = cols.next().unwrap_or("");
        let pos: i64 = cols.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        let depth: u32 = cols.next().and_then(|s| s.parse().ok()).unwrap_or(0);

        // samtools depth positions are 1-based; convert to 0-based for BED.
        let pos0 = pos - 1;

        if depth >= min_depth {
            match &mut current {
                Some(cur) if cur.0 == chrom && cur.2 == pos0 - 1 => {
                    cur.2 = pos0; // extend
                }
                _ => {
                    if let Some(cur) = current.take() {
                        blocks.push(GenomicInterval::new(cur.0, cur.1, cur.2 + 1));
                    }
                    current = Some((chrom.to_string(), pos0, pos0));
                }
            }
        } else if let Some(cur) = current.take() {
            blocks.push(GenomicInterval::new(cur.0, cur.1, cur.2 + 1));
        }
    }
    if let Some(cur) = current {
        blocks.push(GenomicInterval::new(cur.0, cur.1, cur.2 + 1));
    }

    Ok(blocks)
}

// ─────────────────────── target region processing ────────────────────────

/// Process target regions against coverage data, mirroring
/// `slice_bam_by_cov.py::process_target_regions_with_coverage`.
fn process_target_regions(
    target_bed: &Path,
    cov_blocks: &[GenomicInterval],
    chrom_sizes: &ahash::AHashMap<String, i64>,
    delimiter_size: i64,
) -> Result<Vec<GenomicInterval>> {
    let target_regions = sdrecall_io::read_bed(target_bed)?;
    // Merge nearby target regions (within 2 × delimiter).
    let slopped = sdrecall_io::slop(&target_regions, delimiter_size, chrom_sizes)?;
    let merged = sdrecall_io::sort_merge_bed(&slopped, false);

    let min_interval: i64 = 2000;
    let max_interval: i64 = 10000;

    let mut processed = Vec::new();
    let mut small_intervals = Vec::new();

    for region in &merged {
        let size = region.end - region.start;
        if size > max_interval {
            // Large: split by coverage intersection.
            let intersected: Vec<_> = cov_blocks
                .iter()
                .filter(|c| c.chrom == region.chrom && c.start < region.end && c.end > region.start)
                .map(|c| {
                    GenomicInterval::new(&c.chrom, c.start.max(region.start), c.end.min(region.end))
                })
                .collect();
            if intersected.is_empty() {
                log::warn!(
                    "[island] no coverage in large region {}:{}-{}",
                    region.chrom,
                    region.start,
                    region.end
                );
            } else {
                processed.extend(intersected);
            }
        } else if size < min_interval {
            small_intervals.push(region.clone());
        } else {
            // Medium: pad by delimiter_size.
            let start = (region.start - delimiter_size).max(0);
            let end = region.end + delimiter_size;
            processed.push(GenomicInterval::new(&region.chrom, start, end));
        }
    }

    // Small intervals: merge with coverage islands.
    for small in &small_intervals {
        let covering: Vec<_> = cov_blocks
            .iter()
            .filter(|c| c.chrom == small.chrom && c.start < small.end && c.end > small.start)
            .collect();
        if covering.is_empty() {
            log::warn!(
                "[island] no coverage for small region {}:{}-{}",
                small.chrom,
                small.start,
                small.end
            );
            continue;
        }
        // Use the covering island if it's within the size bounds.
        for cov in &covering {
            let cov_size = cov.end - cov.start;
            if cov_size <= max_interval {
                processed.push(GenomicInterval::new(&cov.chrom, cov.start, cov.end));
            } else {
                // Cov island too large; pad the small region by delimiter.
                let start = (small.start - delimiter_size).max(0);
                let end = small.end + delimiter_size;
                processed.push(GenomicInterval::new(&small.chrom, start, end));
            }
        }
    }

    // Final sort + merge to deduplicate overlapping islands.
    Ok(sdrecall_io::sort_merge_bed(&processed, false))
}

// ────────────────────────── actual coverage ──────────────────────────────

/// Get actual read-coverage regions from a BAM using `bedtools bamtobed |
/// merge -d 100`. Mirrors the Python's post-slice step.
fn actual_coverage_bed(bam: &Path, output_bed: &Path, threads: usize) -> Result<()> {
    let script = format!(
        "set -o pipefail; \
         bedtools bamtobed -i {bam} | sort -k1,1 -k2,2n | bedtools merge -d 100 > {out}",
        bam = crate::tools::sq(bam),
        out = crate::tools::sq(output_bed),
    );
    let _ = threads; // bedtools doesn't use threads
    crate::tools::run_bash(&script, "actual_coverage_bed")
}

// ──────────────────────────── helpers ────────────────────────────────────

fn write_intervals(path: &Path, ivs: &[GenomicInterval]) -> Result<()> {
    let mut f = std::fs::File::create(path).map_err(|e| SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    for iv in ivs {
        writeln!(f, "{}\t{}\t{}", iv.chrom, iv.start, iv.end).map_err(|e| SdError::Io {
            path: path.display().to_string(),
            source: e,
        })?;
    }
    Ok(())
}

/// Derive a per-island BAM path: `/path/to/foo.bam` → `/path/to/foo.{id}.bam`.
fn path_with_island(bam: &Path, id: usize) -> PathBuf {
    let stem = bam.file_stem().unwrap_or_default().to_string_lossy();
    let dir = bam.parent().unwrap_or(Path::new("."));
    dir.join(format!("{stem}.{id}.bam"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::header::HeaderRecord;
    use rust_htslib::bam::record::{Cigar, CigarString, Record};
    use std::process::Command;

    fn parse_depth(contents: &str) -> Vec<GenomicInterval> {
        let mut file = tempfile::NamedTempFile::new().expect("create depth fixture");
        file.write_all(contents.as_bytes())
            .expect("write depth fixture");
        extract_depth_blocks(file.path(), 3).expect("parse depth fixture")
    }

    #[test]
    fn sparse_depth_gaps_match_dense_zero_rows() {
        let dense = concat!(
            "chr1\t1\t0\n",
            "chr1\t2\t3\n",
            "chr1\t3\t4\n",
            "chr1\t4\t0\n",
            "chr1\t5\t0\n",
            "chr1\t6\t3\n",
            "chr1\t7\t2\n",
            "chr1\t8\t3\n",
            "chr2\t1\t3\n",
            "chr2\t2\t3\n",
        );
        let sparse = concat!(
            "chr1\t2\t3\n",
            "chr1\t3\t4\n",
            "chr1\t6\t3\n",
            "chr1\t7\t2\n",
            "chr1\t8\t3\n",
            "chr2\t1\t3\n",
            "chr2\t2\t3\n",
        );

        let expected = vec![
            GenomicInterval::new("chr1", 1, 3),
            GenomicInterval::new("chr1", 5, 6),
            GenomicInterval::new("chr1", 7, 8),
            GenomicInterval::new("chr2", 0, 2),
        ];

        assert_eq!(parse_depth(dense), expected);
        assert_eq!(parse_depth(sparse), expected);
    }

    #[test]
    fn routed_intervals_use_half_open_overlap_semantics() {
        let intervals = vec![
            RoutedInterval {
                start: 10,
                end: 20,
                island_index: 0,
            },
            RoutedInterval {
                start: 30,
                end: 40,
                island_index: 1,
            },
        ];
        let routed = |start, end| {
            overlapping_intervals(&intervals, start, end)
                .map(|interval| interval.island_index)
                .collect::<Vec<_>>()
        };

        assert_eq!(routed(0, 10), Vec::<usize>::new());
        assert_eq!(routed(9, 11), vec![0]);
        assert_eq!(routed(20, 30), Vec::<usize>::new());
        assert_eq!(routed(19, 31), vec![0, 1]);
        assert_eq!(routed(40, 41), Vec::<usize>::new());
    }

    #[test]
    fn coverage_merge_matches_bedtools_distance_100() {
        let mut merged = Vec::new();
        merge_coverage_interval(&mut merged, "chr1", 100, 150);
        merge_coverage_interval(&mut merged, "chr1", 250, 260);
        merge_coverage_interval(&mut merged, "chr1", 255, 270);
        merge_coverage_interval(&mut merged, "chr1", 371, 380);

        assert_eq!(
            merged,
            vec![
                GenomicInterval::new("chr1", 100, 270),
                GenomicInterval::new("chr1", 371, 380),
            ]
        );
    }

    #[test]
    fn samtools_pg_links_match_htslib_chain_tail_rules() {
        let mut header = bam::Header::new();
        for (id, previous_id) in [
            ("chain-a-root", None),
            ("chain-a-tail", Some("chain-a-root")),
            ("standalone", None),
            ("chain-b-root", None),
            ("chain-b-middle", Some("chain-b-root")),
            ("chain-b-tail", Some("chain-b-middle")),
            ("samtools", None),
            ("samtools.1", None),
        ] {
            let mut record = HeaderRecord::new(b"PG");
            record.push_tag(b"ID", id);
            if let Some(previous_id) = previous_id {
                record.push_tag(b"PP", previous_id);
            }
            header.push_record(&record);
        }

        let links = next_samtools_pg_links(&bam::HeaderView::from_header(&header));
        assert_eq!(
            links,
            vec![
                SamtoolsPgLink {
                    id: "samtools.2".into(),
                    previous_id: Some("chain-a-tail".into()),
                },
                SamtoolsPgLink {
                    id: "samtools.3".into(),
                    previous_id: Some("chain-b-tail".into()),
                },
            ]
        );
    }

    #[test]
    fn samtools_pg_links_handle_empty_and_singleton_headers() {
        let empty = bam::Header::new();
        assert_eq!(
            next_samtools_pg_links(&bam::HeaderView::from_header(&empty)),
            vec![SamtoolsPgLink {
                id: "samtools".into(),
                previous_id: None,
            }]
        );

        let mut singletons = bam::Header::new();
        for id in ["first", "last"] {
            singletons.push_record(HeaderRecord::new(b"PG").push_tag(b"ID", id));
        }
        assert_eq!(
            next_samtools_pg_links(&bam::HeaderView::from_header(&singletons)),
            vec![SamtoolsPgLink {
                id: "samtools".into(),
                previous_id: Some("last".into()),
            }]
        );
    }

    #[test]
    fn fused_slicing_matches_legacy_records_and_coverage() {
        if !command_available("samtools") || !command_available("bedtools") {
            return;
        }

        let legacy_root = tempfile::tempdir().expect("create legacy fixture root");
        let fused_root = tempfile::tempdir().expect("create fused fixture root");
        let intervals = vec![
            GenomicInterval::new("chr1", 10, 20),
            GenomicInterval::new("chr1", 30, 40),
        ];

        let (legacy_raw, legacy_intrinsic, legacy_plans) =
            prepare_fixture(legacy_root.path(), &intervals);
        let (fused_raw, fused_intrinsic, fused_plans) =
            prepare_fixture(fused_root.path(), &intervals);

        let legacy = slice_bams_legacy(&legacy_raw, &legacy_intrinsic, &legacy_plans, 2, 1)
            .expect("run legacy slicing");
        let fused = slice_bams_fused(&fused_raw, &fused_intrinsic, &fused_plans, 2)
            .expect("run fused slicing");

        assert_eq!(legacy.len(), fused.len());
        for (legacy_paths, fused_paths) in legacy.iter().zip(&fused) {
            assert_eq!(
                normalized_header(&legacy_paths.raw_bam, legacy_root.path()),
                normalized_header(&fused_paths.raw_bam, fused_root.path())
            );
            assert_eq!(
                read_bam_records(&legacy_paths.raw_bam),
                read_bam_records(&fused_paths.raw_bam)
            );
            assert_eq!(
                normalized_header(&legacy_paths.intrinsic_bam, legacy_root.path()),
                normalized_header(&fused_paths.intrinsic_bam, fused_root.path())
            );
            assert_eq!(
                read_bam_records(&legacy_paths.intrinsic_bam),
                read_bam_records(&fused_paths.intrinsic_bam)
            );
            assert_eq!(
                std::fs::read(&legacy_paths.coverage_bed).expect("read legacy coverage"),
                std::fs::read(&fused_paths.coverage_bed).expect("read fused coverage")
            );
            assert!(bam_index_path(&fused_paths.raw_bam).is_file());
            assert!(bam_index_path(&fused_paths.intrinsic_bam).is_file());
        }
    }

    #[test]
    #[ignore = "requires a validated pipeline result fixture"]
    fn fused_slicing_matches_real_pipeline_islands() {
        let baseline_root = PathBuf::from(
            std::env::var("SDRECALL_ISLAND_BASELINE_ROOT")
                .expect("set SDRECALL_ISLAND_BASELINE_ROOT"),
        );
        let output_root = PathBuf::from(
            std::env::var("SDRECALL_ISLAND_BENCH_OUTPUT")
                .expect("set SDRECALL_ISLAND_BENCH_OUTPUT"),
        );
        let threads = std::env::var("SDRECALL_ISLAND_BENCH_THREADS")
            .ok()
            .and_then(|value| value.parse().ok())
            .unwrap_or(25);
        std::fs::create_dir_all(&output_root).expect("create real-data output root");

        let baseline_raw = baseline_root
            .join("recall_results")
            .join("HG002.pooled.raw.deduped.bam");
        let baseline_intrinsic = baseline_root
            .join("intermediates")
            .join("intrinsic.filtered.bam");

        let baseline_islands = baseline_root.join("intermediates").join("islands");
        let mut intervals = Vec::new();
        for id in 1.. {
            let path = baseline_islands.join(format!("island_{id}.bed"));
            if !path.is_file() {
                break;
            }
            let mut rows = sdrecall_io::read_bed(&path).expect("read baseline island BED");
            assert_eq!(rows.len(), 1, "expected one interval in {}", path.display());
            intervals.push(rows.remove(0));
        }
        assert!(!intervals.is_empty(), "no baseline island intervals found");

        let legacy_root = output_root.join("legacy");
        let fused_root = output_root.join("fused");
        let (legacy_raw, legacy_intrinsic, legacy_islands) =
            prepare_real_data_sources(&baseline_raw, &baseline_intrinsic, &legacy_root);
        let (fused_raw, fused_intrinsic, fused_islands) =
            prepare_real_data_sources(&baseline_raw, &baseline_intrinsic, &fused_root);
        let legacy_plans =
            prepare_island_plans(&legacy_raw, &legacy_intrinsic, &intervals, &legacy_islands)
                .expect("prepare legacy real-data island plans");
        let fused_plans =
            prepare_island_plans(&fused_raw, &fused_intrinsic, &intervals, &fused_islands)
                .expect("prepare fused real-data island plans");

        let legacy_started = Instant::now();
        let legacy = slice_bams_legacy(&legacy_raw, &legacy_intrinsic, &legacy_plans, threads, 2)
            .expect("run real-data legacy slicing");
        let legacy_elapsed = legacy_started.elapsed();
        let fused_started = Instant::now();
        let fused = slice_bams_fused(&fused_raw, &fused_intrinsic, &fused_plans, threads)
            .expect("run real-data fused slicing");
        let fused_elapsed = fused_started.elapsed();

        assert_eq!(legacy.len(), fused.len());
        for (legacy_paths, fused_paths) in legacy.iter().zip(&fused) {
            assert_eq!(legacy_paths.id, fused_paths.id);
            assert_eq!(
                normalized_header(&legacy_paths.raw_bam, &legacy_root),
                normalized_header(&fused_paths.raw_bam, &fused_root),
                "raw BAM header differs for island {}",
                fused_paths.id
            );
            assert_bam_records_equal(&legacy_paths.raw_bam, &fused_paths.raw_bam);
            assert_eq!(
                normalized_header(&legacy_paths.intrinsic_bam, &legacy_root),
                normalized_header(&fused_paths.intrinsic_bam, &fused_root),
                "intrinsic BAM header differs for island {}",
                fused_paths.id
            );
            assert_bam_records_equal(&legacy_paths.intrinsic_bam, &fused_paths.intrinsic_bam);
            assert_eq!(
                std::fs::read(&legacy_paths.coverage_bed).expect("read legacy coverage"),
                std::fs::read(&fused_paths.coverage_bed).expect("read fused coverage"),
                "coverage differs for island {}",
                fused_paths.id
            );
        }
        eprintln!(
            "real-data slicing matched {} islands: legacy={:.3}s fused={:.3}s speedup={:.2}x",
            fused.len(),
            legacy_elapsed.as_secs_f64(),
            fused_elapsed.as_secs_f64(),
            legacy_elapsed.as_secs_f64() / fused_elapsed.as_secs_f64(),
        );
    }

    fn command_available(command: &str) -> bool {
        Command::new(command).arg("--version").output().is_ok()
    }

    fn prepare_fixture(
        root: &Path,
        intervals: &[GenomicInterval],
    ) -> (PathBuf, PathBuf, Vec<IslandPlan>) {
        let raw_bam = root.join("raw.bam");
        let intrinsic_bam = root.join("intrinsic.bam");
        write_test_bam(&raw_bam);
        write_test_bam(&intrinsic_bam);
        let tmp_dir = root.join("islands");
        std::fs::create_dir(&tmp_dir).expect("create fixture island directory");
        let plans = prepare_island_plans(&raw_bam, &intrinsic_bam, intervals, &tmp_dir)
            .expect("prepare fixture plans");
        (raw_bam, intrinsic_bam, plans)
    }

    fn write_test_bam(path: &Path) {
        let mut header = bam::Header::new();
        header.push_record(
            HeaderRecord::new(b"HD")
                .push_tag(b"VN", "1.6")
                .push_tag(b"SO", "coordinate"),
        );
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 1_000),
        );
        let mut writer =
            bam::Writer::from_path(path, &header, bam::Format::Bam).expect("create fixture BAM");
        for (name, start, length, flags) in [
            ("left", 5, 10, 0),
            ("placed_unmapped", 12, 10, 4),
            ("both", 15, 20, 0),
            ("boundary", 20, 5, 0),
            ("right", 30, 10, 0),
            ("far", 100, 10, 0),
        ] {
            let mut record = Record::new();
            let sequence = vec![b'A'; length];
            let quality = vec![30; length];
            record.set(
                name.as_bytes(),
                Some(&CigarString(vec![Cigar::Match(length as u32)])),
                &sequence,
                &quality,
            );
            record.set_tid(0);
            record.set_pos(start);
            record.set_mapq(60);
            record.set_flags(flags);
            writer.write(&record).expect("write fixture record");
        }
    }

    fn read_bam_records(path: &Path) -> Vec<Record> {
        bam::Reader::from_path(path)
            .expect("open sliced BAM")
            .records()
            .map(|record| record.expect("read sliced BAM record"))
            .collect()
    }

    fn bam_index_path(path: &Path) -> PathBuf {
        PathBuf::from(format!("{}.bai", path.display()))
    }

    fn normalized_header(path: &Path, root: &Path) -> String {
        let reader = bam::Reader::from_path(path).expect("open BAM header");
        String::from_utf8(reader.header().as_bytes().to_vec())
            .expect("BAM header is UTF-8")
            .replace(root.to_str().expect("fixture root is UTF-8"), "<ROOT>")
    }

    fn link_or_copy(source: &Path, destination: &Path) {
        if destination.exists() {
            return;
        }
        if std::fs::hard_link(source, destination).is_err() {
            std::fs::copy(source, destination).expect("copy real-data source fixture");
        }
    }

    fn prepare_real_data_sources(
        baseline_raw: &Path,
        baseline_intrinsic: &Path,
        output_root: &Path,
    ) -> (PathBuf, PathBuf, PathBuf) {
        std::fs::create_dir_all(output_root).expect("create engine output root");
        let islands = output_root.join("islands");
        std::fs::create_dir_all(&islands).expect("create engine island output");
        let raw = output_root.join("raw.bam");
        let intrinsic = output_root.join("intrinsic.bam");
        link_or_copy(baseline_raw, &raw);
        link_or_copy(baseline_intrinsic, &intrinsic);
        (raw, intrinsic, islands)
    }

    fn assert_bam_records_equal(expected: &Path, observed: &Path) {
        let mut expected_reader = bam::Reader::from_path(expected).expect("open expected BAM");
        let mut observed_reader = bam::Reader::from_path(observed).expect("open observed BAM");
        let mut expected_records = expected_reader.records();
        let mut observed_records = observed_reader.records();
        let mut record_number = 0usize;
        loop {
            match (expected_records.next(), observed_records.next()) {
                (None, None) => break,
                (Some(expected_record), Some(observed_record)) => {
                    record_number += 1;
                    let expected_record = expected_record.expect("read expected BAM record");
                    let observed_record = observed_record.expect("read observed BAM record");
                    assert_record_fields_equal(
                        &expected_record,
                        &observed_record,
                        record_number,
                        expected,
                        observed,
                    );
                }
                _ => panic!(
                    "record count differs between {} and {}",
                    expected.display(),
                    observed.display()
                ),
            }
        }
    }

    fn assert_record_fields_equal(
        expected: &Record,
        observed: &Record,
        record_number: usize,
        expected_path: &Path,
        observed_path: &Path,
    ) {
        let context = || {
            format!(
                "record {record_number} differs: {} vs {}",
                expected_path.display(),
                observed_path.display()
            )
        };
        assert_eq!(expected.qname(), observed.qname(), "{}", context());
        assert_eq!(expected.flags(), observed.flags(), "{}", context());
        assert_eq!(expected.tid(), observed.tid(), "{}", context());
        assert_eq!(expected.pos(), observed.pos(), "{}", context());
        assert_eq!(expected.mapq(), observed.mapq(), "{}", context());
        assert_eq!(expected.mtid(), observed.mtid(), "{}", context());
        assert_eq!(expected.mpos(), observed.mpos(), "{}", context());
        assert_eq!(
            expected.insert_size(),
            observed.insert_size(),
            "{}",
            context()
        );
        assert_eq!(
            expected.cigar().to_string(),
            observed.cigar().to_string(),
            "{}",
            context()
        );
        assert_eq!(
            expected.seq().as_bytes(),
            observed.seq().as_bytes(),
            "{}",
            context()
        );
        assert_eq!(expected.qual(), observed.qual(), "{}", context());

        let mut expected_aux = expected.aux_iter();
        let mut observed_aux = observed.aux_iter();
        loop {
            match (expected_aux.next(), observed_aux.next()) {
                (None, None) => break,
                (
                    Some(Ok((expected_tag, expected_value))),
                    Some(Ok((observed_tag, observed_value))),
                ) => {
                    assert_eq!(expected_tag, observed_tag, "{}", context());
                    assert_eq!(expected_value, observed_value, "{}", context());
                }
                (expected_value, observed_value) => panic!(
                    "{}; auxiliary fields differ: {expected_value:?} vs {observed_value:?}",
                    context()
                ),
            }
        }
    }
}
