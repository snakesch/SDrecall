//! Coverage island detection + BAM slicing — port of
//! `realign_recall/slice_bam_by_cov.py::split_bam_by_cov`.
//!
//! Given a deduplicated raw BAM and the target recall BED, this module:
//! 1. Computes per-base depth via `samtools depth` (subprocess)
//! 2. Extracts continuous coverage blocks (min_depth ≥ 3)
//! 3. Processes target regions: large (>10 kbp) split by coverage, small
//!    (<2 kbp) merged with nearby coverage islands, medium padded
//! 4. Slices the BAM into per-island chunks via `samtools view -L`
//!
//! The result is a `Vec<IslandPaths>` used by the fp-control per-island
//! fan-out in `pipeline.rs`.

use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use rayon::prelude::*;
use sdrecall_utils::{GenomicInterval, Result, SdError};

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
pub fn split_bams_into_islands(
    raw_bam: &Path,
    intrinsic_bam: &Path,
    target_bed: &Path,
    chrom_sizes: &ahash::AHashMap<String, i64>,
    threads: usize,
    tmp_dir: &Path,
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

    // Step 5: slice BAMs per island. The original loop launched thousands of
    // short subprocesses serially. Run two unchanged command sequences at once;
    // keeping the original thread argument preserves samtools headers and BGZF
    // settings for byte-parity testing. IndexedParallelIterator::collect
    // preserves input order, so scheduling cannot reorder the result vector.
    let slice_start = Instant::now();
    let num_jobs = padded.len().clamp(1, 2);
    let threads_per_job = threads;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_jobs)
        .build()
        .map_err(|e| SdError::Compute(e.to_string()))?;
    let slices: Vec<Result<IslandPaths>> = pool.install(|| {
        padded
            .par_iter()
            .enumerate()
            .map(|(i, island)| {
                let id = i + 1;
                let island_bed = tmp_dir.join(format!("island_{id}.bed"));
                write_intervals(&island_bed, std::slice::from_ref(island))?;

                let island_raw = path_with_island(raw_bam, id);
                let island_intrin = path_with_island(intrinsic_bam, id);

                crate::tools::samtools_view_region(
                    raw_bam,
                    &island_bed,
                    &island_raw,
                    threads_per_job,
                )?;
                crate::tools::samtools_view_region(
                    intrinsic_bam,
                    &island_bed,
                    &island_intrin,
                    threads_per_job,
                )?;

                // Get actual coverage from the sliced BAM (captures distant mates).
                let actual_cov = tmp_dir.join(format!("island_{id}.actual_cov.bed"));
                actual_coverage_bed(&island_raw, &actual_cov, threads_per_job)?;

                Ok(IslandPaths {
                    id,
                    raw_bam: island_raw,
                    intrinsic_bam: island_intrin,
                    coverage_bed: actual_cov,
                })
            })
            .collect()
    });
    let mut result = Vec::with_capacity(slices.len());
    for slice in slices {
        result.push(slice?);
    }
    let slice_time = slice_start.elapsed();

    // Sort islands by raw BAM file size (largest first) for load balancing.
    result.sort_by(|a, b| {
        let sa = std::fs::metadata(&a.raw_bam).map(|m| m.len()).unwrap_or(0);
        let sb = std::fs::metadata(&b.raw_bam).map(|m| m.len()).unwrap_or(0);
        sb.cmp(&sa)
    });

    // Filter out islands where BAM is empty/missing.
    result.retain(|ip| {
        std::fs::metadata(&ip.raw_bam)
            .map(|m| m.len() > 1000)
            .unwrap_or(false)
    });

    log::warn!(
        concat!(
            "[island_slice_metrics] islands={} jobs={} threads_per_job={} ",
            "t_depth_s={:.3} t_slice_s={:.3} t_total_s={:.3}"
        ),
        result.len(),
        num_jobs,
        threads_per_job,
        depth_time.as_secs_f64(),
        slice_time.as_secs_f64(),
        total_start.elapsed().as_secs_f64()
    );

    Ok(result)
}

// ─────────────────────────── depth blocks ────────────────────────────────

/// Extract continuous coverage blocks from a `samtools depth -a` file.
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
