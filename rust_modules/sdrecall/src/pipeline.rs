//! Orchestration — the SDrecall pipeline stages wired to their validated
//! crate entry points. See `docs/analysis/tasks/T9_orchestrator.md`.
//!
//! Pipeline spine (mirrors the Python `SDrecall` / `realign_and_recall.py` /
//! `misalignment_elimination.py`):
//!
//! ```text
//! prepare:  sd-prep
//! realign:  region-prep → read-extraction → minimap2 + variant-call
//!           → merge + markdup → slice islands
//!           → fp-control (per island, rayon) → filter BAM + variant-call
//!           → vcf-ops priority-merge → subset to target
//! post:     vcf-ops inhouse-common → vcf-ops merge-with-conventional
//! ```

use std::collections::HashSet;
use std::path::{Path, PathBuf};

use rayon::prelude::*;
use sdrecall_utils::{Result, SdError};

use crate::cli::{PrepareArgs, RealignArgs, RunArgs};
use crate::island::IslandPaths;
use crate::paths::{Paths, RgRef};
use crate::rg_discovery::RgInfo;

/// Resolved thread budget for one pipeline run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ThreadBudget {
    pub total_threads: usize,
    pub num_jobs: usize,
    pub threads_per_job: usize,
}

impl ThreadBudget {
    pub fn new(total_threads: usize, threads_per_job: f64) -> Self {
        let (num_jobs, threads_per_job) =
            sdrecall_utils::configure_parallelism(total_threads, threads_per_job);
        ThreadBudget {
            total_threads,
            num_jobs,
            threads_per_job,
        }
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Top-level subcommand orchestration
// ════════════════════════════════════════════════════════════════════════

pub fn run_full_pipeline(args: &RunArgs, paths: &Paths) -> Result<PathBuf> {
    log::info!(
        "[run] full pipeline for sample={} assembly={} target_tag={}",
        paths.sample_id, paths.assembly, paths.target_tag
    );

    ensure_dirs(paths)?;
    prepare(&args.common, &args.prep, paths)?;
    let sdrecall_vcf = realign_and_recall(args, paths)?;
    post_process_vcf(&sdrecall_vcf, &args.conventional, &args.cohort, paths)
}

pub fn run_preparation_only(args: &PrepareArgs, paths: &Paths) -> Result<()> {
    log::info!("[prepare] preparation-only for sample={}", paths.sample_id);
    ensure_dirs(paths)?;
    prepare(&args.common, &args.prep, paths)
}

pub fn run_realign_only(args: &RealignArgs, paths: &Paths) -> Result<PathBuf> {
    log::info!("[realign] realign-only for sample={}", paths.sample_id);
    ensure_dirs(paths)?;
    let sdrecall_vcf = realign_and_recall_inner(
        paths,
        args.common.threads,
        args.realign.numba_threads,
        args.common.mq_cutoff,
    )?;
    post_process_vcf(&sdrecall_vcf, &args.conventional, &args.cohort, paths)
}

fn ensure_dirs(paths: &Paths) -> Result<()> {
    for dir in [
        &paths.work_dir,
        &paths.recall_results_dir,
        &paths.realign_groups_dir,
        &paths.tmp_dir,
    ] {
        std::fs::create_dir_all(dir).map_err(|e| SdError::Io {
            path: dir.display().to_string(),
            source: e,
        })?;
    }
    Ok(())
}

// ════════════════════════════════════════════════════════════════════════
//  Phase: PREPARE
// ════════════════════════════════════════════════════════════════════════

fn prepare(
    common: &crate::cli::CommonArgs,
    prep: &crate::cli::PreparationArgs,
    paths: &Paths,
) -> Result<()> {
    log::info!("[prepare] building recall regions into {:?}", paths.work_dir);
    let prep_paths = paths.to_prep_paths();
    let prep_params = paths.to_prep_params(common, prep);
    let _prep_result = sd_prep::prepare_recall_regions(&prep_paths, &prep_params)?;
    log::info!(
        "[prepare] Phase 1 complete: {} RGs established",
        _prep_result.rg_outputs.len()
    );
    Ok(())
}

// ════════════════════════════════════════════════════════════════════════
//  Phase: REALIGN + RECALL
// ════════════════════════════════════════════════════════════════════════

fn realign_and_recall(args: &RunArgs, paths: &Paths) -> Result<PathBuf> {
    realign_and_recall_inner(
        paths,
        args.common.threads,
        args.realign.numba_threads,
        args.common.mq_cutoff,
    )
}

fn realign_and_recall_inner(
    paths: &Paths,
    threads: usize,
    numba_threads: usize,
    mq_cutoff: i32,
) -> Result<PathBuf> {
    log::info!(
        "[realign] start for sample={} (threads={threads})",
        paths.sample_id
    );

    // ── Step 1: per-RG region size stats ────────────────────────────────
    let rg_infos = stat_realign_group_regions(paths)?;
    let rg_labels: Vec<&str> = rg_infos.iter().map(|r| r.label.as_str()).collect();
    log::info!("[realign] {} RGs discovered: {:?}", rg_labels.len(), rg_labels);

    // ── Step 2: per-RG masked-align region prep ─────────────────────────
    let prep_budget = ThreadBudget::new(threads, 1.0);
    prepare_masked_align_regions(paths, &rg_infos, prep_budget)?;

    // ── Step 3: per-RG realign + variant call ───────────────────────────
    let realign_budget = ThreadBudget::new(threads, 3.0);
    let (per_rg_bams, per_rg_vcfs) = realign_per_rg(paths, &rg_infos, realign_budget)?;

    // ── Step 4: merge per-RG raw BAMs + dedup ───────────────────────────
    merge_and_markdup_raw_bams(paths, &per_rg_bams, threads)?;

    // ── Step 5: concat per-RG raw VCFs ──────────────────────────────────
    concat_raw_vcfs(paths, &per_rg_vcfs, threads)?;

    // ── Step 6: misalignment elimination ────────────────────────────────
    eliminate_misalignments(paths, threads, numba_threads, mq_cutoff)?;

    // ── Step 7: priority-merge raw vs clean, subset to target ───────────
    merge_and_subset_final_vcf(paths, threads)?;

    Ok(paths.final_recall_vcf_path())
}

// ── Step 1: RG discovery ────────────────────────────────────────────────

fn stat_realign_group_regions(paths: &Paths) -> Result<Vec<RgInfo>> {
    // Discover RGs by listing the realign_groups directory for RG* subdirs.
    let rg_dir = &paths.realign_groups_dir;
    let mut rg_labels = Vec::new();
    let mut bed_paths = Vec::new();

    let entries: Vec<_> = std::fs::read_dir(rg_dir)
        .map_err(|e| SdError::Io {
            path: rg_dir.display().to_string(),
            source: e,
        })?
        .filter_map(|e| e.ok())
        .collect();

    for entry in &entries {
        let name = entry.file_name().to_string_lossy().to_string();
        if name.starts_with("RG") && entry.file_type().map(|t| t.is_dir()).unwrap_or(false) {
            if let Ok(bed) = paths.all_homo_regions_bed_path(RgRef::Label(&name)) {
                if bed.exists() {
                    rg_labels.push(name);
                    bed_paths.push(bed);
                }
            }
        }
    }

    crate::rg_discovery::stat_all_rg_region_size(&rg_labels, &bed_paths)
}

// ── Step 2: region-prep ─────────────────────────────────────────────────

fn prepare_masked_align_regions(
    paths: &Paths,
    rg_infos: &[RgInfo],
    budget: ThreadBudget,
) -> Result<()> {
    log::info!(
        "[region-prep] preparing {} RGs ({} parallel jobs)",
        rg_infos.len(),
        budget.num_jobs
    );

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(budget.num_jobs)
        .build()
        .map_err(|e| SdError::Compute(e.to_string()))?;

    let target_bed = paths.target_bed.as_deref().ok_or_else(|| {
        SdError::Compute("target_bed is required for region-prep".into())
    })?;

    let errors: Vec<_> = pool.install(|| {
        rg_infos
            .par_iter()
            .filter_map(|rg| {
                let rg_ref = RgRef::Label(&rg.label);
                let whole_bed = match paths.all_homo_regions_bed_path(rg_ref) {
                    Ok(p) => p,
                    Err(e) => return Some(e),
                };
                let rg_dir = match paths.rg_dir(rg_ref) {
                    Ok(p) => p,
                    Err(e) => return Some(e),
                };
                let fc_out = rg_dir.join(format!("{}.fc_target.bed", rg.label));

                match region_prep::prepare_masked_align_region_per_rg(
                    &rg.label,
                    &rg.subgroup_ids,
                    target_bed,
                    &whole_bed,
                    &paths.ref_genome,
                    &fc_out,
                    &rg_dir,
                ) {
                    Ok(_records) => {
                        log::info!("[region-prep] {} done", rg.label);
                        None
                    }
                    Err(e) => {
                        log::error!("[region-prep] {} failed: {e}", rg.label);
                        Some(e)
                    }
                }
            })
            .collect()
    });

    if let Some(first_err) = errors.into_iter().next() {
        return Err(first_err);
    }
    Ok(())
}

// ── Step 3: per-RG realign + variant call ───────────────────────────────

fn realign_per_rg(
    paths: &Paths,
    rg_infos: &[RgInfo],
    budget: ThreadBudget,
) -> Result<(Vec<PathBuf>, Vec<PathBuf>)> {
    log::info!(
        "[realign] processing {} RGs ({} parallel, {} threads each)",
        rg_infos.len(),
        budget.num_jobs,
        budget.threads_per_job
    );

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(budget.num_jobs)
        .build()
        .map_err(|e| SdError::Compute(e.to_string()))?;

    let input_bam_str = paths.input_bam.to_string_lossy().to_string();
    let ref_genome = paths.ref_genome.clone();
    let tpj = budget.threads_per_job;

    let results: Vec<Result<(PathBuf, PathBuf)>> = pool.install(|| {
        rg_infos
            .par_iter()
            .map(|rg| {
                let rg_ref = RgRef::Label(&rg.label);
                let rg_dir = paths.rg_dir(rg_ref)?;

                // Merge FC + NFC BEDs for read extraction.
                let query_bed = paths.rg_query_bed_path(rg_ref)?;
                let (r1, r2) = paths.rg_realign_fastqs_path(rg_ref)?;
                let raw_bam = paths.rg_raw_masked_bam_path(rg_ref)?;
                let raw_vcf = raw_bam.with_extension("vcf.gz");
                let masked_genome = paths.masked_genome_path(rg_ref)?;

                // 3a: read extraction (BAM → FASTQ).
                let r1_str = r1.to_string_lossy().to_string();
                let r2_str = r2.to_string_lossy().to_string();
                rust_read_extraction::bam_to_fastq(
                    &input_bam_str,
                    &query_bed.to_string_lossy(),
                    &r1_str,
                    &r2_str,
                    false,
                    tpj,
                )
                .map_err(|e| SdError::Compute(format!("{}: read extraction: {e}", rg.label)))?;

                // Also extract multi-aligned reads from the counterpart regions.
                let counter_bed = paths.rg_counterparts_bed_path(rg_ref)?;
                if counter_bed.exists() {
                    let r1c = rg_dir.join(format!("{}.nfc.r1.fastq", rg.label));
                    let r2c = rg_dir.join(format!("{}.nfc.r2.fastq", rg.label));
                    rust_read_extraction::bam_to_fastq(
                        &input_bam_str,
                        &counter_bed.to_string_lossy(),
                        &r1c.to_string_lossy(),
                        &r2c.to_string_lossy(),
                        true, // multi_aligned
                        tpj,
                    )
                    .map_err(|e| SdError::Compute(format!("{}: NFC extraction: {e}", rg.label)))?;

                    // Concatenate FC + NFC reads.
                    append_file(&r1c, &r1)?;
                    append_file(&r2c, &r2)?;
                }

                // 3b: minimap2 realign.
                crate::tools::minimap2_align(&r1, &r2, &masked_genome, &raw_bam, tpj)?;

                // 3c: variant call.
                crate::tools::bcftools_call(&raw_bam, &ref_genome, &raw_vcf, tpj)?;

                log::info!("[realign] {} done → {:?}", rg.label, raw_bam);
                Ok((raw_bam, raw_vcf))
            })
            .collect()
    });

    let mut bams = Vec::new();
    let mut vcfs = Vec::new();
    for r in results {
        let (b, v) = r?;
        bams.push(b);
        vcfs.push(v);
    }
    Ok((bams, vcfs))
}

fn append_file(src: &Path, dst: &Path) -> Result<()> {
    use std::io::{Read, Write};
    let mut s = std::fs::File::open(src).map_err(|e| SdError::Io {
        path: src.display().to_string(),
        source: e,
    })?;
    let mut d = std::fs::OpenOptions::new()
        .append(true)
        .open(dst)
        .map_err(|e| SdError::Io {
            path: dst.display().to_string(),
            source: e,
        })?;
    let mut buf = Vec::new();
    s.read_to_end(&mut buf).map_err(|e| SdError::Io {
        path: src.display().to_string(),
        source: e,
    })?;
    d.write_all(&buf).map_err(|e| SdError::Io {
        path: dst.display().to_string(),
        source: e,
    })?;
    Ok(())
}

// ── Step 4: merge + markdup ─────────────────────────────────────────────

fn merge_and_markdup_raw_bams(
    paths: &Paths,
    per_rg_bams: &[PathBuf],
    threads: usize,
) -> Result<()> {
    log::info!("[merge] merging {} per-RG BAMs + markdup", per_rg_bams.len());
    let refs: Vec<&Path> = per_rg_bams.iter().map(|p| p.as_path()).collect();

    let pooled = paths.pooled_raw_bam_path();
    crate::tools::samtools_merge(&refs, &pooled, threads)?;

    let deduped = paths.deduped_raw_bam_path();
    crate::tools::samtools_markdup_pipeline(&pooled, &deduped, threads)?;

    log::info!("[merge] markdup done → {:?}", deduped);
    Ok(())
}

// ── Step 5: concat raw VCFs ─────────────────────────────────────────────

fn concat_raw_vcfs(paths: &Paths, per_rg_vcfs: &[PathBuf], threads: usize) -> Result<()> {
    log::info!("[concat] concatenating {} per-RG VCFs", per_rg_vcfs.len());
    let refs: Vec<&Path> = per_rg_vcfs.iter().map(|p| p.as_path()).collect();
    sdrecall_io::concat_sort_vcfs(&refs, &paths.recall_raw_vcf_path(), true, threads as u8)?;
    log::info!("[concat] raw VCF → {:?}", paths.recall_raw_vcf_path());
    Ok(())
}

// ════════════════════════════════════════════════════════════════════════
//  Misalignment elimination (fp-control per-island fan-out)
// ════════════════════════════════════════════════════════════════════════

fn eliminate_misalignments(
    paths: &Paths,
    threads: usize,
    numba_threads: usize,
    mq_cutoff: i32,
) -> Result<()> {
    log::info!(
        "[fp-control] misalignment elimination for {}",
        paths.sample_id
    );

    // Pre-filter intrinsic BAM to target regions.
    let target_bed = paths.target_bed.as_deref().ok_or_else(|| {
        SdError::Compute("target_bed required for fp-control".into())
    })?;
    let filtered_intrin = paths.tmp_dir.join("intrinsic.filtered.bam");
    crate::tools::samtools_view_region(
        &paths.total_intrinsic_bam_path(),
        target_bed,
        &filtered_intrin,
        threads,
    )?;

    // Load chrom sizes from the reference .fai for interval padding.
    let chrom_sizes = load_chrom_sizes(&paths.ref_genome_fai_path())?;

    // Slice into islands.
    let deduped = paths.deduped_raw_bam_path();
    let islands = crate::island::split_bams_into_islands(
        &deduped,
        &filtered_intrin,
        target_bed,
        &chrom_sizes,
        threads,
        &paths.tmp_dir.join("islands"),
    )?;
    log::info!("[fp-control] {} islands to process", islands.len());

    if islands.is_empty() {
        // No islands → copy raw BAM/VCF as filtered.
        std::fs::copy(deduped, paths.pooled_filtered_bam_path()).map_err(|e| SdError::Io {
            path: paths.pooled_filtered_bam_path().display().to_string(),
            source: e,
        })?;
        std::fs::copy(paths.recall_raw_vcf_path(), paths.recall_filtered_vcf_path()).map_err(
            |e| SdError::Io {
                path: paths.recall_filtered_vcf_path().display().to_string(),
                source: e,
            },
        )?;
        return Ok(());
    }

    // NM Poisson cutoff (optional — errors are non-fatal).
    let nm_cutoff = nm_stats::nm_distribution_poisson(&deduped, 0.01, 10000, threads as u8)
        .map(|c| {
            log::info!("[nm-stats] NM cutoff = {} (mean {:.2})", c.cutoff, c.mean);
            c.cutoff
        })
        .unwrap_or_else(|e| {
            log::warn!("[nm-stats] failed ({e}), using default cutoff 0 (no NM filter)");
            0
        });
    let _ = nm_cutoff; // TODO: wire into per-island params when inspect uses it

    // Per-island fp-control (rayon parallel).
    let island_budget = ThreadBudget::new(threads, numba_threads as f64);
    let (clean_bams, clean_vcfs) =
        fp_control_per_island(paths, &islands, island_budget, mq_cutoff)?;

    // Merge per-island outputs.
    merge_island_outputs(paths, &clean_bams, &clean_vcfs, threads)?;

    Ok(())
}

fn fp_control_per_island(
    paths: &Paths,
    islands: &[IslandPaths],
    budget: ThreadBudget,
    mq_cutoff: i32,
) -> Result<(Vec<PathBuf>, Vec<PathBuf>)> {
    log::info!(
        "[fp-control] processing {} islands ({} parallel, {} threads each)",
        islands.len(),
        budget.num_jobs,
        budget.threads_per_job
    );

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(budget.num_jobs)
        .build()
        .map_err(|e| SdError::Compute(e.to_string()))?;

    let ref_genome_str = paths.ref_genome.to_string_lossy().to_string();
    let tpj = budget.threads_per_job;

    let results: Vec<Option<(PathBuf, PathBuf)>> = pool.install(|| {
        islands
            .par_iter()
            .map(|island| {
                // Wrap in catch_unwind per ROB-3.
                let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    process_one_island(island, &ref_genome_str, mq_cutoff as u8, tpj)
                }));

                match result {
                    Ok(Ok(Some(paths))) => Some(paths),
                    Ok(Ok(None)) => {
                        log::debug!("[fp-control] island {} skipped (≤2 haplotypes)", island.id);
                        None
                    }
                    Ok(Err(e)) => {
                        log::error!("[fp-control] island {} failed: {e}", island.id);
                        None
                    }
                    Err(panic) => {
                        let msg = panic
                            .downcast_ref::<&str>()
                            .map(|s| s.to_string())
                            .or_else(|| panic.downcast_ref::<String>().cloned())
                            .unwrap_or_else(|| "unknown panic".to_string());
                        log::error!("[fp-control] island {} panicked: {msg}", island.id);
                        None
                    }
                }
            })
            .collect()
    });

    let mut bams = Vec::new();
    let mut vcfs = Vec::new();
    for opt in results.into_iter().flatten() {
        bams.push(opt.0);
        vcfs.push(opt.1);
    }
    Ok((bams, vcfs))
}

fn process_one_island(
    island: &IslandPaths,
    ref_genome: &str,
    mq_cutoff: u8,
    threads: usize,
) -> Result<Option<(PathBuf, PathBuf)>> {
    let bam_str = island.raw_bam.to_string_lossy().to_string();
    let intrin_str = island.intrinsic_bam.to_string_lossy().to_string();

    let params = fp_control::FpControlParams {
        reference_genome: ref_genome.to_string(),
        mapq_cutoff: mq_cutoff,
        threads: threads as u8,
        ..Default::default()
    };

    let output = fp_control::run_fp_control(&bam_str, &intrin_str, &params)?;

    let Some(output) = output else {
        return Ok(None);
    };

    // Filter BAM to keep only correct reads.
    let correct_set: HashSet<String> = output.correct_qnames.into_iter().collect();
    let clean_bam = island.raw_bam.with_extension("clean.bam");
    crate::bam_filter::filter_bam_by_qnames(
        &island.raw_bam,
        &correct_set,
        &clean_bam,
        threads,
    )?;

    // Variant-call on the clean BAM.
    let clean_vcf = clean_bam.with_extension("vcf.gz");
    crate::tools::bcftools_call(
        &clean_bam,
        Path::new(ref_genome),
        &clean_vcf,
        threads,
    )?;

    log::info!(
        "[fp-control] island {} done: {} correct, {} mismap",
        island.id,
        correct_set.len(),
        output.mismap_qnames.len()
    );

    Ok(Some((clean_bam, clean_vcf)))
}

fn merge_island_outputs(
    paths: &Paths,
    clean_bams: &[PathBuf],
    clean_vcfs: &[PathBuf],
    threads: usize,
) -> Result<()> {
    if clean_bams.is_empty() {
        return Err(SdError::Compute(
            "no clean BAMs from fp-control — all islands failed or skipped".into(),
        ));
    }

    log::info!(
        "[merge-islands] merging {} clean BAMs + {} clean VCFs",
        clean_bams.len(),
        clean_vcfs.len()
    );

    let bam_refs: Vec<&Path> = clean_bams.iter().map(|p| p.as_path()).collect();
    crate::tools::samtools_merge(&bam_refs, &paths.pooled_filtered_bam_path(), threads)?;

    let vcf_refs: Vec<&Path> = clean_vcfs.iter().map(|p| p.as_path()).collect();
    sdrecall_io::concat_sort_vcfs(&vcf_refs, &paths.recall_filtered_vcf_path(), true, threads as u8)?;

    Ok(())
}

// ── Step 7: priority-merge + subset ─────────────────────────────────────

fn merge_and_subset_final_vcf(paths: &Paths, threads: usize) -> Result<()> {
    log::info!("[merge-final] priority merge raw vs clean → final VCF");

    let merged_vcf = paths.merged_recall_vcf_path();

    vcf_ops::merge_with_priority(vcf_ops::MergeParams {
        query_vcf: &paths.recall_raw_vcf_path(),
        reference_vcf: &paths.recall_filtered_vcf_path(),
        output_vcf: &merged_vcf,
        ref_genome: &paths.ref_genome,
        added_filter: Some("MISALIGNED"),
        qv_tag: Some("RAW"),
        rv_tag: Some("CLEAN"),
        modify_gt: false,
        threads: threads as u8,
        tmp_dir: &paths.tmp_dir,
    })?;

    // Subset to target recall regions.
    let target_bed = &paths.total_recall_sd_region_bed_path();
    crate::tools::bcftools_view_regions(
        &merged_vcf,
        target_bed,
        &paths.final_recall_vcf_path(),
        threads,
    )?;

    log::info!("[merge-final] final VCF → {:?}", paths.final_recall_vcf_path());
    Ok(())
}

// ════════════════════════════════════════════════════════════════════════
//  Phase: POST-PROCESS
// ════════════════════════════════════════════════════════════════════════

fn post_process_vcf(
    sdrecall_vcf: &Path,
    conventional: &crate::cli::ConventionalVcfArgs,
    cohort: &crate::cli::CohortArgs,
    paths: &Paths,
) -> Result<PathBuf> {
    let mut final_vcf = sdrecall_vcf.to_path_buf();

    if let Some(cohort_vcf) = &cohort.cohort_vcf {
        log::info!("[post] annotating inhouse-common against {cohort_vcf:?}");
        let annotated = paths
            .tmp_dir
            .join(format!("{}.inhouse_common.vcf.gz", paths.sample_id));

        vcf_ops::annotate_inhouse_common(vcf_ops::InhouseParams {
            query_vcf: &final_vcf,
            cohort_vcf: cohort_vcf.as_path(),
            output_vcf: &annotated,
            ref_genome: &paths.ref_genome,
            added_filter: "INHOUSE_COMMON",
            inhouse_common_cutoff: cohort.inhouse_common_cutoff,
            conf_level: cohort.cohort_conf_level,
            threads: 4, // light task
            tmp_dir: &paths.tmp_dir,
        })?;
        final_vcf = annotated;
    }

    if let Some(conventional_vcf) = &conventional.conventional_vcf {
        log::info!("[post] merging with conventional VCF {conventional_vcf:?}");
        let merged = conventional
            .merged_vcf
            .clone()
            .unwrap_or_else(|| {
                let stem = conventional_vcf
                    .file_stem()
                    .unwrap_or_default()
                    .to_string_lossy();
                conventional_vcf
                    .with_file_name(format!("{stem}.sdrecall_merged.vcf.gz"))
            });

        vcf_ops::merge_with_priority(vcf_ops::MergeParams {
            query_vcf: &final_vcf,
            reference_vcf: conventional_vcf.as_path(),
            output_vcf: &merged,
            ref_genome: &paths.ref_genome,
            added_filter: None,
            qv_tag: Some("SDrecall"),
            rv_tag: Some(&conventional.caller_name),
            modify_gt: false,
            threads: 4,
            tmp_dir: &paths.tmp_dir,
        })?;
        final_vcf = merged;
    }

    Ok(final_vcf)
}

// ════════════════════════════════════════════════════════════════════════
//  Helpers
// ════════════════════════════════════════════════════════════════════════

/// Load chromosome sizes from a `.fai` file.
fn load_chrom_sizes(fai: &Path) -> Result<ahash::AHashMap<String, i64>> {
    let file = std::fs::File::open(fai).map_err(|e| SdError::Io {
        path: fai.display().to_string(),
        source: e,
    })?;
    let reader = std::io::BufReader::new(file);
    let mut sizes = ahash::AHashMap::new();
    for line in std::io::BufRead::lines(reader) {
        let line = line.map_err(|e| SdError::Io {
            path: fai.display().to_string(),
            source: e,
        })?;
        let mut cols = line.split('\t');
        if let (Some(name), Some(len)) = (cols.next(), cols.next()) {
            if let Ok(l) = len.parse::<i64>() {
                sizes.insert(name.to_string(), l);
            }
        }
    }
    Ok(sizes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thread_budget_mirrors_python_rule() {
        let b = ThreadBudget::new(12, 1.0);
        assert_eq!(b.num_jobs, 12);
        assert_eq!(b.threads_per_job, 1);

        let b = ThreadBudget::new(12, 3.0);
        assert_eq!(b.num_jobs, 4);
        assert_eq!(b.threads_per_job, 3);

        let b = ThreadBudget::new(10, 2.0);
        assert_eq!((b.num_jobs, b.threads_per_job), (5, 2));

        let b = ThreadBudget::new(10, 3.0);
        assert_eq!(b.num_jobs, 4);
        assert_eq!(b.threads_per_job, 3);
    }

    #[test]
    fn rgref_normalization_smoke() {
        assert_eq!(Paths::normalize_rg_label(RgRef::from(0u32)).unwrap(), "RG0");
    }
}
