//! Orchestration skeleton — the SDrecall pipeline stage SEQUENCE as a
//! documented call graph.
//!
//! This is the **scaffold** of T9. Each stage below is a STUB function with a
//! precise `// TODO(T9)` naming the exact validated entry point it will invoke
//! in the integration pass. The stage-crate dependencies are deliberately NOT
//! wired yet (the in-process threading + the external minimap2/bcftools calls
//! are the full T9 job, and a concurrent agent is editing `sd-prep`). The stubs
//! exist so the call graph, the thread budget plumbing, and the per-island
//! parallelism axis are all visible and testable now.
//!
//! Pipeline spine (mirrors `docs/analysis/tasks/T9_orchestrator.md` §"Data flow"
//! and the Python `SDrecall` / `realign_and_recall.py` / `misalignment_elimination.py`):
//!
//! ```text
//! prepare:  sd-prep
//! realign:  region-prep → read-extraction → minimap2 + variant-call
//!           → merge + markdup → slice islands
//!           → fp-control (per island, rayon) → filter BAM + variant-call
//!           → vcf-ops priority-merge → subset to target
//! post:     vcf-ops inhouse-common → vcf-ops merge-with-conventional
//! ```

use sdrecall_utils::{configure_parallelism, Result};

use crate::cli::{PrepareArgs, RealignArgs, RunArgs};
use crate::paths::Paths;

/// Resolved thread budget for one pipeline run.
///
/// Mirrors the Python `configure_parallelism(threads, per_job)` split applied at
/// several points in `realign_and_recall.py` (per-RG prep uses `per_job=1`,
/// realign uses `per_job=3`) and `misalignment_elimination.py` (per-island
/// filtering uses `per_job=numba_threads`). At the orchestrator the `num_jobs`
/// axis becomes the rayon pool width and `threads_per_job` is forwarded to each
/// stage's own htslib/aligner thread pool — the project's
/// `job_num × threads_per_job = total_threads` rule, with no nested
/// oversubscription (T0 DESIGN §5.9).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ThreadBudget {
    pub total_threads: usize,
    pub num_jobs: usize,
    pub threads_per_job: usize,
}

impl ThreadBudget {
    /// Build the budget for a given `threads_per_job` policy using the one
    /// budget unit, [`sdrecall_utils::configure_parallelism`].
    pub fn new(total_threads: usize, threads_per_job: f64) -> Self {
        let (num_jobs, threads_per_job) = configure_parallelism(total_threads, threads_per_job);
        ThreadBudget {
            total_threads,
            num_jobs,
            threads_per_job,
        }
    }
}

// ════════════════════════════════════════════════════════════════════════
//  Top-level subcommand orchestration (Python run_full_pipeline /
//  run_preparation_only / run_realign_only, SDrecall:117-299)
// ════════════════════════════════════════════════════════════════════════

/// `run`: full pipeline — preparation then realignment+recall, then post-process
/// (Python `run_full_pipeline`, `SDrecall:117-183`).
///
/// Returns the path to the final recalled VCF (possibly merged with a
/// conventional VCF).
pub fn run_full_pipeline(args: &RunArgs, paths: &Paths) -> Result<std::path::PathBuf> {
    log::info!(
        "[run] full pipeline for sample={} assembly={} target_tag={}",
        paths.sample_id,
        paths.assembly,
        paths.target_tag
    );

    // Python: skip preparation if check_preparation_validity() (SDrecall:151).
    // TODO(T9): freshness check via sdrecall_io::is_preparation_valid(paths).
    prepare(&args.common, &args.prep, paths)?;

    // Python: skip realign if check_final_vcf_validity() (SDrecall:166).
    // TODO(T9): freshness check via sdrecall_io::is_final_vcf_valid(paths).
    let sdrecall_vcf = realign_and_recall(args, paths)?;

    post_process_vcf(&sdrecall_vcf, &args.conventional, &args.cohort, paths)
}

/// `prepare`: preparation phase only (Python `run_preparation_only`,
/// `SDrecall:186-241`).
pub fn run_preparation_only(args: &PrepareArgs, paths: &Paths) -> Result<()> {
    log::info!("[prepare] preparation-only for sample={}", paths.sample_id);
    prepare(&args.common, &args.prep, paths)
}

/// `realign`: realignment + recall only, assuming preparation already ran
/// (Python `run_realign_only`, `SDrecall:244-299`).
pub fn run_realign_only(args: &RealignArgs, paths: &Paths) -> Result<std::path::PathBuf> {
    log::info!("[realign] realign-only for sample={}", paths.sample_id);
    // Python verifies multiplex_graph_path() exists, else errors (SDrecall:270).
    // TODO(T9): assert paths.multiplex_graph_path() exists before proceeding.

    // `realign`'s flatten layout differs from `run`'s, so adapt the call.
    let sdrecall_vcf = realign_and_recall_inner(
        paths,
        args.common.threads,
        args.realign.numba_threads,
        args.common.mq_cutoff,
    )?;
    post_process_vcf(&sdrecall_vcf, &args.conventional, &args.cohort, paths)
}

// ════════════════════════════════════════════════════════════════════════
//  Phase: PREPARE (Python prepare_recall_regions.py → future sd-prep crate)
// ════════════════════════════════════════════════════════════════════════

/// Preparation phase: build the multi-align BED, the multiplexed SD graph, the
/// per-RG masked genomes and the intrinsic-alignment BAM (Python
/// `prepare_recall_regions`, `prepare_recall_regions.py:57-243`).
fn prepare(
    common: &crate::cli::CommonArgs,
    _prep: &crate::cli::PreparationArgs,
    paths: &Paths,
) -> Result<()> {
    log::info!("[prepare] building recall regions into {:?}", paths.work_dir);
    // TODO(T9): call sd-prep crate (T8) — port of prepare_recall_regions.py.
    //   Entry point TBD when sd-prep stabilises (a concurrent agent owns it);
    //   it produces paths.multi_align_bed_path(), paths.multiplex_graph_path(),
    //   paths.annotated_graph_path(), per-RG masked genomes + minimap indices,
    //   and paths.total_intrinsic_bam_path().
    //   The thread budget here mirrors Python configure_parallelism(threads, 1).
    let _budget = ThreadBudget::new(common.threads, 1.0);
    stub("prepare/sd-prep")
}

// ════════════════════════════════════════════════════════════════════════
//  Phase: REALIGN + RECALL (Python SDrecall_per_sample, realign_and_recall.py)
// ════════════════════════════════════════════════════════════════════════

/// Adapter from the `run` subcommand's argument layout.
fn realign_and_recall(args: &RunArgs, paths: &Paths) -> Result<std::path::PathBuf> {
    realign_and_recall_inner(
        paths,
        args.common.threads,
        args.realign.numba_threads,
        args.common.mq_cutoff,
    )
}

/// The realignment + recall phase proper (Python `SDrecall_per_sample`,
/// `realign_and_recall.py:31-237`).
///
/// Stage sequence (each a stub below):
///   1. stat per-RG region sizes (load balancing)
///   2. per-RG region prep (region-prep, parallel over RGs)
///   3. per-RG read-extraction → minimap2 realign → variant call (parallel)
///   4. merge per-RG raw BAMs → collate/fixmate/sort/markdup
///   5. concat per-RG raw VCFs
///   6. misalignment elimination (the fp-control fan-out — see below)
///   7. priority-merge raw vs clean VCFs, then subset to the target recall BED
fn realign_and_recall_inner(
    paths: &Paths,
    threads: usize,
    numba_threads: usize,
    mq_cutoff: i32,
) -> Result<std::path::PathBuf> {
    log::info!(
        "[realign] start for sample={} (threads={threads})",
        paths.sample_id
    );

    // ── Step 1: per-RG region size stats for load balancing ──────────────
    // Python: stat_all_RG_region_size(sdrecall_paths) (realign_and_recall.py:44).
    // TODO(T9): absorb realign_recall/stat_realign_group_regions.py
    //   (stat_all_RG_region_size) — likely a region-prep helper over the per-RG
    //   homo-region BEDs at paths.all_homo_regions_bed_path(rg).
    stat_realign_group_regions(paths)?;

    // ── Step 2: per-RG masked-align region prep (parallel over RGs) ───────
    // Python: pool.imap_unordered(imap_prepare_masked_align_region_per_RG, ...)
    //   with configure_parallelism(threads, 1) (realign_and_recall.py:58-65).
    let prep_budget = ThreadBudget::new(threads, 1.0);
    log::debug!("[realign] region-prep budget {prep_budget:?}");
    prepare_masked_align_regions(paths, prep_budget)?;

    // ── Step 3: per-RG realign + variant call (parallel over RGs) ─────────
    // Python: pool.imap_unordered(imap_process_masked_bam, ...) with
    //   configure_parallelism(threads, 3) (realign_and_recall.py:100-114).
    let realign_budget = ThreadBudget::new(threads, 3.0);
    log::debug!("[realign] per-RG realign budget {realign_budget:?}");
    realign_per_rg(paths, realign_budget)?;

    // ── Step 4: merge per-RG raw BAMs + dedup (markdup) ──────────────────
    // Python: merge_bams(...) then samtools collate|fixmate|sort|markdup
    //   (realign_and_recall.py:139-153).
    merge_and_markdup_raw_bams(paths, threads)?;

    // ── Step 5: concat per-RG raw VCFs into the pooled raw VCF ────────────
    // Python: bcftools_concatvcfs over the per-RG raw VCFs
    //   (realign_and_recall.py:156-173).
    concat_raw_vcfs(paths, threads)?;

    // ── Step 6: misalignment elimination (the fp-control fan-out) ─────────
    // Python: filter_redundant_intrinsic_origins(...) then
    //   eliminate_misalignments(...) (realign_and_recall.py:178-211).
    eliminate_misalignments(paths, threads, numba_threads, mq_cutoff)?;

    // ── Step 7: priority-merge raw vs clean, then subset to target ────────
    // Python: merge_with_priority(RAW vs CLEAN) then
    //   bcftools view -R total_recall_SD_region_bed (realign_and_recall.py:213-231).
    merge_and_subset_final_vcf(paths, threads)?;

    Ok(paths.final_recall_vcf_path())
}

/// Step 1 stub: per-RG region size stats (Python `stat_all_RG_region_size`).
fn stat_realign_group_regions(_paths: &Paths) -> Result<()> {
    // TODO(T9): absorb realign_recall/stat_realign_group_regions.py
    //   (the thin wrapper this orchestrator subsumes). Reads
    //   paths.all_homo_regions_bed_path(rg) per RG and ranks RGs by NFC size.
    stub("stat_realign_group_regions")
}

/// Step 2 stub: per-RG masked-align region prep, parallel over RGs (Python
/// `imap_prepare_masked_align_region_per_RG`).
fn prepare_masked_align_regions(_paths: &Paths, _budget: ThreadBudget) -> Result<()> {
    // TODO(T9): for each RG (rayon over RGs, _budget.num_jobs workers) call
    //   region_prep::prepare_masked_align_region_per_rg(...)
    //   (rust_modules/region-prep/src/per_rg.rs:338 — validated). Forward
    //   _budget.threads_per_job to its inner htslib pool.
    stub("region-prep/prepare_masked_align_region_per_rg")
}

/// Step 3 stub: per-RG read-extraction → minimap2 realign → variant call,
/// parallel over RGs (Python `imap_process_masked_bam`,
/// `realign_recall/realign_per_RG.py:10-130`).
fn realign_per_rg(_paths: &Paths, _budget: ThreadBudget) -> Result<()> {
    // For each RG (rayon over RGs, _budget.num_jobs workers):
    //   TODO(T9): read-extraction — call the read_extraction crate
    //     (rust_read_extraction; BAM→FASTQ over the FC + NFC BEDs, validated),
    //     replacing realign_per_RG.py's bam_to_fastq_biobambam + seqkit rmdup/pair.
    //   TODO(T9): realign — minimap2 via minimap2-rs FFI against
    //     paths.masked_genome_path(rg) / paths.minimap_index_path(rg),
    //     replacing shell_utils.sh::independent_minimap2_masked, writing
    //     paths.rg_raw_masked_bam_path(rg).
    //   TODO(T9): variant call — bcftools call (FFI or isolated leaf
    //     `bcftools call` subprocess), replacing shell_utils.sh::bcftools_call_per_RG,
    //     writing the per-RG raw VCF (raw_masked_bam.replace(".bam", ".vcf.gz")).
    //   Forward _budget.threads_per_job to minimap2 + htslib.
    stub("realign_per_rg (read-extraction + minimap2 + variant-call)")
}

/// Step 4 stub: merge per-RG raw BAMs + dedup (Python `merge_bams` +
/// collate/fixmate/sort/markdup, `realign_and_recall.py:139-153`).
fn merge_and_markdup_raw_bams(_paths: &Paths, _threads: usize) -> Result<()> {
    // TODO(T9): merge per-RG raw BAMs (paths.rg_raw_masked_bam_path per RG) into
    //   paths.pooled_raw_bam_path() via the sdrecall-io BAM merge (replaces
    //   src/utils.py::merge_bams), then collate→fixmate→sort→markdup -r into the
    //   ".deduped.bam" sibling. markdup is reimplemented in rust-htslib or run as
    //   an isolated leaf samtools step (T9 policy decision; see T0 DESIGN risk).
    stub("merge_and_markdup_raw_bams")
}

/// Step 5 stub: concat per-RG raw VCFs (Python `bcftools_concatvcfs`,
/// `realign_and_recall.py:156-173`).
fn concat_raw_vcfs(_paths: &Paths, _threads: usize) -> Result<()> {
    // TODO(T9): call vcf_ops concat/sort (the sdrecall-io concat_sort_vcfs k-way
    //   merge that replaces combine_vcfs / bcftools_concatvcfs) over the per-RG
    //   raw VCFs, writing paths.recall_raw_vcf_path().
    stub("concat_raw_vcfs (vcf-ops concat)")
}

// ════════════════════════════════════════════════════════════════════════
//  Sub-phase: MISALIGNMENT ELIMINATION (Python eliminate_misalignments,
//  misalignment_elimination.py:49-249 — the fp-control per-island fan-out)
// ════════════════════════════════════════════════════════════════════════

/// Misalignment elimination: slice the pooled BAM into coverage islands, then
/// run fp-control per island in parallel (Python `eliminate_misalignments`,
/// `misalignment_elimination.py:49-249`).
fn eliminate_misalignments(
    paths: &Paths,
    threads: usize,
    numba_threads: usize,
    mq_cutoff: i32,
) -> Result<()> {
    log::info!("[fp-control] misalignment elimination for {}", paths.sample_id);

    // Python: filter_redundant_intrinsic_origins(total_intrinsic_bam)
    //   (realign_and_recall.py:183-184) — done once before slicing.
    // TODO(T9): call fp_control::filter_redundant_intrinsic_origins(...) on
    //   paths.total_intrinsic_bam_path() (or its sdrecall-io equivalent).
    filter_redundant_intrinsic_origins(paths)?;

    // Python: split_bam_by_cov(...) on the deduped raw BAM + on the intrinsic BAM,
    //   slop+merge the actual-coverage BEDs by 1000bp (misalignment_elimination.py:75-108).
    slice_islands(paths)?;

    // Python: optional NM Poisson cutoff (cal_edge_NM_values).
    // TODO(T9): per-island (or pooled) call nm_stats::nm_distribution_poisson(bam,
    //   conf_level, sample_size, threads) (rust_modules/nm-stats/src/lib.rs:122 —
    //   validated) to derive the edit-distance cutoff, replacing
    //   realign_recall/cal_edge_NM_values.py::calculate_NM_distribution_poisson.
    nm_cutoff(paths)?;

    // Python: configure_parallelism(threads, numba_threads) → job_num islands at
    //   once, then pool.imap_unordered(imap_filter_out, ...)
    //   (misalignment_elimination.py:131-158). In Rust this is rayon over islands.
    let island_budget = ThreadBudget::new(threads, numba_threads as f64);
    log::debug!("[fp-control] per-island budget {island_budget:?}");
    fp_control_per_island(paths, island_budget, mq_cutoff)?;

    // Python: merge_bams over per-island clean + raw BAMs, then bcftools_concatvcfs
    //   over per-island clean VCFs (misalignment_elimination.py:218-244).
    merge_island_outputs(paths, threads)?;

    Ok(())
}

/// Stub: filter redundant intrinsic origins once before slicing (Python
/// `filter_redundant_intrinsic_origins`).
fn filter_redundant_intrinsic_origins(_paths: &Paths) -> Result<()> {
    // TODO(T9): call fp_control::filter_redundant_intrinsic_origins(...) (the
    //   validated fp-control crate; identify_misaligned_haps origin filter).
    stub("filter_redundant_intrinsic_origins")
}

/// Stub: slice the pooled BAM into coverage islands (Python `split_bam_by_cov`,
/// the absorbed `realign_recall/slice_bam_by_cov.py`).
fn slice_islands(_paths: &Paths) -> Result<()> {
    // TODO(T9): absorb realign_recall/slice_bam_by_cov.py::split_bam_by_cov as an
    //   sdrecall-io coverage-island slicer (rust-htslib pileup + bedrs slop/merge,
    //   delimiter = ceil(avg_frag_size * 1.5)); produces per-island raw BAMs +
    //   padded coverage BEDs, and the matching per-island intrinsic BAM slices.
    stub("slice_islands (slice_bam_by_cov)")
}

/// Stub: NM (edit-distance) Poisson cutoff (Python `calculate_NM_distribution_poisson`).
fn nm_cutoff(_paths: &Paths) -> Result<()> {
    // TODO(T9): nm_stats::nm_distribution_poisson(bam, conf_level, sample_size,
    //   threads) (rust_modules/nm-stats/src/lib.rs:122 — validated).
    stub("nm_stats/nm_distribution_poisson")
}

/// Stub: the per-island fp-control fan-out — the headline parallel hot path
/// (Python `imap_filter_out` → `fp_control.realign_filter_per_cov`,
/// `misalignment_elimination.py:138-158`).
///
/// This is where rayon parallelises over coverage islands (the independent
/// axis), each island getting `budget.threads_per_job` inner threads. Per
/// ROB-3, each island task must be wrapped (`std::panic::catch_unwind` /
/// rayon propagation) so one island's M-CIGAR guard panic is isolated as a
/// failed island, not a whole-run abort.
fn fp_control_per_island(_paths: &Paths, _budget: ThreadBudget, _mq_cutoff: i32) -> Result<()> {
    // For each island (rayon par-iter, _budget.num_jobs at a time):
    //   TODO(T9): call fp_control::run_fp_control(bam, intrinsic_bam, &params)
    //     (rust_modules/fp-control/src/lib.rs:156 — validated; this fuses
    //     build_phasing_graph + phasing + haplotype_inspection in-process, the
    //     2.5-3.7× hot path). Then:
    //   TODO(T9): filter the island BAM by the returned correct/mismap qname sets
    //     (rust-htslib write, replacing realign_filter_per_cov.py's BAM rewrite),
    //     writing the per-island clean BAM.
    //   TODO(T9): variant call on the clean BAM (bcftools, as in step 3) →
    //     per-island clean VCF.
    //   Wrap the whole island body in catch_unwind so a panic → skipped island.
    stub("fp_control/run_fp_control (per-island, rayon)")
}

/// Stub: merge per-island clean + raw BAMs and concat per-island clean VCFs
/// (Python `merge_bams` ×2 + `bcftools_concatvcfs`,
/// `misalignment_elimination.py:218-244`).
fn merge_island_outputs(_paths: &Paths, _threads: usize) -> Result<()> {
    // TODO(T9): sdrecall-io BAM merge of per-island clean BAMs →
    //   paths.pooled_filtered_bam_path(), and of per-island raw BAMs back into
    //   the deduped raw BAM; vcf_ops concat of per-island clean VCFs →
    //   paths.recall_filtered_vcf_path().
    stub("merge_island_outputs")
}

/// Step 7 stub: priority-merge raw vs clean VCFs, then subset to the target
/// recall BED (Python `merge_with_priority` + `bcftools view -R`,
/// `realign_and_recall.py:213-231`).
fn merge_and_subset_final_vcf(_paths: &Paths, _threads: usize) -> Result<()> {
    // TODO(T9): vcf_ops::merge_with_priority(MergeParams{ query=recall_raw_vcf,
    //   reference=recall_filtered_vcf, added_filter="MISALIGNED", qv_tag="RAW",
    //   rv_tag="CLEAN", ... }) (rust_modules/vcf-ops/src/priority_merge.rs:110 —
    //   validated), then subset to paths.total_recall_sd_region_bed_path() via the
    //   sdrecall-io VCF region filter, writing paths.final_recall_vcf_path().
    stub("vcf-ops/merge_with_priority + subset")
}

// ════════════════════════════════════════════════════════════════════════
//  Phase: POST-PROCESS (Python post_process_vcf, SDrecall:41-114)
// ════════════════════════════════════════════════════════════════════════

/// Post-process the recall VCF: optional inhouse-common annotation against a
/// cohort VCF, then optional merge with a conventional caller's VCF (Python
/// `post_process_vcf`, `SDrecall:41-114`).
fn post_process_vcf(
    sdrecall_vcf: &std::path::Path,
    conventional: &crate::cli::ConventionalVcfArgs,
    cohort: &crate::cli::CohortArgs,
    _paths: &Paths,
) -> Result<std::path::PathBuf> {
    let mut final_vcf = sdrecall_vcf.to_path_buf();

    // Python: annotate_inhouse_common(...) if args.cohort_vcf (SDrecall:55-76).
    if let Some(cohort_vcf) = &cohort.cohort_vcf {
        log::info!("[post] annotating inhouse-common against {cohort_vcf:?}");
        // TODO(T9): vcf_ops::annotate_inhouse_common(InhouseParams{ query=final_vcf,
        //   cohort=cohort_vcf, cutoff=cohort.inhouse_common_cutoff,
        //   conf_level=cohort.cohort_conf_level, ... })
        //   (rust_modules/vcf-ops/src/inhouse_common.rs:90 — validated).
        final_vcf = annotate_inhouse_common(&final_vcf, cohort_vcf)?;
    }

    // Python: merge_with_priority(...) if args.conventional_vcf (SDrecall:79-104).
    if let Some(conventional_vcf) = &conventional.conventional_vcf {
        log::info!("[post] merging with conventional VCF {conventional_vcf:?}");
        // TODO(T9): vcf_ops::merge_with_priority(MergeParams{ query=final_vcf,
        //   reference=conventional_vcf, qv_tag="SDrecall", rv_tag=caller_name, ... })
        //   (rust_modules/vcf-ops/src/priority_merge.rs:110 — validated). Output is
        //   conventional.merged_vcf or derived from the conventional path.
        final_vcf = merge_with_conventional(&final_vcf, conventional)?;
    }

    // TODO(T9): shutil.rmtree(paths.tmp_dir) cleanup (SDrecall:106-112).
    Ok(final_vcf)
}

/// Stub: inhouse-common annotation (Python `annotate_inhouse_common`).
fn annotate_inhouse_common(
    query_vcf: &std::path::Path,
    _cohort_vcf: &std::path::Path,
) -> Result<std::path::PathBuf> {
    // TODO(T9): vcf_ops::annotate_inhouse_common (inhouse_common.rs:90).
    let _ = stub::<()>("vcf-ops/annotate_inhouse_common");
    Ok(query_vcf.to_path_buf())
}

/// Stub: merge with a conventional caller's VCF (Python `merge_with_priority`).
fn merge_with_conventional(
    query_vcf: &std::path::Path,
    _conventional: &crate::cli::ConventionalVcfArgs,
) -> Result<std::path::PathBuf> {
    // TODO(T9): vcf_ops::merge_with_priority (priority_merge.rs:110).
    let _ = stub::<()>("vcf-ops/merge_with_priority (conventional)");
    Ok(query_vcf.to_path_buf())
}

// ════════════════════════════════════════════════════════════════════════

/// Scaffold marker: a stage that is not yet wired. Logs which entry point the
/// integration pass must call and returns `Ok` so the call graph is walkable
/// end-to-end without doing any work (the scaffold pass deliberately wires no
/// stage crate).
fn stub<T: Default>(stage: &str) -> Result<T> {
    log::warn!("[stub] stage '{stage}' not yet wired (T9 integration pass)");
    Ok(T::default())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thread_budget_mirrors_python_rule() {
        // realign per-RG prep policy: configure_parallelism(threads, 1).
        let b = ThreadBudget::new(12, 1.0);
        assert_eq!(b.num_jobs, 12);
        assert_eq!(b.threads_per_job, 1);
        assert_eq!(b.total_threads, 12);

        // realign policy: configure_parallelism(threads, 3) → ceil(12/3)=4 jobs.
        let b = ThreadBudget::new(12, 3.0);
        assert_eq!(b.num_jobs, 4);
        assert_eq!(b.threads_per_job, 3);
        assert_eq!(b.num_jobs * b.threads_per_job, b.total_threads);

        // fp-control per-island policy with numba_threads=2 over 10 threads:
        // ceil(10/2)=5 islands at once, 2 threads each.
        let b = ThreadBudget::new(10, 2.0);
        assert_eq!((b.num_jobs, b.threads_per_job), (5, 2));

        // non-divisible total rounds the job count UP (ceil), like np.ceil.
        let b = ThreadBudget::new(10, 3.0);
        assert_eq!(b.num_jobs, 4); // ceil(10/3) = 4
        assert_eq!(b.threads_per_job, 3);
    }

    #[test]
    fn rgref_normalization_smoke() {
        // The pipeline addresses RGs through RgRef; confirm the conversion the
        // stage stubs rely on.
        use crate::paths::RgRef;
        assert_eq!(
            Paths::normalize_rg_label(RgRef::from(0u32)).unwrap(),
            "RG0"
        );
    }
}
