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

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

use rayon::prelude::*;
use rust_htslib::{bam, bam::Read};
use sdrecall_utils::{clamp_threads_u8, Result, SdError};

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
        paths.sample_id,
        paths.assembly,
        paths.target_tag
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
        args.realign.strict_islands,
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
    log::info!(
        "[prepare] building recall regions into {:?}",
        paths.work_dir
    );
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
        args.realign.strict_islands,
    )
}

fn realign_and_recall_inner(
    paths: &Paths,
    threads: usize,
    numba_threads: usize,
    mq_cutoff: i32,
    strict_islands: bool,
) -> Result<PathBuf> {
    log::info!(
        "[realign] start for sample={} (threads={threads})",
        paths.sample_id
    );

    // ── Step 1: per-RG region size stats ────────────────────────────────
    let rg_infos = stat_realign_group_regions(paths)?;
    let rg_labels: Vec<&str> = rg_infos.iter().map(|r| r.label.as_str()).collect();
    log::info!(
        "[realign] {} RGs discovered: {:?}",
        rg_labels.len(),
        rg_labels
    );

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
    eliminate_misalignments(paths, threads, numba_threads, mq_cutoff, strict_islands)?;

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

    let target_bed = paths
        .target_bed
        .as_deref()
        .ok_or_else(|| SdError::Compute("target_bed is required for region-prep".into()))?;

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
                dedup_and_pair_fastqs(&r1, &r2)?;

                // 3b: minimap2 realign onto the masked genome (contigs `{chrom}:{start}`)
                // → an intermediate LOCAL-coordinate BAM.
                let local_bam = raw_bam.with_extension("local.bam");
                crate::tools::minimap2_align(&r1, &r2, &masked_genome, &local_bam, tpj)?;

                // 3b': remap masked(local) → ORIGINAL-genome coordinates (port of
                // shell_utils.sh independent_minimap2_masked's modify_bam_sq_lines +
                // modify_masked_genome_coords). This makes `raw_bam` genomic so the
                // later merge of all per-RG BAMs shares the original-ref header (Bug B)
                // and bcftools mpileup -f ref_genome is coordinate-consistent.
                sdrecall_io::remap_masked_bam_to_genomic(
                    &local_bam,
                    &ref_genome,
                    &rg.label,
                    &raw_bam,
                )?;
                let _ = std::fs::remove_file(&local_bam);
                let _ = std::fs::remove_file(format!("{}.bai", local_bam.display()));

                // 3c: variant call on the GENOMIC BAM.
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

#[derive(Clone, Debug)]
struct FastqRecord {
    header: String,
    seq: String,
    plus: String,
    qual: String,
}

fn dedup_and_pair_fastqs(r1: &Path, r2: &Path) -> Result<()> {
    let (r1_order, r1_records) = read_fastq_dedup_by_name(r1)?;
    let (_r2_order, r2_records) = read_fastq_dedup_by_name(r2)?;

    let tmp_r1 = fastq_tmp_path(r1, "rmdup.paired");
    let tmp_r2 = fastq_tmp_path(r2, "rmdup.paired");
    {
        let mut w1 =
            std::io::BufWriter::new(std::fs::File::create(&tmp_r1).map_err(|e| SdError::Io {
                path: tmp_r1.display().to_string(),
                source: e,
            })?);
        let mut w2 =
            std::io::BufWriter::new(std::fs::File::create(&tmp_r2).map_err(|e| SdError::Io {
                path: tmp_r2.display().to_string(),
                source: e,
            })?);
        for name in &r1_order {
            if let (Some(rec1), Some(rec2)) = (r1_records.get(name), r2_records.get(name)) {
                write_fastq_record(&mut w1, rec1, &tmp_r1)?;
                write_fastq_record(&mut w2, rec2, &tmp_r2)?;
            }
        }
        w1.flush().map_err(|e| SdError::Io {
            path: tmp_r1.display().to_string(),
            source: e,
        })?;
        w2.flush().map_err(|e| SdError::Io {
            path: tmp_r2.display().to_string(),
            source: e,
        })?;
    }

    std::fs::rename(&tmp_r1, r1).map_err(|e| SdError::Io {
        path: r1.display().to_string(),
        source: e,
    })?;
    std::fs::rename(&tmp_r2, r2).map_err(|e| SdError::Io {
        path: r2.display().to_string(),
        source: e,
    })?;
    Ok(())
}

fn read_fastq_dedup_by_name(path: &Path) -> Result<(Vec<String>, HashMap<String, FastqRecord>)> {
    let file = std::fs::File::open(path).map_err(|e| SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    let mut lines = std::io::BufReader::new(file).lines();
    let mut order = Vec::new();
    let mut records = HashMap::new();

    loop {
        let Some(header) = next_fastq_line(&mut lines, path)? else {
            break;
        };
        let seq = required_fastq_line(&mut lines, path, "sequence")?;
        let plus = required_fastq_line(&mut lines, path, "plus")?;
        let qual = required_fastq_line(&mut lines, path, "quality")?;
        let name = fastq_name(&header)?;
        if !records.contains_key(&name) {
            order.push(name.clone());
            records.insert(
                name,
                FastqRecord {
                    header,
                    seq,
                    plus,
                    qual,
                },
            );
        }
    }

    Ok((order, records))
}

fn next_fastq_line(
    lines: &mut impl Iterator<Item = std::io::Result<String>>,
    path: &Path,
) -> Result<Option<String>> {
    match lines.next() {
        Some(Ok(line)) => Ok(Some(line)),
        Some(Err(e)) => Err(SdError::Io {
            path: path.display().to_string(),
            source: e,
        }),
        None => Ok(None),
    }
}

fn required_fastq_line(
    lines: &mut impl Iterator<Item = std::io::Result<String>>,
    path: &Path,
    field: &str,
) -> Result<String> {
    next_fastq_line(lines, path)?.ok_or_else(|| {
        SdError::Compute(format!(
            "truncated FASTQ {} while reading {field}",
            path.display()
        ))
    })
}

fn fastq_name(header: &str) -> Result<String> {
    let Some(rest) = header.strip_prefix('@') else {
        return Err(SdError::Compute(format!(
            "invalid FASTQ header without '@': {header}"
        )));
    };
    Ok(rest.split_whitespace().next().unwrap_or(rest).to_string())
}

fn write_fastq_record(writer: &mut impl Write, rec: &FastqRecord, path: &Path) -> Result<()> {
    writeln!(writer, "{}", rec.header).map_err(|e| SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    writeln!(writer, "{}", rec.seq).map_err(|e| SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    writeln!(writer, "{}", rec.plus).map_err(|e| SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    writeln!(writer, "{}", rec.qual).map_err(|e| SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    Ok(())
}

fn fastq_tmp_path(path: &Path, suffix: &str) -> PathBuf {
    let name = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("reads.fastq");
    path.with_file_name(format!("{}.{}.{}.fastq", name, std::process::id(), suffix))
}

// ── Step 4: merge + markdup ─────────────────────────────────────────────

fn merge_and_markdup_raw_bams(
    paths: &Paths,
    per_rg_bams: &[PathBuf],
    threads: usize,
) -> Result<()> {
    log::info!(
        "[merge] merging {} per-RG BAMs + markdup",
        per_rg_bams.len()
    );
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
    sdrecall_io::concat_sort_vcfs(
        &refs,
        &paths.recall_raw_vcf_path(),
        true,
        clamp_threads_u8(threads),
    )?;
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
    strict_islands: bool,
) -> Result<()> {
    log::info!(
        "[fp-control] misalignment elimination for {}",
        paths.sample_id
    );

    // Python feeds fp-control with the multi-align BED, not the user target BED.
    let target_bed = paths.multi_align_bed_path();
    let total_intrinsic_bam =
        filter_redundant_intrinsic_origins(&paths.total_intrinsic_bam_path(), threads)?;
    let filtered_intrin = paths.tmp_dir.join("intrinsic.filtered.bam");
    crate::tools::samtools_view_region(
        &total_intrinsic_bam,
        &target_bed,
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
        &target_bed,
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
        std::fs::copy(
            paths.recall_raw_vcf_path(),
            paths.recall_filtered_vcf_path(),
        )
        .map_err(|e| SdError::Io {
            path: paths.recall_filtered_vcf_path().display().to_string(),
            source: e,
        })?;
        return Ok(());
    }

    // NM Poisson cutoff: NOT YET WIRED (T9). The Python pipeline applies an NM
    // cutoff during per-island inspection; the Rust inspect path does not consume
    // one yet. Rather than run a 10k-read scan whose result is discarded —
    // pretending to apply a filter we don't — the scan is omitted until `nm_cutoff`
    // is threaded into `fp_control::FpControlParams` and the per-island filter.
    // `nm_stats::nm_distribution_poisson` stays available and unit-tested for that
    // wiring.

    // Per-island fp-control (rayon parallel).
    let island_budget = ThreadBudget::new(threads, numba_threads as f64);
    let (clean_bams, clean_vcfs) =
        fp_control_per_island(paths, &islands, island_budget, mq_cutoff, strict_islands)?;

    // Merge per-island outputs.
    merge_island_outputs(paths, &clean_bams, &clean_vcfs, threads)?;

    Ok(())
}

#[derive(Clone, Debug)]
struct IntrinsicOrigin {
    chrom: String,
    start: i64,
    end: i64,
}

fn filter_redundant_intrinsic_origins(intrinsic_bam: &Path, threads: usize) -> Result<PathBuf> {
    let mut reader = bam::Reader::from_path(intrinsic_bam).map_err(hts_err)?;
    if threads > 1 {
        reader.set_threads(threads - 1).map_err(hts_err)?;
    }

    let mut origins: HashMap<String, IntrinsicOrigin> = HashMap::new();
    for result in reader.records() {
        let record = result.map_err(hts_err)?;
        let qname = String::from_utf8_lossy(record.qname()).to_string();
        if origins.contains_key(&qname) {
            continue;
        }
        if let Some(origin) = parse_intrinsic_origin(&qname) {
            origins.insert(qname, origin);
        }
    }

    if origins.is_empty() {
        log::warn!(
            "[fp-control] no intrinsic origin intervals parsed from {}",
            intrinsic_bam.display()
        );
        return Ok(intrinsic_bam.to_path_buf());
    }

    let removable = redundant_origin_qnames(&origins);
    if removable.is_empty() {
        log::info!(
            "[fp-control] intrinsic origin filtering: {} unique origins, none redundant",
            origins.len()
        );
        return Ok(intrinsic_bam.to_path_buf());
    }

    let filtered_bam = intrinsic_bam.with_extension("filtered.bam");
    let unsorted_bam = filtered_bam.with_extension("unsorted.bam");

    let mut reader = bam::Reader::from_path(intrinsic_bam).map_err(hts_err)?;
    if threads > 1 {
        reader.set_threads(threads - 1).map_err(hts_err)?;
    }
    let header = bam::Header::from_template(reader.header());
    let mut writer =
        bam::Writer::from_path(&unsorted_bam, &header, bam::Format::Bam).map_err(hts_err)?;
    if threads > 1 {
        writer.set_threads(threads - 1).map_err(hts_err)?;
    }

    let mut kept = 0usize;
    let mut removed = 0usize;
    for result in reader.records() {
        let record = result.map_err(hts_err)?;
        let qname = String::from_utf8_lossy(record.qname());
        if removable.contains(qname.as_ref()) {
            removed += 1;
        } else {
            writer.write(&record).map_err(hts_err)?;
            kept += 1;
        }
    }
    drop(writer);

    crate::tools::samtools_sort_index(&unsorted_bam, &filtered_bam, threads)?;
    let _ = std::fs::remove_file(&unsorted_bam);
    log::info!(
        "[fp-control] intrinsic origin filtering: {} unique origins -> {} kept, {} removed; alignments kept={}, removed={}",
        origins.len(),
        origins.len() - removable.len(),
        removable.len(),
        kept,
        removed
    );

    Ok(filtered_bam)
}

fn parse_intrinsic_origin(qname: &str) -> Option<IntrinsicOrigin> {
    let (chrom, rest) = qname.split_once(':')?;
    if chrom.is_empty()
        || !chrom
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'_')
    {
        return None;
    }
    let (start, end) = rest.split_once('-')?;
    let start = start.parse::<i64>().ok()?;
    let end = end.parse::<i64>().ok()?;
    Some(IntrinsicOrigin {
        chrom: chrom.to_string(),
        start,
        end,
    })
}

fn redundant_origin_qnames(origins: &HashMap<String, IntrinsicOrigin>) -> HashSet<String> {
    let mut by_chrom: HashMap<&str, Vec<(i64, i64, &str)>> = HashMap::new();
    for (qname, origin) in origins {
        by_chrom
            .entry(&origin.chrom)
            .or_default()
            .push((origin.start, origin.end, qname));
    }

    let mut removable = HashSet::new();
    for intervals in by_chrom.values_mut() {
        intervals.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| b.1.cmp(&a.1)));
        let mut max_end = -1i64;
        let mut max_qname: Option<&str> = None;
        for &(start, end, qname) in intervals.iter() {
            if end <= start {
                continue;
            }
            if let Some(container) = max_qname {
                if end <= max_end && container != qname {
                    removable.insert(qname.to_string());
                }
            }
            if end > max_end {
                max_end = end;
                max_qname = Some(qname);
            }
        }
    }
    removable
}

/// One island's terminal outcome in the rayon fan-out.
enum IslandOutcome {
    Done(PathBuf, PathBuf),
    /// `≤2` haplotypes — a legitimate non-result (not a failure), stays non-fatal.
    Skipped,
    /// Errored or panicked — recorded so the policy below can surface it.
    Failed {
        id: String,
        reason: String,
    },
}

fn fp_control_per_island(
    paths: &Paths,
    islands: &[IslandPaths],
    budget: ThreadBudget,
    mq_cutoff: i32,
    strict_islands: bool,
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

    let outcomes: Vec<IslandOutcome> = pool.install(|| {
        islands
            .par_iter()
            .map(|island| {
                // catch_unwind per ROB-3: a panic in one island must NOT abort the
                // whole batch. The outcome is recorded and the failure policy is
                // applied after the join.
                let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    process_one_island(island, &ref_genome_str, mq_cutoff as u8, tpj)
                }));

                match result {
                    Ok(Ok(Some((bam, vcf)))) => IslandOutcome::Done(bam, vcf),
                    Ok(Ok(None)) => {
                        log::debug!("[fp-control] island {} skipped (≤2 haplotypes)", island.id);
                        IslandOutcome::Skipped
                    }
                    Ok(Err(e)) => {
                        log::error!("[fp-control] island {} FAILED: {e}", island.id);
                        IslandOutcome::Failed {
                            id: island.id.to_string(),
                            reason: format!("error: {e}"),
                        }
                    }
                    Err(panic) => {
                        let msg = panic
                            .downcast_ref::<&str>()
                            .map(|s| s.to_string())
                            .or_else(|| panic.downcast_ref::<String>().cloned())
                            .unwrap_or_else(|| "unknown panic".to_string());
                        log::error!("[fp-control] island {} PANICKED: {msg}", island.id);
                        IslandOutcome::Failed {
                            id: island.id.to_string(),
                            reason: format!("panic: {msg}"),
                        }
                    }
                }
            })
            .collect()
    });

    let mut bams = Vec::new();
    let mut vcfs = Vec::new();
    let mut failures: Vec<(String, String)> = Vec::new();
    let mut skipped = 0usize;
    for outcome in outcomes {
        match outcome {
            IslandOutcome::Done(bam, vcf) => {
                bams.push(bam);
                vcfs.push(vcf);
            }
            IslandOutcome::Skipped => skipped += 1,
            IslandOutcome::Failed { id, reason } => failures.push((id, reason)),
        }
    }

    log::info!(
        "[fp-control] islands: {} succeeded, {} skipped (≤2 haps), {} failed",
        bams.len(),
        skipped,
        failures.len()
    );

    // Failure policy: tolerate by default (Python parity — failures yield a NaN
    // sentinel and the run continues) but never silently. Always ERROR-log (above)
    // + write a manifest; only abort when the operator opts into `--strict_islands`.
    if !failures.is_empty() {
        let manifest = write_failed_island_manifest(paths, &failures)?;
        log::warn!(
            "[fp-control] {} island(s) failed; their variants are NOT in the output. Manifest: {}",
            failures.len(),
            manifest.display()
        );
        if strict_islands {
            return Err(SdError::Compute(format!(
                "{} island(s) failed and --strict_islands is set; see {}",
                failures.len(),
                manifest.display()
            )));
        }
    }

    Ok((bams, vcfs))
}

/// Write a TSV of failed islands (`island_id\treason`) into the recall-results
/// dir so a tolerated (non-strict) failure is auditable, never silent. Returns
/// the manifest path.
fn write_failed_island_manifest(paths: &Paths, failures: &[(String, String)]) -> Result<PathBuf> {
    let manifest = paths
        .recall_results_dir
        .join(format!("{}.failed_islands.tsv", paths.sample_id));
    let mut body = String::from("island_id\treason\n");
    for (id, reason) in failures {
        body.push_str(id);
        body.push('\t');
        body.push_str(reason);
        body.push('\n');
    }
    std::fs::write(&manifest, body).map_err(|e| SdError::Io {
        path: manifest.display().to_string(),
        source: e,
    })?;
    Ok(manifest)
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
        basequal_median_cutoff: 15,
        threads: clamp_threads_u8(threads),
        ..Default::default()
    };

    let output = fp_control::run_fp_control(&bam_str, &intrin_str, &params)?;

    let Some(output) = output else {
        return Ok(None);
    };

    let correct_set: HashSet<String> = output.correct_qnames.iter().cloned().collect();
    let mismap_set: HashSet<String> = output.mismap_qnames.iter().cloned().collect();
    let lowqual_set: HashSet<String> = output.lowqual_qnames.iter().cloned().collect();
    let qname_hap = output.qname_hap.clone();

    // Replace the island raw BAM with HP-tagged primary alignments for parity
    // with Python's visualization path, then derive the clean BAM from it.
    crate::bam_filter::annotate_hp_tags(
        &island.raw_bam,
        &correct_set,
        &mismap_set,
        &lowqual_set,
        &qname_hap,
        &island.raw_bam,
        island.id,
        threads,
    )?;

    // Filter BAM to keep only correct, primary, non-duplicate, non-QC-fail reads.
    let clean_bam = island.raw_bam.with_extension("clean.bam");
    let clean_alignments = crate::bam_filter::filter_bam_by_qnames(
        &island.raw_bam,
        &correct_set,
        &lowqual_set,
        &qname_hap,
        &clean_bam,
        island.id,
        threads,
    )?;
    if clean_alignments == 0 {
        log::warn!(
            "[fp-control] island {} produced an empty clean BAM after Python-policy filtering",
            island.id
        );
        return Ok(None);
    }

    // Variant-call on the clean BAM, then annotate clean variants with HP support.
    let clean_vcf = clean_bam.with_extension("vcf.gz");
    let unannotated_vcf = clean_bam.with_extension("unannotated.vcf.gz");
    crate::tools::bcftools_call(&clean_bam, Path::new(ref_genome), &unannotated_vcf, threads)?;
    crate::vcf_hp::annotate_vcf_hp(&unannotated_vcf, &clean_vcf, &clean_bam, threads)?;
    let _ = std::fs::remove_file(&unannotated_vcf);
    let _ = std::fs::remove_file(format!("{}.csi", unannotated_vcf.display()));
    let _ = std::fs::remove_file(format!("{}.tbi", unannotated_vcf.display()));

    log::info!(
        "[fp-control] island {} done: {} correct, {} mismap",
        island.id,
        correct_set.len(),
        mismap_set.len()
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
    sdrecall_io::concat_sort_vcfs(
        &vcf_refs,
        &paths.recall_filtered_vcf_path(),
        true,
        clamp_threads_u8(threads),
    )?;

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
        threads: clamp_threads_u8(threads),
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

    log::info!(
        "[merge-final] final VCF → {:?}",
        paths.final_recall_vcf_path()
    );
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
        let merged = conventional.merged_vcf.clone().unwrap_or_else(|| {
            let stem = conventional_vcf
                .file_stem()
                .unwrap_or_default()
                .to_string_lossy();
            conventional_vcf.with_file_name(format!("{stem}.sdrecall_merged.vcf.gz"))
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

fn hts_err(e: impl std::fmt::Display) -> SdError {
    SdError::Htslib(e.to_string())
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
