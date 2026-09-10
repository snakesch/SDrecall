//! `sdrecall` (SDrecall T9) — top-level orchestrator binary.
//!
//! Threads the stage crates in-process (no Python interpreter, no PyO3, no
//! subprocess CLI-to-CLI on the hot path) behind the `run` / `prepare` /
//! `realign` subcommands, mirroring the Python `SDrecall` executable.

mod bam_filter;
mod cli;
mod island;
mod paths;
mod pipeline;
mod rg_discovery;
mod tools;
mod vcf_hp;

use std::process::ExitCode;
use std::time::Instant;

use clap::Parser;

use cli::{Cli, Command};
use paths::Paths;

fn main() -> ExitCode {
    let cli = Cli::parse();
    sdrecall_utils::init_console_logger(cli.verbose.to_level_filter());

    match dispatch(cli) {
        Ok(Some(final_vcf)) => {
            log::info!("Final output: {}", final_vcf.display());
            ExitCode::SUCCESS
        }
        Ok(None) => ExitCode::SUCCESS,
        Err(msg) => {
            eprintln!("sdrecall error: {msg}");
            ExitCode::FAILURE
        }
    }
}

fn dispatch(cli: Cli) -> Result<Option<std::path::PathBuf>, String> {
    match cli.command {
        Command::Run(args) => {
            let command_start = Instant::now();
            args.common.validate()?;
            let mut paths = derive_paths(&args.common)?;
            compute_frag_stats(&mut paths);
            let setup_time = command_start.elapsed();
            let vcf = pipeline::run_full_pipeline(&args, &paths).map_err(|e| e.to_string())?;
            log::warn!(
                "[command_stage_metrics] command=run t_setup_s={:.3} t_total_s={:.3}",
                setup_time.as_secs_f64(),
                command_start.elapsed().as_secs_f64()
            );
            Ok(Some(vcf))
        }
        Command::Prepare(args) => {
            args.common.validate()?;
            let mut paths = derive_paths(&args.common)?;
            compute_frag_stats(&mut paths);
            pipeline::run_preparation_only(&args, &paths).map_err(|e| e.to_string())?;
            Ok(None)
        }
        Command::Realign(args) => {
            args.common.validate()?;
            let mut paths = derive_paths(&args.common)?;
            compute_frag_stats(&mut paths);
            let vcf = pipeline::run_realign_only(&args, &paths).map_err(|e| e.to_string())?;
            Ok(Some(vcf))
        }
    }
}

fn derive_paths(common: &cli::CommonArgs) -> Result<Paths, String> {
    Paths::derive(
        &common.ref_genome,
        &common.input_bam,
        &common.reference_sd_map,
        &common.outdir,
        Some(&common.target_bed),
        common.sample_id.as_deref(),
        Some(&common.target_tag),
        Some(&common.ref_genome_tag),
    )
    .map_err(|e| e.to_string())
}

/// Compute fragment-size stats from the input BAM and set them on `Paths`.
fn compute_frag_stats(paths: &mut Paths) {
    match sdrecall_io::get_insert_size_distribution(&paths.input_bam) {
        Ok(Some(stats)) => {
            log::info!(
                "Fragment size: mean={:.1} median={:.1} std={:.1}",
                stats.mean,
                stats.median,
                stats.std
            );
            paths.avg_frag_size = Some(stats.mean);
            paths.median_frag_size = Some(stats.median);
            paths.frag_size_std = Some(stats.std);
        }
        Ok(None) => {
            log::warn!("Could not derive fragment-size distribution; using defaults");
        }
        Err(e) => {
            log::warn!("Fragment-size estimation failed: {e}; using defaults");
        }
    }
}
