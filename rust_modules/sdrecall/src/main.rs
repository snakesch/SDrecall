//! `sdrecall` (SDrecall T9) — top-level orchestrator binary.
//!
//! Threads the stage crates in-process (no Python interpreter, no PyO3, no
//! subprocess CLI-to-CLI on the hot path) behind the `run` / `prepare` /
//! `realign` subcommands, mirroring the Python `SDrecall` executable.
//!
//! This is the **scaffold** of T9: the clap CLI, the `Paths` port, the
//! thread-budget plumbing and the orchestration call graph are in place, but the
//! stage crates are deliberately NOT wired yet (each pipeline stage is a stub
//! with a precise `// TODO(T9)` naming the validated entry point it will call).
//! See `docs/analysis/tasks/T9_orchestrator.md`.

// Scaffold pass: the full `Paths` surface (all getters + the input/frag-size
// fields) and the per-subcommand orchestration call graph are deliberately
// defined ahead of their wiring, so the T9 integration pass has the complete,
// tested API to call into. They are exercised by the unit tests but not yet by
// the (stubbed) pipeline body, so silence the unused-surface warnings until
// integration consumes them.
#![allow(dead_code)]

mod cli;
mod paths;
mod pipeline;

use std::path::Path;
use std::process::ExitCode;

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
            // Errors to stderr (CLI convention); non-zero exit for scripts.
            eprintln!("sdrecall error: {msg}");
            ExitCode::FAILURE
        }
    }
}

/// Validate the common args, derive [`Paths`], and dispatch to the per-subcommand
/// orchestration entry point. Returns the final VCF path for the phases that
/// produce one (`run`, `realign`); `prepare` returns `None`.
fn dispatch(cli: Cli) -> Result<Option<std::path::PathBuf>, String> {
    match cli.command {
        Command::Run(args) => {
            args.common.validate()?;
            let paths = derive_paths(&args.common)?;
            let vcf = pipeline::run_full_pipeline(&args, &paths).map_err(|e| e.to_string())?;
            Ok(Some(vcf))
        }
        Command::Prepare(args) => {
            args.common.validate()?;
            let paths = derive_paths(&args.common)?;
            pipeline::run_preparation_only(&args, &paths).map_err(|e| e.to_string())?;
            Ok(None)
        }
        Command::Realign(args) => {
            args.common.validate()?;
            let paths = derive_paths(&args.common)?;
            let vcf = pipeline::run_realign_only(&args, &paths).map_err(|e| e.to_string())?;
            Ok(Some(vcf))
        }
    }
}

/// Build the [`Paths`] layout from the common CLI args (pure derivation).
///
/// `repo_dir` is the directory containing the `data/` tree of intrinsic-align
/// BAMs. Python derives it from the source file location (`SDrecallPaths`
/// `os.path.dirname(os.path.dirname(__file__))`); here it comes from the
/// `SDRECALL_REPO_DIR` env var, falling back to the binary's parent dir.
/// TODO(T9-integration): resolve repo_dir the way the Python install does
/// (relative to the package root) once packaging is decided.
fn derive_paths(common: &cli::CommonArgs) -> Result<Paths, String> {
    let repo_dir = std::env::var("SDRECALL_REPO_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| {
            std::env::current_exe()
                .ok()
                .and_then(|p| p.parent().map(Path::to_path_buf))
                .unwrap_or_else(|| std::path::PathBuf::from("."))
        });

    Paths::derive(
        &common.ref_genome,
        &common.input_bam,
        &common.reference_sd_map,
        &common.outdir,
        Some(&common.target_bed),
        common.sample_id.as_deref(),
        Some(&common.target_tag),
        Some(&common.ref_genome_tag),
        &repo_dir,
    )
    .map_err(|e| e.to_string())
}
