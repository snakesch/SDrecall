//! `nm-stats` CLI — file in → scalars out.
//!
//! Computes the NM (edit-distance) Poisson cutoff for a BAM and prints
//! `cutoff\tmean` to stdout, the same shape the Python differential dump emits.

use clap::Parser;
use nm_stats::nm_distribution_poisson;
use std::path::PathBuf;
use std::process::ExitCode;

#[derive(Parser)]
#[command(
    name = "nm-stats",
    about = "Compute the NM (edit-distance) Poisson cutoff for a BAM (SDrecall T6)"
)]
struct Cli {
    /// Input BAM (coordinate-sorted).
    #[arg(long)]
    bam: PathBuf,

    /// Confidence level (tail probability); the cutoff reaches the 1-conf_level percentile.
    #[arg(long, default_value_t = 0.01)]
    conf_level: f64,

    /// Maximum number of usable reads to sample.
    #[arg(long, default_value_t = 3_000_000)]
    sample_size: usize,

    /// htslib BGZF decompression threads.
    #[arg(long, default_value_t = 4)]
    threads: u8,

    /// Log level: error | warn | info | debug | trace.
    #[arg(long, default_value = "info")]
    log_level: String,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let level = cli
        .log_level
        .parse()
        .unwrap_or(log::LevelFilter::Info);
    sdrecall_utils::init_console_logger(level);

    match nm_distribution_poisson(&cli.bam, cli.conf_level, cli.sample_size, cli.threads) {
        Ok(nm) => {
            // Match the Python dump format: cutoff<TAB>mean.
            println!("{}\t{}", nm.cutoff, nm.mean);
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("nm-stats error: {e}");
            ExitCode::FAILURE
        }
    }
}
