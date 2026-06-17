//! Standalone BAM phaser CLI.
//!
//! Given a BAM file, builds a phasing graph from read-pair overlaps, clusters
//! reads into haplotypes, and writes an HP-tagged output BAM.
//!
//! Usage:
//!   phasing --bam <INPUT.BAM> --reference <REF.FA> --output <OUTPUT.BAM>

use std::time::Instant;

use anyhow::{Context, Result};
use clap::Parser;

use phasing::{phase_bam_with_intrinsic, PhaserParams};

#[derive(Parser, Debug)]
#[command(about = "Standalone BAM phaser for paired short-read data")]
struct Args {
    /// Input BAM file (paired-end short reads).
    #[arg(long)]
    bam: String,

    /// Reference genome FASTA (with .fai index).
    #[arg(long)]
    reference: String,

    /// Output HP-tagged BAM path.
    #[arg(long)]
    output: String,

    /// Optional intrinsic BAM (paralog reference alignments). When supplied, its
    /// allele-depth pileup drives paralogous-sequence-variant (PSV) detection in
    /// the edge-weight formula (Python `intrinsic_ad_dict`); omit it to disable PSV
    /// deductions (the standalone-phaser default).
    #[arg(long)]
    intrinsic_bam: Option<String>,

    /// Edge-weight cutoff for phasing clique rounds.
    #[arg(long, default_value_t = 0.301)]
    edge_weight_cutoff: f32,

    /// Mean read length (drives haplotype scoring).
    #[arg(long, default_value_t = 148.0)]
    mean_read_length: f32,

    /// Minimum MAPQ for graph + inspection.
    #[arg(long, default_value_t = 10)]
    mapq_cutoff: u8,

    /// Minimum median base quality.
    #[arg(long, default_value_t = 15)]
    basequal_median_cutoff: u8,

    /// htslib / samtools thread budget.
    #[arg(long, default_value_t = 4)]
    threads: u8,
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Args::parse();

    let params = PhaserParams {
        edge_weight_cutoff: args.edge_weight_cutoff,
        mean_read_length: args.mean_read_length,
        mapq_cutoff: args.mapq_cutoff,
        basequal_median_cutoff: args.basequal_median_cutoff,
        threads: args.threads,
    };

    let t0 = Instant::now();
    let out = phase_bam_with_intrinsic(
        &args.bam,
        &args.reference,
        &args.output,
        args.intrinsic_bam.as_deref(),
        &params,
    )
    .context("phase_bam failed")?;

    let elapsed = t0.elapsed();
    println!(
        "Done in {:.1}s: {} reads, {} haplotypes → {}",
        elapsed.as_secs_f64(),
        out.n_reads,
        out.n_haplotypes,
        out.output_bam.display()
    );

    Ok(())
}
