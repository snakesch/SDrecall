//! `fp-control` CLI — files in → qname sets out.
//!
//! Runs the fused Phase-2c FP-control core
//! (`build_phasing_graph → phasing → haplotype_inspection`) on one island BAM
//! and writes the two qname columns to a TSV: `qname<TAB>label` where label is
//! `correct` or `mismap`. This is the independently-runnable face of the lib;
//! the in-process `sdrecall` orchestrator (T9) calls the library directly.

use clap::Parser;
use fp_control::{run_fp_control, FpControlParams, PairingEngine};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process::ExitCode;

#[derive(Parser)]
#[command(
    name = "fp-control",
    about = "Fused Phase-2c FP control: BAM -> graph -> phasing -> inspection (SDrecall T4)"
)]
struct Cli {
    /// Island BAM (the realigned per-chunk BAM).
    #[arg(long)]
    bam: String,

    /// Intrinsic (reference-sequence) BAM for the same chunk.
    #[arg(long)]
    intrinsic: String,

    /// Reference genome FASTA (with `.fai`); needed by the graph's allele-depth mpileup.
    #[arg(long)]
    reference: String,

    /// Region label (informational; the BAMs are already sliced per chunk). Reserved.
    #[arg(long, default_value = "")]
    region: String,

    /// Output TSV: `qname<TAB>{correct|mismap}`.
    #[arg(short = 'o', long)]
    output: PathBuf,

    /// Phasing edge-weight cutoff.
    #[arg(long, default_value_t = 0.301)]
    edge_weight_cutoff: f32,

    /// Mean read length.
    #[arg(long, default_value_t = 148.0)]
    mean_read_length: f32,

    /// MAPQ floor (recall_mq_cutoff).
    #[arg(long, default_value_t = 10)]
    mapq_cutoff: u8,

    /// Median base-quality floor.
    #[arg(long, default_value_t = 15)]
    basequal_median_cutoff: u8,

    /// htslib / collate threads.
    #[arg(long, default_value_t = 4)]
    threads: u8,

    /// Read-pair grouping engine.
    #[arg(long, value_enum, default_value_t = PairingEngine::SamtoolsPipe)]
    pairing_engine: PairingEngine,

    /// Optional haplotype-comparison meta TSV path (inspect's debug dump). Empty = skip.
    #[arg(long, default_value = "")]
    compare_meta: String,

    /// Log level: error | warn | info | debug | trace.
    #[arg(long, default_value = "info")]
    log_level: String,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let level = cli.log_level.parse().unwrap_or(log::LevelFilter::Info);
    sdrecall_utils::init_console_logger(level);

    let params = FpControlParams {
        reference_genome: cli.reference,
        edge_weight_cutoff: cli.edge_weight_cutoff,
        mean_read_length: cli.mean_read_length,
        mapq_cutoff: cli.mapq_cutoff,
        basequal_median_cutoff: cli.basequal_median_cutoff,
        threads: cli.threads,
        pairing_engine: cli.pairing_engine,
        resources: None,
        compare_haplotype_meta_tab: cli.compare_meta,
    };

    match run_fp_control(&cli.bam, &cli.intrinsic, &params) {
        Ok(Some(out)) => match write_output(&cli.output, &out) {
            Ok(()) => {
                log::info!(
                    "wrote {} correct + {} mismap qnames to {}",
                    out.correct_qnames.len(),
                    out.mismap_qnames.len(),
                    cli.output.display()
                );
                ExitCode::SUCCESS
            }
            Err(e) => {
                eprintln!("fp-control: writing {}: {e}", cli.output.display());
                ExitCode::FAILURE
            }
        },
        Ok(None) => {
            // Island skipped (≤2 vertices / empty matrix); emit an empty TSV so
            // downstream tooling sees a present-but-empty result, not a missing file.
            log::warn!("island {} skipped (no classification produced)", cli.bam);
            match write_output(&cli.output, &Default::default()) {
                Ok(()) => ExitCode::SUCCESS,
                Err(e) => {
                    eprintln!("fp-control: writing {}: {e}", cli.output.display());
                    ExitCode::FAILURE
                }
            }
        }
        Err(e) => {
            eprintln!("fp-control error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn write_output(path: &PathBuf, out: &fp_control::FpControlOutput) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "qname\tlabel")?;
    for q in &out.correct_qnames {
        writeln!(w, "{q}\tcorrect")?;
    }
    for q in &out.mismap_qnames {
        writeln!(w, "{q}\tmismap")?;
    }
    w.flush()
}
