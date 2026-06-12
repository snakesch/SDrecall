//! `region-prep` CLI — files in → fc/nfc BEDs out for one RG.
//!
//! Reads the 7-col `{rg}_related_homo_regions.bed`, the sample target BED and the
//! reference (for `<ref>.fasta.fai` contig sizes), then writes the shared FC
//! target BED and one nfc BED per subgroup, printing the 6-field CSV records the
//! Python differential dump emits (`rg,sub,fc_bed,nfc_bed,fc_size,nfc_size`).

use clap::Parser;
use region_prep::prepare_masked_align_region_per_rg;
use std::path::PathBuf;
use std::process::ExitCode;

#[derive(Parser)]
#[command(
    name = "region-prep",
    about = "Per-RG fc/nfc realignment-region projection (SDrecall T7)"
)]
struct Cli {
    /// RG label (re-derived from --whole-region-bed basename if they disagree, R5).
    #[arg(long)]
    rg_label: String,

    /// Comma-separated subgroup ids to process (e.g. `0,1,2`). If omitted, every
    /// subgroup with an `FC:{rg}_{sub}` tag in the whole-region BED is processed.
    #[arg(long, value_delimiter = ',')]
    subgroups: Option<Vec<String>>,

    /// The 7-col `{rg}_related_homo_regions.bed` (whole_region_bed).
    #[arg(long)]
    whole_region_bed: PathBuf,

    /// The sample target BED (first 3 columns used).
    #[arg(long)]
    target_bed: PathBuf,

    /// Reference genome; its `<ref>.fasta.fai` supplies contig sizes for slop (R7).
    #[arg(long)]
    ref_genome: PathBuf,

    /// Output path for the shared per-RG FC target BED.
    #[arg(long)]
    fc_target_out: PathBuf,

    /// Output directory for the per-subgroup nfc BEDs (`{rg}_{sub}.nfc.bed`).
    #[arg(long)]
    nfc_out_dir: PathBuf,

    /// Log level: error | warn | info | debug | trace.
    #[arg(long, default_value = "info")]
    log_level: String,
}

/// Resolve the subgroup id list: explicit `--subgroups`, else every `FC:` tag.
fn resolve_subgroups(cli: &Cli) -> sdrecall_utils::Result<Vec<String>> {
    if let Some(subs) = &cli.subgroups {
        return Ok(subs.clone());
    }
    let rows = region_prep::read_all_region_bed(&cli.whole_region_bed)?;
    let mut subs: Vec<String> = region_prep::fc_rows(&rows)
        .iter()
        .filter_map(|r| match &r.tag {
            region_prep::RgTag::Fc { sub, .. } => Some(sub.clone()),
            _ => None,
        })
        .collect();
    subs.sort();
    subs.dedup();
    Ok(subs)
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let level = cli.log_level.parse().unwrap_or(log::LevelFilter::Info);
    sdrecall_utils::init_console_logger(level);

    let subgroups = match resolve_subgroups(&cli) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("region-prep error: {e}");
            return ExitCode::FAILURE;
        }
    };

    match prepare_masked_align_region_per_rg(
        &cli.rg_label,
        &subgroups,
        &cli.target_bed,
        &cli.whole_region_bed,
        &cli.ref_genome,
        &cli.fc_target_out,
        &cli.nfc_out_dir,
    ) {
        Ok(records) => {
            for rec in &records {
                println!("{}", rec.to_csv());
            }
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("region-prep error: {e}");
            ExitCode::FAILURE
        }
    }
}
