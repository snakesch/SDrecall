//! `vcf-ops` CLI — files in → files out.
//!
//! `merge`           — priority merge of a query VCF against a reference VCF
//!                     (ports `merge_variants_with_priority.py`).
//! `inhouse-common`  — binomial inhouse-common annotation against a cohort VCF
//!                     (ports `identify_common_vars.py`).

use clap::{Parser, Subcommand};
use std::path::PathBuf;
use std::process::ExitCode;
use vcf_ops::{
    annotate_inhouse_common, merge_with_priority, InhouseParams, MergeParams,
};

#[derive(Parser)]
#[command(
    name = "vcf-ops",
    about = "Sorted-VCF priority merge + inhouse-common binomial annotation (SDrecall T5)"
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,

    /// Log level: error | warn | info | debug | trace.
    #[arg(long, global = true, default_value = "info")]
    log_level: String,
}

#[derive(Subcommand)]
enum Cmd {
    /// Priority-merge a query VCF against a reference VCF.
    Merge(MergeArgs),
    /// Annotate inhouse-common variants against a cohort VCF.
    InhouseCommon(InhouseArgs),
}

#[derive(Parser)]
struct MergeArgs {
    #[arg(long)]
    query_vcf: PathBuf,
    #[arg(long)]
    reference_vcf: PathBuf,
    #[arg(long)]
    output_vcf: PathBuf,
    #[arg(long)]
    ref_genome: PathBuf,
    /// Filter tag added to query-only records (Python "MISALIGNED").
    #[arg(long)]
    added_filter: Option<String>,
    /// Source tag for query records (Python "RAW").
    #[arg(long, default_value = "RAW")]
    qv_tag: String,
    /// Source tag for reference records (Python "CLEAN").
    #[arg(long, default_value = "CLEAN")]
    rv_tag: String,
    /// Apply GT correction from AD/GQ/HPSUP (Python production call passes false).
    #[arg(long, default_value_t = false)]
    modify_gt: bool,
    #[arg(long, default_value_t = 4)]
    threads: u8,
    #[arg(long, default_value = "/tmp")]
    tmp_dir: PathBuf,
}

#[derive(Parser)]
struct InhouseArgs {
    #[arg(long)]
    query_vcf: PathBuf,
    #[arg(long)]
    cohort_vcf: PathBuf,
    #[arg(long)]
    output_vcf: PathBuf,
    #[arg(long)]
    ref_genome: PathBuf,
    #[arg(long, default_value = "INHOUSE_COMMON")]
    added_filter: String,
    #[arg(long, default_value_t = 0.01)]
    inhouse_common_cutoff: f64,
    #[arg(long, default_value_t = 0.999)]
    conf_level: f64,
    #[arg(long, default_value_t = 4)]
    threads: u8,
    #[arg(long, default_value = "/tmp")]
    tmp_dir: PathBuf,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let level = cli.log_level.parse().unwrap_or(log::LevelFilter::Info);
    sdrecall_utils::init_console_logger(level);

    let result = match &cli.cmd {
        Cmd::Merge(a) => merge_with_priority(MergeParams {
            query_vcf: &a.query_vcf,
            reference_vcf: &a.reference_vcf,
            output_vcf: &a.output_vcf,
            ref_genome: &a.ref_genome,
            added_filter: a.added_filter.as_deref(),
            qv_tag: Some(a.qv_tag.as_str()),
            rv_tag: Some(a.rv_tag.as_str()),
            modify_gt: a.modify_gt,
            threads: a.threads,
            tmp_dir: &a.tmp_dir,
        }),
        Cmd::InhouseCommon(a) => annotate_inhouse_common(InhouseParams {
            query_vcf: &a.query_vcf,
            cohort_vcf: &a.cohort_vcf,
            output_vcf: &a.output_vcf,
            ref_genome: &a.ref_genome,
            added_filter: &a.added_filter,
            inhouse_common_cutoff: a.inhouse_common_cutoff,
            conf_level: a.conf_level,
            threads: a.threads,
            tmp_dir: &a.tmp_dir,
        }),
    };

    match result {
        Ok(()) => {
            let out = match &cli.cmd {
                Cmd::Merge(a) => &a.output_vcf,
                Cmd::InhouseCommon(a) => &a.output_vcf,
            };
            println!("{}", out.display());
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("vcf-ops error: {e}");
            ExitCode::FAILURE
        }
    }
}
