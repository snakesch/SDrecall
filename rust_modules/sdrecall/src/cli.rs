//! clap CLI — mirrors the Python `SDrecall` argparse (the `SDrecall` executable,
//! `main()` + `_add_*_args` at `SDrecall:302-423`).
//!
//! Three subcommands: `run` (full pipeline), `prepare` (preparation only),
//! `realign` (realignment + recall only). Argument names, short flags, defaults
//! and the per-subcommand argument groups match the Python 1:1 so existing
//! invocations keep working against the Rust binary.

use std::path::PathBuf;

use clap::{Parser, Subcommand, ValueEnum};

/// SDrecall — complement SNV and small-indel detection within Segmental
/// Duplications from NGS data.
#[derive(Parser, Debug)]
#[command(
    name = "sdrecall",
    about = "Complement SNV and small Indel detection within Segmental Duplication based on NGS data",
    version
)]
pub struct Cli {
    /// Logging verbosity level (Python `-v/--verbose`, default INFO).
    #[arg(short, long, default_value = "INFO", value_enum, global = true)]
    pub verbose: Verbosity,

    #[command(subcommand)]
    pub command: Command,
}

/// Logging verbosity choices — matches the Python `choices=[...]` set
/// (uppercase: DEBUG/INFO/WARNING/ERROR/CRITICAL). The lowercase forms are
/// accepted too for ergonomics.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum Verbosity {
    #[value(name = "DEBUG", alias = "debug")]
    Debug,
    #[value(name = "INFO", alias = "info")]
    Info,
    #[value(name = "WARNING", alias = "warning")]
    Warning,
    #[value(name = "ERROR", alias = "error")]
    Error,
    #[value(name = "CRITICAL", alias = "critical")]
    Critical,
}

impl Verbosity {
    /// Map to the `log` crate's level filter.
    pub fn to_level_filter(self) -> log::LevelFilter {
        match self {
            Verbosity::Debug => log::LevelFilter::Debug,
            Verbosity::Info => log::LevelFilter::Info,
            Verbosity::Warning => log::LevelFilter::Warn,
            Verbosity::Error => log::LevelFilter::Error,
            Verbosity::Critical => log::LevelFilter::Error,
        }
    }
}

/// The three operation modes (Python subparsers `run` / `prepare` / `realign`).
#[derive(Subcommand, Debug)]
pub enum Command {
    /// Run the complete SDrecall pipeline (preparation → realignment + recall).
    Run(RunArgs),
    /// Run only the preparation phase.
    Prepare(PrepareArgs),
    /// Run only the realignment and recall phase.
    Realign(RealignArgs),
}

/// Arguments shared by every subcommand (Python `_add_common_args`,
/// `SDrecall:366-387`).
#[derive(clap::Args, Debug, Clone)]
pub struct CommonArgs {
    /// Reference genome FASTA (hg19 or hg38). File-name suffix must be `.fasta`.
    #[arg(short = 'r', long = "ref_genome")]
    pub ref_genome: PathBuf,

    /// Base directory for output files.
    #[arg(short = 'o', long = "outdir")]
    pub outdir: PathBuf,

    /// Input BAM file.
    #[arg(short = 'i', long = "input_bam")]
    pub input_bam: PathBuf,

    /// Reference SD-map file.
    #[arg(short = 'm', long = "reference_sd_map")]
    pub reference_sd_map: PathBuf,

    /// Target BED file.
    #[arg(short = 'b', long = "target_bed")]
    pub target_bed: PathBuf,

    /// Sample ID (extracted from the BAM filename if not provided).
    #[arg(short = 's', long = "sample_id")]
    pub sample_id: Option<String>,

    /// Target region tag.
    #[arg(long = "target_tag", default_value = "exome")]
    pub target_tag: String,

    /// Reference genome tag (short-circuits assembly detection).
    #[arg(long = "ref_genome_tag", default_value = "hg38")]
    pub ref_genome_tag: String,

    /// Number of threads to use.
    #[arg(short = 't', long = "threads", default_value_t = 10)]
    pub threads: usize,

    /// Mapping-quality cutoff.
    #[arg(long = "mq_cutoff", default_value_t = 41)]
    pub mq_cutoff: i32,
}

/// Preparation-phase arguments (Python `_add_preparation_args`,
/// `SDrecall:390-397`).
#[derive(clap::Args, Debug, Clone)]
pub struct PreparationArgs {
    /// High-quality depth cutoff.
    #[arg(long = "high_quality_depth", default_value_t = 10)]
    pub high_quality_depth: i32,

    /// Minimum depth cutoff.
    #[arg(long = "minimum_depth", default_value_t = 5)]
    pub minimum_depth: i32,

    /// Multi-align fraction cutoff.
    #[arg(long = "multialign_frac", default_value_t = 0.5)]
    pub multialign_frac: f64,
}

/// Realignment-phase arguments (Python `_add_realignment_args`,
/// `SDrecall:400-403`).
#[derive(clap::Args, Debug, Clone)]
pub struct RealignmentArgs {
    /// Number of threads for numba acceleration (kept for CLI parity; the Rust
    /// path uses rayon, so this becomes the per-island inner-thread hint).
    #[arg(long = "numba_threads", default_value_t = 2)]
    pub numba_threads: usize,

    /// Abort the run if any island's false-positive control fails or panics.
    /// The default (`false`) mirrors the Python pipeline — per-island failures are
    /// tolerated and the run continues — but unlike Python they are now
    /// ERROR-logged and written to a `<sample>.failed_islands.tsv` manifest, so a
    /// dropped island's variants are never lost silently.
    #[arg(long = "strict_islands", default_value_t = false)]
    pub strict_islands: bool,
}

/// Conventional-VCF merge arguments (Python `_add_conventional_vcf_args`,
/// `SDrecall:406-413`).
#[derive(clap::Args, Debug, Clone)]
pub struct ConventionalVcfArgs {
    /// VCF from a conventional caller (e.g. GATK, DeepVariant) to merge with.
    #[arg(long = "conventional_vcf")]
    pub conventional_vcf: Option<PathBuf>,

    /// Name of the conventional caller.
    #[arg(long = "caller_name", default_value = "conventional")]
    pub caller_name: String,

    /// Output path for the merged VCF (default derived from `conventional_vcf`).
    #[arg(long = "merged_vcf")]
    pub merged_vcf: Option<PathBuf>,
}

/// Cohort-VCF annotation arguments (Python `_add_cohort_args`,
/// `SDrecall:416-423`).
#[derive(clap::Args, Debug, Clone)]
pub struct CohortArgs {
    /// Cohort-level VCF with control samples.
    #[arg(long = "cohort_vcf")]
    pub cohort_vcf: Option<PathBuf>,

    /// Frequency cutoff for common variants in the cohort.
    #[arg(long = "inhouse_common_cutoff", default_value_t = 0.01)]
    pub inhouse_common_cutoff: f64,

    /// Confidence-level threshold for common-variant determination.
    #[arg(long = "cohort_conf_level", default_value_t = 0.999)]
    pub cohort_conf_level: f64,
}

/// `run`: full pipeline — common + preparation + realignment + conventional +
/// cohort args (Python `full_parser`, `SDrecall:320-326`).
#[derive(clap::Args, Debug)]
pub struct RunArgs {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub prep: PreparationArgs,
    #[command(flatten)]
    pub realign: RealignmentArgs,
    #[command(flatten)]
    pub conventional: ConventionalVcfArgs,
    #[command(flatten)]
    pub cohort: CohortArgs,
}

/// `prepare`: common + preparation args (Python `prep_parser`,
/// `SDrecall:329-332`).
#[derive(clap::Args, Debug)]
pub struct PrepareArgs {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub prep: PreparationArgs,
}

/// `realign`: common + realignment + conventional + cohort args (Python
/// `realign_parser`, `SDrecall:335-339`).
#[derive(clap::Args, Debug)]
pub struct RealignArgs {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub realign: RealignmentArgs,
    #[command(flatten)]
    pub conventional: ConventionalVcfArgs,
    #[command(flatten)]
    pub cohort: CohortArgs,
}

impl CommonArgs {
    /// Validate the inputs that the Python `main()` checks up-front
    /// (`SDrecall:346-348`): the reference genome must end in `.fasta`.
    ///
    /// File-existence checks are I/O and deferred to T9 integration (the stage
    /// crates fail with typed errors when a file is missing); the scaffold
    /// validates only the cheap, pure preconditions.
    pub fn validate(&self) -> Result<(), String> {
        let ref_str = self.ref_genome.to_string_lossy();
        if !ref_str.ends_with(".fasta") {
            return Err(format!(
                "Reference genome file name suffix has to be .fasta, but {ref_str} is not."
            ));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_definition_is_valid() {
        // clap's own consistency assertions (duplicate args, bad defaults, …).
        Cli::command().debug_assert();
    }

    #[test]
    fn parses_run_with_short_flags() {
        let cli = Cli::try_parse_from([
            "sdrecall",
            "run",
            "-r",
            "/refs/g.fasta",
            "-o",
            "/tmp/out",
            "-i",
            "/data/s.bam",
            "-m",
            "/refs/map.tsv",
            "-b",
            "/refs/t.bed",
            "-t",
            "8",
        ])
        .expect("run subcommand should parse");
        match cli.command {
            Command::Run(args) => {
                assert_eq!(args.common.threads, 8);
                assert_eq!(args.common.target_tag, "exome");
                assert_eq!(args.common.mq_cutoff, 41);
                assert_eq!(args.prep.minimum_depth, 5);
                assert_eq!(args.realign.numba_threads, 2);
                assert!(args.conventional.conventional_vcf.is_none());
            }
            _ => panic!("expected Run subcommand"),
        }
    }

    #[test]
    fn validate_rejects_non_fasta_reference() {
        let cli = Cli::try_parse_from([
            "sdrecall",
            "prepare",
            "-r",
            "/refs/genome.fa",
            "-o",
            "/tmp/out",
            "-i",
            "/data/s.bam",
            "-m",
            "/refs/map.tsv",
            "-b",
            "/refs/t.bed",
        ])
        .unwrap();
        let common = match cli.command {
            Command::Prepare(a) => a.common,
            _ => unreachable!(),
        };
        assert!(common.validate().is_err());
    }

    #[test]
    fn verbosity_maps_to_log_filter() {
        assert_eq!(Verbosity::Debug.to_level_filter(), log::LevelFilter::Debug);
        assert_eq!(
            Verbosity::Critical.to_level_filter(),
            log::LevelFilter::Error
        );
    }
}
