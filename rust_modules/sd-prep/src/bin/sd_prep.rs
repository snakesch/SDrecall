//! `sd-prep` CLI — exposes the IMPLEMENTED Phase-1 graph units (the deliverable).
//!
//! The full Phase-1 driver (depth → SD pairs → graph → traversal → grouping →
//! masked genomes) is staged; the traversal + I/O parts are STUBBED (see the lib
//! module docs). This CLI wires the parts that ARE done so they can be run /
//! differentially checked against the Python outputs:
//!
//!   sd-prep graph --sd-map filtered_SD_binary_map.tsv
//!     → builds the multiplex graph, prints node/edge/SD-edge counts + the
//!       all-edge component partition sizes (vs `gt.label_components`).
//!
//! The `traverse` / `prepare` subcommands are intentionally absent until the
//! traversal route-walk + minimap2 FFI land (TODO(T8)).

use clap::{Parser, Subcommand};
use sd_prep::graph_build::{build_multiplex_graph, NodeKey, SdPairRow};
use sd_prep::graph_core::component_labels;
use sdrecall_utils::Strand;
use std::path::PathBuf;
use std::process::ExitCode;

#[derive(Parser)]
#[command(name = "sd-prep", about = "SDrecall Phase-1 graph units (T8, partial)")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,

    /// Log level: error | warn | info | debug | trace.
    #[arg(long, default_value = "info", global = true)]
    log_level: String,
}

#[derive(Subcommand)]
enum Cmd {
    /// Build the multiplex SD+PO graph from a `filtered_SD_binary_map.tsv` and
    /// report node/edge/SD-edge counts + component partition sizes.
    Graph {
        /// The 9-col filtered SD binary map (chr_1,start_1,end_1,strand1,chr_2,
        /// start_2,end_2,strand2,mismatch_rate) — Python's `filtered_SD_binary_map.tsv`.
        #[arg(long)]
        sd_map: PathBuf,
        /// Threads (forwarded to the per-chr PO build; currently sequential).
        #[arg(long, default_value_t = 1)]
        threads: usize,
    },
    /// Build a masked genome from a query BED (isolates the masking logic for the
    /// md5 differential against Python's `RG<n>.masked.fasta`).
    Mask {
        /// Query BED (BED3/BED6) — the per-RG query intervals.
        #[arg(long)]
        query_bed: PathBuf,
        /// Reference genome FASTA (with `.fai`).
        #[arg(long)]
        ref_fa: PathBuf,
        /// Output masked FASTA path.
        #[arg(long)]
        out: PathBuf,
        /// Average fragment size (median insert size).
        #[arg(long, default_value_t = 569.9)]
        avg_frag: f64,
        /// Fragment-size standard deviation.
        #[arg(long, default_value_t = 151.9)]
        std_frag: f64,
    },
}

/// Parse a strand token (`+`/`-`) into [`Strand`].
fn parse_strand(s: &str) -> Strand {
    match s {
        "+" => Strand::Forward,
        "-" => Strand::Reverse,
        _ => Strand::Unknown,
    }
}

/// Read the filtered SD binary map TSV (header row skipped) into `SdPairRow`s.
fn read_sd_map(path: &PathBuf) -> sdrecall_utils::Result<Vec<SdPairRow>> {
    let text = std::fs::read_to_string(path).map_err(|e| sdrecall_utils::SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    let mut rows = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if i == 0 && line.starts_with("chr_1") {
            continue; // header
        }
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 {
            return Err(sdrecall_utils::SdError::BedParse {
                line: i + 1,
                msg: format!("expected ≥9 columns, got {}", f.len()),
            });
        }
        let parse_i64 = |s: &str, col: &str| -> sdrecall_utils::Result<i64> {
            s.parse::<i64>()
                .map_err(|_| sdrecall_utils::SdError::BedParse {
                    line: i + 1,
                    msg: format!("non-integer {col}: {s:?}"),
                })
        };
        let a = NodeKey::new(
            f[0],
            parse_i64(f[1], "start_1")?,
            parse_i64(f[2], "end_1")?,
            parse_strand(f[3]),
        );
        let b = NodeKey::new(
            f[4],
            parse_i64(f[5], "start_2")?,
            parse_i64(f[6], "end_2")?,
            parse_strand(f[7]),
        );
        let mismatch_rate = f[8].parse::<f64>().unwrap_or(0.0);
        rows.push(SdPairRow {
            a,
            b,
            mismatch_rate,
        });
    }
    Ok(rows)
}

fn run_graph(sd_map: &PathBuf, threads: usize) -> sdrecall_utils::Result<()> {
    let rows = read_sd_map(sd_map)?;
    log::info!("Read {} SD pair rows from {}", rows.len(), sd_map.display());
    let g = build_multiplex_graph(&rows, threads);
    let sd_edges = g.g.edge_weights().filter(|a| a.is_sd).count();
    let po_edges = g.edge_count() - sd_edges;

    let labels = component_labels(&g, |_| true);
    let n_components = labels.iter().copied().max().map(|m| m + 1).unwrap_or(0);
    // Component size histogram.
    let mut sizes: rustc_hash::FxHashMap<u32, usize> = rustc_hash::FxHashMap::default();
    for &l in &labels {
        *sizes.entry(l).or_insert(0) += 1;
    }
    let mut size_vec: Vec<usize> = sizes.into_values().collect();
    size_vec.sort_unstable_by(|a, b| b.cmp(a));

    println!("nodes\t{}", g.node_count());
    println!("edges\t{}", g.edge_count());
    println!("sd_edges\t{sd_edges}");
    println!("po_edges\t{po_edges}");
    println!("components\t{n_components}");
    println!(
        "component_sizes\t{}",
        size_vec
            .iter()
            .map(|s| s.to_string())
            .collect::<Vec<_>>()
            .join(",")
    );
    Ok(())
}

/// Build a masked genome from a query BED and print its md5 (isolates the masking
/// logic for the differential against Python's `RG<n>.masked.fasta`).
fn run_mask(
    query_bed: &std::path::Path,
    ref_fa: &std::path::Path,
    out: &std::path::Path,
    avg_frag: f64,
    std_frag: f64,
) -> sdrecall_utils::Result<()> {
    let query = sdrecall_io::read_bed(query_bed)?;
    log::info!(
        "Masking {} query intervals from {}",
        query.len(),
        query_bed.display()
    );
    sd_prep::mask_genome(&query, ref_fa, out, avg_frag, std_frag, 1000)?;
    // print md5 + contig count.
    use std::io::Read;
    let mut f = std::fs::File::open(out).map_err(|e| sdrecall_utils::SdError::Io {
        path: out.display().to_string(),
        source: e,
    })?;
    let mut ctx = md5::Context::new();
    let mut buf = [0u8; 65536];
    loop {
        let n = f.read(&mut buf).map_err(|e| sdrecall_utils::SdError::Io {
            path: out.display().to_string(),
            source: e,
        })?;
        if n == 0 {
            break;
        }
        ctx.consume(&buf[..n]);
    }
    println!("masked_fasta\t{}", out.display());
    println!("md5\t{:x}", ctx.compute());
    Ok(())
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let level = cli.log_level.parse().unwrap_or(log::LevelFilter::Info);
    sdrecall_utils::init_console_logger(level);

    let res = match &cli.cmd {
        Cmd::Graph { sd_map, threads } => run_graph(sd_map, *threads),
        Cmd::Mask {
            query_bed,
            ref_fa,
            out,
            avg_frag,
            std_frag,
        } => run_mask(query_bed, ref_fa, out, *avg_frag, *std_frag),
    };
    match res {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("sd-prep error: {e}");
            ExitCode::FAILURE
        }
    }
}
