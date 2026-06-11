//! `phasing` CLI / differential harness.
//!
//! Reads the per-island dumps written by `fp_control/diff_dump.py::dump_phasing`
//! (weight matrix + graph edges + read vectors + the Python partition), runs the Rust
//! phasing on identical inputs, and asserts the qname **partition** matches up to relabeling.
//!
//! Usage:
//!   phasing --path <DUMP_ROOT>            # compare every island under a dump root
//!   phasing --path <ISLAND_DIR> --single  # compare a single island directory

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use ndarray::Array2;
use serde::Deserialize;

use phasing::{phase, qname_partition, PhasingInput};

#[derive(Parser, Debug)]
#[command(about = "Differential harness for the Rust phasing port (T3)")]
struct Args {
    /// Dump root (one sub-dir per island) or a single island directory with `--single`.
    #[arg(long)]
    path: PathBuf,
    /// Treat `--path` as one island directory rather than a root of islands.
    #[arg(long)]
    single: bool,
}

#[derive(Deserialize)]
struct Meta {
    edge_weight_cutoff: f32,
}

fn read_json<T: for<'de> Deserialize<'de>>(p: &Path) -> Result<T> {
    let s = std::fs::read_to_string(p).with_context(|| format!("reading {}", p.display()))?;
    serde_json::from_str(&s).with_context(|| format!("parsing {}", p.display()))
}

struct Island {
    input: PhasingInput,
    vertex_qname: Vec<String>,
    py_partition: HashSet<Vec<String>>,
}

fn load_island(dir: &Path) -> Result<Island> {
    let wm: Array2<f32> =
        ndarray_npy::read_npy(dir.join("weight_matrix.npy")).context("reading weight_matrix.npy")?;
    let n = wm.nrows();

    let meta: Meta = read_json(&dir.join("meta.json"))?;

    let edges_raw: Vec<[i32; 2]> = read_json(&dir.join("edges.json"))?;
    let edges: Vec<(i32, i32)> = edges_raw.into_iter().map(|e| (e[0], e[1])).collect();

    let nri_map: HashMap<String, Vec<String>> = read_json(&dir.join("node_read_ids.json"))?;
    let mut node_read_ids = vec![Vec::new(); n];
    for (k, v) in nri_map {
        let idx: usize = k.parse().with_context(|| format!("vertex key {k}"))?;
        if idx < n {
            node_read_ids[idx] = v;
        }
    }

    let vq_map: HashMap<String, String> = read_json(&dir.join("vertex_qname.json"))?;
    let mut vertex_qname = vec![String::new(); n];
    for (k, v) in vq_map {
        let idx: usize = k.parse().with_context(|| format!("vertex key {k}"))?;
        if idx < n {
            vertex_qname[idx] = v;
        }
    }

    let read_hap: HashMap<String, Vec<i16>> =
        read_json(&dir.join("read_hap.json")).unwrap_or_default();
    let read_err: HashMap<String, Vec<f32>> =
        read_json(&dir.join("read_err.json")).unwrap_or_default();

    let py_hap: HashMap<String, Vec<String>> = read_json(&dir.join("phasing_hap_qname_info.json"))?;
    let py_partition: HashSet<Vec<String>> = py_hap
        .into_values()
        .map(|mut v| {
            v.sort();
            v
        })
        .collect();

    Ok(Island {
        input: PhasingInput {
            weight_matrix: wm,
            edges,
            edge_weight_cutoff: meta.edge_weight_cutoff,
            node_read_ids,
            read_hap,
            read_err,
        },
        vertex_qname,
        py_partition,
    })
}

fn compare_island(dir: &Path) -> Result<bool> {
    let island = load_island(dir)?;
    let vh = phase(&island.input);
    let rust_partition = qname_partition(&vh, &island.vertex_qname);

    let matched = rust_partition == island.py_partition;
    let name = dir.file_name().and_then(|s| s.to_str()).unwrap_or("?");
    if matched {
        println!(
            "  MATCH    {name}  ({} haplotypes, {} qnames)",
            rust_partition.len(),
            rust_partition.iter().map(|s| s.len()).sum::<usize>()
        );
    } else {
        println!(
            "  MISMATCH {name}  rust={} haps py={} haps",
            rust_partition.len(),
            island.py_partition.len()
        );
        // surface a couple of differing groups for triage
        for g in rust_partition.difference(&island.py_partition).take(2) {
            println!("    only-rust group (n={}): {:?}", g.len(), &g[..g.len().min(4)]);
        }
        for g in island.py_partition.difference(&rust_partition).take(2) {
            println!("    only-py   group (n={}): {:?}", g.len(), &g[..g.len().min(4)]);
        }
    }
    Ok(matched)
}

fn island_dirs(root: &Path) -> Result<Vec<PathBuf>> {
    let mut dirs = Vec::new();
    for entry in std::fs::read_dir(root).with_context(|| format!("listing {}", root.display()))? {
        let entry = entry?;
        let p = entry.path();
        if p.is_dir() && p.join("weight_matrix.npy").exists() {
            dirs.push(p);
        }
    }
    dirs.sort();
    Ok(dirs)
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Args::parse();

    let dirs = if args.single {
        vec![args.path.clone()]
    } else {
        island_dirs(&args.path)?
    };

    if dirs.is_empty() {
        println!("No island dumps found under {}", args.path.display());
        return Ok(());
    }

    println!("Comparing {} island(s):", dirs.len());
    let mut total = 0usize;
    let mut matched = 0usize;
    let mut errored = 0usize;
    for dir in &dirs {
        total += 1;
        match compare_island(dir) {
            Ok(true) => matched += 1,
            Ok(false) => {}
            Err(e) => {
                errored += 1;
                println!("  ERROR    {}: {e:#}", dir.display());
            }
        }
    }

    println!("\nResult: {matched}/{total} islands match");
    if errored > 0 {
        println!("{errored} island(s) errored while loading");
    }
    if matched == total {
        Ok(())
    } else {
        std::process::exit(1);
    }
}
