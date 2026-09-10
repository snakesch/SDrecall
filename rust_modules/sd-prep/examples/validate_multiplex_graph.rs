//! Differential harness: build the multiplex SD+PO graph from the Python
//! `filtered_SD_binary_map.tsv` and compare structural counts against the Python
//! `*_multiplexed_SDs.graphml` ground truth.
//!
//! Pass criteria for the structural part of the T8 differential:
//!   - node_count   == Python multiplex graph `<node>` count
//!   - sd_edge_count == filtered SD map row count == Python `segmental_duplication` count
//!   - total_edge_count == Python multiplex graph `<edge>` count
//!
//! Run:
//!   cargo run -p sd-prep --example validate_multiplex_graph -- \
//!     <filtered_SD_binary_map.tsv> <multiplexed_SDs.graphml>
//!
//! The graphml is parsed only to count `<node>` / `<edge>` / `segmental_duplication`
//! lines (a grep, not a full parse) so this harness has no extra deps.

use sd_prep::graph_build::{build_multiplex_graph, NodeKey, SdPairRow};
use sd_prep::graph_core::component_labels;
use sdrecall_utils::Strand;
use std::path::PathBuf;

fn parse_strand(s: &str) -> Strand {
    match s {
        "+" => Strand::Forward,
        "-" => Strand::Reverse,
        _ => Strand::Unknown,
    }
}

fn read_sd_map(path: &PathBuf) -> Vec<SdPairRow> {
    let text = std::fs::read_to_string(path).expect("read sd map");
    let mut rows = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if (i == 0 && line.starts_with("chr_1")) || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        assert!(f.len() >= 9, "row {i} has {} cols", f.len());
        let a = NodeKey::new(
            f[0],
            f[1].parse().unwrap(),
            f[2].parse().unwrap(),
            parse_strand(f[3]),
        );
        let b = NodeKey::new(
            f[4],
            f[5].parse().unwrap(),
            f[6].parse().unwrap(),
            parse_strand(f[7]),
        );
        rows.push(SdPairRow {
            a,
            b,
            mismatch_rate: f[8].parse().unwrap_or(0.0),
        });
    }
    rows
}

/// Count `<node>`, `<edge>`, `segmental_duplication` occurrences in a graphml.
fn graphml_counts(path: &PathBuf) -> (usize, usize, usize) {
    let text = std::fs::read_to_string(path).expect("read graphml");
    let nodes = text.matches("<node ").count();
    let edges = text.matches("<edge ").count();
    let sd = text.matches("segmental_duplication").count();
    (nodes, edges, sd)
}

fn main() {
    sdrecall_utils::init_console_logger(log::LevelFilter::Info);
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "usage: {} <filtered_SD_binary_map.tsv> <multiplexed_SDs.graphml>",
            args[0]
        );
        std::process::exit(2);
    }
    let sd_map = PathBuf::from(&args[1]);
    let graphml = PathBuf::from(&args[2]);

    let rows = read_sd_map(&sd_map);
    println!("Read {} SD pair rows", rows.len());

    let g = build_multiplex_graph(&rows, 1);
    let sd_edges = g.g.edge_weights().filter(|a| a.is_sd).count();
    let total_edges = g.edge_count();
    let nodes = g.node_count();

    let labels = component_labels(&g, |_| true);
    let n_comp = labels.iter().copied().max().map(|m| m + 1).unwrap_or(0);

    let (py_nodes, py_edges, py_sd) = graphml_counts(&graphml);

    println!("─────────────────────────────────────────────");
    println!("metric            rust      python");
    println!("nodes             {nodes:<9} {py_nodes}");
    println!("total_edges       {total_edges:<9} {py_edges}");
    println!("sd_edges          {sd_edges:<9} {py_sd}");
    println!("components        {n_comp}");
    println!("─────────────────────────────────────────────");

    let mut ok = true;
    if nodes != py_nodes {
        eprintln!("MISMATCH nodes: rust {nodes} != python {py_nodes}");
        ok = false;
    }
    if sd_edges != py_sd {
        eprintln!("MISMATCH sd_edges: rust {sd_edges} != python {py_sd}");
        ok = false;
    }
    if total_edges != py_edges {
        eprintln!(
            "DIFF total_edges: rust {total_edges} != python {py_edges} (PO-edge dedup / tie-break differences are possible; investigate if large)"
        );
        // total-edge equality is the strictest check; treat a difference as a
        // soft fail to surface but not abort, since PO direction tie-breaks on
        // size ties can legitimately differ in count near book-ended overlaps.
        ok = false;
    }
    if ok {
        println!("PASS: structural counts match Python multiplex graph");
    } else {
        println!("FAIL: see mismatches above");
        std::process::exit(1);
    }
}
