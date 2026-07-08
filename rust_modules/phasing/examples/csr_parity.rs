//! Minimal dense-vs-CSR parity harness for the phasing matrix/GCE refactor.
//!
//! Run from `rust_modules/`:
//!
//! ```text
//! cargo run -p phasing --example csr_parity
//! ```
//!
//! The dense path mirrors the current source shape: a square `Array2<f32>` with
//! diagonal `1.0`, zero/no-overlap entries, and positive same-haplotype weights.
//! The CSR path builds the equivalent sparse representation directly from entries,
//! then runs the CSR-native GCE and phasing wrappers.

use std::collections::{HashMap, HashSet};

use ndarray::Array2;
use phasing::gce::{gce_algorithm, gce_algorithm_csr};
use phasing::kernels::Csr;
use phasing::{phase, phase_sparse, qname_partition, PhasingInput, SparsePhasingInput};

fn canonical_cliques(cliques: &[HashSet<i32>]) -> HashSet<Vec<i32>> {
    cliques
        .iter()
        .map(|clique| {
            let mut vertices: Vec<i32> = clique.iter().copied().collect();
            vertices.sort_unstable();
            vertices
        })
        .collect()
}

fn build_dense_current_source_graph() -> (Array2<f32>, Vec<(i32, i32)>) {
    let n = 14usize;
    let mut weights = Array2::<f32>::zeros((n, n));
    let mut edges = Vec::new();

    for i in 0..n {
        weights[[i, i]] = 1.0;
    }

    for a in 0..n {
        for b in (a + 1)..n {
            let same_component = (a < 7 && b < 7) || (a >= 7 && b >= 7);
            if !same_component {
                continue; // zero/no-overlap entry: present in dense, absent from CSR
            }
            let weight = 0.8;
            weights[[a, b]] = weight;
            weights[[b, a]] = weight;
            edges.push((a as i32, b as i32));
        }
    }

    (weights, edges)
}

fn build_refactored_csr_graph(size: usize) -> Csr {
    let mut entries = Vec::new();

    for i in 0..size {
        entries.push((i, i, 1.0));
    }

    for a in 0..size {
        for b in (a + 1)..size {
            let same_component = (a < 7 && b < 7) || (a >= 7 && b >= 7);
            if !same_component {
                continue; // zero/no-overlap entry: present in dense, absent from CSR
            }
            let weight = 0.8;
            entries.push((a, b, weight));
            entries.push((b, a, weight));
        }
    }

    Csr::from_entries(size, entries)
}

fn empty_read_evidence() -> (
    Vec<Vec<String>>,
    HashMap<String, Vec<i16>>,
    HashMap<String, Vec<f32>>,
) {
    (Vec::new(), HashMap::new(), HashMap::new())
}

fn main() {
    let (dense_weights, edges) = build_dense_current_source_graph();
    let csr_weights = build_refactored_csr_graph(dense_weights.nrows());
    let selected_indices: Vec<i32> = (0..dense_weights.nrows() as i32).collect();
    let cutoff = 0.301;

    let dense_cliques = gce_algorithm(&selected_indices, &dense_weights, cutoff);
    let sparse_cliques = gce_algorithm_csr(&selected_indices, csr_weights.clone(), cutoff);
    assert_eq!(
        canonical_cliques(&sparse_cliques),
        canonical_cliques(&dense_cliques)
    );

    let (node_read_ids, read_hap, read_err) = empty_read_evidence();
    let dense_input = PhasingInput {
        weight_matrix: dense_weights.clone(),
        edges: edges.clone(),
        edge_weight_cutoff: cutoff,
        node_read_ids: node_read_ids.clone(),
        read_hap: read_hap.clone(),
        read_err: read_err.clone(),
    };
    let sparse_input = SparsePhasingInput {
        weights: csr_weights.clone(),
        edges,
        edge_weight_cutoff: cutoff,
        node_read_ids,
        read_hap,
        read_err,
    };

    let vertex_qname: Vec<String> = (0..dense_weights.nrows())
        .map(|i| format!("read_pair_{i}"))
        .collect();
    let dense_partition = qname_partition(&phase(&dense_input), &vertex_qname);
    let sparse_partition = qname_partition(&phase_sparse(&sparse_input), &vertex_qname);
    assert_eq!(sparse_partition, dense_partition);

    let dense_nonzero = dense_weights.iter().filter(|&&value| value != 0.0).count();
    println!("CSR parity OK");
    println!("  GCE cliques: {}", dense_cliques.len());
    println!("  final haplotypes: {}", dense_partition.len());
    println!("  dense cells: {}", dense_weights.len());
    println!("  dense non-zero entries: {dense_nonzero}");
    println!("  CSR stored entries: {}", csr_weights.data.len());
}
