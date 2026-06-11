//! Greedy-Clique-Expansion core, ported from `fp_control/gce_algorithm.py`.
//!
//! `heuristic_largest_clique` finds one clique by seed-and-expand on a CSR weight
//! sub-matrix; `gce_algorithm` repeatedly carves cliques out of a component until no
//! vertices remain. The port is 1:1 with the numba version, including:
//!   * `>=` tie-breaking in the seed/argmax (ties resolve to the highest index),
//!   * the 9-member lookback that lets an earlier clique member contribute a stronger
//!     extension edge, gated by a `1e-4` "the candidate doesn't prefer someone else" test,
//!   * the cutoff backtrack that scans all prior members when the greedy step stalls.

use std::collections::HashSet;

use ndarray::Array2;

use crate::kernels::{
    and_masks, count_true, efficient_mask, efficient_row_max, max_idx_mem, reverse_boolean_mask,
    row_wise_max_with_mask, Csr,
};

/// Find the single largest-edge-weight clique reachable from the best seed within
/// `initial_index_mask`. Returns `(clique_local_indices, remaining_mask)` where indices
/// are local to this CSR and `remaining_mask` marks the in-scope vertices NOT taken.
pub fn heuristic_largest_clique(
    csr: &Csr,
    initial_index_mask: &[bool],
    cutoff: f32,
) -> (Vec<i32>, Vec<bool>) {
    let size = csr.size;

    // Seed: the vertex carrying the single highest valid edge.
    let row_wise = row_wise_max_with_mask(csr, initial_index_mask, -1.0);
    let (mut max_row_ind, _seed_val) = max_idx_mem(&row_wise, Some(initial_index_mask));

    let mut select_indices: Vec<i32> = Vec::with_capacity(size);
    select_indices.push(max_row_ind);

    let mut index_mask = initial_index_mask.to_vec();
    {
        let (rd, rc) = csr.row(max_row_ind as usize);
        let mut neg1 = efficient_mask(rd, rc, size, -1.0);
        neg1[max_row_ind as usize] = false; // exclude the seed itself
        index_mask = and_masks(&index_mask, &neg1);
    }

    let mut i = 1usize;
    while index_mask.iter().any(|&b| b) {
        // Greedy best extension from the most-recently-added member.
        let (mut next_max_ind, mut next_max_value) = {
            let (rd, rc) = csr.row(max_row_ind as usize);
            efficient_row_max(rd, rc, &index_mask)
        };

        // Lookback over the last 9 members: an earlier member may offer a better edge,
        // but only if the candidate it points to doesn't itself prefer a different vertex.
        let lookback_start = i.saturating_sub(9);
        for &member in &select_indices[lookback_start..i] {
            let trial_member = member as usize;
            let (trial_max_ind, trial_max_value) = {
                let (td, tc) = csr.row(trial_member);
                efficient_row_max(td, tc, &index_mask)
            };
            if trial_max_value <= cutoff {
                continue;
            }
            let (_tm_ind, tm_value) = {
                let (bd, bc) = csr.row(trial_max_ind as usize);
                efficient_row_max(bd, bc, &index_mask)
            };
            if tm_value - trial_max_value > 1e-4 {
                continue;
            }
            if trial_max_value > next_max_value {
                next_max_ind = trial_max_ind;
                next_max_value = trial_max_value;
            }
        }

        if next_max_value <= cutoff {
            // Greedy stalled: scan every prior member for an above-cutoff extension.
            let mut trial = 0usize;
            while next_max_value <= cutoff && trial < i {
                let sm = select_indices[trial] as usize;
                let (sd, sc) = csr.row(sm);
                let (ni, nv) = efficient_row_max(sd, sc, &index_mask);
                next_max_ind = ni;
                next_max_value = nv;
                trial += 1;
            }
            if next_max_value <= cutoff && trial >= i {
                break; // cannot extend the clique any further
            } else {
                max_row_ind = next_max_ind;
            }
        } else {
            max_row_ind = next_max_ind;
        }

        select_indices.push(max_row_ind);
        i += 1;
        if i >= size {
            break;
        }

        let (rd, rc) = csr.row(max_row_ind as usize);
        let mut neg1 = efficient_mask(rd, rc, size, -1.0);
        neg1[max_row_ind as usize] = false;
        index_mask = and_masks(&index_mask, &neg1);
    }

    let remain = reverse_boolean_mask(size, &select_indices);
    let remain = and_masks(&remain, initial_index_mask);
    (select_indices, remain)
}

/// Iteratively carve cliques out of one component's weight sub-matrix.
///
/// `selected_indices[k]` maps local row `k` of `weight_matrix` back to its index in the
/// original graph. Returns each clique as a set of **original-graph** vertex indices,
/// in discovery order. Mirrors the `gce_algorithm` generator.
pub fn gce_algorithm(
    selected_indices: &[i32],
    weight_matrix: &Array2<f32>,
    cutoff: f32,
) -> Vec<HashSet<i32>> {
    let mut cliques: Vec<HashSet<i32>> = Vec::new();
    let mut sel = selected_indices.to_vec();
    let mut csr = Csr::from_dense(weight_matrix);
    let mut size = csr.size;
    assert_eq!(sel.len(), size, "selected_indices must match the sub-matrix size");

    loop {
        if size == 0 {
            break;
        }
        let index_mask = vec![true; size];
        let (clique_local, drop_mask) = heuristic_largest_clique(&csr, &index_mask, cutoff);

        let raw: Vec<i32> = clique_local.iter().map(|&c| sel[c as usize]).collect();
        let remain_count = count_true(&drop_mask);

        if !raw.is_empty() {
            cliques.push(raw.into_iter().collect());
        } else if remain_count > 0 {
            log::error!(
                "GCE produced an empty clique while {remain_count} vertices remain; aborting component"
            );
            break;
        }

        if remain_count == 0 {
            break;
        }

        // Shrink to the vertices not yet assigned to a clique, reindexing in lockstep.
        sel = drop_mask
            .iter()
            .enumerate()
            .filter_map(|(k, &keep)| if keep { Some(sel[k]) } else { None })
            .collect();
        csr = csr.select(&drop_mask);
        size = csr.size;
    }

    cliques
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    fn clique_set(cliques: &[HashSet<i32>]) -> HashSet<Vec<i32>> {
        cliques
            .iter()
            .map(|c| {
                let mut v: Vec<i32> = c.iter().copied().collect();
                v.sort_unstable();
                v
            })
            .collect()
    }

    #[test]
    fn fully_connected_yields_single_clique() {
        let m = array![
            [0.0, 0.8, 0.8, 0.8],
            [0.8, 0.0, 0.8, 0.8],
            [0.8, 0.8, 0.0, 0.8],
            [0.8, 0.8, 0.8, 0.0]
        ];
        let sel: Vec<i32> = vec![0, 1, 2, 3];
        let cliques = gce_algorithm(&sel, &m, 0.2);
        assert_eq!(cliques.len(), 1);
        let mut v: Vec<i32> = cliques[0].iter().copied().collect();
        v.sort_unstable();
        assert_eq!(v, vec![0, 1, 2, 3]);
    }

    #[test]
    fn two_triangles_split_on_incompatibility() {
        // {0,1,2} and {3,4,5} are internally cohesive (0.8); every cross pair is -1.
        let n1 = -1.0f32;
        let m = array![
            [0.0, 0.8, 0.8, n1, n1, n1],
            [0.8, 0.0, 0.8, n1, n1, n1],
            [0.8, 0.8, 0.0, n1, n1, n1],
            [n1, n1, n1, 0.0, 0.8, 0.8],
            [n1, n1, n1, 0.8, 0.0, 0.8],
            [n1, n1, n1, 0.8, 0.8, 0.0]
        ];
        let sel: Vec<i32> = (0..6).collect();
        let cliques = gce_algorithm(&sel, &m, 0.2);
        let got = clique_set(&cliques);
        let mut want = HashSet::new();
        want.insert(vec![0, 1, 2]);
        want.insert(vec![3, 4, 5]);
        assert_eq!(got, want);
    }

    #[test]
    fn selected_indices_remap_to_original_graph() {
        // Same two triangles, but the sub-matrix rows map to original ids 10..16.
        let n1 = -1.0f32;
        let m = array![
            [0.0, 0.8, 0.8, n1, n1, n1],
            [0.8, 0.0, 0.8, n1, n1, n1],
            [0.8, 0.8, 0.0, n1, n1, n1],
            [n1, n1, n1, 0.0, 0.8, 0.8],
            [n1, n1, n1, 0.8, 0.0, 0.8],
            [n1, n1, n1, 0.8, 0.8, 0.0]
        ];
        let sel: Vec<i32> = vec![10, 11, 12, 13, 14, 15];
        let cliques = gce_algorithm(&sel, &m, 0.2);
        let got = clique_set(&cliques);
        let mut want = HashSet::new();
        want.insert(vec![10, 11, 12]);
        want.insert(vec![13, 14, 15]);
        assert_eq!(got, want);
    }
}
