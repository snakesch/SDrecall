//! Orchestration ported from `fp_control/phasing.py`.
//!
//! Pipeline:
//!   1. connected components of the phasing graph (replaces `gt.label_components`) via union-find;
//!   2. three rounds of seed-and-expand clique finding at decreasing edge-weight cutoffs
//!      (`find_cliques_in_components`), the later rounds rescuing fragmented low-variant read-pairs;
//!   3. split each clique into haplotypes by the connectivity of its *induced* sub-graph plus
//!      supplemented weak edges (`find_components_inside_filtered_cliques`).

use std::collections::{HashMap, HashSet};

use ndarray::Array2;

use crate::gce::{
    gce_algorithm, gce_algorithm_csr_subset_with_context, structural_transpose, GceContext,
};
use crate::kernels::{and_masks, apply_index_mask, dense_submatrix, isin_arange, Csr};

/// Which clique-finding round produced a clique (kept for parity/debugging; does not affect the
/// final partition).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Round {
    Main,
    Second,
    Third,
}

/// Everything the phasing needs, mirroring the arguments of `phasing_realigned_reads`.
pub struct PhasingInput {
    /// N×N weights: `-1` incompatible, `0` no overlap, `0..1` edge weight.
    pub weight_matrix: Array2<f32>,
    /// Phasing-graph adjacency (positive-weight overlaps); defines connected components.
    pub edges: Vec<(i32, i32)>,
    pub edge_weight_cutoff: f32,
    /// `vertex -> [read_id]` (1–2 mate read ids per read-pair vertex).
    pub node_read_ids: Vec<Vec<String>>,
    /// `read_id -> per-base haplotype vector` (1 == matches reference).
    pub read_hap: HashMap<String, Vec<i16>>,
    /// `read_id -> per-base error probability` (≤ 0.03 == high quality).
    pub read_err: HashMap<String, Vec<f32>>,
}

/// Sparse equivalent of [`PhasingInput`] for the CSR refactor.
pub struct SparsePhasingInput {
    /// CSR weights: `-1` incompatible, absent/`0` no overlap, `0..1` edge weight.
    pub weights: Csr,
    /// Phasing-graph adjacency (positive-weight overlaps); defines connected components.
    pub edges: Vec<(i32, i32)>,
    pub edge_weight_cutoff: f32,
    /// `vertex -> [read_id]` (1–2 mate read ids per read-pair vertex).
    pub node_read_ids: Vec<Vec<String>>,
    /// `read_id -> per-base haplotype vector` (1 == matches reference).
    pub read_hap: HashMap<String, Vec<i16>>,
    /// `read_id -> per-base error probability` (≤ 0.03 == high quality).
    pub read_err: HashMap<String, Vec<f32>>,
}

impl From<&PhasingInput> for SparsePhasingInput {
    fn from(input: &PhasingInput) -> Self {
        Self {
            weights: Csr::from_dense(&input.weight_matrix),
            edges: input.edges.clone(),
            edge_weight_cutoff: input.edge_weight_cutoff,
            node_read_ids: input.node_read_ids.clone(),
            read_hap: input.read_hap.clone(),
            read_err: input.read_err.clone(),
        }
    }
}

/// Iterative union-find with path compression and union-by-rank.
struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    fn find(&mut self, x: usize) -> usize {
        let mut root = x;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        // path compression
        let mut cur = x;
        while self.parent[cur] != root {
            let next = self.parent[cur];
            self.parent[cur] = root;
            cur = next;
        }
        root
    }

    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        match self.rank[ra].cmp(&self.rank[rb]) {
            std::cmp::Ordering::Less => self.parent[ra] = rb,
            std::cmp::Ordering::Greater => self.parent[rb] = ra,
            std::cmp::Ordering::Equal => {
                self.parent[rb] = ra;
                self.rank[ra] += 1;
            }
        }
    }
}

/// `vertex -> component id` over `size` vertices. Isolated vertices form singleton components.
/// Component ids are assigned in order of first vertex appearance (matching graph-tool's
/// insertion-ordered labeling), so downstream ordering matches Python.
fn connected_components(size: usize, edges: &[(i32, i32)]) -> Vec<i32> {
    let mut uf = UnionFind::new(size);
    for &(u, v) in edges {
        uf.union(u as usize, v as usize);
    }
    let roots: Vec<usize> = (0..size).map(|v| uf.find(v)).collect();
    let mut comp = vec![0i32; size];
    let mut root_id: HashMap<usize, i32> = HashMap::new();
    let mut next = 0i32;
    for (slot, &r) in comp.iter_mut().zip(roots.iter()) {
        *slot = *root_id.entry(r).or_insert_with(|| {
            let id = next;
            next += 1;
            id
        });
    }
    comp
}

/// `np.all(weight_matrix <= 0.1, axis=1)`: a row is "small" when it has no entry above 0.1.
fn small_row_mask(m: &Array2<f32>) -> Vec<bool> {
    let n = m.nrows();
    (0..n)
        .map(|i| (0..m.ncols()).all(|j| m[[i, j]] <= 0.1))
        .collect()
}

/// Does any member read of `clique` carry a high-quality non-reference base?
/// Mirrors the round-2 no-variant gate in `phasing.py` (an empty high-quality set counts as
/// "no variant", matching numpy's `array([]).all() == True`).
fn clique_has_variant(clique: &HashSet<i32>, input: &PhasingInput) -> bool {
    for &qid in clique {
        let rids = match input.node_read_ids.get(qid as usize) {
            Some(r) => r,
            None => continue,
        };
        for rid in rids {
            let (hap, err) = match (input.read_hap.get(rid), input.read_err.get(rid)) {
                (Some(h), Some(e)) => (h, e),
                _ => continue,
            };
            let mut any_valid = false;
            let mut all_one = true;
            for (h, e) in hap.iter().zip(err.iter()) {
                if *e <= 0.03 {
                    any_valid = true;
                    if *h != 1 {
                        all_one = false;
                        break;
                    }
                }
            }
            if any_valid && !all_one {
                return true;
            }
        }
    }
    false
}

fn clique_has_variant_sparse(clique: &HashSet<i32>, input: &SparsePhasingInput) -> bool {
    for &qid in clique {
        let rids = match input.node_read_ids.get(qid as usize) {
            Some(r) => r,
            None => continue,
        };
        for rid in rids {
            let (hap, err) = match (input.read_hap.get(rid), input.read_err.get(rid)) {
                (Some(h), Some(e)) => (h, e),
                _ => continue,
            };
            let mut any_valid = false;
            let mut all_one = true;
            for (h, e) in hap.iter().zip(err.iter()) {
                if *e <= 0.03 {
                    any_valid = true;
                    if *h != 1 {
                        all_one = false;
                        break;
                    }
                }
            }
            if any_valid && !all_one {
                return true;
            }
        }
    }
    false
}

/// Three-round seed-and-expand clique finding. Returns `(round, clique)` pairs where each clique
/// is a set of original-graph vertex indices.
fn find_cliques_in_components(input: &PhasingInput) -> Vec<(Round, HashSet<i32>)> {
    let m = &input.weight_matrix;
    let size = m.nrows();
    let cutoff = input.edge_weight_cutoff;

    let comp = connected_components(size, &input.edges);
    let mut component_dict: HashMap<i32, Vec<i32>> = HashMap::new();
    for (v, &cid) in comp.iter().enumerate() {
        component_dict.entry(cid).or_default().push(v as i32);
    }

    let small_mask = small_row_mask(m);
    let big_mask: Vec<bool> = small_mask.iter().map(|&b| !b).collect();
    let mut small_row_indices: HashSet<i32> = small_mask
        .iter()
        .enumerate()
        .filter_map(|(i, &b)| if b { Some(i as i32) } else { None })
        .collect();

    let mut result: Vec<(Round, HashSet<i32>)> = Vec::new();

    // Iterate components in ascending id order (== order of first vertex appearance).
    let mut comp_ids: Vec<i32> = component_dict.keys().copied().collect();
    comp_ids.sort_unstable();

    // ---- Round 1: largest cliques per component (cutoff) ----
    for cid in comp_ids {
        let comp_verts = &component_dict[&cid];
        if comp_verts.len() <= 5 {
            small_row_indices.extend(comp_verts.iter().copied());
            continue;
        }
        let mut comp_index_mask = vec![false; size];
        for &v in comp_verts {
            comp_index_mask[v as usize] = true;
        }
        let comp_index_mask = and_masks(&comp_index_mask, &big_mask);
        let selected = apply_index_mask(&comp_index_mask);
        if selected.is_empty() {
            continue; // all of this component's rows are "small"; handled in round 2
        }
        let big_wm = dense_submatrix(m, &comp_index_mask);
        for clique in gce_algorithm(&selected, &big_wm, cutoff) {
            if clique.len() <= 5 {
                small_row_indices.extend(clique);
            } else {
                result.push((Round::Main, clique));
            }
        }
    }

    // ---- Round 2: rescue fragmented read-pairs (cutoff * 2/3) ----
    if !small_row_indices.is_empty() {
        let mask = isin_arange(size, &small_row_indices);
        let selected = apply_index_mask(&mask);
        let small_wm = dense_submatrix(m, &mask);
        small_row_indices = HashSet::new();
        for clique in gce_algorithm(&selected, &small_wm, cutoff * 2.0 / 3.0) {
            if clique.len() > 5 || clique_has_variant(&clique, input) {
                result.push((Round::Second, clique));
            } else {
                small_row_indices.extend(clique);
            }
        }
    }

    // ---- Round 3: assemble whatever remains (cutoff / 3) ----
    if !small_row_indices.is_empty() {
        let mask = isin_arange(size, &small_row_indices);
        let selected = apply_index_mask(&mask);
        let small_wm = dense_submatrix(m, &mask);
        for clique in gce_algorithm(&selected, &small_wm, cutoff / 3.0) {
            result.push((Round::Third, clique));
        }
    }

    result
}

/// CSR-native version of [`find_cliques_in_components`].
fn find_cliques_in_components_sparse(
    input: &SparsePhasingInput,
    threads: usize,
) -> Vec<(Round, HashSet<i32>)> {
    let weights = &input.weights;
    let size = weights.size;
    let cutoff = input.edge_weight_cutoff;
    let gce_ctx = GceContext::with_threads(weights, threads);
    let refinement_transpose = structural_transpose(weights);

    let comp = connected_components(size, &input.edges);
    let mut component_dict: HashMap<i32, Vec<i32>> = HashMap::new();
    for (v, &cid) in comp.iter().enumerate() {
        component_dict.entry(cid).or_default().push(v as i32);
    }

    let small_mask = weights.small_row_mask(0.1);
    let big_mask: Vec<bool> = small_mask.iter().map(|&b| !b).collect();
    let mut small_row_indices: HashSet<i32> = small_mask
        .iter()
        .enumerate()
        .filter_map(|(i, &b)| if b { Some(i as i32) } else { None })
        .collect();

    let mut result: Vec<(Round, HashSet<i32>)> = Vec::new();

    let mut comp_ids: Vec<i32> = component_dict.keys().copied().collect();
    comp_ids.sort_unstable();

    // ---- Round 1: largest cliques per component (cutoff) ----
    for cid in comp_ids {
        let comp_verts = &component_dict[&cid];
        if comp_verts.len() <= 5 {
            small_row_indices.extend(comp_verts.iter().copied());
            continue;
        }
        let mut comp_index_mask = vec![false; size];
        for &v in comp_verts {
            comp_index_mask[v as usize] = true;
        }
        let comp_index_mask = and_masks(&comp_index_mask, &big_mask);
        let selected = apply_index_mask(&comp_index_mask);
        if selected.is_empty() {
            continue;
        }
        for clique in gce_algorithm_csr_subset_with_context(
            &selected,
            weights,
            &gce_ctx,
            Some(&refinement_transpose),
            cutoff,
        ) {
            if clique.len() <= 5 {
                small_row_indices.extend(clique);
            } else {
                result.push((Round::Main, clique));
            }
        }
    }

    // ---- Round 2: rescue fragmented read-pairs (cutoff * 2/3) ----
    if !small_row_indices.is_empty() {
        let mask = isin_arange(size, &small_row_indices);
        let selected = apply_index_mask(&mask);
        small_row_indices = HashSet::new();
        for clique in gce_algorithm_csr_subset_with_context(
            &selected,
            weights,
            &gce_ctx,
            Some(&refinement_transpose),
            cutoff * 2.0 / 3.0,
        ) {
            if clique.len() > 5 || clique_has_variant_sparse(&clique, input) {
                result.push((Round::Second, clique));
            } else {
                small_row_indices.extend(clique);
            }
        }
    }

    // ---- Round 3: assemble whatever remains (cutoff / 3) ----
    if !small_row_indices.is_empty() {
        let mask = isin_arange(size, &small_row_indices);
        let selected = apply_index_mask(&mask);
        for clique in gce_algorithm_csr_subset_with_context(
            &selected,
            weights,
            &gce_ctx,
            Some(&refinement_transpose),
            cutoff / 3.0,
        ) {
            result.push((Round::Third, clique));
        }
    }

    result
}

/// Split each clique into haplotypes by the connectivity of its induced sub-graph (original edges
/// among clique members) plus supplemented weak edges (`0.1 < w <= cutoff`). Returns `vertex -> hap_id`.
fn find_components_inside_cliques(
    cliques: &[(Round, HashSet<i32>)],
    input: &PhasingInput,
) -> HashMap<i32, i32> {
    let m = &input.weight_matrix;
    let cutoff = input.edge_weight_cutoff;

    // Global undirected adjacency of the phasing graph for fast in-clique edge lookup.
    let mut adj: HashMap<i32, HashSet<i32>> = HashMap::new();
    for &(u, v) in &input.edges {
        adj.entry(u).or_default().insert(v);
        adj.entry(v).or_default().insert(u);
    }

    let mut vertex_hap: HashMap<i32, i32> = HashMap::new();
    let mut haplotype_idx = 0i32;

    for (_round, clique) in cliques {
        let mut members: Vec<i32> = clique.iter().copied().collect();
        members.sort_unstable();
        let idx_of: HashMap<i32, usize> =
            members.iter().enumerate().map(|(i, &v)| (v, i)).collect();
        let mut uf = UnionFind::new(members.len());

        // original graph edges among clique members
        for (li, &u) in members.iter().enumerate() {
            if let Some(neigh) = adj.get(&u) {
                for &w in neigh {
                    if let Some(&lj) = idx_of.get(&w) {
                        uf.union(li, lj);
                    }
                }
            }
        }
        // supplemented weak edges in the (0.1, cutoff] band
        for a in 0..members.len() {
            for b in (a + 1)..members.len() {
                let w = m[[members[a] as usize, members[b] as usize]];
                if w > 0.1 && w <= cutoff {
                    uf.union(a, b);
                }
            }
        }

        let mut root_hap: HashMap<usize, i32> = HashMap::new();
        for (li, &v) in members.iter().enumerate() {
            let r = uf.find(li);
            let hid = *root_hap.entry(r).or_insert_with(|| {
                let h = haplotype_idx;
                haplotype_idx += 1;
                h
            });
            vertex_hap.insert(v, hid);
        }
    }

    vertex_hap
}

/// CSR-native version of [`find_components_inside_cliques`].
fn find_components_inside_cliques_sparse(
    cliques: &[(Round, HashSet<i32>)],
    input: &SparsePhasingInput,
) -> HashMap<i32, i32> {
    let weights = &input.weights;
    let cutoff = input.edge_weight_cutoff;

    let member_lists: Vec<Vec<i32>> = cliques
        .iter()
        .map(|(_, clique)| {
            let mut members: Vec<i32> = clique.iter().copied().collect();
            members.sort_unstable();
            members.dedup();
            members
        })
        .collect();

    let mut clique_of = vec![usize::MAX; weights.size];
    let mut local_of = vec![usize::MAX; weights.size];
    for (clique_id, members) in member_lists.iter().enumerate() {
        for (local, &node) in members.iter().enumerate() {
            let node = node as usize;
            if node >= weights.size {
                continue;
            }
            if clique_of[node] != usize::MAX {
                return find_components_inside_cliques_sparse_legacy(cliques, input);
            }
            clique_of[node] = clique_id;
            local_of[node] = local;
        }
    }

    let mut finders: Vec<UnionFind> = member_lists
        .iter()
        .map(|members| UnionFind::new(members.len()))
        .collect();

    for &(u, v) in &input.edges {
        let (left, right) = (u as usize, v as usize);
        if left >= weights.size || right >= weights.size {
            continue;
        }
        let clique_id = clique_of[left];
        if clique_id != usize::MAX && clique_id == clique_of[right] {
            finders[clique_id].union(local_of[left], local_of[right]);
        }
    }

    for members in &member_lists {
        for &left_i32 in members {
            let left = left_i32 as usize;
            if left >= weights.size {
                continue;
            }
            let clique_id = clique_of[left];
            let (row_data, row_cols) = weights.row(left);
            for (&weight, &right_i32) in row_data.iter().zip(row_cols.iter()) {
                if right_i32 > left_i32
                    && weight > 0.1
                    && weight <= cutoff
                    && clique_of[right_i32 as usize] == clique_id
                {
                    finders[clique_id].union(local_of[left], local_of[right_i32 as usize]);
                }
            }
        }
    }

    let mut vertex_hap: HashMap<i32, i32> = HashMap::new();
    let mut haplotype_idx = 0i32;
    for (clique_id, members) in member_lists.iter().enumerate() {
        let finder = &mut finders[clique_id];
        let mut root_hap: HashMap<usize, i32> = HashMap::new();
        for (li, &v) in members.iter().enumerate() {
            let r = finder.find(li);
            let hid = *root_hap.entry(r).or_insert_with(|| {
                let h = haplotype_idx;
                haplotype_idx += 1;
                h
            });
            vertex_hap.insert(v, hid);
        }
    }

    vertex_hap
}

fn find_components_inside_cliques_sparse_legacy(
    cliques: &[(Round, HashSet<i32>)],
    input: &SparsePhasingInput,
) -> HashMap<i32, i32> {
    let weights = &input.weights;
    let cutoff = input.edge_weight_cutoff;

    let mut adj: HashMap<i32, HashSet<i32>> = HashMap::new();
    for &(u, v) in &input.edges {
        adj.entry(u).or_default().insert(v);
        adj.entry(v).or_default().insert(u);
    }

    let mut vertex_hap: HashMap<i32, i32> = HashMap::new();
    let mut haplotype_idx = 0i32;
    for (_round, clique) in cliques {
        let mut members: Vec<i32> = clique.iter().copied().collect();
        members.sort_unstable();
        let idx_of: HashMap<i32, usize> =
            members.iter().enumerate().map(|(i, &v)| (v, i)).collect();
        let mut uf = UnionFind::new(members.len());

        for (li, &u) in members.iter().enumerate() {
            if let Some(neigh) = adj.get(&u) {
                for &w in neigh {
                    if let Some(&lj) = idx_of.get(&w) {
                        uf.union(li, lj);
                    }
                }
            }
        }

        for a in 0..members.len() {
            for b in (a + 1)..members.len() {
                let w = weights.get(members[a] as usize, members[b] as usize);
                if w > 0.1 && w <= cutoff {
                    uf.union(a, b);
                }
            }
        }

        let mut root_hap: HashMap<usize, i32> = HashMap::new();
        for (li, &v) in members.iter().enumerate() {
            let r = uf.find(li);
            let hid = *root_hap.entry(r).or_insert_with(|| {
                let h = haplotype_idx;
                haplotype_idx += 1;
                h
            });
            vertex_hap.insert(v, hid);
        }
    }

    vertex_hap
}

/// Run the full phasing: returns `vertex -> hap_id` (the `qname_hap_info` map, keyed by vertex).
pub fn phase(input: &PhasingInput) -> HashMap<i32, i32> {
    let cliques = find_cliques_in_components(input);
    find_components_inside_cliques(&cliques, input)
}

/// CSR-native phasing prototype. Intended to prove parity before replacing the
/// production dense `Array2<f32>` graph build.
pub fn phase_sparse(input: &SparsePhasingInput) -> HashMap<i32, i32> {
    phase_sparse_with_threads(input, 1)
}

/// Sparse phasing with an explicit GCE worker budget.
pub fn phase_sparse_with_threads(input: &SparsePhasingInput, threads: usize) -> HashMap<i32, i32> {
    let cliques = find_cliques_in_components_sparse(input, threads.max(1));
    find_components_inside_cliques_sparse(&cliques, input)
}

/// Collapse a `vertex -> hap_id` map into a qname **partition** (set of qname sets), for
/// relabeling-invariant comparison against the Python `hap_qname_info`.
pub fn qname_partition(
    vertex_hap: &HashMap<i32, i32>,
    vertex_qname: &[String],
) -> HashSet<Vec<String>> {
    let mut hap_qnames: HashMap<i32, HashSet<String>> = HashMap::new();
    for (&v, &h) in vertex_hap {
        if let Some(qname) = vertex_qname.get(v as usize) {
            hap_qnames.entry(h).or_default().insert(qname.clone());
        }
    }
    hap_qnames
        .into_values()
        .map(|s| {
            let mut v: Vec<String> = s.into_iter().collect();
            v.sort();
            v
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    type ReadEvidence = (
        Vec<Vec<String>>,
        HashMap<String, Vec<i16>>,
        HashMap<String, Vec<f32>>,
    );

    fn empty_reads() -> ReadEvidence {
        (Vec::new(), HashMap::new(), HashMap::new())
    }

    fn vertex_partition(vertex_hap: &HashMap<i32, i32>) -> HashSet<Vec<i32>> {
        let mut groups: HashMap<i32, Vec<i32>> = HashMap::new();
        for (&vertex, &hap) in vertex_hap {
            groups.entry(hap).or_default().push(vertex);
        }
        groups
            .into_values()
            .map(|mut group| {
                group.sort_unstable();
                group
            })
            .collect()
    }

    fn dense_sparse_pair() -> (PhasingInput, SparsePhasingInput) {
        let n = 14usize;
        let mut m = Array2::<f32>::zeros((n, n));
        let mut entries = Vec::new();
        let mut edges = Vec::new();

        for i in 0..n {
            m[[i, i]] = 1.0;
            entries.push((i, i, 1.0));
        }

        for a in 0..n {
            for b in (a + 1)..n {
                let same = (a < 7) == (b < 7);
                let w = if same { 0.8 } else { -1.0 };
                m[[a, b]] = w;
                m[[b, a]] = w;
                entries.push((a, b, w));
                entries.push((b, a, w));
                if w > 0.0 {
                    edges.push((a as i32, b as i32));
                }
            }
        }

        let (node_read_ids, read_hap, read_err) = empty_reads();
        let dense = PhasingInput {
            weight_matrix: m,
            edges: edges.clone(),
            edge_weight_cutoff: 0.301,
            node_read_ids: node_read_ids.clone(),
            read_hap: read_hap.clone(),
            read_err: read_err.clone(),
        };
        let sparse = SparsePhasingInput {
            weights: Csr::from_entries(n, entries),
            edges,
            edge_weight_cutoff: 0.301,
            node_read_ids,
            read_hap,
            read_err,
        };
        (dense, sparse)
    }

    #[test]
    fn connected_components_splits_two_groups() {
        // 0-1-2 connected; 3-4 connected; 5 isolated.
        let edges = vec![(0, 1), (1, 2), (3, 4)];
        let comp = connected_components(6, &edges);
        assert_eq!(comp[0], comp[1]);
        assert_eq!(comp[1], comp[2]);
        assert_eq!(comp[3], comp[4]);
        assert_ne!(comp[0], comp[3]);
        assert_ne!(comp[5], comp[0]);
        assert_ne!(comp[5], comp[3]);
        // ids assigned in first-appearance order
        assert_eq!(comp[0], 0);
        assert_eq!(comp[3], 1);
        assert_eq!(comp[5], 2);
    }

    #[test]
    fn small_row_mask_detects_no_strong_edge() {
        let m = array![[0.0, 0.8, 0.0], [0.8, 0.0, 0.05], [0.0, 0.05, 0.0]];
        // row 2 has max 0.05 (<= 0.1) -> small; rows 0,1 have 0.8 -> big
        assert_eq!(small_row_mask(&m), vec![false, false, true]);
    }

    #[test]
    fn two_haplotypes_partition_by_incompatibility() {
        // Two cohesive groups of 6 read-pairs each (cliques must exceed size 5 to survive round 1).
        let n = 12usize;
        let mut m = Array2::<f32>::zeros((n, n));
        let mut edges = Vec::new();
        for a in 0..n {
            for b in (a + 1)..n {
                let same = (a < 6) == (b < 6);
                let w = if same { 0.8 } else { -1.0 };
                m[[a, b]] = w;
                m[[b, a]] = w;
                if w > 0.0 {
                    edges.push((a as i32, b as i32));
                }
            }
        }
        let (nri, hap, err) = empty_reads();
        let input = PhasingInput {
            weight_matrix: m,
            edges,
            edge_weight_cutoff: 0.201,
            node_read_ids: nri,
            read_hap: hap,
            read_err: err,
        };
        let vh = phase(&input);
        // every vertex assigned
        assert_eq!(vh.len(), n);
        // exactly two haplotypes, splitting {0..5} and {6..11}
        let group_a: HashSet<i32> = (0..6).map(|v| vh[&v]).collect();
        let group_b: HashSet<i32> = (6..12).map(|v| vh[&v]).collect();
        assert_eq!(group_a.len(), 1, "first 6 vertices share one haplotype");
        assert_eq!(group_b.len(), 1, "last 6 vertices share one haplotype");
        assert_ne!(
            group_a.iter().next().unwrap(),
            group_b.iter().next().unwrap(),
            "the two groups are different haplotypes"
        );
    }

    #[test]
    fn sparse_phase_matches_dense_phase_partition() {
        let (dense, sparse) = dense_sparse_pair();
        let dense_hap = phase(&dense);
        let sparse_hap = phase_sparse(&sparse);
        assert_eq!(vertex_partition(&sparse_hap), vertex_partition(&dense_hap));
    }

    #[test]
    fn sparse_phase_is_thread_budget_invariant() {
        let (_, sparse) = dense_sparse_pair();
        let one = phase_sparse_with_threads(&sparse, 1);
        let two = phase_sparse_with_threads(&sparse, 2);
        let four = phase_sparse_with_threads(&sparse, 4);
        assert_eq!(vertex_partition(&two), vertex_partition(&one));
        assert_eq!(vertex_partition(&four), vertex_partition(&one));
    }

    #[test]
    fn sparse_final_split_matches_dense_weak_edge_supplement() {
        let mut m = Array2::<f32>::zeros((3, 3));
        for i in 0..3 {
            m[[i, i]] = 1.0;
        }
        m[[0, 1]] = 0.8;
        m[[1, 0]] = 0.8;
        m[[1, 2]] = 0.2;
        m[[2, 1]] = 0.2;

        let cliques = vec![(Round::Main, HashSet::from([0, 1, 2]))];
        let (node_read_ids, read_hap, read_err) = empty_reads();
        let dense = PhasingInput {
            weight_matrix: m,
            edges: vec![(0, 1)], // weak 1-2 is intentionally only in the matrix
            edge_weight_cutoff: 0.301,
            node_read_ids: node_read_ids.clone(),
            read_hap: read_hap.clone(),
            read_err: read_err.clone(),
        };
        let sparse = SparsePhasingInput {
            weights: Csr::from_entries(
                3,
                [
                    (0, 0, 1.0),
                    (1, 1, 1.0),
                    (2, 2, 1.0),
                    (0, 1, 0.8),
                    (1, 0, 0.8),
                    (1, 2, 0.2),
                    (2, 1, 0.2),
                ],
            ),
            edges: vec![(0, 1)],
            edge_weight_cutoff: 0.301,
            node_read_ids,
            read_hap,
            read_err,
        };

        let dense_hap = find_components_inside_cliques(&cliques, &dense);
        let sparse_hap = find_components_inside_cliques_sparse(&cliques, &sparse);
        assert_eq!(vertex_partition(&sparse_hap), vertex_partition(&dense_hap));
        assert_eq!(
            vertex_partition(&sparse_hap),
            HashSet::from([vec![0, 1, 2]])
        );
    }

    #[test]
    fn qname_partition_is_relabeling_invariant() {
        let mut vh = HashMap::new();
        vh.insert(0, 7);
        vh.insert(1, 7);
        vh.insert(2, 3);
        let names = vec!["qA".to_string(), "qB".to_string(), "qC".to_string()];
        let part = qname_partition(&vh, &names);
        let mut want = HashSet::new();
        want.insert(vec!["qA".to_string(), "qB".to_string()]);
        want.insert(vec!["qC".to_string()]);
        assert_eq!(part, want);
    }
}
