//! Greedy-Clique-Expansion core, ported from `fp_control/gce_algorithm.py`.
//!
//! `heuristic_largest_clique` finds one clique by seed-and-expand on a CSR weight
//! sub-matrix; `gce_algorithm` repeatedly carves cliques out of a component until no
//! vertices remain. The port is 1:1 with the numba version, including:
//!   * `>=` tie-breaking in the seed/argmax (ties resolve to the highest index),
//!   * the 9-member lookback that lets an earlier clique member contribute a stronger
//!     extension edge, gated by a `1e-4` "the candidate doesn't prefer someone else" test,
//!   * the cutoff backtrack that scans all prior members when the greedy step stalls.

use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashSet};
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};

use ndarray::Array2;

use crate::kernels::{
    and_masks, efficient_mask, efficient_row_max, max_idx_mem, reverse_boolean_mask,
    row_wise_max_with_mask, Csr,
};

const REFINE_SWEEP_CEILING: usize = 1000;
const REFINE_EPS: f32 = 1e-6;

struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    fn find(&mut self, x: usize) -> usize {
        let mut root = x;
        while self.parent[root] != root {
            root = self.parent[root];
        }
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
            Ordering::Less => self.parent[ra] = rb,
            Ordering::Greater => self.parent[rb] = ra,
            Ordering::Equal => {
                self.parent[rb] = ra;
                self.rank[ra] += 1;
            }
        }
    }
}

fn par_chunk_reduce<T, M, R>(n: usize, threads: usize, map: M, mut reduce: R, init: T) -> T
where
    T: Send,
    M: Fn(std::ops::Range<usize>) -> T + Sync,
    R: FnMut(T, T) -> T,
{
    let threads = threads.clamp(1, n.max(1));
    if threads == 1 || n == 0 {
        return reduce(init, map(0..n));
    }
    let chunk = n.div_ceil(threads);
    let results = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|t| {
                let start = t * chunk;
                let end = ((t + 1) * chunk).min(n);
                let map = &map;
                scope.spawn(move || map(start..end))
            })
            .collect();
        handles
            .into_iter()
            .map(|handle| handle.join().expect("GCE worker thread panicked"))
            .collect::<Vec<T>>()
    });
    let mut acc = init;
    for value in results {
        acc = reduce(acc, value);
    }
    acc
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct HeapWeight(f32);

impl Eq for HeapWeight {}

impl PartialOrd for HeapWeight {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for HeapWeight {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

struct FastGraph<'a> {
    n: usize,
    indptr: &'a [i32],
    cols: &'a [i32],
    vals: &'a [f32],
    safe_mode: bool,
    t_indptr: Vec<usize>,
    t_rows: Vec<u32>,
}

impl<'a> FastGraph<'a> {
    fn from_csr(csr: &'a Csr) -> Self {
        assert!(
            u32::try_from(csr.size).is_ok() && csr.size < u32::MAX as usize,
            "fast GCE supports fewer than u32::MAX nodes"
        );

        let mut safe_mode = false;
        for &value in &csr.data {
            assert!(
                value.is_finite(),
                "fast GCE requires finite edge weights, found {value}"
            );
            if value < 0.0 && value != -1.0 {
                safe_mode = true;
            }
        }

        let (t_indptr, t_rows) = if safe_mode {
            build_structural_transpose(csr.size, &csr.indptr, &csr.indices)
        } else {
            (Vec::new(), Vec::new())
        };

        Self {
            n: csr.size,
            indptr: &csr.indptr,
            cols: &csr.indices,
            vals: &csr.data,
            safe_mode,
            t_indptr,
            t_rows,
        }
    }

    #[inline]
    fn row(&self, row: usize) -> (&[f32], &[i32]) {
        let start = self.indptr[row] as usize;
        let end = self.indptr[row + 1] as usize;
        (&self.vals[start..end], &self.cols[start..end])
    }
}

/// Reusable GCE engine state for one CSR matrix.
///
/// The production phasing pipeline calls GCE for many component/round subsets
/// of the same graph. This context pays the finite-weight/safe-mode scan once,
/// then partitions every subset as a masked view of the original CSR.
pub struct GceContext<'a> {
    graph: FastGraph<'a>,
    threads: usize,
}

impl<'a> GceContext<'a> {
    pub fn new(csr: &'a Csr) -> Self {
        Self::with_threads(csr, 1)
    }

    pub fn with_threads(csr: &'a Csr, threads: usize) -> Self {
        Self {
            graph: FastGraph::from_csr(csr),
            threads: threads.max(1),
        }
    }

    fn partition(&self, cutoff: f32) -> Vec<Vec<usize>> {
        let n = self.graph.n;
        partition_on_mask(&self.graph, vec![true; n], n, cutoff, self.threads)
    }

    pub fn partition_subset(&self, nodes: &[usize], cutoff: f32) -> Vec<Vec<usize>> {
        if nodes.is_empty() {
            return Vec::new();
        }
        assert!(
            nodes.windows(2).all(|pair| pair[0] < pair[1]),
            "subset nodes must be strictly ascending"
        );
        assert!(
            *nodes.last().expect("nonempty node subset") < self.graph.n,
            "selected node index out of bounds"
        );
        let mut active = vec![false; self.graph.n];
        for &node in nodes {
            active[node] = true;
        }
        partition_on_mask(&self.graph, active, nodes.len(), cutoff, self.threads)
    }
}

fn build_structural_transpose(n: usize, indptr: &[i32], cols: &[i32]) -> (Vec<usize>, Vec<u32>) {
    let mut t_indptr = vec![0usize; n + 1];
    for &col in cols {
        t_indptr[col as usize + 1] += 1;
    }
    for i in 0..n {
        t_indptr[i + 1] += t_indptr[i];
    }
    let mut cursor = t_indptr[..n].to_vec();
    let mut t_rows = vec![0u32; cols.len()];
    for row in 0..n {
        #[allow(clippy::cast_possible_truncation)]
        let row_u32 = row as u32;
        for &col in &cols[indptr[row] as usize..indptr[row + 1] as usize] {
            let col = col as usize;
            t_rows[cursor[col]] = row_u32;
            cursor[col] += 1;
        }
    }
    (t_indptr, t_rows)
}

pub(crate) fn structural_transpose(csr: &Csr) -> (Vec<usize>, Vec<u32>) {
    build_structural_transpose(csr.size, &csr.indptr, &csr.indices)
}

#[inline]
fn effective_best(graph: &FastGraph<'_>, node: usize, active: &[bool]) -> f32 {
    let (vals, cols) = graph.row(node);
    let mut best = f32::NEG_INFINITY;
    let mut found = false;
    for (&value, &col) in vals.iter().zip(cols.iter()) {
        let col = col as usize;
        if col != node && value != -1.0 && active[col] {
            if value > best {
                best = value;
            }
            found = true;
        }
    }
    let effective = if found { best } else { 0.0 };
    if effective == 0.0 {
        0.0
    } else {
        effective
    }
}

#[inline]
fn row_best_candidate(vals: &[f32], cols: &[i32], cand: &[bool]) -> (usize, f32) {
    let mut best_col = usize::MAX;
    let mut best_value = -2.0f32;
    let mut found = false;
    for (&value, &col) in vals.iter().zip(cols.iter()) {
        let col = col as usize;
        if value != -1.0 && cand[col] && (!found || value >= best_value) {
            best_col = col;
            best_value = value;
            found = true;
        }
    }
    (best_col, best_value)
}

#[inline]
fn apply_member_to_candidates(
    graph: &FastGraph<'_>,
    member: usize,
    cand: &mut [bool],
    cand_count: &mut usize,
) {
    if cand[member] {
        cand[member] = false;
        *cand_count -= 1;
    }
    let (vals, cols) = graph.row(member);
    for (&value, &col) in vals.iter().zip(cols.iter()) {
        if value == -1.0 {
            let col = col as usize;
            if cand[col] {
                cand[col] = false;
                *cand_count -= 1;
            }
        }
    }
}

fn grow_clique(
    graph: &FastGraph<'_>,
    seed: usize,
    cutoff: f32,
    active: &[bool],
    active_count: usize,
    cand: &mut [bool],
) -> Vec<usize> {
    cand.copy_from_slice(active);
    let mut cand_count = active_count;
    apply_member_to_candidates(graph, seed, cand, &mut cand_count);

    let mut select_indices = vec![seed];
    let mut member_count = 1usize;
    let mut max_row_ind = seed;

    while cand_count > 0 {
        let (mut next_max_ind, mut next_max_value) = {
            let (vals, cols) = graph.row(max_row_ind);
            row_best_candidate(vals, cols, cand)
        };

        let lookback_start = member_count.saturating_sub(9);
        for &member in &select_indices[lookback_start..member_count] {
            if member == max_row_ind {
                continue;
            }
            let (trial_max_ind, trial_max_value) = {
                let (vals, cols) = graph.row(member);
                row_best_candidate(vals, cols, cand)
            };
            if trial_max_value <= cutoff || trial_max_ind == usize::MAX {
                continue;
            }
            let (_best_neighbor, best_value) = {
                let (vals, cols) = graph.row(trial_max_ind);
                row_best_candidate(vals, cols, cand)
            };
            if best_value - trial_max_value > 1e-4 {
                continue;
            }
            if trial_max_value > next_max_value {
                next_max_ind = trial_max_ind;
                next_max_value = trial_max_value;
            }
        }

        if next_max_value <= cutoff {
            let mut trial = 0usize;
            while next_max_value <= cutoff && trial < member_count {
                let member = select_indices[trial];
                let (vals, cols) = graph.row(member);
                let (candidate, value) = row_best_candidate(vals, cols, cand);
                next_max_ind = candidate;
                next_max_value = value;
                trial += 1;
            }
            if next_max_value <= cutoff && trial >= member_count {
                break;
            }
            max_row_ind = next_max_ind;
        } else {
            max_row_ind = next_max_ind;
        }

        if max_row_ind == usize::MAX {
            break;
        }

        select_indices.push(max_row_ind);
        member_count += 1;
        if member_count >= active_count {
            break;
        }

        apply_member_to_candidates(graph, max_row_ind, cand, &mut cand_count);
    }

    select_indices
}

fn push_affected_rows(
    graph: &FastGraph<'_>,
    removed: &[usize],
    active: &[bool],
    heap: &mut BinaryHeap<(HeapWeight, u32)>,
    touched_epoch: &mut [u32],
    epoch: u32,
) {
    for &member in removed {
        let start = graph.t_indptr[member];
        let end = graph.t_indptr[member + 1];
        for &row in &graph.t_rows[start..end] {
            let row = row as usize;
            if active[row] && touched_epoch[row] != epoch {
                touched_epoch[row] = epoch;
                let value = effective_best(graph, row, active);
                #[allow(clippy::cast_possible_truncation)]
                heap.push((HeapWeight(value), row as u32));
            }
        }
    }
}

fn partition_on_mask_sequential(
    graph: &FastGraph<'_>,
    mut active: Vec<bool>,
    mut active_count: usize,
    cutoff: f32,
    threads: usize,
) -> Vec<Vec<usize>> {
    let n = graph.n;
    if n == 0 || active_count == 0 {
        return Vec::new();
    }

    let mut cand = vec![false; n];
    let seed_threads = if active_count < 20_000 { 1 } else { threads };
    let entries = par_chunk_reduce(
        n,
        seed_threads,
        |range| {
            let mut chunk = Vec::new();
            for node in range {
                if !active[node] {
                    continue;
                }
                let value = effective_best(graph, node, &active);
                #[allow(clippy::cast_possible_truncation)]
                chunk.push((HeapWeight(value), node as u32));
            }
            chunk
        },
        |mut acc, mut chunk| {
            acc.append(&mut chunk);
            acc
        },
        Vec::with_capacity(active_count),
    );
    let mut heap: BinaryHeap<(HeapWeight, u32)> = BinaryHeap::from(entries);

    let mut touched_epoch = vec![u32::MAX; if graph.safe_mode { n } else { 0 }];
    let mut epoch = 0u32;
    let mut cliques = Vec::new();

    while active_count > 0 {
        let (HeapWeight(cached), node_u32) = heap
            .pop()
            .expect("every active node keeps at least one heap entry");
        let node = node_u32 as usize;
        if !active[node] {
            continue;
        }
        let current = effective_best(graph, node, &active);
        if current != cached {
            heap.push((HeapWeight(current), node_u32));
            continue;
        }

        if cutoff >= 0.0 && current <= cutoff {
            active[node] = false;
            active_count -= 1;
            if graph.safe_mode {
                push_affected_rows(
                    graph,
                    &[node],
                    &active,
                    &mut heap,
                    &mut touched_epoch,
                    epoch,
                );
                epoch += 1;
            }
            cliques.push(vec![node]);
            continue;
        }

        let mut clique = grow_clique(graph, node, cutoff, &active, active_count, &mut cand);
        for &member in &clique {
            active[member] = false;
        }
        active_count -= clique.len();
        if graph.safe_mode {
            push_affected_rows(
                graph,
                &clique,
                &active,
                &mut heap,
                &mut touched_epoch,
                epoch,
            );
            epoch += 1;
        }
        clique.sort_unstable();
        cliques.push(clique);
    }

    cliques
}

type CliqueStream = Vec<(HeapWeight, u32, Vec<usize>)>;

fn run_component_phase(
    graph: &FastGraph<'_>,
    comp: &[usize],
    cutoff: f32,
    local_active: &mut [bool],
    cand: &mut [bool],
) -> CliqueStream {
    for &node in comp {
        local_active[node] = true;
    }
    let mut local_count = comp.len();

    let mut entries = Vec::with_capacity(comp.len());
    for &node in comp {
        let value = effective_best(graph, node, local_active);
        #[allow(clippy::cast_possible_truncation)]
        entries.push((HeapWeight(value), node as u32));
    }
    let mut heap = BinaryHeap::from(entries);

    let mut stream = Vec::new();
    while local_count > 0 {
        let Some((HeapWeight(cached), node_u32)) = heap.pop() else {
            break;
        };
        let node = node_u32 as usize;
        if !local_active[node] {
            continue;
        }
        let current = effective_best(graph, node, local_active);
        if current != cached {
            heap.push((HeapWeight(current), node_u32));
            continue;
        }
        if current <= cutoff {
            break;
        }
        let mut clique = grow_clique(graph, node, cutoff, local_active, local_count, cand);
        for &member in &clique {
            local_active[member] = false;
        }
        local_count -= clique.len();
        clique.sort_unstable();
        stream.push((HeapWeight(current), node_u32, clique));
    }

    for &node in comp {
        local_active[node] = false;
    }
    stream
}

fn partition_on_mask(
    graph: &FastGraph<'_>,
    mut active: Vec<bool>,
    mut active_count: usize,
    cutoff: f32,
    threads: usize,
) -> Vec<Vec<usize>> {
    let n = graph.n;
    if n == 0 || active_count == 0 {
        return Vec::new();
    }
    if cutoff < 0.0 {
        return partition_on_mask_sequential(graph, active, active_count, cutoff, threads);
    }

    let mut finder = UnionFind::new(n);
    for row in 0..n {
        if !active[row] {
            continue;
        }
        let (vals, cols) = graph.row(row);
        for (&value, &col) in vals.iter().zip(cols.iter()) {
            let col = col as usize;
            if value > cutoff && col != row && active[col] {
                finder.union(row, col);
            }
        }
    }

    let mut comp_size = vec![0u32; n];
    for (node, &is_active) in active.iter().enumerate() {
        if is_active {
            comp_size[finder.find(node)] += 1;
        }
    }
    let mut comp_index = vec![usize::MAX; n];
    let mut comps: Vec<Vec<usize>> = Vec::new();
    for (node, &is_active) in active.iter().enumerate() {
        if !is_active {
            continue;
        }
        let root = finder.find(node);
        if comp_size[root] < 2 {
            continue;
        }
        let idx = if comp_index[root] == usize::MAX {
            comp_index[root] = comps.len();
            comps.push(Vec::with_capacity(comp_size[root] as usize));
            comps.len() - 1
        } else {
            comp_index[root]
        };
        comps[idx].push(node);
    }

    if comps.is_empty() {
        return partition_on_mask_sequential(graph, active, active_count, cutoff, threads);
    }

    let mut order: Vec<usize> = (0..comps.len()).collect();
    order.sort_unstable_by_key(|&idx| std::cmp::Reverse(comps[idx].len()));
    let workers = threads.max(1).min(comps.len());

    let mut streams: Vec<CliqueStream> = Vec::with_capacity(comps.len());
    if workers <= 1 {
        let mut local_active = vec![false; n];
        let mut cand = vec![false; n];
        streams.resize_with(comps.len(), Vec::new);
        for &idx in &order {
            streams[idx] =
                run_component_phase(graph, &comps[idx], cutoff, &mut local_active, &mut cand);
        }
    } else {
        let next = AtomicUsize::new(0);
        let mut collected: Vec<Vec<(usize, CliqueStream)>> = std::thread::scope(|scope| {
            let handles: Vec<_> = (0..workers)
                .map(|_| {
                    let order = &order;
                    let comps = &comps;
                    let next = &next;
                    scope.spawn(move || {
                        let mut local_active = vec![false; n];
                        let mut cand = vec![false; n];
                        let mut produced = Vec::new();
                        loop {
                            let slot = next.fetch_add(1, AtomicOrdering::Relaxed);
                            if slot >= order.len() {
                                break;
                            }
                            let idx = order[slot];
                            let stream = run_component_phase(
                                graph,
                                &comps[idx],
                                cutoff,
                                &mut local_active,
                                &mut cand,
                            );
                            produced.push((idx, stream));
                        }
                        produced
                    })
                })
                .collect();
            handles
                .into_iter()
                .map(|handle| handle.join().expect("GCE component worker panicked"))
                .collect()
        });
        streams.resize_with(comps.len(), Vec::new);
        for batch in &mut collected {
            for (idx, stream) in batch.drain(..) {
                streams[idx] = stream;
            }
        }
    }

    let mut cliques = Vec::new();
    let mut merge: BinaryHeap<(HeapWeight, u32, usize)> = BinaryHeap::new();
    let mut cursor = vec![0usize; streams.len()];
    for (idx, stream) in streams.iter().enumerate() {
        if let Some(&(key, seed, _)) = stream.first() {
            merge.push((key, seed, idx));
        }
    }
    while let Some((_, _, idx)) = merge.pop() {
        let position = cursor[idx];
        cursor[idx] += 1;
        let clique = std::mem::take(&mut streams[idx][position].2);
        for &member in &clique {
            active[member] = false;
        }
        active_count -= clique.len();
        cliques.push(clique);
        if let Some(&(key, seed, _)) = streams[idx].get(cursor[idx]) {
            merge.push((key, seed, idx));
        }
    }

    let mut tail = partition_on_mask_sequential(graph, active, active_count, cutoff, threads);
    cliques.append(&mut tail);
    cliques
}

fn partition_disjoint_cliques(csr: &Csr, cutoff: f32) -> Vec<Vec<usize>> {
    GceContext::new(csr).partition(cutoff)
}

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

#[cfg(test)]
fn partition_disjoint_cliques_legacy(csr: &Csr, cutoff: f32) -> Vec<Vec<usize>> {
    let mut cliques = Vec::new();
    let mut selected_indices: Vec<usize> = (0..csr.size).collect();
    let mut current_csr = csr.clone();
    let mut size = current_csr.size;

    loop {
        if size == 0 {
            break;
        }

        let index_mask = vec![true; size];
        let (clique_local, drop_mask) = heuristic_largest_clique(&current_csr, &index_mask, cutoff);
        let mut clique: Vec<usize> = clique_local
            .iter()
            .map(|&local_idx| selected_indices[local_idx as usize])
            .collect();
        clique.sort_unstable();

        let remain_count = crate::kernels::count_true(&drop_mask);
        if !clique.is_empty() {
            cliques.push(clique);
        } else if remain_count > 0 {
            break;
        }

        if remain_count == 0 {
            break;
        }

        selected_indices = drop_mask
            .iter()
            .enumerate()
            .filter_map(|(idx, &keep)| keep.then_some(selected_indices[idx]))
            .collect();
        current_csr = current_csr.select(&drop_mask);
        size = current_csr.size;
    }

    cliques
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct RefineStats {
    rounds_run: usize,
    total_moves: usize,
    converged: bool,
}

struct CliqueScore {
    max_w: f32,
    sum_w: f32,
    blocked: bool,
}

fn refine_partition(
    csr: &Csr,
    cliques: &[Vec<usize>],
    cutoff: f32,
    max_rounds: usize,
) -> (Vec<Vec<usize>>, RefineStats) {
    let nodes: Vec<usize> = (0..csr.size).collect();
    refine_partition_in_subset(csr, &nodes, cliques, cutoff, max_rounds)
}

fn refine_partition_in_subset(
    csr: &Csr,
    nodes: &[usize],
    cliques: &[Vec<usize>],
    cutoff: f32,
    max_rounds: usize,
) -> (Vec<Vec<usize>>, RefineStats) {
    refine_partition_in_subset_with(csr, nodes, cliques, cutoff, max_rounds, None)
}

fn refine_partition_in_subset_with(
    csr: &Csr,
    nodes: &[usize],
    cliques: &[Vec<usize>],
    cutoff: f32,
    max_rounds: usize,
    shared_transpose: Option<&(Vec<usize>, Vec<u32>)>,
) -> (Vec<Vec<usize>>, RefineStats) {
    let n = csr.size;
    let clique_count = cliques.len();
    let mut assignment = vec![usize::MAX; n];
    let mut members = cliques.to_vec();
    for (clique_id, clique) in members.iter().enumerate() {
        for &node in clique {
            if node < n {
                assignment[node] = clique_id;
            }
        }
    }

    let mut stamp = vec![0u64; clique_count];
    let mut scores: Vec<CliqueScore> = (0..clique_count)
        .map(|_| CliqueScore {
            max_w: 0.0,
            sum_w: 0.0,
            blocked: false,
        })
        .collect();
    let mut touched = Vec::new();
    let mut epoch = 0u64;

    let mut total_moves = 0usize;
    let mut rounds_run = 0usize;
    let mut converged = false;

    let owned_transpose;
    let transpose = if let Some(shared) = shared_transpose {
        (max_rounds > 1).then_some(shared)
    } else {
        owned_transpose = (max_rounds > 1).then(|| structural_transpose(csr));
        owned_transpose.as_ref()
    };
    let mut dirty = vec![false; n];
    for &node in nodes {
        dirty[node] = true;
    }

    for _ in 0..max_rounds {
        let mut moves_this_round = 0usize;
        for &node in nodes {
            if !dirty[node] {
                continue;
            }
            let current = assignment[node];
            if current == usize::MAX {
                dirty[node] = false;
                continue;
            }

            epoch += 1;
            touched.clear();
            let (row_data, row_cols) = csr.row(node);
            for (&weight, &neighbor) in row_data.iter().zip(row_cols.iter()) {
                let neighbor = neighbor as usize;
                if neighbor == node {
                    continue;
                }
                let clique_id = assignment[neighbor];
                if clique_id == usize::MAX {
                    continue;
                }
                if stamp[clique_id] != epoch {
                    stamp[clique_id] = epoch;
                    scores[clique_id] = CliqueScore {
                        max_w: 0.0,
                        sum_w: 0.0,
                        blocked: false,
                    };
                    touched.push(clique_id);
                }
                if weight == -1.0 {
                    scores[clique_id].blocked = true;
                } else if weight > 0.0 {
                    let score = &mut scores[clique_id];
                    if weight > score.max_w {
                        score.max_w = weight;
                    }
                    score.sum_w += weight;
                }
            }

            let (current_max, current_sum) = if stamp[current] == epoch {
                (scores[current].max_w, scores[current].sum_w)
            } else {
                (0.0, 0.0)
            };

            let mut candidates: Vec<usize> = touched
                .iter()
                .copied()
                .filter(|&clique_id| {
                    clique_id != current
                        && !scores[clique_id].blocked
                        && scores[clique_id].max_w > cutoff
                })
                .collect();
            candidates.sort_unstable_by(|&a, &b| {
                scores[b]
                    .max_w
                    .total_cmp(&scores[a].max_w)
                    .then(scores[b].sum_w.total_cmp(&scores[a].sum_w))
                    .then(a.cmp(&b))
            });

            let mut moved = false;
            for &target in &candidates {
                if scores[target].max_w < current_max - REFINE_EPS {
                    break;
                }
                let improves_max = scores[target].max_w > current_max + REFINE_EPS;
                let ties_max = (scores[target].max_w - current_max).abs() <= REFINE_EPS;
                let improves_sum = scores[target].sum_w > current_sum + REFINE_EPS;
                if !(improves_max || (ties_max && improves_sum)) {
                    continue;
                }

                let reverse_blocked = members[target]
                    .iter()
                    .any(|&member| csr.get(member, node) == -1.0);
                if reverse_blocked {
                    continue;
                }

                let position = members[current]
                    .iter()
                    .position(|&member| member == node)
                    .expect("GCE refinement assignment and member lists stay in sync");
                members[current].swap_remove(position);
                members[target].push(node);
                assignment[node] = target;
                moves_this_round += 1;
                moved = true;
                break;
            }

            if moved {
                if let Some((t_indptr, t_rows)) = &transpose {
                    let (_, row_cols) = csr.row(node);
                    for &neighbor in row_cols {
                        dirty[neighbor as usize] = true;
                    }
                    for &row in &t_rows[t_indptr[node]..t_indptr[node + 1]] {
                        dirty[row as usize] = true;
                    }
                }
                dirty[node] = true;
            } else {
                dirty[node] = false;
            }
        }

        rounds_run += 1;
        total_moves += moves_this_round;
        if moves_this_round == 0 {
            converged = true;
            break;
        }
    }

    let refined: Vec<Vec<usize>> = members
        .into_iter()
        .filter(|clique| !clique.is_empty())
        .map(|mut clique| {
            clique.sort_unstable();
            clique
        })
        .collect();

    (
        refined,
        RefineStats {
            rounds_run,
            total_moves,
            converged,
        },
    )
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
    gce_algorithm_csr(selected_indices, Csr::from_dense(weight_matrix), cutoff)
}

/// CSR-native Greedy-Clique-Expansion entry point.
///
/// `selected_indices[k]` maps local CSR row `k` back to its index in the original
/// graph. This is behavior-equivalent to [`gce_algorithm`] but avoids the dense
/// matrix allocation and the dense-to-CSR conversion.
pub fn gce_algorithm_csr(
    selected_indices: &[i32],
    weight_csr: Csr,
    cutoff: f32,
) -> Vec<HashSet<i32>> {
    assert_eq!(
        selected_indices.len(),
        weight_csr.size,
        "selected_indices must match the sub-matrix size"
    );

    let cliques = partition_disjoint_cliques(&weight_csr, cutoff);
    let (cliques, refine_stats) =
        refine_partition(&weight_csr, &cliques, cutoff, REFINE_SWEEP_CEILING);
    if !refine_stats.converged {
        log::warn!(
            "GCE refinement was still moving nodes after {REFINE_SWEEP_CEILING} sweeps; using final sweep state"
        );
    }

    cliques
        .into_iter()
        .map(|clique| {
            clique
                .into_iter()
                .map(|local_idx| selected_indices[local_idx])
                .collect()
        })
        .collect()
}

/// CSR entry point for an ordered subset of the original graph.
///
/// `selected_indices` are original CSR row ids. The order defines compacted
/// tie-breaking, matching `weight_csr.select(mask)` when the ids are ascending.
pub fn gce_algorithm_csr_subset(
    selected_indices: &[i32],
    weight_csr: &Csr,
    cutoff: f32,
) -> Vec<HashSet<i32>> {
    let ctx = GceContext::new(weight_csr);
    gce_algorithm_csr_subset_with_context(selected_indices, weight_csr, &ctx, None, cutoff)
}

pub(crate) fn gce_algorithm_csr_subset_with_context(
    selected_indices: &[i32],
    weight_csr: &Csr,
    ctx: &GceContext<'_>,
    shared_transpose: Option<&(Vec<usize>, Vec<u32>)>,
    cutoff: f32,
) -> Vec<HashSet<i32>> {
    let nodes: Vec<usize> = selected_indices
        .iter()
        .map(|&idx| {
            usize::try_from(idx).expect("selected_indices must be non-negative CSR row ids")
        })
        .collect();
    let cliques = ctx.partition_subset(&nodes, cutoff);
    let (cliques, refine_stats) = refine_partition_in_subset_with(
        weight_csr,
        &nodes,
        &cliques,
        cutoff,
        REFINE_SWEEP_CEILING,
        shared_transpose,
    );
    if !refine_stats.converged {
        log::warn!(
            "GCE refinement was still moving nodes after {REFINE_SWEEP_CEILING} sweeps; using final sweep state"
        );
    }

    cliques
        .into_iter()
        .map(|clique| clique.into_iter().map(|node| node as i32).collect())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    fn canonical_local_cliques(cliques: &[Vec<usize>]) -> HashSet<Vec<usize>> {
        cliques
            .iter()
            .map(|clique| {
                let mut v = clique.clone();
                v.sort_unstable();
                v
            })
            .collect()
    }

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

    #[test]
    fn csr_entry_point_matches_dense_entry_point() {
        let n1 = -1.0f32;
        let m = array![
            [1.0, 0.8, 0.8, n1, n1, n1],
            [0.8, 1.0, 0.8, n1, n1, n1],
            [0.8, 0.8, 1.0, n1, n1, n1],
            [n1, n1, n1, 1.0, 0.8, 0.8],
            [n1, n1, n1, 0.8, 1.0, 0.8],
            [n1, n1, n1, 0.8, 0.8, 1.0]
        ];
        let sel: Vec<i32> = (0..6).collect();
        let dense = gce_algorithm(&sel, &m, 0.2);
        let sparse = gce_algorithm_csr(&sel, Csr::from_dense(&m), 0.2);
        assert_eq!(clique_set(&sparse), clique_set(&dense));
    }

    #[test]
    fn subset_entry_point_matches_compacted_csr() {
        let n1 = -1.0f32;
        let csr = Csr::from_entries(
            8,
            [
                (0, 0, 1.0),
                (1, 1, 1.0),
                (2, 2, 1.0),
                (3, 3, 1.0),
                (4, 4, 1.0),
                (5, 5, 1.0),
                (6, 6, 1.0),
                (7, 7, 1.0),
                (1, 2, 0.9),
                (2, 1, 0.9),
                (1, 4, 0.7),
                (4, 1, 0.7),
                (2, 4, 0.8),
                (4, 2, 0.8),
                (3, 5, 0.85),
                (5, 3, 0.85),
                (4, 5, n1),
                (5, 4, n1),
                (6, 7, 0.4),
                (7, 6, 0.4),
            ],
        );
        let selected = vec![1, 2, 3, 4, 5, 7];
        let mut mask = vec![false; csr.size];
        for &node in &selected {
            mask[node as usize] = true;
        }
        let compacted = csr.select(&mask);
        let compacted_cliques = gce_algorithm_csr(&selected, compacted, 0.301);
        let subset_cliques = gce_algorithm_csr_subset(&selected, &csr, 0.301);
        assert_eq!(clique_set(&subset_cliques), clique_set(&compacted_cliques));
    }

    #[test]
    fn fast_partition_matches_legacy_before_refinement() {
        let n1 = -1.0f32;
        let csr = Csr::from_entries(
            8,
            [
                (0, 0, 1.0),
                (1, 1, 1.0),
                (2, 2, 1.0),
                (3, 3, 1.0),
                (4, 4, 1.0),
                (5, 5, 1.0),
                (6, 6, 1.0),
                (7, 7, 1.0),
                (0, 1, 0.9),
                (1, 0, 0.9),
                (0, 2, 0.8),
                (2, 0, 0.8),
                (1, 2, 0.7),
                (2, 1, 0.7),
                (3, 4, 0.85),
                (4, 3, 0.85),
                (3, 5, 0.6),
                (5, 3, 0.6),
                (4, 5, 0.55),
                (5, 4, 0.55),
                (0, 3, n1),
                (3, 0, n1),
                (1, 4, n1),
                (4, 1, n1),
                (6, 7, 0.05),
                (7, 6, 0.05),
            ],
        );

        let fast = partition_disjoint_cliques(&csr, 0.301);
        let legacy = partition_disjoint_cliques_legacy(&csr, 0.301);
        assert_eq!(fast, legacy);
    }

    #[test]
    fn fast_partition_handles_negative_nonblocked_safe_mode() {
        let csr = Csr::from_entries(
            5,
            [
                (0, 0, 1.0),
                (1, 1, 1.0),
                (2, 2, 1.0),
                (3, 3, 1.0),
                (4, 4, 1.0),
                (0, 1, -0.2),
                (1, 0, -0.2),
                (0, 2, 0.8),
                (2, 0, 0.8),
                (1, 3, 0.8),
                (3, 1, 0.8),
                (2, 4, -1.0),
                (4, 2, -1.0),
            ],
        );

        let fast = partition_disjoint_cliques(&csr, 0.0);
        let legacy = partition_disjoint_cliques_legacy(&csr, 0.0);
        assert_eq!(
            canonical_local_cliques(&fast),
            canonical_local_cliques(&legacy)
        );
    }

    #[test]
    fn explicit_thread_budgets_produce_identical_partitions() {
        let mut entries = Vec::new();
        for component in 0..4 {
            let start = component * 3;
            for node in start..(start + 3) {
                entries.push((node, node, 1.0));
            }
            for left in start..(start + 3) {
                for right in start..(start + 3) {
                    if left != right {
                        entries.push((left, right, 0.8));
                    }
                }
            }
        }
        let csr = Csr::from_entries(12, entries);
        let one = GceContext::with_threads(&csr, 1).partition(0.3);
        let two = GceContext::with_threads(&csr, 2).partition(0.3);
        let four = GceContext::with_threads(&csr, 4).partition(0.3);

        assert_eq!(canonical_local_cliques(&two), canonical_local_cliques(&one));
        assert_eq!(
            canonical_local_cliques(&four),
            canonical_local_cliques(&one)
        );
    }

    #[test]
    fn refinement_moves_node_to_strongest_edge_clique() {
        let entries = [
            (0, 0, 1.0),
            (1, 1, 1.0),
            (2, 2, 1.0),
            (3, 3, 1.0),
            (4, 4, 1.0),
            (0, 1, 0.95),
            (1, 0, 0.95),
            (0, 4, 0.60),
            (4, 0, 0.60),
            (1, 4, 0.55),
            (4, 1, 0.55),
            (2, 3, 0.50),
            (3, 2, 0.50),
            (2, 4, 0.80),
            (4, 2, 0.80),
            (3, 4, 0.45),
            (4, 3, 0.45),
            (0, 2, -1.0),
            (2, 0, -1.0),
            (0, 3, -1.0),
            (3, 0, -1.0),
            (1, 2, -1.0),
            (2, 1, -1.0),
            (1, 3, -1.0),
            (3, 1, -1.0),
        ];
        let csr = Csr::from_entries(5, entries);
        let cutoff = 0.30;
        let cliques = partition_disjoint_cliques(&csr, cutoff);
        let (refined, stats) = refine_partition(&csr, &cliques, cutoff, 10);

        assert!(stats.total_moves >= 1);
        assert!(stats.converged);
        assert!(
            refined
                .iter()
                .any(|clique| clique.contains(&2) && clique.contains(&4)),
            "node 4 should join the clique containing node 2: {refined:?}"
        );
    }

    #[test]
    fn refinement_respects_reverse_blocked_pairs() {
        let csr = Csr::from_entries(
            5,
            [
                (0, 0, 1.0),
                (1, 1, 1.0),
                (2, 2, 1.0),
                (3, 3, 1.0),
                (4, 4, 1.0),
                (0, 1, 0.95),
                (1, 0, 0.95),
                (0, 4, 0.60),
                (4, 0, 0.60),
                (1, 4, 0.55),
                (4, 1, 0.55),
                (2, 3, 0.50),
                (3, 2, 0.50),
                (2, 4, 0.80),
                (4, 2, 0.80),
                (4, 3, 0.45),
                (3, 4, -1.0),
                (0, 2, -1.0),
                (2, 0, -1.0),
                (0, 3, -1.0),
                (3, 0, -1.0),
                (1, 2, -1.0),
                (2, 1, -1.0),
                (1, 3, -1.0),
                (3, 1, -1.0),
            ],
        );
        let cutoff = 0.30;
        let cliques = partition_disjoint_cliques(&csr, cutoff);
        let (refined, stats) = refine_partition(&csr, &cliques, cutoff, 10);

        assert_eq!(stats.total_moves, 0);
        assert!(
            !refined
                .iter()
                .any(|clique| clique.contains(&3) && clique.contains(&4)),
            "blocked pair should not be created: {refined:?}"
        );
    }
}
