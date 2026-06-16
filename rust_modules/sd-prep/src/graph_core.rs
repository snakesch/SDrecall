//! ★ THE graph-tool replacement — the deliverable of T8.
//!
//! Three pure-combinatorics units that replace graph-tool's `label_components`,
//! `shortest_path` and `sequential_vertex_coloring`. They are decoupled from
//! `petgraph`'s graph type wherever the algorithm is pure combinatorics (coloring
//! takes a CSR-ish adjacency), and take an explicit edge-filter predicate so ONE
//! `component_labels` serves both the all-edge partition and the overlap-only
//! sub-partition (Python's `gt.GraphView(efilt=...)` + `label_components`).
//!
//! ## Parity notes (the highest-risk part of the migration)
//!
//! - **`component_labels`** — labels are canonicalized to "smallest member node
//!   index ⇒ component id", then densified to `0..k` in ascending first-appearance
//!   order, so the labeling is deterministic and order-stable. The PASS CRITERION
//!   compares components as a **set partition**, not by label id, so exact gt label
//!   ids do not matter — only the grouping.
//! - **`greedy_vertex_coloring`** — reproduces `gt.sequential_vertex_coloring`
//!   EXACTLY: vertices are colored in **ascending vertex-index order**, each given
//!   the **smallest color not used by any already-colored neighbor** (graph-tool's
//!   default greedy in index order). petgraph's `dsatur_coloring` is a *different*
//!   heuristic and must NOT be used for parity. ⚠ The vertex-index order is set by
//!   the caller's node-insertion order into the `connected_qnodes` graph
//!   (`grouping`/`traversal`), which must reproduce the Python insertion order —
//!   see `grouping.rs` for that contract.
//! - **`dijkstra_route`** — a predecessor-recording binary-heap Dijkstra (petgraph's
//!   `dijkstra` returns only costs). Equal-cost tie-breaking may differ from
//!   `gt.shortest_path`; this is a documented parity risk (T8 §6). We break ties by
//!   the **lower target node index** when popping equal-cost frontier entries, a
//!   deterministic rule; if a differential surfaces a divergence the tie-rule is
//!   the first place to look.

use crate::graph_build::{EdgeAttr, SdGraph};
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use rustc_hash::FxHashMap;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// Component label per node index, treating the (directed) `SdGraph` as
/// **undirected** and keeping only edges for which `edge_filter` returns true.
///
/// Replaces `gt.label_components` (which returns per-vertex labels, not a count)
/// AND the overlap-only sub-partition (`gt.GraphView(efilt=overlap)` +
/// `label_components`) — ONE versatile unit, the predicate selects which.
///
/// Returns a `Vec<u32>` indexed by `NodeIndex::index()`; isolated nodes get their
/// own singleton component. Labels are densified `0..k` in ascending
/// first-appearance (node-index) order so the result is deterministic. Compare as
/// a partition (the pass criterion), not by raw label value.
///
/// `&SdGraph` + `Fn(&EdgeAttr)->bool`: borrow-only; the predicate lets the caller
/// pick all-edges (`|_| true`) or overlap-only (`|a| a.is_overlap`).
pub fn component_labels(g: &SdGraph, edge_filter: impl Fn(&EdgeAttr) -> bool) -> Vec<u32> {
    let n = g.g.node_count();
    let mut uf = UnionFind::new(n);
    for e in g.g.edge_references() {
        if edge_filter(e.weight()) {
            uf.union(e.source().index(), e.target().index());
        }
    }
    // Densify roots to 0..k in ascending node-index order.
    let mut root_to_label: FxHashMap<usize, u32> = FxHashMap::default();
    let mut next: u32 = 0;
    let mut labels = vec![0u32; n];
    for (i, slot) in labels.iter_mut().enumerate() {
        let r = uf.find(i);
        let label = *root_to_label.entry(r).or_insert_with(|| {
            let l = next;
            next += 1;
            l
        });
        *slot = label;
    }
    labels
}

/// Greedy vertex coloring in **ascending vertex-index order** — exact parity with
/// `gt.sequential_vertex_coloring`.
///
/// `adj[i]` is the list of neighbor indices of vertex `i` (an undirected
/// adjacency; both directions present). Each vertex `i = 0,1,2,...` is assigned
/// the smallest color (`u32`, from 0) not used by any neighbor `j < ... ` that is
/// already colored. Returns the color per vertex.
///
/// `&[Vec<u32>]`: pure combinatorics on the connected-qnodes adjacency, decoupled
/// from petgraph (cache-friendly, independently unit-testable). The caller builds
/// the adjacency in the node-insertion order that fixes the coloring index order.
pub fn greedy_vertex_coloring(adj: &[Vec<u32>]) -> Vec<u32> {
    let n = adj.len();
    let mut colors = vec![u32::MAX; n]; // MAX = uncolored
    // Scratch set of forbidden colors, reused per vertex.
    let mut forbidden: Vec<bool> = Vec::new();
    for v in 0..n {
        // Mark colors used by already-colored neighbors.
        forbidden.clear();
        for &nb in &adj[v] {
            let c = colors[nb as usize];
            if c != u32::MAX {
                let ci = c as usize;
                if ci >= forbidden.len() {
                    forbidden.resize(ci + 1, false);
                }
                forbidden[ci] = true;
            }
        }
        // Smallest free color.
        let mut chosen = 0u32;
        while (chosen as usize) < forbidden.len() && forbidden[chosen as usize] {
            chosen += 1;
        }
        colors[v] = chosen;
    }
    colors
}

/// Weighted single-source → single-target shortest path returning the **edge
/// route** (petgraph's `dijkstra` returns only distances). Treats the directed
/// `SdGraph` as **undirected** (Python runs `shortest_path` on an undirected
/// graph). Edge cost = `EdgeAttr::weight()` (Python `graph.ep["weight"]`).
///
/// Returns `Some((total_cost, route))` where `route` is the ordered list of edge
/// indices from `src` to `tgt`, or `None` if `tgt` is unreachable. The route is an
/// owned `Vec<EdgeIndex>`; the caller resolves each `EdgeIndex` to its `EdgeAttr`
/// + endpoint nodes (the route is consumed by `inspect_cnode_along_route`).
///
/// Tie-break: the heap pops the lowest-cost entry, and on equal cost the lowest
/// target-node index (deterministic). See the module doc parity note.
pub fn dijkstra_route(
    g: &SdGraph,
    src: NodeIndex,
    tgt: NodeIndex,
) -> Option<(f64, Vec<EdgeIndex>)> {
    let n = g.g.node_count();
    let mut dist = vec![f64::INFINITY; n];
    // came_from[v] = (predecessor node, edge used to reach v)
    let mut came_from: Vec<Option<(NodeIndex, EdgeIndex)>> = vec![None; n];
    let mut visited = vec![false; n];

    dist[src.index()] = 0.0;
    let mut heap = BinaryHeap::new();
    heap.push(HeapEntry {
        cost: 0.0,
        node: src.index(),
    });

    while let Some(HeapEntry { cost, node }) = heap.pop() {
        if visited[node] {
            continue;
        }
        visited[node] = true;
        if node == tgt.index() {
            break;
        }
        let nidx = NodeIndex::new(node);
        // Undirected: iterate both outgoing and incoming edges.
        for e in g.g.edges(nidx).chain(g.g.edges_directed(nidx, petgraph::Direction::Incoming)) {
            // `e.weight()` is petgraph's edge payload (&EdgeAttr); `.weight()` is
            // Python's `graph.ep["weight"]` (PO weight wins on a combined edge).
            let w = e.weight().weight();
            // The neighbor is the endpoint that is not `node`.
            let other = if e.source() == nidx { e.target() } else { e.source() };
            let oi = other.index();
            if visited[oi] {
                continue;
            }
            let nd = cost + w;
            if nd < dist[oi] {
                dist[oi] = nd;
                came_from[oi] = Some((nidx, e.id()));
                heap.push(HeapEntry { cost: nd, node: oi });
            }
        }
    }

    if !visited[tgt.index()] && came_from[tgt.index()].is_none() && src != tgt {
        return None;
    }
    if src == tgt {
        return Some((0.0, Vec::new()));
    }
    if dist[tgt.index()].is_infinite() {
        return None;
    }

    // Reconstruct edge route src -> tgt.
    let mut route: Vec<EdgeIndex> = Vec::new();
    let mut cur = tgt;
    while cur != src {
        let (prev, eidx) = came_from[cur.index()]?;
        route.push(eidx);
        cur = prev;
    }
    route.reverse();
    Some((dist[tgt.index()], route))
}

// ── internal: a small, dependency-free union-find ────────────────────────────

/// Union-find with path compression + union by rank. Equivalent to
/// `petgraph::unionfind::UnionFind` but kept local so `component_labels` controls
/// the densification order; petgraph's `into_labeling` would also work but yields
/// labels keyed on its internal representative choice.
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
        // Iterative path-halving.
        let mut x = x;
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }

    fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
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

// ── internal: min-heap entry for Dijkstra ────────────────────────────────────

/// A `BinaryHeap` is a max-heap; we invert the `Ord` on `cost` to get a min-heap.
/// Equal costs break ties on the **lower node index** (deterministic). NaN costs
/// are not expected (weights are finite mismatch rates / overlap fractions); if
/// one appears it sorts as "greater" (popped last), which is safe.
struct HeapEntry {
    cost: f64,
    node: usize,
}

impl PartialEq for HeapEntry {
    fn eq(&self, other: &Self) -> bool {
        self.cost == other.cost && self.node == other.node
    }
}
impl Eq for HeapEntry {}

impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        // Min-heap: lower cost = "greater" so it is popped first. On equal cost,
        // lower node index = "greater" (popped first) for a deterministic order.
        match other
            .cost
            .partial_cmp(&self.cost)
            .unwrap_or(Ordering::Equal)
        {
            Ordering::Equal => other.node.cmp(&self.node),
            ord => ord,
        }
    }
}
impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph_build::{EdgeKind, NodeKey, SdGraph};
    use sdrecall_utils::Strand;

    fn nk(i: i64) -> NodeKey {
        NodeKey::new("chr1", i * 1000, i * 1000 + 500, Strand::Forward)
    }

    /// Add an edge with a given kind/weight to an SdGraph (interning nodes by an
    /// integer label so tests read cleanly). `weight` is the SD `mismatch_rate` for
    /// an SD edge, or the PO weight `1/overlap_frac` for an overlap edge.
    fn add_edge(g: &mut SdGraph, a: i64, b: i64, kind: EdgeKind, weight: f64) {
        let u = g.node(nk(a));
        let v = g.node(nk(b));
        let attr = match kind {
            EdgeKind::SegmentalDuplication => EdgeAttr::sd(weight),
            EdgeKind::Overlap => EdgeAttr::po(weight, 0),
        };
        g.g.add_edge(u, v, attr);
    }

    /// Partition (set of frozensets) induced by a labeling, for order-independent
    /// comparison.
    fn partition(labels: &[u32]) -> std::collections::BTreeSet<std::collections::BTreeSet<usize>> {
        let mut by_label: FxHashMap<u32, std::collections::BTreeSet<usize>> = FxHashMap::default();
        for (i, &l) in labels.iter().enumerate() {
            by_label.entry(l).or_default().insert(i);
        }
        by_label.into_values().collect()
    }

    // ── component_labels ─────────────────────────────────────────────────────

    #[test]
    fn components_single_connected() {
        // 0-1-2-3 path, all SD edges → one component.
        let mut g = SdGraph::new();
        add_edge(&mut g, 0, 1, EdgeKind::SegmentalDuplication, 0.1);
        add_edge(&mut g, 1, 2, EdgeKind::SegmentalDuplication, 0.1);
        add_edge(&mut g, 2, 3, EdgeKind::SegmentalDuplication, 0.1);
        let labels = component_labels(&g, |_| true);
        assert_eq!(labels.iter().collect::<std::collections::BTreeSet<_>>().len(), 1);
    }

    #[test]
    fn components_two_disjoint() {
        // {0-1} and {2-3} → two components.
        let mut g = SdGraph::new();
        add_edge(&mut g, 0, 1, EdgeKind::SegmentalDuplication, 0.1);
        add_edge(&mut g, 2, 3, EdgeKind::SegmentalDuplication, 0.1);
        let labels = component_labels(&g, |_| true);
        let part = partition(&labels);
        assert_eq!(part.len(), 2);
        // Each component has exactly 2 members.
        assert!(part.iter().all(|c| c.len() == 2));
    }

    #[test]
    fn components_overlap_only_predicate_splits_off_po_bridge() {
        // 0 -SD- 1 -PO- 2 -SD- 3.
        // All-edges: one component {0,1,2,3}.
        // Overlap-only: only the PO edge 1-2 survives → components {1,2}, plus
        // singletons {0} and {3}.
        let mut g = SdGraph::new();
        add_edge(&mut g, 0, 1, EdgeKind::SegmentalDuplication, 0.1);
        add_edge(&mut g, 1, 2, EdgeKind::Overlap, 0.5);
        add_edge(&mut g, 2, 3, EdgeKind::SegmentalDuplication, 0.1);

        let all = component_labels(&g, |_| true);
        assert_eq!(partition(&all).len(), 1);

        let ov = component_labels(&g, |a| a.is_overlap);
        let part = partition(&ov);
        // {0}, {1,2}, {3}
        assert_eq!(part.len(), 3);
        assert!(part.contains(&[1usize, 2].into_iter().collect()));
        assert!(part.contains(&[0usize].into_iter().collect()));
        assert!(part.contains(&[3usize].into_iter().collect()));
    }

    #[test]
    fn components_isolated_node_is_own_component() {
        let mut g = SdGraph::new();
        g.node(nk(0)); // isolated
        add_edge(&mut g, 1, 2, EdgeKind::SegmentalDuplication, 0.1);
        let labels = component_labels(&g, |_| true);
        assert_eq!(partition(&labels).len(), 2);
    }

    // ── greedy_vertex_coloring ───────────────────────────────────────────────

    /// Group vertices by color into a partition for parity comparison.
    fn color_partition(
        colors: &[u32],
    ) -> std::collections::BTreeSet<std::collections::BTreeSet<usize>> {
        let mut by_color: FxHashMap<u32, std::collections::BTreeSet<usize>> = FxHashMap::default();
        for (i, &c) in colors.iter().enumerate() {
            by_color.entry(c).or_default().insert(i);
        }
        by_color.into_values().collect()
    }

    /// Check the coloring is proper (no edge joins same-color vertices).
    fn is_proper(adj: &[Vec<u32>], colors: &[u32]) -> bool {
        for (v, nbs) in adj.iter().enumerate() {
            for &nb in nbs {
                if colors[v] == colors[nb as usize] {
                    return false;
                }
            }
        }
        true
    }

    #[test]
    fn coloring_path_p4() {
        // P4: 0-1-2-3. gt sequential coloring in index order:
        // v0=0; v1: nb {0}->1; v2: nb {1}->0; v3: nb {2(=0)}->1.
        // colors = [0,1,0,1]; partition {{0,2},{1,3}}.
        let adj = vec![vec![1], vec![0, 2], vec![1, 3], vec![2]];
        let colors = greedy_vertex_coloring(&adj);
        assert_eq!(colors, vec![0, 1, 0, 1]);
        assert!(is_proper(&adj, &colors));
        let part = color_partition(&colors);
        assert_eq!(
            part,
            [[0usize, 2].into_iter().collect(), [1usize, 3].into_iter().collect()]
                .into_iter()
                .collect()
        );
    }

    #[test]
    fn coloring_five_cycle_needs_three_colors() {
        // C5: 0-1-2-3-4-0. Greedy index-order:
        // v0=0; v1 nb{0}->1; v2 nb{1}->0; v3 nb{2}->1; v4 nb{3,0}->{1,0}->2.
        // colors=[0,1,0,1,2]; 3 colors, proper.
        let adj = vec![vec![1, 4], vec![0, 2], vec![1, 3], vec![2, 4], vec![3, 0]];
        let colors = greedy_vertex_coloring(&adj);
        assert_eq!(colors, vec![0, 1, 0, 1, 2]);
        assert!(is_proper(&adj, &colors));
    }

    #[test]
    fn coloring_k3_plus_isolated() {
        // Triangle 0-1-2 (K3) + isolated vertex 3.
        // v0=0; v1 nb{0}->1; v2 nb{0,1}->2; v3 nb{} ->0.
        let adj = vec![vec![1, 2], vec![0, 2], vec![0, 1], vec![]];
        let colors = greedy_vertex_coloring(&adj);
        assert_eq!(colors, vec![0, 1, 2, 0]);
        assert!(is_proper(&adj, &colors));
    }

    #[test]
    fn coloring_empty_graph_no_edges() {
        // 4 isolated vertices → all color 0 (one group).
        let adj = vec![vec![], vec![], vec![], vec![]];
        let colors = greedy_vertex_coloring(&adj);
        assert_eq!(colors, vec![0, 0, 0, 0]);
    }

    #[test]
    fn coloring_is_always_proper_on_random_graphs() {
        // Property: the coloring is always proper. Deterministic LCG, several graphs.
        let mut seed: u64 = 0x1234_5678;
        let mut rng = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (seed >> 33) as u32
        };
        for _ in 0..50 {
            let n = 4 + (rng() % 12) as usize;
            let mut adj = vec![Vec::new(); n];
            for i in 0..n {
                for j in (i + 1)..n {
                    if rng() % 3 == 0 {
                        adj[i].push(j as u32);
                        adj[j].push(i as u32);
                    }
                }
            }
            let colors = greedy_vertex_coloring(&adj);
            assert!(is_proper(&adj, &colors), "improper coloring on random graph");
        }
    }

    // ── dijkstra_route ───────────────────────────────────────────────────────

    #[test]
    fn dijkstra_simple_two_edge_route() {
        // 0 -(0.1)- 1 -(0.2)- 2.  src=0 tgt=2 → route [e01, e12], cost 0.3.
        let mut g = SdGraph::new();
        add_edge(&mut g, 0, 1, EdgeKind::SegmentalDuplication, 0.1);
        add_edge(&mut g, 1, 2, EdgeKind::SegmentalDuplication, 0.2);
        let src = g.index.get(&nk(0)).copied().unwrap();
        let tgt = g.index.get(&nk(2)).copied().unwrap();
        let (cost, route) = dijkstra_route(&g, src, tgt).unwrap();
        assert!((cost - 0.3).abs() < 1e-9);
        assert_eq!(route.len(), 2);
        // Walk the route and confirm it connects src->tgt.
        let mut cur = src;
        for &e in &route {
            let (a, b) = g.g.edge_endpoints(e).unwrap();
            cur = if a == cur { b } else { a };
        }
        assert_eq!(cur, tgt);
    }

    #[test]
    fn dijkstra_prefers_lower_cost_path() {
        // 0-1 (cost 1.0), 0-2 (0.1), 2-1 (0.1). Shortest 0->1 is via 2 (0.2).
        let mut g = SdGraph::new();
        add_edge(&mut g, 0, 1, EdgeKind::SegmentalDuplication, 1.0);
        add_edge(&mut g, 0, 2, EdgeKind::SegmentalDuplication, 0.1);
        add_edge(&mut g, 2, 1, EdgeKind::SegmentalDuplication, 0.1);
        let src = g.index.get(&nk(0)).copied().unwrap();
        let tgt = g.index.get(&nk(1)).copied().unwrap();
        let (cost, route) = dijkstra_route(&g, src, tgt).unwrap();
        assert!((cost - 0.2).abs() < 1e-9, "cost {cost}");
        assert_eq!(route.len(), 2, "two-hop route via node 2");
    }

    #[test]
    fn dijkstra_unreachable_is_none() {
        let mut g = SdGraph::new();
        add_edge(&mut g, 0, 1, EdgeKind::SegmentalDuplication, 0.1);
        g.node(nk(5)); // disconnected
        let src = g.index.get(&nk(0)).copied().unwrap();
        let tgt = g.index.get(&nk(5)).copied().unwrap();
        assert!(dijkstra_route(&g, src, tgt).is_none());
    }

    #[test]
    fn dijkstra_src_eq_tgt_empty_route() {
        let mut g = SdGraph::new();
        add_edge(&mut g, 0, 1, EdgeKind::SegmentalDuplication, 0.1);
        let src = g.index.get(&nk(0)).copied().unwrap();
        let (cost, route) = dijkstra_route(&g, src, src).unwrap();
        assert_eq!(cost, 0.0);
        assert!(route.is_empty());
    }
}
