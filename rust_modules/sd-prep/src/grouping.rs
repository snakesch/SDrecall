//! RG node grouping via greedy vertex coloring — the `optimal_node_grouping`
//! fast path.
//!
//! Ports `preparation/graph_query.py::optimal_node_grouping` (l.337-428), but
//! ONLY the production-reachable **`min_distance==2, max_similarity is None`**
//! fast path (l.355-367): color the connected-qnodes graph with
//! `gt.sequential_vertex_coloring` and group vertices by color. The
//! general/conflict-graph branch (l.369-428) is **DEAD for parity** (the
//! production caller `query_connected_nodes` never passes overrides) and is NOT
//! ported — porting it would create the near-duplicate the #1 rule forbids.
//!
//! ## The vertex-insertion-order contract (the parity hazard)
//!
//! `gt.sequential_vertex_coloring` colors vertices in **ascending internal index
//! order**, and that index order is the **node-insertion order** into the
//! `connected_qnodes_gt` graph (graph_query.py l.216-273):
//!   1. all qnodes first, in `sorted_query_nodes` order (degree×component-size
//!      descending, l.231-234);
//!   2. then counter-qnodes (cnodes that are themselves query nodes), in the
//!      order they were yielded by the per-qnode traversal results (l.264-270).
//!
//! [`ConnectedQnodes`] is the order-preserving builder that reproduces this: nodes
//! are interned in insertion order, so vertex index == insertion rank, feeding
//! [`crate::graph_core::greedy_vertex_coloring`] in the exact order Python does.
//! The production traversal drives insertion; this module owns the
//! ordering contract + the color→group collapse so the contract lives in one place.

use crate::graph_core::greedy_vertex_coloring;
use rustc_hash::FxHashMap;

/// An order-preserving undirected graph over query/counter-query node keys, built
/// in the Python `connected_qnodes_gt` insertion order. Generic over the node key
/// `K` (the pipeline uses `NodeKey`; tests use small integer keys).
///
/// Vertex index == insertion order, which is what `greedy_vertex_coloring`
/// consumes — so the caller MUST insert nodes in the Python order (qnodes by
/// `sorted_query_nodes`, then cnodes in traversal-result order). This type does
/// not sort; it preserves the order it is given.
pub struct ConnectedQnodes<K: Clone + std::hash::Hash + Eq> {
    keys: Vec<K>,
    index: FxHashMap<K, usize>,
    adj: Vec<Vec<u32>>,
    // Dedup edges so a repeated qnode↔cnode insertion does not double the adjacency
    // (Python uses a graph where add_edge twice would create a multigraph, but the
    // coloring only cares about the neighbor SET; we keep the set to be safe).
    edge_seen: ahash::AHashSet<(usize, usize)>,
}

impl<K: Clone + std::hash::Hash + Eq> ConnectedQnodes<K> {
    pub fn new() -> Self {
        Self {
            keys: Vec::new(),
            index: FxHashMap::default(),
            adj: Vec::new(),
            edge_seen: ahash::AHashSet::new(),
        }
    }

    /// Intern a node by key in insertion order, returning its vertex index. If the
    /// node already exists its index is returned unchanged (insertion order is set
    /// by first appearance — matching Python's "add vertex if not in
    /// data_to_vertex").
    pub fn node(&mut self, key: K) -> usize {
        if let Some(&i) = self.index.get(&key) {
            return i;
        }
        let i = self.keys.len();
        self.keys.push(key.clone());
        self.index.insert(key, i);
        self.adj.push(Vec::new());
        i
    }

    /// Add an undirected edge between two interned nodes (deduped).
    pub fn edge(&mut self, a: usize, b: usize) {
        if a == b {
            return;
        }
        let key = if a < b { (a, b) } else { (b, a) };
        if self.edge_seen.insert(key) {
            self.adj[a].push(b as u32);
            self.adj[b].push(a as u32);
        }
    }

    pub fn len(&self) -> usize {
        self.keys.len()
    }
    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    /// Number of distinct (deduped, undirected) edges wired so far — for
    /// diagnostics and tests (e.g. asserting the FIX-#9 gating wires the right
    /// number of grouping edges).
    pub fn edge_count(&self) -> usize {
        self.edge_seen.len()
    }

    /// The node key at a vertex index (for mapping groups back to keys).
    pub fn key_at(&self, idx: usize) -> &K {
        &self.keys[idx]
    }

    /// Run the fast-path coloring and group vertex indices by color. Groups are
    /// returned in ascending color order; within a group, vertices are in
    /// ascending index (insertion) order. Mirrors `optimal_node_grouping`'s
    /// `defaultdict(list)` grouped-by-color, whose iteration order over colors is
    /// first-seen — but ascending color order is equivalent for the fast path
    /// since colors are assigned `0,1,2,...` in increasing order, and the PASS
    /// CRITERION compares groups as a set partition (color order does not matter).
    pub fn color_groups(&self) -> Vec<Vec<usize>> {
        let colors = greedy_vertex_coloring(&self.adj);
        let max_color = colors.iter().copied().max().unwrap_or(0);
        let mut groups: Vec<Vec<usize>> = vec![Vec::new(); (max_color as usize) + 1];
        for (v, &c) in colors.iter().enumerate() {
            groups[c as usize].push(v);
        }
        groups.into_iter().filter(|g| !g.is_empty()).collect()
    }

    /// Like [`color_groups`] but maps each group to the node keys.
    pub fn color_groups_keys(&self) -> Vec<Vec<K>> {
        self.color_groups()
            .into_iter()
            .map(|g| g.into_iter().map(|i| self.keys[i].clone()).collect())
            .collect()
    }
}

impl<K: Clone + std::hash::Hash + Eq> Default for ConnectedQnodes<K> {
    fn default() -> Self {
        Self::new()
    }
}

/// Convenience: color-group a prebuilt adjacency directly (the fast-path body of
/// `optimal_node_grouping`). Exposed so callers that already hold a
/// connected-qnodes adjacency (e.g. read back from the qnode-grouping GraphML for
/// a differential) can group it without rebuilding a `ConnectedQnodes`.
pub fn optimal_node_grouping(adj: &[Vec<u32>]) -> Vec<Vec<usize>> {
    let colors = greedy_vertex_coloring(adj);
    let max_color = colors.iter().copied().max().unwrap_or(0);
    let mut groups: Vec<Vec<usize>> = vec![Vec::new(); (max_color as usize) + 1];
    for (v, &c) in colors.iter().enumerate() {
        groups[c as usize].push(v);
    }
    groups.into_iter().filter(|g| !g.is_empty()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn insertion_order_is_vertex_index() {
        let mut cq: ConnectedQnodes<&str> = ConnectedQnodes::new();
        let a = cq.node("qA");
        let b = cq.node("qB");
        let c = cq.node("cC");
        assert_eq!((a, b, c), (0, 1, 2));
        // re-inserting returns the same index (first appearance wins).
        assert_eq!(cq.node("qA"), 0);
        assert_eq!(cq.len(), 3);
    }

    #[test]
    fn color_groups_path_partition() {
        // Build P4: qA-qB-qC-qD in insertion order 0,1,2,3.
        let mut cq: ConnectedQnodes<&str> = ConnectedQnodes::new();
        let a = cq.node("qA");
        let b = cq.node("qB");
        let c = cq.node("qC");
        let d = cq.node("qD");
        cq.edge(a, b);
        cq.edge(b, c);
        cq.edge(c, d);
        // colors [0,1,0,1] → groups {0,2},{1,3}.
        let groups = cq.color_groups();
        assert_eq!(groups, vec![vec![0, 2], vec![1, 3]]);
        let key_groups = cq.color_groups_keys();
        assert_eq!(key_groups, vec![vec!["qA", "qC"], vec!["qB", "qD"]]);
    }

    #[test]
    fn color_groups_independent_set_one_group() {
        // No edges → all one color → one group with everyone.
        let mut cq: ConnectedQnodes<i32> = ConnectedQnodes::new();
        for k in 0..5 {
            cq.node(k);
        }
        let groups = cq.color_groups();
        assert_eq!(groups, vec![vec![0, 1, 2, 3, 4]]);
    }

    #[test]
    fn color_groups_triangle_three_groups() {
        // K3: each vertex its own color → 3 singleton groups.
        let mut cq: ConnectedQnodes<i32> = ConnectedQnodes::new();
        let a = cq.node(0);
        let b = cq.node(1);
        let c = cq.node(2);
        cq.edge(a, b);
        cq.edge(b, c);
        cq.edge(a, c);
        let groups = cq.color_groups();
        assert_eq!(groups, vec![vec![0], vec![1], vec![2]]);
    }

    #[test]
    fn dedup_edge_does_not_double_adjacency() {
        let mut cq: ConnectedQnodes<i32> = ConnectedQnodes::new();
        let a = cq.node(0);
        let b = cq.node(1);
        cq.edge(a, b);
        cq.edge(b, a); // same undirected edge
                       // adjacency should list b once for a.
        let groups = cq.color_groups();
        assert_eq!(groups, vec![vec![0], vec![1]]);
    }

    #[test]
    fn free_function_matches_method() {
        let adj = vec![vec![1], vec![0, 2], vec![1, 3], vec![2]];
        let groups = optimal_node_grouping(&adj);
        assert_eq!(groups, vec![vec![0, 2], vec![1, 3]]);
    }
}
