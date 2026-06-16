//! SD-graph traversal → SD paralog pairs (the Phase-1 hotspot).
//!
//! Ports `preparation/graph_query.py::extract_SD_paralog_pairs_from_graph`
//! (l.133-283) + `graph_traversal.py` (the per-qnode shortest-path walk +
//! minimap2 candidate validation).
//!
//! ## Port status (implemented)
//!
//! - [`sort_query_nodes`] — the load-balance sort (`degree × component_size`
//!   descending) that ALSO fixes the coloring index order (graph_query.py
//!   l.103-126). Pure given the graph + component labels; the parity-critical
//!   input to [`crate::grouping`]. UNIT-TESTED.
//! - [`prune_graph`] — the two graph-pruning passes (graph_query.py l.156-185):
//!   drop small-SD nodes ([`is_small_sd`]: `size <= max(mean_read_len, avg_frag -
//!   std)`), then drop self-loops and weak PO edges ([`should_prune_po_edge`]:
//!   `overlap_size < cutoff && 1/weight < 0.5`). Applied at the top of
//!   [`extract_sd_paralog_pairs`] so labeling/sorting see the pruned graph, as in
//!   Python. The predicates are pure + UNIT-TESTED; `prune_graph` rebuilds a fresh
//!   graph (the `DiGraph` is not stable, so `remove_node` is avoided).
//! - [`extract_sd_paralog_pairs`] / `traverse_qnode` — the per-qnode component
//!   walk + `summarize_shortest_paths_per_subgraph` (route walk via
//!   [`inspect_cnode_along_route`] + minimap2 similarity through
//!   [`crate::minimap::align_similarity`]). Implemented on the graph primitives
//!   (`component_labels` all-edge + overlap-only, `dijkstra_route`,
//!   `HomoseqRegion`); the connected-qnodes graph is built in Python insertion
//!   order for coloring parity. Validated end-to-end by the
//!   `validate_phase1_e2e` differential harness.

use crate::graph_build::{EdgeAttr, EdgeKind, NodeKey, SdGraph};
use crate::graph_core::{component_labels, dijkstra_route};
use crate::homoseq::{HomoseqRegion, RouteStep};
use crate::minimap::{align_similarity, Preset, SIMILARITY_THRESHOLD};
use bio::alphabets::dna::revcomp;
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use rustc_hash::FxHashMap;
use sdrecall_utils::{fatal_invariant, Result, SdError};
use std::path::Path;

/// The fragment-size + read-length scalars threaded through the route walk — one
/// shared bundle (so the route-walk + traversal fns take a single `&FragParams`
/// instead of three scalars each; reused by `inspect_cnode_along_route`,
/// `traverse_qnode` and `extract_sd_paralog_pairs`).
#[derive(Clone, Copy, Debug)]
pub struct FragParams {
    pub avg_frag: f64,
    pub std_frag: f64,
    pub mean_read_length: f64,
}

impl FragParams {
    /// The minimum-region floor `max(mean_read_length, avg_frag - 2*std_frag)`
    /// (graph_traversal.py l.53, l.67).
    fn floor(&self) -> f64 {
        self.mean_read_length.max(self.avg_frag - 2.0 * self.std_frag)
    }
}

/// Sort query-node vertex indices for load balancing AND coloring-order parity
/// (`sort_query_nodes`, graph_query.py l.103-126): key = `degree(node) ×
/// component_size(node)`, **descending**.
///
/// `component_labels` is the all-edge component labeling (from
/// [`crate::graph_core::component_labels`] on the undirected graph). Ties: Python's
/// `sorted(..., reverse=True)` is a stable sort, so equal keys keep their input
/// order — here the input order is the `query_nodes` order the driver passes
/// (`total_bin_sd_df` row order). We reproduce a **stable descending** sort.
///
/// `degree` is the undirected degree (in + out, since SD/PO edges are stored
/// directed but treated undirected). Returns the query nodes reordered.
pub fn sort_query_nodes(
    query_nodes: &[NodeIndex],
    g: &SdGraph,
    component_labels: &[u32],
) -> Vec<NodeIndex> {
    // component size = count of nodes sharing each label.
    let mut comp_size: rustc_hash::FxHashMap<u32, usize> = rustc_hash::FxHashMap::default();
    for &l in component_labels {
        *comp_size.entry(l).or_insert(0) += 1;
    }

    // Stable sort descending by key. We attach the original position to keep the
    // sort stable under descending order (sort_by is stable, but we make the
    // ordering total so equal keys preserve input order explicitly).
    let mut indexed: Vec<(usize, NodeIndex, i64)> = query_nodes
        .iter()
        .enumerate()
        .map(|(pos, &v)| {
            let deg = undirected_degree(g, v) as i64;
            let csize = component_labels
                .get(v.index())
                .and_then(|l| comp_size.get(l))
                .copied()
                .unwrap_or(0) as i64;
            (pos, v, deg * csize)
        })
        .collect();
    // Descending key; stable on ties via original position ascending.
    indexed.sort_by(|a, b| b.2.cmp(&a.2).then(a.0.cmp(&b.0)));
    indexed.into_iter().map(|(_, v, _)| v).collect()
}

/// Undirected degree of a node: count incident edges in both directions.
fn undirected_degree(g: &SdGraph, v: NodeIndex) -> usize {
    g.g.edges(v).count() + g.g.edges_directed(v, petgraph::Direction::Incoming).count()
}

/// Predicate: a node is "small" and should be pruned iff its size is `<= cutoff`,
/// where `cutoff = max(mean_read_length, avg_frag - std)` (graph_query.py
/// l.157-167). Pure; the driver removes matching nodes before traversal.
pub fn is_small_sd(size: i64, mean_read_length: f64, avg_frag: f64, std_frag: f64) -> bool {
    let cutoff = mean_read_length.max(avg_frag - std_frag);
    (size as f64) <= cutoff
}

/// Predicate: a PO edge should be pruned iff `overlap_size < cutoff && 1/weight <
/// 0.5` (graph_query.py l.170-184), where `overlap_size = smaller_node_size *
/// (1/weight)` and `weight` is the stored edge weight (`1/overlap_frac`). Self
/// loops are pruned separately by the caller.
///
/// `smaller_node_size` is the size of the edge's TARGET node (Python uses
/// `edge[1]` — the small node, since PO edges point large→small).
pub fn should_prune_po_edge(attr: &EdgeAttr, smaller_node_size: i64, cutoff: f64) -> bool {
    if !attr.is_overlap() {
        return false;
    }
    // PO weight is `1/overlap_fraction ∈ [1, ∞)`; a non-finite or non-positive
    // weight is degenerate (never produced by the builder, where Python would hit
    // a ZeroDivisionError). Guard the division so it yields a definite `false`
    // (keep the edge) instead of an `inf`/NaN that would make the predicate
    // unpredictable once a NaN-policy is applied elsewhere.
    if !attr.weight.is_finite() || attr.weight <= 0.0 {
        return false;
    }
    let inv_weight = 1.0 / attr.weight; // = overlap fraction
    let overlap_size = smaller_node_size as f64 * inv_weight;
    overlap_size < cutoff && inv_weight < 0.5
}

/// Prune the multiplex graph before traversal — port of `graph_query.py`
/// l.156-185 (`extract_SD_paralog_pairs_from_graph`):
///   1. drop **small-SD nodes** (`size <= cutoff`, [`is_small_sd`]) where
///      `cutoff = max(mean_read_length, avg_frag - std)`;
///   2. drop **self-loops** (`edge[0] == edge[1]`) and **weak PO edges**
///      ([`should_prune_po_edge`]: `overlap_size < cutoff && 1/weight < 0.5`).
///
/// This MUST run before component labeling + `sort_query_nodes`, because the
/// component sizes and degrees those use (and hence the coloring order) are
/// computed on the **pruned** graph — exactly as Python prunes `directed_graph`
/// before `to_undirected()` / `label_components` / `sort_query_nodes`.
///
/// `SdGraph.g` is a (non-stable) `DiGraph`, so `remove_node` would invalidate
/// every later `NodeIndex` and desync the intern map. We therefore REBUILD a
/// fresh graph from the surviving nodes + edges rather than removing in place.
pub fn prune_graph(g: &SdGraph, frag: &FragParams) -> SdGraph {
    let cutoff = frag.mean_read_length.max(frag.avg_frag - frag.std_frag);
    let mut pruned = SdGraph::new();

    // (1) carry over nodes that are NOT small SDs.
    for v in g.g.node_indices() {
        let key = &g.g[v];
        if !is_small_sd(key.size(), frag.mean_read_length, frag.avg_frag, frag.std_frag) {
            pruned.node(key.clone());
        }
    }

    // (2) carry over edges whose BOTH endpoints survived, dropping self-loops and
    // weak PO edges. PO edges point large→small, so the smaller node is the edge
    // TARGET (`kb`), matching Python's `edge[1]`.
    for e in g.g.edge_indices() {
        let (a, b) = match g.g.edge_endpoints(e) {
            Some(ends) => ends,
            None => continue,
        };
        if a == b {
            continue; // self-loop
        }
        let ka = &g.g[a];
        let kb = &g.g[b];
        if !pruned.has_node(ka) || !pruned.has_node(kb) {
            continue; // an endpoint was removed as a small SD
        }
        let attr = g.g[e];
        if should_prune_po_edge(&attr, kb.size(), cutoff) {
            continue; // weak PO edge
        }
        let u = pruned.node(ka.clone());
        let w = pruned.node(kb.clone());
        pruned.g.add_edge(u, w, attr);
    }

    log::info!(
        "Pruned multiplex graph (cutoff {:.1}): {} → {} nodes, {} → {} edges",
        cutoff,
        g.node_count(),
        pruned.node_count(),
        g.edge_count(),
        pruned.edge_count(),
    );
    pruned
}

/// Reconstruct the **vertex sequence** of a dijkstra route from its edge indices.
/// `dijkstra_route` returns ordered `EdgeIndex`es; this walks them from `src` to
/// derive `[src, v1, v2, ..., tgt]` by following the endpoint that is not the
/// current node at each step (the graph is treated undirected).
///
/// # Invariant
///
/// Every `EdgeIndex` in `edges` was produced by `dijkstra_route` on the same
/// immutable `&SdGraph`. The graph cannot have been mutated between the dijkstra
/// call and this one (both take `&SdGraph`), so `edge_endpoints` returning `None`
/// is structurally impossible — it would require either a petgraph soundness bug
/// or memory corruption. We abort on violation rather than propagating a
/// recoverable error that would be silently swallowed upstream.
fn route_vertices(g: &SdGraph, src: NodeIndex, edges: &[petgraph::graph::EdgeIndex]) -> Vec<NodeIndex> {
    let mut verts = Vec::with_capacity(edges.len() + 1);
    verts.push(src);
    let mut cur = src;
    for &e in edges {
        let (a, b) = g.g.edge_endpoints(e).unwrap_or_else(|| {
            fatal_invariant!(
                "dijkstra route references edge {:?} absent from the graph \
                 (graph has {} nodes, {} edges) — this is a bug in graph_core or memory corruption",
                e,
                g.node_count(),
                g.edge_count()
            )
        });
        cur = if a == cur { b } else { a };
        verts.push(cur);
    }
    verts
}

/// The route-walk that derives the counterpart node's relative window + its
/// `traverse_route` — port of `inspect_cnode_along_route` (graph_traversal.py
/// l.13-75). Returns the final cnode as a [`HomoseqRegion`] (with `rela_start/end`
/// plus the route of upstream `(NodeKey, EdgeKind)` snapshots), or `None` when an
/// overlap window collapses below the read-length floor (the Python early returns).
///
/// `verts` is `[qnode, v1, ..., cnode]` (the dijkstra route vertices); `edges` is
/// the parallel edge route. The floor is `max(mean_read_length, avg_frag - 2*std)`.
fn inspect_cnode_along_route(
    g: &SdGraph,
    verts: &[NodeIndex],
    edges: &[petgraph::graph::EdgeIndex],
    frag: &FragParams,
) -> Option<HomoseqRegion> {
    let floor = frag.floor();

    // qnode = HOMOSEQ_REGION(verts[0]) — full window. We track its mutable
    // rela/ups/down windows in a small working struct.
    let mut qnode = HomoseqRegion::new(g.g[verts[0]].clone(), verts[0]);
    let mut cnode: Option<HomoseqRegion> = None;

    for i in 0..edges.len() {
        let v = verts[i + 1];
        let edge_attr = g.g[edges[i]];
        let mut c = HomoseqRegion::new(g.g[v].clone(), v);

        match edge_attr.kind {
            EdgeKind::SegmentalDuplication => {
                if qnode.key.strand == c.key.strand {
                    c.ups_rela_start = qnode.rela_start.max(0);
                    c.ups_rela_end = qnode.rela_end.min(c.size);
                } else {
                    c.ups_rela_end = (qnode.size - qnode.rela_start).min(c.size);
                    c.ups_rela_start = (qnode.size - qnode.rela_end).max(0);
                }
                c.rela_start = c.ups_rela_start;
                c.rela_end = c.ups_rela_end;
                if c.rela_end <= c.rela_start {
                    // Python sys.exit(1) here; we treat as a dropped route.
                    log::error!(
                        "SD route step produced empty window for {:?} (rela {}..{})",
                        c.key, c.rela_start, c.rela_end
                    );
                    return None;
                }
                let mut route = qnode.route.clone();
                route.push(RouteStep {
                    node: qnode.key.clone(),
                    edge_kind: EdgeKind::SegmentalDuplication,
                });
                c.route = route;
            }
            EdgeKind::Overlap => {
                // qnode downstream window from cnode's absolute position.
                let mut qnode_rela_start = c.key.start - qnode.key.start;
                let qnode_rela_end = (qnode_rela_start + c.size).min(qnode.size);
                qnode_rela_start = qnode_rela_start.max(0);

                qnode.down_rela_start = qnode_rela_start;
                qnode.down_rela_end = qnode_rela_end;
                let q_overlap_start = qnode.down_rela_start.max(qnode.ups_rela_start);
                let q_overlap_end = qnode.down_rela_end.min(qnode.ups_rela_end);
                let q_overlap_size = q_overlap_end - q_overlap_start;
                qnode.rela_start = q_overlap_start;
                qnode.rela_end = q_overlap_end;

                if (q_overlap_size as f64) <= floor {
                    return None;
                }
                let q_overlap_abs_start = qnode.key.start + q_overlap_start;
                let mut c_rela_start = q_overlap_abs_start - c.key.start;
                let c_rela_end = c.size.min(c_rela_start + q_overlap_size);
                c_rela_start = c_rela_start.max(0);
                c.ups_rela_start = c_rela_start;
                c.ups_rela_end = c_rela_end;
                c.rela_start = c_rela_start;
                c.rela_end = c_rela_end;

                if (c_rela_end - c_rela_start) as f64 <= floor {
                    return None;
                }
                let mut route = qnode.route.clone();
                route.push(RouteStep {
                    node: qnode.key.clone(),
                    edge_kind: EdgeKind::Overlap,
                });
                c.route = route;
            }
        }
        qnode = c.clone();
        cnode = Some(c);
    }
    cnode
}

/// Fetch a reference subsequence `[start, stop)` on `chrom`, reverse-complementing
/// when `revcomp_it` (the Python `samtools faidx -i` for opposite-strand cnodes).
fn fetch_seq(
    reader: &mut bio::io::fasta::IndexedReader<std::fs::File>,
    chrom: &str,
    start: i64,
    stop: i64,
    revcomp_it: bool,
) -> Result<Vec<u8>> {
    reader
        .fetch(chrom, start.max(0) as u64, stop.max(0) as u64)
        .map_err(|e| SdError::Io {
            path: format!("{chrom}:{start}-{stop}"),
            source: e,
        })?;
    let mut buf = Vec::new();
    reader.read(&mut buf).map_err(|e| SdError::Io {
        path: format!("{chrom}:{start}-{stop}"),
        source: e,
    })?;
    Ok(if revcomp_it { revcomp(&buf) } else { buf })
}

/// Compare the homologous sequences of the original qnode and a candidate cnode
/// (port of `compare_homologous_sequences`, graph_traversal.py l.78-175). Returns
/// the similarity (`match_len/block_len`), or `0.0` when the back-projection
/// collapses (`qnode_relative_region` → NaN) or there is no alignment.
///
/// The qnode region is the cnode's window back-projected onto `ori_qnode`; the
/// cnode region is its own `fix_coord` window. Opposite-strand cnodes are
/// reverse-complemented (the Python `-i` flag). minimap2 `asm10` similarity.
fn compare_homologous_sequences(
    ori_qnode: &NodeKey,
    cnode: &HomoseqRegion,
    reader: &mut bio::io::fasta::IndexedReader<std::fs::File>,
) -> Result<f64> {
    let (q_rela_start, q_rela_end) = match cnode.qnode_relative_region(ori_qnode) {
        Some(w) => w,
        None => return Ok(0.0), // Python "NaN" → (False, 0.0)
    };
    let q_start = ori_qnode.start + q_rela_start;
    let q_end = ori_qnode.start + q_rela_end;
    if q_end <= q_start {
        return Ok(0.0);
    }
    // cnode region = its fix_coord window.
    let (c_chrom, c_start, c_end, c_strand) = cnode.fix_coord();
    if c_end <= c_start {
        return Ok(0.0);
    }
    // Query (qnode) is always forward (Python extracts it without -i).
    let q_seq = fetch_seq(reader, &ori_qnode.chrom, q_start, q_end, false)?;
    // cnode is revcomp'd iff its strand differs from the qnode's.
    let revcomp_it = c_strand != ori_qnode.strand;
    let c_seq = fetch_seq(reader, &c_chrom, c_start, c_end, revcomp_it)?;
    // Python target = qnode (first positional), query = cnode (second).
    align_similarity(&c_seq, &q_seq, Preset::Asm10)
}

/// A counterpart found for a qnode: the cnode region + its similarity, and whether
/// the cnode is itself a query node (a "counter-qnode" feeding the coloring graph).
pub struct Counterpart {
    pub cnode: HomoseqRegion,
    pub similarity: f64,
    /// Counter-qnode keys to wire into the connected-qnodes graph (the cnode itself
    /// if it is a qnode, plus its overlap-neighbour qnodes).
    pub counter_qnodes: Vec<NodeKey>,
}

/// The output of [`extract_sd_paralog_pairs`].
pub struct TraversalResult {
    /// `qnode → its accepted counterpart cnodes` (kept at `similarity >= 0.95`).
    pub sd_paralog_pairs: FxHashMap<NodeKey, Vec<HomoseqRegion>>,
    /// The connected-qnodes adjacency, built in the Python insertion order
    /// (qnodes first in `sort_query_nodes` order, then counter-qnodes in
    /// traversal-result order). Colored by [`crate::grouping`].
    pub connected: crate::grouping::ConnectedQnodes<NodeKey>,
}

/// Traverse one qnode's component to find its homology counterparts — port of
/// `traverse_network_to_get_homology_counterparts` + `summarize_shortest_paths_per_subgraph`
/// (graph_traversal.py l.209-362).
///
/// Within the qnode's connected component: sub-partition by overlap-only edges,
/// then for every other vertex run a weighted shortest path from the qnode, reject
/// routes ending on a PO edge / with ≥2 adjacent PO edges / with SD-product
/// `∏(1-weight) <= 0.8`, walk the route to get the cnode window, and validate by
/// minimap2 similarity. Counterparts kept at `>= 0.95`; counter-qnodes kept at any
/// positive similarity (the looser coloring threshold, l.350).
fn traverse_qnode(
    qnode: &NodeKey,
    qnode_v: NodeIndex,
    g: &SdGraph,
    comp_labels: &[u32],
    qnode_vertices: &std::collections::HashSet<NodeIndex>,
    reader: &mut bio::io::fasta::IndexedReader<std::fs::File>,
    frag: &FragParams,
) -> Result<(Vec<HomoseqRegion>, Vec<NodeKey>)> {
    let my_comp = comp_labels[qnode_v.index()];
    // Overlap-only sub-partition (graph_traversal.py l.328-329).
    let overlap_labels = component_labels(g, |a| a.is_overlap());

    let mut counterparts: Vec<HomoseqRegion> = Vec::new();
    let mut counter_qnodes: Vec<NodeKey> = Vec::new();

    // Iterate vertices in the qnode's component (the Python `subgraph.vertices()`
    // over each overlap-subgraph; equivalent to iterating the whole component once,
    // since each vertex lives in exactly one overlap-subgraph).
    for v in g.g.node_indices() {
        if v == qnode_v || comp_labels[v.index()] != my_comp {
            continue;
        }
        // Both endpoints must share the same overlap-subgraph as... actually Python
        // runs shortest_path on the full component graph for every vertex in the
        // component (the overlap-subgraph only scopes which vertices are iterated;
        // all component vertices are covered across subgraphs). So we probe every
        // component vertex once. (overlap_labels is computed for parity of the
        // iteration grouping; unused for path scoping.)
        let _ = &overlap_labels;

        let (_, edges) = match dijkstra_route(g, qnode_v, v) {
            Some(r) => r,
            None => continue,
        };
        // Reject empty routes, and routes whose last edge is a PO edge.
        let Some(&last_edge) = edges.last() else {
            continue;
        };
        if g.g[last_edge].is_overlap() {
            continue;
        }
        // Reject: ≥2 adjacent PO-edge pairs (Python `len([...]) > 1`).
        let mut adjacent_po = 0;
        for w in edges.windows(2) {
            if g.g[w[0]].is_overlap() && g.g[w[1]].is_overlap() {
                adjacent_po += 1;
            }
        }
        if adjacent_po > 1 {
            continue;
        }
        // Reject: SD-similarity product ∏(1-weight) <= 0.8 over SD edges.
        let sd_product: f64 = edges
            .iter()
            .map(|&e| g.g[e])
            .filter(|a| a.is_sd())
            .map(|a| 1.0 - a.weight)
            .product();
        if sd_product <= 0.8 {
            continue;
        }

        // Walk the route → cnode window + route.
        let verts = route_vertices(g, qnode_v, &edges);
        let cnode = match inspect_cnode_along_route(g, &verts, &edges, frag) {
            Some(c) => c,
            None => continue,
        };
        let similarity = compare_homologous_sequences(qnode, &cnode, reader)?;
        if similarity > 0.9 {
            // Python: keep cnodes at >= 0.95 (committed below by caller threshold).
            let is_qnode = qnode_vertices.contains(&cnode.vertex);
            if similarity >= SIMILARITY_THRESHOLD {
                counterparts.push(cnode.clone());
            }
            // Counter-qnodes use the looser threshold (any kept cnode that is a
            // qnode), to improve coloring partitioning (l.350).
            if is_qnode {
                counter_qnodes.push(cnode.key.clone());
                // overlap-neighbour qnodes of the cnode (get_overlap_neighbors over
                // the query-node overlap view).
                for nb in overlap_neighbor_qnodes(g, cnode.vertex, qnode_vertices) {
                    counter_qnodes.push(nb);
                }
            }
        }
    }
    Ok((counterparts, counter_qnodes))
}

/// The query-node overlap neighbours of `v` — port of `get_overlap_neighbors`
/// restricted to the query-node overlap view (graph_traversal.py l.179-205, called
/// on `query_nodes_view` which is vfilt=query AND efilt=overlap). Returns the
/// NodeKeys of query-node neighbours reachable from `v` via an overlap edge.
fn overlap_neighbor_qnodes(
    g: &SdGraph,
    v: NodeIndex,
    qnode_vertices: &std::collections::HashSet<NodeIndex>,
) -> Vec<NodeKey> {
    let mut out = Vec::new();
    if !qnode_vertices.contains(&v) {
        return out;
    }
    for e in g
        .g
        .edges(v)
        .chain(g.g.edges_directed(v, petgraph::Direction::Incoming))
    {
        if !e.weight().is_overlap() {
            continue;
        }
        let other = if e.source() == v { e.target() } else { e.source() };
        if qnode_vertices.contains(&other) {
            out.push(g.g[other].clone());
        }
    }
    out
}

/// Drive the per-qnode traversal over the pruned multiplex graph — port of
/// `extract_SD_paralog_pairs_from_graph` (graph_query.py l.133-283), the Phase-1
/// hotspot. `query_nodes` are the target-overlapping SD query nodes (already
/// interned in `g`); `ref_fa` is the reference FASTA (with a `.fai` alongside).
///
/// Returns the SD paralog pairs (`qnode → cnodes`) and the connected-qnodes graph
/// built in the **exact Python insertion order** (qnodes first in
/// `sort_query_nodes` order, then counter-qnodes in traversal-result order) — the
/// coloring-parity contract. Single-threaded over qnodes for now (the FASTA reader
/// plus the minimap2 index are per-call; rayon parallelism needs a per-thread
/// reader and is a follow-up — see the report). `g` is shared read-only.
pub fn extract_sd_paralog_pairs(
    query_nodes: &[NodeKey],
    g: &SdGraph,
    ref_fa: &Path,
    avg_frag: f64,
    std_frag: f64,
    mean_read_length: f64,
) -> Result<TraversalResult> {
    let frag = FragParams {
        avg_frag,
        std_frag,
        mean_read_length,
    };

    // Prune small SDs + weak PO edges + self-loops BEFORE labeling/sorting. The
    // component sizes and degrees that drive `sort_query_nodes` (and hence the
    // coloring order) must be computed on the pruned graph, matching Python's
    // prune-then-`to_undirected()`-then-`label_components` order. Small query
    // nodes that get pruned simply have no vertex below (the `g.index.get` lookup
    // returns `None`), as in Python.
    let pruned = prune_graph(g, &frag);
    let g = &pruned;

    // All-edge component labels (gt.label_components on the undirected graph).
    let comp_labels = component_labels(g, |_| true);

    // qnode vertex indices, then sort by degree×component_size (load balance AND
    // the coloring index order — the parity contract).
    let mut qnode_vs: Vec<NodeIndex> = Vec::new();
    let mut qnode_set: std::collections::HashSet<NodeIndex> = std::collections::HashSet::new();
    for k in query_nodes {
        if let Some(&v) = g.index.get(k) {
            qnode_vs.push(v);
            qnode_set.insert(v);
        }
    }
    let sorted_qnode_vs = sort_query_nodes(&qnode_vs, g, &comp_labels);

    // Build the connected-qnodes graph in Python insertion order: all qnodes first.
    let mut connected: crate::grouping::ConnectedQnodes<NodeKey> =
        crate::grouping::ConnectedQnodes::new();
    for &v in &sorted_qnode_vs {
        connected.node(g.g[v].clone());
    }

    // Per-qnode traversal, PARALLEL across cores — this is the Phase-1 hotspot
    // (dijkstra walk + minimap2 similarity per candidate). Each rayon worker opens
    // its OWN faidx reader (the reader is stateful and not shareable) via
    // `map_init`. Results are collected in qnode order (rayon's indexed `collect`
    // preserves it), then assembled SEQUENTIALLY below so `sd_paralog_pairs` and
    // the `connected` graph come out byte-identical to the serial walk — the
    // coloring-order parity contract. `g`/`comp_labels`/`qnode_set` are shared
    // read-only; each `traverse_qnode` is deterministic and side-effect-free, so
    // the parallel result equals the serial one.
    use rayon::prelude::*;
    let per_qnode: Vec<Result<(Vec<HomoseqRegion>, Vec<NodeKey>)>> = sorted_qnode_vs
        .par_iter()
        .map_init(
            || bio::io::fasta::IndexedReader::from_file(&ref_fa),
            |reader, &qv| {
                let reader = match reader {
                    Ok(r) => r,
                    Err(e) => {
                        return Err(SdError::Io {
                            path: ref_fa.display().to_string(),
                            source: std::io::Error::other(e.to_string()),
                        })
                    }
                };
                let qkey = g.g[qv].clone();
                traverse_qnode(&qkey, qv, g, &comp_labels, &qnode_set, reader, &frag)
            },
        )
        .collect();

    let mut sd_paralog_pairs: FxHashMap<NodeKey, Vec<HomoseqRegion>> = FxHashMap::default();

    // Sequential, order-preserving assembly (identical to the serial walk).
    for (&qv, res) in sorted_qnode_vs.iter().zip(per_qnode) {
        let (counterparts, counter_qnodes) = res?;
        let qkey = g.g[qv].clone();
        if !counterparts.is_empty() {
            sd_paralog_pairs.insert(qkey.clone(), counterparts);
        }
        // Wire qnode↔counter-qnode edges (insertion order: cnodes appended in the
        // order produced above).
        let qi = connected.node(qkey);
        for ck in &counter_qnodes {
            let ci = connected.node(ck.clone());
            connected.edge(qi, ci);
        }
    }

    log::info!(
        "Traversal: {} qnodes → {} with counterparts; connected-qnodes graph has {} nodes",
        sorted_qnode_vs.len(),
        sd_paralog_pairs.len(),
        connected.len()
    );
    Ok(TraversalResult {
        sd_paralog_pairs,
        connected,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph_build::{EdgeKind, NodeKey};
    use crate::graph_core::component_labels;
    use sdrecall_utils::Strand;

    /// Standard frag params for the route-walk tests (avg 500, std 150, read 147).
    fn frag() -> FragParams {
        FragParams {
            avg_frag: 500.0,
            std_frag: 150.0,
            mean_read_length: 147.0,
        }
    }

    fn nk(i: i64) -> NodeKey {
        NodeKey::new("chr1", i * 1000, i * 1000 + 500, Strand::Forward)
    }
    fn add_edge(g: &mut SdGraph, a: i64, b: i64, kind: EdgeKind) {
        let u = g.node(nk(a));
        let v = g.node(nk(b));
        g.g.add_edge(
            u,
            v,
            EdgeAttr {
                kind,
                weight: 0.1,
                nonoverlap_smaller: 0,
            },
        );
    }

    #[test]
    fn sort_query_nodes_by_degree_times_component_size_desc() {
        // Component A: star 0-1, 0-2, 0-3 (node 0 deg 3, comp size 4).
        // Component B: edge 4-5 (deg 1, comp size 2).
        let mut g = SdGraph::new();
        add_edge(&mut g, 0, 1, EdgeKind::SegmentalDuplication);
        add_edge(&mut g, 0, 2, EdgeKind::SegmentalDuplication);
        add_edge(&mut g, 0, 3, EdgeKind::SegmentalDuplication);
        add_edge(&mut g, 4, 5, EdgeKind::SegmentalDuplication);
        let labels = component_labels(&g, |_| true);

        let v0 = g.index.get(&nk(0)).copied().unwrap();
        let v1 = g.index.get(&nk(1)).copied().unwrap();
        let v4 = g.index.get(&nk(4)).copied().unwrap();
        // query nodes: [v1, v4, v0] → keys: v0 deg3*size4=12; v1 deg1*size4=4;
        // v4 deg1*size2=2. Descending → [v0, v1, v4].
        let sorted = sort_query_nodes(&[v1, v4, v0], &g, &labels);
        assert_eq!(sorted, vec![v0, v1, v4]);
    }

    #[test]
    fn sort_query_nodes_stable_on_ties() {
        // Two separate edges → all nodes deg1, comp size2 → key all equal (2).
        // Stable sort must preserve input order.
        let mut g = SdGraph::new();
        add_edge(&mut g, 0, 1, EdgeKind::SegmentalDuplication);
        add_edge(&mut g, 2, 3, EdgeKind::SegmentalDuplication);
        let labels = component_labels(&g, |_| true);
        let v0 = g.index.get(&nk(0)).copied().unwrap();
        let v2 = g.index.get(&nk(2)).copied().unwrap();
        let v3 = g.index.get(&nk(3)).copied().unwrap();
        let input = vec![v3, v0, v2];
        let sorted = sort_query_nodes(&input, &g, &labels);
        assert_eq!(sorted, input, "ties preserve input order (stable)");
    }

    #[test]
    fn small_sd_predicate() {
        // cutoff = max(147, 500-150)=350. size 300 <= 350 → small; 400 > 350 → keep.
        assert!(is_small_sd(300, 147.0, 500.0, 150.0));
        assert!(!is_small_sd(400, 147.0, 500.0, 150.0));
        // mean_read_length dominates when frag spread is small.
        assert!(is_small_sd(100, 147.0, 200.0, 100.0)); // cutoff max(147,100)=147; 100<=147
    }

    #[test]
    fn prune_po_edge_predicate() {
        // weight stored = 1/overlap_frac. overlap_frac 0.4 → weight 2.5.
        // smaller_node_size 300, overlap_size = 300*0.4 = 120 < cutoff 350, and
        // inv_weight 0.4 < 0.5 → prune.
        let attr = EdgeAttr {
            kind: EdgeKind::Overlap,
            weight: 2.5,
            nonoverlap_smaller: 0,
        };
        assert!(should_prune_po_edge(&attr, 300, 350.0));
        // overlap_frac 0.6 → inv_weight 0.6 >= 0.5 → keep regardless.
        let attr2 = EdgeAttr {
            kind: EdgeKind::Overlap,
            weight: 1.0 / 0.6,
            nonoverlap_smaller: 0,
        };
        assert!(!should_prune_po_edge(&attr2, 300, 350.0));
        // SD edges are never PO-pruned.
        let sd = EdgeAttr {
            kind: EdgeKind::SegmentalDuplication,
            weight: 0.05,
            nonoverlap_smaller: 0,
        };
        assert!(!should_prune_po_edge(&sd, 300, 350.0));
    }

    #[test]
    fn should_prune_po_edge_guards_degenerate_weight() {
        // Zero / non-finite weight must not panic or prune via inf/NaN arithmetic.
        let zero = EdgeAttr { kind: EdgeKind::Overlap, weight: 0.0, nonoverlap_smaller: 0 };
        assert!(!should_prune_po_edge(&zero, 100, 350.0));
        let inf = EdgeAttr {
            kind: EdgeKind::Overlap,
            weight: f64::INFINITY,
            nonoverlap_smaller: 0,
        };
        assert!(!should_prune_po_edge(&inf, 100, 350.0));
    }

    #[test]
    fn prune_graph_drops_small_nodes_self_loops_and_weak_po() {
        // cutoff = max(147, 500-150) = 350.
        let big_a = NodeKey::new("chr1", 0, 1000, Strand::Forward); // size 1000 (keep)
        let big_b = NodeKey::new("chr1", 2000, 3000, Strand::Forward); // size 1000 (keep)
        let small = NodeKey::new("chr1", 5000, 5200, Strand::Forward); // size 200 (drop: small)
        let po_tgt = NodeKey::new("chr1", 4000, 4400, Strand::Forward); // size 400 (keep as node)

        let mut g = SdGraph::new();
        let a = g.node(big_a.clone());
        let b = g.node(big_b.clone());
        let s = g.node(small.clone());
        let t = g.node(po_tgt.clone());

        // SD edge between survivors → kept.
        g.g.add_edge(a, b, EdgeAttr { kind: EdgeKind::SegmentalDuplication, weight: 0.1, nonoverlap_smaller: 0 });
        // SD edge into the small node → dropped (endpoint pruned).
        g.g.add_edge(a, s, EdgeAttr { kind: EdgeKind::SegmentalDuplication, weight: 0.1, nonoverlap_smaller: 0 });
        // Self-loop → dropped.
        g.g.add_edge(a, a, EdgeAttr { kind: EdgeKind::SegmentalDuplication, weight: 0.1, nonoverlap_smaller: 0 });
        // Weak PO edge a→t: weight 3.0 → inv 0.333 < 0.5, overlap_size 400*0.333 ≈ 133 < 350 → dropped.
        g.g.add_edge(a, t, EdgeAttr { kind: EdgeKind::Overlap, weight: 3.0, nonoverlap_smaller: 0 });

        let pruned = prune_graph(&g, &frag());

        assert!(!pruned.has_node(&small), "small SD node should be dropped");
        assert!(pruned.has_node(&big_a) && pruned.has_node(&big_b) && pruned.has_node(&po_tgt));
        assert_eq!(pruned.node_count(), 3);

        // Only the surviving SD edge a-b remains (self-loop, small-endpoint, weak-PO dropped).
        assert_eq!(pruned.edge_count(), 1);
        let au = pruned.index.get(&big_a).copied().unwrap();
        let bu = pruned.index.get(&big_b).copied().unwrap();
        assert!(pruned.g.find_edge(au, bu).is_some());
    }

    #[test]
    fn prune_graph_keeps_strong_po_edge() {
        // A PO edge with inv_weight 0.667 (≥ 0.5) survives even with a small overlap.
        let a = NodeKey::new("chr1", 0, 1000, Strand::Forward);
        let b = NodeKey::new("chr1", 2000, 3000, Strand::Forward);
        let mut g = SdGraph::new();
        let au = g.node(a.clone());
        let bu = g.node(b.clone());
        g.g.add_edge(au, bu, EdgeAttr { kind: EdgeKind::Overlap, weight: 1.5, nonoverlap_smaller: 0 });
        let pruned = prune_graph(&g, &frag());
        assert_eq!(pruned.node_count(), 2);
        assert_eq!(pruned.edge_count(), 1);
    }

    // ── route walk (inspect_cnode_along_route) ───────────────────────────────

    fn nk_full(chrom: &str, start: i64, end: i64, s: Strand) -> NodeKey {
        NodeKey::new(chrom, start, end, s)
    }
    /// Add an edge between two arbitrary keys, interning them.
    fn add_edge_keys(g: &mut SdGraph, a: NodeKey, b: NodeKey, kind: EdgeKind, weight: f64) {
        let u = g.node(a);
        let v = g.node(b);
        g.g.add_edge(
            u,
            v,
            EdgeAttr {
                kind,
                weight,
                nonoverlap_smaller: 0,
            },
        );
    }

    #[test]
    fn route_vertices_follows_undirected_route() {
        // qnode(0) -SD- v1 -SD- v2. route edges -> verts [0,1,2].
        let mut g = SdGraph::new();
        add_edge(&mut g, 0, 1, EdgeKind::SegmentalDuplication);
        add_edge(&mut g, 1, 2, EdgeKind::SegmentalDuplication);
        let v0 = g.index.get(&nk(0)).copied().unwrap();
        let v2 = g.index.get(&nk(2)).copied().unwrap();
        let (_, edges) = dijkstra_route(&g, v0, v2).unwrap();
        let verts = route_vertices(&g, v0, &edges);
        assert_eq!(verts.len(), 3);
        assert_eq!(verts[0], v0);
        assert_eq!(*verts.last().unwrap(), v2);
    }

    #[test]
    fn inspect_single_sd_same_strand_full_window() {
        // qnode [10000,12000) +, cnode [20000,22000) + via ONE SD edge.
        // Same strand → cnode window = full qnode window clipped to cnode size:
        // ups_rela_start = max(0,0)=0, ups_rela_end = min(2000, 2000)=2000.
        let mut g = SdGraph::new();
        let q = nk_full("chr1", 10000, 12000, Strand::Forward);
        let c = nk_full("chr2", 20000, 22000, Strand::Forward);
        add_edge_keys(&mut g, q.clone(), c.clone(), EdgeKind::SegmentalDuplication, 0.02);
        let qv = g.index.get(&q).copied().unwrap();
        let cv = g.index.get(&c).copied().unwrap();
        let (_, edges) = dijkstra_route(&g, qv, cv).unwrap();
        let verts = route_vertices(&g, qv, &edges);
        let cnode = inspect_cnode_along_route(&g, &verts, &edges, &frag()).unwrap();
        assert_eq!(cnode.key, c);
        assert_eq!((cnode.rela_start, cnode.rela_end), (0, 2000));
        // route records the upstream qnode + SD edge kind.
        assert_eq!(cnode.route.len(), 1);
        assert_eq!(cnode.route[0].node, q);
        assert_eq!(cnode.route[0].edge_kind, EdgeKind::SegmentalDuplication);
        // back-project to the qnode → full window.
        assert_eq!(cnode.qnode_relative_region(&q), Some((0, 2000)));
    }

    #[test]
    fn inspect_single_sd_opposite_strand_flips_window() {
        // qnode [10000,12000) +, cnode [20000,22000) - via ONE SD edge.
        // Opposite strand: ups_rela_end = min(qsize - rela_start, csize) =
        //   min(2000-0, 2000)=2000; ups_rela_start = max(qsize - rela_end,0)=
        //   max(2000-2000,0)=0 → window [0,2000).
        let mut g = SdGraph::new();
        let q = nk_full("chr1", 10000, 12000, Strand::Forward);
        let c = nk_full("chr2", 20000, 22000, Strand::Reverse);
        add_edge_keys(&mut g, q.clone(), c.clone(), EdgeKind::SegmentalDuplication, 0.02);
        let qv = g.index.get(&q).copied().unwrap();
        let cv = g.index.get(&c).copied().unwrap();
        let (_, edges) = dijkstra_route(&g, qv, cv).unwrap();
        let verts = route_vertices(&g, qv, &edges);
        let cnode = inspect_cnode_along_route(&g, &verts, &edges, &frag()).unwrap();
        assert_eq!((cnode.rela_start, cnode.rela_end), (0, 2000));
    }

    #[test]
    fn inspect_overlap_too_small_returns_none() {
        // qnode [10000,12000) +, then a PO overlap to a node sharing only a tiny
        // span < floor → inspect returns None. Build qnode -PO- small node where the
        // overlap window is below floor=max(147, 500-2*150=200)=200.
        let mut g = SdGraph::new();
        let q = nk_full("chr1", 10000, 12000, Strand::Forward);
        // cnode starts near the end of qnode so the overlap span is ~100bp (<200).
        let c = nk_full("chr1", 11900, 13000, Strand::Forward);
        add_edge_keys(&mut g, q.clone(), c.clone(), EdgeKind::Overlap, 2.0);
        let qv = g.index.get(&q).copied().unwrap();
        let cv = g.index.get(&c).copied().unwrap();
        let (_, edges) = dijkstra_route(&g, qv, cv).unwrap();
        let verts = route_vertices(&g, qv, &edges);
        let cnode = inspect_cnode_along_route(&g, &verts, &edges, &frag());
        assert!(cnode.is_none(), "tiny overlap window should drop the route");
    }
}
