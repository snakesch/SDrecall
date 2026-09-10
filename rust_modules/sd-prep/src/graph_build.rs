//! Multiplex SD + physical-overlap (PO) graph — the `NodeKey`/`EdgeAttr`/`SdGraph`
//! types and `build_multiplex_graph`.
//!
//! Port of `preparation/graph_build.py::create_multiplex_graph` (l.83-144) +
//! `compose_PO_graph_per_chr` (l.35-81).
//!
//! ## What this builds (matching the Python control flow)
//!
//! 1. An **SD-edge** set, one undirected edge per filtered SD pair, weight =
//!    `mismatch_rate`, `is_sd = true` (Python l.102-106).
//! 2. A per-chromosome **directed PO graph**: for every pair of SD intervals that
//!    physically overlap on the same chrom, a PO edge **large → small**, weight
//!    `1/max(span/size_a, span/size_b)`, `is_overlap = true` (Python l.52-76, via
//!    an `IntervalTree`; here an insertion-order overlap sweep).
//! 3. The per-chr PO graphs are merged (union of nodes + edges) and the SD edges
//!    are overlaid onto them: an SD pair that already has a PO edge in the SAME
//!    direction KEEPS that edge's PO weight + `is_overlap` flag and is *also* marked
//!    `is_sd` (a "combined" edge); otherwise a fresh SD edge is added — even when
//!    only the *reverse* PO edge exists, in which case both edges coexist (Python
//!    l.127-135; see [`build_multiplex_graph`] for the directional overlay).
//!
//! ## Representation
//!
//! `petgraph::DiGraph<NodeKey, EdgeAttr>` with an `FxHashMap<NodeKey, NodeIndex>`
//! intern map (the Python tuple-node identity). The graph **owns** fresh node
//! copies; `build_multiplex_graph` borrows the SD-pair table read-only.
//!
//! An edge carries its physical-overlap and SD natures **independently**
//! ([`EdgeAttr::is_overlap`] / [`EdgeAttr::is_sd`]) so a single edge can be both —
//! matching the Python multiplex graph where the overlay leaves an upgraded PO
//! edge `overlap="True"` AND `type="segmental_duplication"` with the PO weight.
//!
//! Directedness matches Python: PO edges are directed large→small; SD edges are
//! logically undirected but stored on the same DiGraph (the traversal core treats
//! the graph as undirected via `component_labels`, exactly as Python does its
//! `to_undirected()` before `label_components`).

use petgraph::graph::{DiGraph, NodeIndex};
use rustc_hash::FxHashMap;
use sdrecall_utils::Strand;

/// Node identity — the Python 4-tuple `(chrom, start, end, strand)`. Owned
/// `String` chrom so a node outlives any single table row and hashes by value.
///
/// `Strand` here is the SD-map strand (`+`/`-`); `Unknown` should not occur for a
/// real SD node but is allowed so the type composes with `sdrecall-utils`.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct NodeKey {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
}

impl NodeKey {
    pub fn new(chrom: impl Into<String>, start: i64, end: i64, strand: Strand) -> Self {
        Self {
            chrom: chrom.into(),
            start,
            end,
            strand,
        }
    }

    /// Node span `end - start` (the Python `size` node attribute).
    pub fn size(&self) -> i64 {
        self.end - self.start
    }
}

/// Route-step kind — which branch the route walk took when crossing an edge
/// (`inspect_cnode_along_route`): an SD homology step or a physical-overlap step.
///
/// This is **no longer an `EdgeAttr` field**. Python dispatches the route walk by
/// `if type == "segmental_duplication" … elif overlap == "True" …` (so SD wins on a
/// combined edge) and records the branch taken on each `RouteStep`, so the
/// back-projection (`qnode_relative_region`) can replay it. An edge's *nature* is
/// the independent [`EdgeAttr::is_overlap`] / [`EdgeAttr::is_sd`] flags; this enum
/// only records the route-walk branch a `RouteStep` was reached by.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum EdgeKind {
    SegmentalDuplication,
    Overlap,
}

/// Edge attributes — an edge carries its physical-overlap (PO) and segmental-
/// duplication (SD) natures **independently**, matching the Python multiplex graph
/// where a single edge can be *both* a PO edge (`overlap == "True"`) **and** an SD
/// edge (`type == "segmental_duplication"`).
///
/// ## Why two flags + two weights (not one `kind` enum)
///
/// The Python SD overlay (`create_multiplex_graph` l.127-135) overlays SD edges
/// onto the merged PO graph: when an SD pair already has a PO edge it sets ONLY
/// `type = "segmental_duplication"` and **keeps the PO edge's `weight`
/// (`1/overlap_frac`) and `overlap="True"`**. So a "combined" edge ends up
/// `overlap="True"` + `type="segmental_duplication"` carrying the **PO** weight. A
/// single `kind` enum cannot represent that — it conflates the two natures and
/// loses either the overlap flag or the PO weight (the bug this struct fixes).
///
/// - [`is_overlap`](Self::is_overlap): Python `overlap == "True"`.
/// - [`is_sd`](Self::is_sd): Python `type == "segmental_duplication"`.
/// - [`po_weight`](Self::po_weight): the PO weight `1/overlap_frac` (Python stores
///   `weight = 1/max(span/size)`); `NaN` when `!is_overlap`.
/// - [`sd_weight`](Self::sd_weight): the SD `mismatch_rate` (kept available even on
///   a combined edge, where Python discards it); `NaN` when `!is_sd`.
/// - [`weight()`](Self::weight) reproduces Python's `ep["weight"]`: the PO weight
///   wins on a combined edge (because the overlay keeps it), else the SD weight.
/// - `nonoverlap_smaller`: Python `nonoverlapping_smaller_interval` (PO edges only;
///   `0` for a pure SD edge).
#[derive(Clone, Copy, Debug)]
pub struct EdgeAttr {
    pub is_overlap: bool,
    pub is_sd: bool,
    pub po_weight: f64,
    pub sd_weight: f64,
    pub nonoverlap_smaller: i64,
}

impl EdgeAttr {
    /// A fresh **PO** (physical-overlap) edge: `overlap="True"`, weight
    /// `1/overlap_frac`, not (yet) SD.
    pub fn po(po_weight: f64, nonoverlap_smaller: i64) -> Self {
        Self {
            is_overlap: true,
            is_sd: false,
            po_weight,
            sd_weight: f64::NAN,
            nonoverlap_smaller,
        }
    }

    /// A fresh **SD** (segmental-duplication) edge: `type="segmental_duplication"`,
    /// weight = `mismatch_rate`, not a PO edge.
    pub fn sd(mismatch_rate: f64) -> Self {
        Self {
            is_overlap: false,
            is_sd: true,
            po_weight: f64::NAN,
            sd_weight: mismatch_rate,
            nonoverlap_smaller: 0,
        }
    }

    /// Python `ep["weight"]` — the single weight every weight-consuming site uses
    /// (dijkstra cost, the SD-product route filter, the PO-edge prune). The **PO
    /// weight wins on a combined edge** (Python's overlay keeps the PO weight when
    /// it upgrades a PO edge to SD); otherwise the SD weight. Never reads a `NaN`
    /// slot: a PO/combined edge (`is_overlap`) reads `po_weight`, a pure-SD edge
    /// reads `sd_weight`.
    pub fn weight(&self) -> f64 {
        if self.is_overlap {
            self.po_weight
        } else {
            self.sd_weight
        }
    }
}

/// The multiplex graph: a `DiGraph` plus the node-intern map.
pub struct SdGraph {
    pub g: DiGraph<NodeKey, EdgeAttr>,
    pub index: FxHashMap<NodeKey, NodeIndex>,
}

impl SdGraph {
    pub fn new() -> Self {
        Self {
            g: DiGraph::new(),
            index: FxHashMap::default(),
        }
    }

    /// Intern a node by key, returning its (existing or freshly created) index.
    pub fn node(&mut self, key: NodeKey) -> NodeIndex {
        if let Some(&idx) = self.index.get(&key) {
            return idx;
        }
        let idx = self.g.add_node(key.clone());
        self.index.insert(key, idx);
        idx
    }

    /// True iff `key` is already a node in the graph.
    pub fn has_node(&self, key: &NodeKey) -> bool {
        self.index.contains_key(key)
    }

    pub fn node_count(&self) -> usize {
        self.g.node_count()
    }

    pub fn edge_count(&self) -> usize {
        self.g.edge_count()
    }
}

impl Default for SdGraph {
    fn default() -> Self {
        Self::new()
    }
}

/// One row of the filtered SD-pair table (`filtered_SD_binary_map.tsv`): the two
/// SD intervals plus the pair's `mismatch_rate`. Mirrors the columns the Python
/// driver keeps after umbrella filtering (`prepare_recall_regions.py` l.168-170).
#[derive(Clone, Debug)]
pub struct SdPairRow {
    pub a: NodeKey,
    pub b: NodeKey,
    pub mismatch_rate: f64,
}

/// Build the multiplex graph (SD edges + per-chr PO edges + overlay), matching
/// `create_multiplex_graph`.
///
/// `&[SdPairRow]` read-only over the filtered SD table; the graph owns fresh node
/// copies. The PO build is **per chromosome** (Python parallelizes this over a
/// Pool; here it is a sequential per-chrom sweep — parallelizing it with rayon is
/// straightforward but the chrom count is tiny vs the traversal hotspot, so it is
/// left sequential for determinism and simplicity; see T8 perf notes).
///
/// Parity-critical points reproduced from `compose_PO_graph_per_chr`:
/// - PO edges are drawn between every overlapping pair of SD intervals **on the
///   same chrom**, irrespective of overlap extent (l.50, l.63-76).
/// - The edge points from the **larger** interval to the **smaller** (l.72-76);
///   on a size tie (`interval_size >= o_size`), the *current* interval is treated
///   as the larger (Python's `>=`), so direction follows insertion order on ties.
/// - Weight = `1 / max(span/o_size, span/interval_size)` (l.70, l.73).
/// - `nonoverlapping_smaller = min(o_size, interval_size) - span` (l.71).
/// - An interval is only **added to the tree** when `chrom == interval.chrom`
///   (l.59-61), but overlap is queried against *both* intervals of every pair
///   (the second interval may be on another chrom and still probe the tree).
pub fn build_multiplex_graph(sd_pairs: &[SdPairRow], _threads: usize) -> SdGraph {
    let mut sd = SdGraph::new();

    // ---- collect the set of chromosomes touched by any SD interval ----
    let mut chroms: Vec<String> = Vec::new();
    {
        let mut seen = ahash::AHashSet::new();
        for row in sd_pairs {
            for c in [&row.a.chrom, &row.b.chrom] {
                if seen.insert(c.clone()) {
                    chroms.push(c.clone());
                }
            }
        }
    }

    // ---- build the merged PO graph (over all chroms) into `sd` ----
    // We dedup PO edges by an ordered (large_node, small_node) key so the same
    // directed overlap is not added twice when two chroms both visit a pair.
    let mut po_edge_seen: ahash::AHashSet<(NodeKey, NodeKey)> = ahash::AHashSet::new();
    for chrom in &chroms {
        compose_po_per_chr(chrom, sd_pairs, &mut sd, &mut po_edge_seen);
    }

    // ---- overlay SD edges onto the merged PO graph (Python l.127-135) ----
    //
    //     for u, v, d in G.edges(data=True):          # G = undirected SD graph
    //         if not PO.has_edge(u, v): PO.add_edge(u, v, **d)               # fresh SD
    //         else:                     PO[u][v]["type"]="segmental_duplication"  # upgrade
    //
    // Two parity-critical details, both previously WRONG in Rust:
    //   1. DIRECTIONALITY — `has_edge(u, v)` probes the SPECIFIC (u, v) direction.
    //      If only the *reverse* PO edge (v, u) exists, Python ADDS a new SD edge
    //      (u, v) and leaves the reverse PO edge intact (the two edges coexist).
    //      So we check `find_edge(u, v)` ONLY — NOT `find_edge(u,v).or(find_edge(v,u))`.
    //   2. PRESERVATION — upgrading an existing PO edge sets ONLY `type`; it KEEPS
    //      the PO `weight` (1/overlap_frac) + `overlap="True"` + nonoverlap. So we
    //      set `is_sd = true` and DO NOT touch `po_weight`/`is_overlap`/nonoverlap
    //      (the old code overwrote the weight + collapsed both directions).
    //
    // (u, v) ORIENTATION: Python yields each SD edge from the undirected `G` in
    // (earlier-inserted, later-inserted) node order (networkx adjacency iteration);
    // `G`'s nodes are inserted by `G.add_edge(row.a, row.b)` in row order (a before
    // b). We reproduce that orientation via each node's first-appearance index so
    // `find_edge` probes the SAME direction Python's `has_edge` does — that is what
    // decides "combine onto the PO edge" vs "add a reciprocal SD edge".
    let mut first_seen: ahash::AHashMap<NodeKey, usize> = ahash::AHashMap::new();
    {
        let mut order = 0usize;
        for row in sd_pairs {
            for k in [&row.a, &row.b] {
                if !first_seen.contains_key(k) {
                    first_seen.insert(k.clone(), order);
                    order += 1;
                }
            }
        }
    }

    // Dedup each unordered SD pair once (Python's undirected `G` has one edge per
    // pair); the sorted (min,max) key is canonical regardless of orientation.
    let mut sd_edge_seen: ahash::AHashSet<(NodeKey, NodeKey)> = ahash::AHashSet::new();
    for row in sd_pairs {
        let (lo, hi) = ordered_pair(&row.a, &row.b);
        if !sd_edge_seen.insert((lo.clone(), hi.clone())) {
            continue;
        }
        // Orient (u_key → v_key) by first-appearance order (networkx `G` edge order).
        let (u_key, v_key) = if first_seen[&row.a] <= first_seen[&row.b] {
            (row.a.clone(), row.b.clone())
        } else {
            (row.b.clone(), row.a.clone())
        };
        let u = sd.node(u_key);
        let v = sd.node(v_key);
        // SPECIFIC (u, v) direction only (Python `has_edge(u, v)`).
        match sd.g.find_edge(u, v) {
            Some(e) => {
                // A PO edge already runs u→v: upgrade it to ALSO be SD, KEEPING its
                // PO weight + overlap flag + nonoverlap (Python sets only `type`).
                // We additionally stash the SD `mismatch_rate` in `sd_weight` so it
                // stays available, though `weight()` (= ep["weight"]) keeps using
                // the PO weight on this combined edge, exactly as Python does.
                let attr = sd.g.edge_weight_mut(e).expect("edge exists");
                attr.is_sd = true;
                attr.sd_weight = row.mismatch_rate;
            }
            None => {
                // No u→v edge: add a fresh SD edge (Python `add_edge(u, v, **d)`).
                // A reverse PO edge (v, u), if any, is intentionally left intact.
                sd.g.add_edge(u, v, EdgeAttr::sd(row.mismatch_rate));
            }
        }
    }

    log::info!(
        "Multiplex graph: {} nodes, {} edges ({} SD edges)",
        sd.node_count(),
        sd.edge_count(),
        sd.g.edge_weights().filter(|a| a.is_sd).count()
    );
    sd
}

/// Sorted `(min, max)` key for unordered-pair dedup (the Rust analog of the
/// Python `frozenset`). Compares on `(chrom, start, end, strand-rank)`.
fn ordered_pair(a: &NodeKey, b: &NodeKey) -> (NodeKey, NodeKey) {
    if node_rank(a) <= node_rank(b) {
        (a.clone(), b.clone())
    } else {
        (b.clone(), a.clone())
    }
}

fn node_rank(n: &NodeKey) -> (&str, i64, i64, u8) {
    let s = match n.strand {
        Strand::Forward => 1,
        Strand::Reverse => 2,
        Strand::Unknown => 0,
    };
    (n.chrom.as_str(), n.start, n.end, s)
}

/// Build PO edges for one chromosome and fold them into `sd`. Ports
/// `compose_PO_graph_per_chr` (graph_build.py l.35-81) using an incremental
/// overlap sweep in place of the Python `IntervalTree`.
///
/// Intervals are processed in Python row order and added only when they lie on
/// `chrom`, so each probe sees only intervals inserted earlier. This preserves
/// Python's tie-break and edge-direction semantics. The direct scan is linear per
/// probe; per-chromosome SD counts are modest and incremental visibility matters.
fn compose_po_per_chr(
    chrom: &str,
    sd_pairs: &[SdPairRow],
    sd: &mut SdGraph,
    po_edge_seen: &mut ahash::AHashSet<(NodeKey, NodeKey)>,
) {
    // Reproduce Python's iteration: for each row touching this chrom, process
    // interval_1 then interval_2; an interval is only tree-eligible when it lies
    // on `chrom`. We accumulate (insertion_order, NodeKey) for on-chrom intervals
    // and, for every probe interval (whether on-chrom or not), query the
    // already-inserted on-chrom intervals for overlap.
    let mut inserted: Vec<NodeKey> = Vec::new();

    let mut probe = |iv: &NodeKey, sd: &mut SdGraph, inserted: &mut Vec<NodeKey>| {
        // Query overlaps among already-inserted on-chrom intervals.
        for o in inserted.iter() {
            // Same chrom guaranteed (inserted are all on `chrom`); skip identical.
            if o.chrom != iv.chrom
                || (o.start == iv.start && o.end == iv.end && o.strand == iv.strand)
            {
                continue;
            }
            // Half-open overlap.
            if iv.start < o.end && o.start < iv.end {
                let span = iv.end.min(o.end) - iv.start.max(o.start);
                if span <= 0 {
                    continue;
                }
                let o_size = (o.end - o.start) as f64;
                let iv_size = (iv.end - iv.start) as f64;
                let span_f = span as f64;
                let weight = (span_f / o_size).max(span_f / iv_size);
                let nonoverlap = (o_size.min(iv_size) as i64) - span;
                // Direction: large -> small. Python: if interval_size >= o_size,
                // edge interval -> o_interval (iv is the larger / equal).
                let (large, small) = if iv_size >= o_size {
                    (iv.clone(), o.clone())
                } else {
                    (o.clone(), iv.clone())
                };
                if po_edge_seen.insert((large.clone(), small.clone())) {
                    let lu = sd.node(large);
                    let sv = sd.node(small);
                    sd.g.add_edge(lu, sv, EdgeAttr::po(1.0 / weight, nonoverlap));
                }
            }
        }
        // Make this interval visible to later probes iff it lies on `chrom`.
        if iv.chrom == chrom {
            inserted.push(iv.clone());
        }
    };

    for row in sd_pairs {
        if row.a.chrom != chrom && row.b.chrom != chrom {
            continue;
        }
        // Ensure both nodes exist as graph nodes (Python adds both nodes per row).
        sd.node(row.a.clone());
        sd.node(row.b.clone());
        probe(&row.a, sd, &mut inserted);
        probe(&row.b, sd, &mut inserted);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nk(chrom: &str, start: i64, end: i64, s: Strand) -> NodeKey {
        NodeKey::new(chrom, start, end, s)
    }

    #[test]
    fn intern_returns_same_index() {
        let mut g = SdGraph::new();
        let a = g.node(nk("chr1", 10, 20, Strand::Forward));
        let b = g.node(nk("chr1", 10, 20, Strand::Forward));
        assert_eq!(a, b);
        assert_eq!(g.node_count(), 1);
    }

    #[test]
    fn sd_edge_added_for_pair() {
        let rows = vec![SdPairRow {
            a: nk("chr1", 100, 1100, Strand::Forward),
            b: nk("chr2", 200, 1200, Strand::Reverse),
            mismatch_rate: 0.03,
        }];
        let g = build_multiplex_graph(&rows, 1);
        assert_eq!(g.node_count(), 2);
        assert_eq!(g.edge_count(), 1);
        let sd_edges = g.g.edge_weights().filter(|a| a.is_sd).count();
        assert_eq!(sd_edges, 1);
        // A pure inter-chrom SD edge: is_sd, NOT overlap, weight() == mismatch_rate.
        let attr = g.g.edge_weights().next().unwrap();
        assert!(attr.is_sd && !attr.is_overlap);
        assert!((attr.weight() - 0.03).abs() < 1e-9);
        assert!((attr.sd_weight - 0.03).abs() < 1e-9);
        assert!(attr.po_weight.is_nan(), "pure SD edge has no PO weight");
    }

    #[test]
    fn po_edge_large_to_small_on_same_chrom() {
        // Two pairs both placing overlapping intervals on chr1.
        // Interval L = [100, 1100) (size 1000), S = [600, 1000) (size 400), overlap [600,1000)=400.
        // L is larger → PO edge L -> S, weight = 1/max(400/400, 400/1000) = 1/1.0 = 1.0.
        let rows = vec![
            SdPairRow {
                a: nk("chr1", 100, 1100, Strand::Forward),
                b: nk("chrX", 5000, 6000, Strand::Forward),
                mismatch_rate: 0.05,
            },
            SdPairRow {
                a: nk("chr1", 600, 1000, Strand::Forward),
                b: nk("chrY", 7000, 7400, Strand::Forward),
                mismatch_rate: 0.05,
            },
        ];
        let g = build_multiplex_graph(&rows, 1);
        // Find the PO edge among the chr1 nodes.
        let l = g
            .index
            .get(&nk("chr1", 100, 1100, Strand::Forward))
            .copied()
            .unwrap();
        let s = g
            .index
            .get(&nk("chr1", 600, 1000, Strand::Forward))
            .copied()
            .unwrap();
        let e = g.g.find_edge(l, s).expect("PO edge L->S exists");
        let attr = g.g.edge_weight(e).unwrap();
        assert!(attr.is_overlap);
        assert!(!attr.is_sd, "pure PO edge is not (yet) an SD edge");
        assert!(
            (attr.po_weight - 1.0).abs() < 1e-9,
            "po_weight {}",
            attr.po_weight
        );
        // ep["weight"] == PO weight for a pure PO edge.
        assert!((attr.weight() - 1.0).abs() < 1e-9);
        assert!(attr.sd_weight.is_nan(), "pure PO edge has no SD weight");
        // No reverse PO edge.
        assert!(g.g.find_edge(s, l).is_none());
    }

    // ── FIX #8: SD overlay onto a PO edge → combined edge (PO weight + flag kept) ──

    #[test]
    fn sd_overlay_combines_with_existing_po_edge_keeping_po_weight() {
        // L = [100, 1100) (size 1000) and S = [600, 1000) (size 400) overlap on chr1
        // (PO edge L->S, po_weight = 1/max(400/400,400/1000) = 1.0) AND are an SD
        // pair (mismatch_rate 0.02). first_seen: L (row.a) before S (row.b), so the
        // overlay probes find_edge(L, S) — the PO edge's OWN direction → COMBINE.
        let l = nk("chr1", 100, 1100, Strand::Forward);
        let s = nk("chr1", 600, 1000, Strand::Forward);
        let rows = vec![SdPairRow {
            a: l.clone(),
            b: s.clone(),
            mismatch_rate: 0.02,
        }];
        let g = build_multiplex_graph(&rows, 1);

        // Still exactly ONE edge between them (overlay upgraded, did not duplicate).
        assert_eq!(g.edge_count(), 1, "combine keeps a single edge");
        let lu = g.index.get(&l).copied().unwrap();
        let sv = g.index.get(&s).copied().unwrap();
        let e = g.g.find_edge(lu, sv).expect("combined edge L->S");
        let attr = g.g.edge_weight(e).unwrap();

        // The required end state for an overlapping SD pair:
        assert!(attr.is_overlap, "overlap flag KEPT");
        assert!(attr.is_sd, "ALSO marked SD");
        assert!(
            (attr.po_weight - 1.0).abs() < 1e-9,
            "PO weight 1/overlap_frac preserved"
        );
        assert!(
            (attr.sd_weight - 0.02).abs() < 1e-9,
            "SD mismatch_rate kept available"
        );
        // ep["weight"] = PO weight on a combined edge (Python keeps the PO weight).
        assert!(
            (attr.weight() - 1.0).abs() < 1e-9,
            "ep[weight] is the PO weight, not 0.02"
        );
        // No reverse edge was created.
        assert!(g.g.find_edge(sv, lu).is_none());
    }

    #[test]
    fn sd_overlay_adds_reciprocal_edge_when_only_reverse_po_exists() {
        // Same overlapping intervals, but the SD row is oriented SMALL→LARGE
        // (a = S, b = L). The PO edge is still LARGE→SMALL (L->S), so the overlay
        // probes find_edge(S, L) — the *reverse* direction → MISS → a fresh SD edge
        // S->L is added and the reverse PO edge L->S is left intact (two coexist),
        // matching Python's directional `has_edge(u, v)` overlay.
        let l = nk("chr1", 100, 1100, Strand::Forward); // size 1000 (larger)
        let s = nk("chr1", 600, 1000, Strand::Forward); // size 400  (smaller)
        let rows = vec![SdPairRow {
            a: s.clone(), // SD-edge orientation small→large (first-seen: s before l)
            b: l.clone(),
            mismatch_rate: 0.02,
        }];
        let g = build_multiplex_graph(&rows, 1);

        let lu = g.index.get(&l).copied().unwrap();
        let sv = g.index.get(&s).copied().unwrap();

        // TWO edges coexist: PO L->S and a reciprocal SD S->L.
        assert_eq!(
            g.edge_count(),
            2,
            "reciprocal SD edge coexists with the reverse PO edge"
        );

        // PO edge L->S is intact (overlap, NOT SD, PO weight 1.0).
        let po =
            g.g.edge_weight(g.g.find_edge(lu, sv).expect("PO L->S"))
                .unwrap();
        assert!(po.is_overlap && !po.is_sd, "reverse PO edge untouched");
        assert!((po.po_weight - 1.0).abs() < 1e-9);

        // Fresh SD edge S->L (SD, NOT overlap, weight() == mismatch_rate).
        let sd_e =
            g.g.edge_weight(g.g.find_edge(sv, lu).expect("SD S->L"))
                .unwrap();
        assert!(sd_e.is_sd && !sd_e.is_overlap, "fresh reciprocal SD edge");
        assert!((sd_e.weight() - 0.02).abs() < 1e-9);
    }
}
