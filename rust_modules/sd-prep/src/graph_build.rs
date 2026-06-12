//! Multiplex SD + physical-overlap (PO) graph — the `NodeKey`/`EdgeAttr`/`SdGraph`
//! types and `build_multiplex_graph`.
//!
//! Port of `preparation/graph_build.py::create_multiplex_graph` (l.83-144) +
//! `compose_PO_graph_per_chr` (l.35-81).
//!
//! ## What this builds (matching the Python control flow)
//!
//! 1. An **SD-edge** set, one undirected edge per filtered SD pair, weight =
//!    `mismatch_rate`, `kind = SegmentalDuplication` (Python l.102-106).
//! 2. A per-chromosome **directed PO graph**: for every pair of SD intervals that
//!    physically overlap on the same chrom, a PO edge **large → small**, weight
//!    `1/max(span/size_a, span/size_b)`, `overlap = true` (Python l.52-76, via an
//!    `IntervalTree`; here a `rust_lapper::Lapper`).
//! 3. The per-chr PO graphs are merged (union of nodes + edges) and the SD edges
//!    are overlaid: if a PO edge already connects the two SD nodes its `kind` is
//!    upgraded to `SegmentalDuplication`, else a fresh SD edge is added (l.127-135).
//!
//! ## Representation
//!
//! `petgraph::DiGraph<NodeKey, EdgeAttr>` with an `FxHashMap<NodeKey, NodeIndex>`
//! intern map (the Python tuple-node identity). The graph **owns** fresh node
//! copies; `build_multiplex_graph` borrows the SD-pair table read-only.
//!
//! Directedness matches Python: PO edges are directed large→small; SD edges are
//! logically undirected but stored on the same DiGraph (the traversal core treats
//! the graph as undirected via `component_labels`, exactly as Python does its
//! `to_undirected()` before `label_components`).

use petgraph::graph::{DiGraph, NodeIndex};
use rust_lapper::{Interval, Lapper};
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

/// Edge kind: an SD homology edge or a physical-overlap edge (Python `type ==
/// "segmental_duplication"` vs `overlap == "True"`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum EdgeKind {
    SegmentalDuplication,
    Overlap,
}

/// Edge attributes. `weight` is the SD `mismatch_rate` for SD edges, or `1/overlap_frac`
/// for PO edges (Python stores `weight=1/weight` where `weight=max(span/size)`).
/// `nonoverlap_smaller` is only meaningful for PO edges (Python
/// `nonoverlapping_smaller_interval`); `0` for SD edges.
#[derive(Clone, Copy, Debug)]
pub struct EdgeAttr {
    pub kind: EdgeKind,
    pub weight: f64,
    pub nonoverlap_smaller: i64,
}

impl EdgeAttr {
    pub fn is_overlap(&self) -> bool {
        matches!(self.kind, EdgeKind::Overlap)
    }
    pub fn is_sd(&self) -> bool {
        matches!(self.kind, EdgeKind::SegmentalDuplication)
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

    // ---- overlay SD edges (dedup by unordered pair key) ----
    // Python uses an undirected SD graph G then overlays onto the merged PO
    // DiGraph. We dedup SD pairs by the sorted (min,max) tuple key so the same SD
    // pair in either order produces one SD edge (matches the frozenset dedup the
    // driver does at l.175-176 before graph build — but we dedup here defensively
    // too, since build_multiplex_graph may receive an undeduped table).
    let mut sd_edge_seen: ahash::AHashSet<(NodeKey, NodeKey)> = ahash::AHashSet::new();
    for row in sd_pairs {
        let (lo, hi) = ordered_pair(&row.a, &row.b);
        if !sd_edge_seen.insert((lo.clone(), hi.clone())) {
            continue;
        }
        let u = sd.node(row.a.clone());
        let v = sd.node(row.b.clone());
        // If a PO edge already connects u and v (either direction), upgrade its
        // kind to SegmentalDuplication; else add a fresh SD edge u->v.
        let existing = sd
            .g
            .find_edge(u, v)
            .or_else(|| sd.g.find_edge(v, u));
        match existing {
            Some(e) => {
                let attr = sd.g.edge_weight_mut(e).expect("edge exists");
                attr.kind = EdgeKind::SegmentalDuplication;
                attr.weight = row.mismatch_rate;
            }
            None => {
                sd.g.add_edge(
                    u,
                    v,
                    EdgeAttr {
                        kind: EdgeKind::SegmentalDuplication,
                        weight: row.mismatch_rate,
                        nonoverlap_smaller: 0,
                    },
                );
            }
        }
    }

    log::info!(
        "Multiplex graph: {} nodes, {} edges ({} SD edges)",
        sd.node_count(),
        sd.edge_count(),
        sd.g
            .edge_weights()
            .filter(|a| a.is_sd())
            .count()
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
/// `compose_PO_graph_per_chr` (graph_build.py l.35-81) using a `Lapper` in place
/// of the Python `IntervalTree`.
///
/// The tree is grown incrementally exactly as Python does (an interval is only
/// added when it lies on `chrom`), so each interval only overlaps intervals seen
/// *earlier* in iteration — preserving the Python tie-break / direction semantics.
/// Because `Lapper` is a static structure (built once from a Vec), we replicate
/// the incremental behavior by collecting the on-chrom intervals first, then for
/// the i-th interval querying only intervals `< i` in insertion order. The
/// insertion order is the SD-pair-row order, matching Python's `iterrows()`.
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
    //
    // `inserted`: the on-chrom intervals added to the "tree" so far, with a
    // Lapper rebuilt lazily. For modest per-chrom counts a fresh small Lapper per
    // probe is fine; to stay O(n log n) we instead keep a growing Vec and do the
    // overlap query against it via a Lapper rebuilt only when needed. Simpler and
    // still correct: keep a Vec and binary-search by building the Lapper once at
    // the end is NOT possible (incremental visibility matters). So we keep the
    // growing Vec and query it directly with a linear scan guarded by start/end —
    // per-chrom SD counts are small (tens–hundreds), so this is not a hotspot.
    let mut inserted: Vec<NodeKey> = Vec::new();

    let mut probe = |iv: &NodeKey, sd: &mut SdGraph, inserted: &mut Vec<NodeKey>| {
        // Query overlaps among already-inserted on-chrom intervals.
        for o in inserted.iter() {
            // Same chrom guaranteed (inserted are all on `chrom`); skip identical.
            if o.chrom != iv.chrom || (o.start == iv.start && o.end == iv.end && o.strand == iv.strand) {
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
                    sd.g.add_edge(
                        lu,
                        sv,
                        EdgeAttr {
                            kind: EdgeKind::Overlap,
                            weight: 1.0 / weight,
                            nonoverlap_smaller: nonoverlap,
                        },
                    );
                }
            }
        }
        // Insert iv into the tree iff it lies on `chrom`.
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

/// Build a `Lapper<u32, NodeKey>` from a set of on-chrom intervals. Kept as a
/// thin helper so the per-chr PO build can switch to the static-index path for
/// large chromosomes (currently the incremental linear probe is used; see
/// `compose_po_per_chr`). Exposed (pub(crate)) so a future optimization or test
/// can build the same index the design specifies.
#[allow(dead_code)]
pub(crate) fn build_chr_lapper(intervals: &[NodeKey]) -> Lapper<u32, NodeKey> {
    let ivs: Vec<Interval<u32, NodeKey>> = intervals
        .iter()
        .map(|k| Interval {
            start: k.start.max(0) as u32,
            stop: k.end.max(0) as u32,
            val: k.clone(),
        })
        .collect();
    Lapper::new(ivs)
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
        let sd_edges = g.g.edge_weights().filter(|a| a.is_sd()).count();
        assert_eq!(sd_edges, 1);
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
        let l = g.index.get(&nk("chr1", 100, 1100, Strand::Forward)).copied().unwrap();
        let s = g.index.get(&nk("chr1", 600, 1000, Strand::Forward)).copied().unwrap();
        let e = g.g.find_edge(l, s).expect("PO edge L->S exists");
        let attr = g.g.edge_weight(e).unwrap();
        assert!(attr.is_overlap());
        assert!((attr.weight - 1.0).abs() < 1e-9, "weight {}", attr.weight);
        // No reverse PO edge.
        assert!(g.g.find_edge(s, l).is_none());
    }

    #[test]
    fn sd_overlay_upgrades_existing_po_edge() {
        // Make chr1 intervals overlap (creating a PO edge) AND be an SD pair.
        let l = nk("chr1", 100, 1100, Strand::Forward);
        let s = nk("chr1", 600, 1000, Strand::Forward);
        let rows = vec![SdPairRow {
            a: l.clone(),
            b: s.clone(),
            mismatch_rate: 0.02,
        }];
        let g = build_multiplex_graph(&rows, 1);
        // Only one edge between them, upgraded to SD.
        let lu = g.index.get(&l).copied().unwrap();
        let sv = g.index.get(&s).copied().unwrap();
        let e = g.g.find_edge(lu, sv).or_else(|| g.g.find_edge(sv, lu)).unwrap();
        let attr = g.g.edge_weight(e).unwrap();
        assert!(attr.is_sd(), "edge should be upgraded to SD");
        assert!((attr.weight - 0.02).abs() < 1e-9);
        assert_eq!(g.edge_count(), 1, "PO + SD on same pair collapse to one edge");
    }
}
