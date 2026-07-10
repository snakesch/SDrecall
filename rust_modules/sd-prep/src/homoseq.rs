//! `HomoseqRegion` — the route-carrying coordinate object (Python
//! `preparation/homoseq_region.py::HOMOSEQ_REGION`).
//!
//! Each `HomoseqRegion` carries a node's identity plus the relative-coordinate
//! window derived while walking a shortest path from a query node (qnode) to this
//! counterpart node (cnode), and the `route` of `(node, edge_kind)` steps that
//! produced it. `qnode_relative_region` back-projects this cnode's window onto a
//! given upstream qnode by replaying the route — the load-bearing coordinate math
//! the BED/masked-genome outputs depend on.
//!
//! ## Port status
//!
//! The struct + `qnode_relative_region` (the pure coordinate math) are ported and
//! UNIT-TESTED. The route is built by [`crate::traversal`] (STUBBED) — see that
//! module. The `traverse_route` here is `Vec<RouteStep>` where each step owns a
//! snapshot of the upstream node it was reached from (matching Python's
//! `(HOMOSEQ_REGION, edge_type)` tuples, which the back-projection reads for
//! `start`/`end`/`strand`).

use crate::graph_build::{EdgeKind, NodeKey};
use petgraph::graph::NodeIndex;
use sdrecall_utils::Strand;

/// One step of a traverse route: the upstream node's coordinate snapshot and the
/// edge kind crossed to leave it. Mirrors Python's `(HOMOSEQ_REGION, edge_type)`
/// tuple — the back-projection only needs the upstream node's `start`/`end`/
/// `strand` (its absolute coordinates) plus the edge kind, so we store a light
/// snapshot rather than a full recursive `HomoseqRegion`.
#[derive(Clone, Debug)]
pub struct RouteStep {
    pub node: NodeKey,
    pub edge_kind: EdgeKind,
}

/// The route-carrying coordinate object. Ports `HOMOSEQ_REGION`'s fields:
/// `chrom/start/end/strand/size`, the three relative-window pairs
/// (`ups_*`, `down_*`, `rela_*`), the graph `vertex` index, and the owned `route`.
///
/// Owned `route` Vec (the design's choice): each accepted counterpart needs an
/// independent route snapshot; borrowing would tie its lifetime to the traversal
/// scratch state, which is mutated per shortest-path probe.
#[derive(Clone, Debug)]
pub struct HomoseqRegion {
    pub key: NodeKey,
    pub size: i64,
    pub ups_rela_start: i64,
    pub ups_rela_end: i64,
    pub down_rela_start: i64,
    pub down_rela_end: i64,
    pub rela_start: i64,
    pub rela_end: i64,
    pub vertex: NodeIndex,
    pub route: Vec<RouteStep>,
}

impl HomoseqRegion {
    /// Construct from a node key + graph vertex, with the initial relative window
    /// equal to the whole node (`HOMOSEQ_REGION.__init__`, homoseq_region.py
    /// l.7-28): `ups_rela_start=0`, `ups_rela_end=size`, `down=0`, `rela=ups`.
    pub fn new(key: NodeKey, vertex: NodeIndex) -> Self {
        let size = key.end - key.start;
        Self {
            key,
            size,
            ups_rela_start: 0,
            ups_rela_end: size,
            down_rela_start: 0,
            down_rela_end: 0,
            rela_start: 0,
            rela_end: size,
            vertex,
            route: Vec::new(),
        }
    }

    pub fn chrom(&self) -> &str {
        &self.key.chrom
    }
    pub fn start(&self) -> i64 {
        self.key.start
    }
    pub fn end(&self) -> i64 {
        self.key.end
    }
    pub fn strand(&self) -> Strand {
        self.key.strand
    }

    /// The fixed-coordinate 4-tuple the BED writer emits via Python's `__iter__`:
    /// `(chrom, start + ups_rela_start, start + ups_rela_end, strand)`
    /// (homoseq_region.py l.39-41, l.53-54).
    pub fn fix_coord(&self) -> (String, i64, i64, Strand) {
        (
            self.key.chrom.clone(),
            self.key.start + self.ups_rela_start,
            self.key.start + self.ups_rela_end,
            self.key.strand,
        )
    }

    /// Back-project this cnode's relative window onto an upstream qnode by
    /// replaying the route from that qnode forward (`qnode_relative_region`,
    /// homoseq_region.py l.56-104).
    ///
    /// Returns `Some((rela_start, rela_end))` giving the qnode-relative window, or
    /// `None` when the projection collapses to a non-positive span (Python returns
    /// `("NaN","NaN")`).
    ///
    /// The Python algorithm: find `qnode_tuple` in `self.route`, take the
    /// **suffix** of the route from that qnode onward, then iterate that suffix in
    /// **reverse** (`.pop()` from the end), folding the current relative window
    /// back through each upstream node. SD steps clip/flip by strand; overlap
    /// steps re-anchor by absolute coordinate. We reproduce this with the route
    /// suffix popped from the back.
    pub fn qnode_relative_region(&self, qnode: &NodeKey) -> Option<(i64, i64)> {
        // Find the qnode position in the route (by node identity).
        let pos = self.route.iter().position(|step| &step.node == qnode)?;
        // Suffix from qnode onward; iterate in reverse (Python pop() from end).
        let suffix = &self.route[pos..];

        // current_node starts as self; its window is rela_start/rela_end.
        // We track current_node's (start, end, strand) explicitly.
        let mut cur_start = self.key.start;
        let mut cur_end = self.key.end;
        let mut cur_strand = self.key.strand;
        let mut rela_start = self.rela_start;
        let mut rela_end = self.rela_end;
        if rela_end <= rela_start {
            // Python asserts rela_end > rela_start before the loop.
            return None;
        }

        for step in suffix.iter().rev() {
            let up = &step.node;
            let up_size = up.end - up.start;
            match step.edge_kind {
                EdgeKind::SegmentalDuplication => {
                    if up.strand == cur_strand {
                        rela_start = rela_start.max(0);
                        rela_end = rela_end.min(up_size);
                    } else {
                        let qsize = rela_end - rela_start;
                        let cur_size = cur_end - cur_start;
                        rela_end = (cur_size - rela_start).min(up_size);
                        rela_start = (rela_end - qsize).max(0);
                    }
                }
                EdgeKind::Overlap => {
                    let qnode_abs_start = rela_start + cur_start;
                    let qsize = rela_end - rela_start;
                    rela_start = qnode_abs_start - up.start;
                    rela_end = (rela_start + qsize).min(up_size);
                    rela_start = rela_start.max(0);
                }
            }
            if rela_end - rela_start <= 0 {
                return None;
            }
            cur_start = up.start;
            cur_end = up.end;
            cur_strand = up.strand;
        }
        Some((rela_start, rela_end))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nk(chrom: &str, start: i64, end: i64, s: Strand) -> NodeKey {
        NodeKey::new(chrom, start, end, s)
    }

    #[test]
    fn new_initializes_full_window() {
        let k = nk("chr1", 1000, 2000, Strand::Forward);
        let h = HomoseqRegion::new(k.clone(), NodeIndex::new(0));
        assert_eq!(h.size, 1000);
        assert_eq!((h.ups_rela_start, h.ups_rela_end), (0, 1000));
        assert_eq!((h.rela_start, h.rela_end), (0, 1000));
    }

    #[test]
    fn fix_coord_applies_ups_window() {
        let k = nk("chr1", 1000, 2000, Strand::Forward);
        let mut h = HomoseqRegion::new(k, NodeIndex::new(0));
        h.ups_rela_start = 100;
        h.ups_rela_end = 400;
        let (c, s, e, st) = h.fix_coord();
        assert_eq!(
            (c.as_str(), s, e, st),
            ("chr1", 1100, 1400, Strand::Forward)
        );
    }

    #[test]
    fn back_project_single_sd_same_strand_clips() {
        // qnode [10000,11000) +, cnode [20000,21000) + reached via ONE SD edge.
        // cnode rela window [200,500); back-project to qnode: same strand →
        // rela_start=max(200,0)=200, rela_end=min(500, qsize=1000)=500.
        let qnode = nk("chr1", 10000, 11000, Strand::Forward);
        let cnode = nk("chr2", 20000, 21000, Strand::Forward);
        let mut h = HomoseqRegion::new(cnode, NodeIndex::new(1));
        h.rela_start = 200;
        h.rela_end = 500;
        h.route = vec![RouteStep {
            node: qnode.clone(),
            edge_kind: EdgeKind::SegmentalDuplication,
        }];
        assert_eq!(h.qnode_relative_region(&qnode), Some((200, 500)));
    }

    #[test]
    fn back_project_single_sd_opposite_strand_flips() {
        // cnode size 1000, rela [200,500) (qsize 300). qnode size 1000, opp strand.
        // rela_end = min(cur_size - rela_start, up_size) = min(1000-200,1000)=800.
        // rela_start = max(rela_end - qsize, 0) = max(800-300,0)=500.
        let qnode = nk("chr1", 10000, 11000, Strand::Reverse);
        let cnode = nk("chr2", 20000, 21000, Strand::Forward);
        let mut h = HomoseqRegion::new(cnode, NodeIndex::new(1));
        h.rela_start = 200;
        h.rela_end = 500;
        h.route = vec![RouteStep {
            node: qnode.clone(),
            edge_kind: EdgeKind::SegmentalDuplication,
        }];
        assert_eq!(h.qnode_relative_region(&qnode), Some((500, 800)));
    }

    #[test]
    fn back_project_qnode_not_in_route_is_none() {
        let qnode = nk("chr1", 10000, 11000, Strand::Forward);
        let other = nk("chrZ", 1, 2, Strand::Forward);
        let cnode = nk("chr2", 20000, 21000, Strand::Forward);
        let mut h = HomoseqRegion::new(cnode, NodeIndex::new(1));
        h.rela_start = 200;
        h.rela_end = 500;
        h.route = vec![RouteStep {
            node: qnode,
            edge_kind: EdgeKind::SegmentalDuplication,
        }];
        assert!(h.qnode_relative_region(&other).is_none());
    }
}
