//! `sd-prep` (SDrecall T8) — Phase-1 preparation (the graph-tool replacement).
//!
//! Replaces `prepare_recall_regions.py` + the `preparation/` package: SD-graph
//! traversal, RG grouping via vertex coloring, multiplex graph build, umbrella
//! SD-pair filtering, multi-align depth, masked genomes, intrinsic alignment.
//!
//! ## Build order / status (see `docs/analysis/tasks/T8_sd_prep.md`)
//!
//! Code-complete: the graph core (the deliverable), the interval logic, the
//! route-walk traversal + minimap2 FFI, the masked-genome + intrinsic I/O, and the
//! end-to-end `driver`. Validated on HG002 (4 RG groupings matching Python, all 4
//! masked FASTAs md5-identical, paralog-pair set 96% concordant). Per-module status
//! is documented at the top of each module.
//!
//! ## Units (one versatile unit per job)
//!
//! - [`graph_core`] — ★ the graph-tool replacement: `component_labels`
//!   (UnionFind, all-edge OR overlap-only via one predicate), `dijkstra_route`
//!   (predecessor-recording weighted Dijkstra) and `greedy_vertex_coloring`
//!   (gt.sequential_vertex_coloring index-order parity). UNIT-TESTED.
//! - [`graph_build`] — `NodeKey`/`EdgeAttr`/`SdGraph` + `build_multiplex_graph`
//!   (per-chr PO via Lapper + SD overlay + dedup). UNIT-TESTED.
//! - [`sd_pairs`] — `Pair` + the ONE `umbrella_to_remove` sweep (raw & granular).
//!   UNIT-TESTED.
//! - [`multialign`] — the ONE per-base depth sweep kernel + 4-way AND filter.
//!   UNIT-TESTED (sweep + filter); BAM read STUBBED.
//! - [`homoseq`] — `HomoseqRegion` route-carrying coordinate object +
//!   `qnode_relative_region` back-projection. UNIT-TESTED.
//! - [`grouping`] — `optimal_node_grouping` fast path (color → RG clusters).
//!   UNIT-TESTED.
//! - [`traversal`] — the route walk (`inspect_cnode_along_route`) + minimap2
//!   similarity validation + `ConnectedQnodes` insertion-order wiring. UNIT-TESTED.
//! - [`minimap`] — the ONE `align_similarity` FFI wrapper (the only external
//!   algorithm); see its dependency-availability + version-skew note.
//! - [`masking`] — `mask_genome` (faidx fetch + N-mask + N-bridge merge + wrapped
//!   FASTA). UNIT-TESTED; md5-identical to Python on HG002.
//! - [`intrinsic`] — counterpart→masked-genome minimap2 alignment → filtered
//!   intrinsic BAM. UNIT-TESTED.
//! - [`driver`] — the 7-step `prepare_recall_regions` orchestration. UNIT-TESTED;
//!   HG002 end-to-end differential in `examples/validate_phase1_e2e.rs`.

pub mod driver;
pub mod graph_build;
pub mod graph_core;
pub mod grouping;
pub mod homoseq;
pub mod intrinsic;
pub mod masking;
pub mod minimap;
pub mod multialign;
pub mod sd_pairs;
pub mod traversal;

pub use driver::{prepare_recall_regions, PrepParams, PrepPaths, PrepResult, RgOutputs};
pub use graph_build::{build_multiplex_graph, EdgeAttr, EdgeKind, NodeKey, SdGraph};
pub use graph_core::{component_labels, dijkstra_route, greedy_vertex_coloring};
pub use grouping::{optimal_node_grouping, ConnectedQnodes};
pub use homoseq::HomoseqRegion;
pub use intrinsic::{intrinsic_bam, merge_total_intrinsic_bam};
pub use masking::mask_genome;
pub use minimap::{align_similarity, Preset, SIMILARITY_THRESHOLD};
pub use multialign::{
    depth_sweep, inferred_coverage, multialign_filter_mask, pick_multialigned_regions, DepthPass,
};
pub use sd_pairs::{umbrella_to_remove, Pair};
pub use traversal::{extract_sd_paralog_pairs, Counterpart, TraversalResult};
