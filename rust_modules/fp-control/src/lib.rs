//! `fp-control` (SDrecall T4) — fused Phase-2c FP control.
//!
//! Threads `build_phasing_graph` → `phasing` → `haplotype_inspection` into one
//! in-process Rust call
//! `(bam, intrinsic_bam, region/params) → (correct_qnames, mismap_qnames)`,
//! ending the current Rust→Python→Rust round-trip.
//!
//! ## What this replaces
//!
//! The orchestration in `fp_control/realign_filter_per_cov.py:220-346`, which
//! today: calls the Rust `build_phasing_graph` (ships a graph + dense matrix +
//! 7 dicts out to Python), runs Python `phasing_realigned_reads`, then calls the
//! Rust `inspect_haplotypes_rust` (which re-opens the BAM). The fuse removes the
//! two PyO3 boundaries and the Python phasing hop — the three stage libraries
//! become one in-process pipeline passing Rust structs.
//!
//! ## Data flow (study `realign_filter_per_cov.py` for the exact wiring)
//!
//! 1. `build_phasing_graph` (graph_builder) → `PhasingGraphResult`:
//!    `weight_matrix` (N×N `Array2<f32>`), graph edges, `node_read_ids`
//!    (per-vertex read ids, keyed `"{qname}:{flag}"`), `read_hap_vectors` /
//!    `read_error_vectors` (keyed by the same read id), `lowqual_qnames`.
//!    Vertex index == `qname_idx` (petgraph's sequential-assignment guarantee).
//! 2. `phasing::phase` → `vertex_hap: HashMap<i32, i32>` (== Python's
//!    `qname_hap_info`, keyed by vertex index). From it we derive the two maps
//!    Python materialised: `hap_qname_info` (hap_id → [qname]) and
//!    `qname_to_node` (qname → vertex index).
//! 3. `haplotype_inspection::inspect_haplotypes` → `(correct, mismap)`.
//!
//! ## Scope (T4 first pass — deferred items)
//!
//! * **PERF pass (PERF-1/2/3) deferred** — sequential inspect kept; no
//!   parallel/`Arc` rewrite yet. Correctness first.
//! * **Golden-encoding switch (DivA/DivB) deferred** — crates used as-is.
//! * **Single-shared-BAM-read deferred** — each stage opens the BAM itself, as
//!   today. `inspect_haplotypes` takes `bam_path`/`intrinsic_bam_path` and
//!   re-opens internally; collapsing to one read needs a new shared-Lapper entry
//!   on `haplotype_inspection` (a blocker — see report), so it is not attempted
//!   here. BAM-open count is therefore 2 (graph + inspect) plus htslib's own
//!   index/header touches, not 1.

use std::collections::{HashMap, HashSet};

use phasing::bam_reading::{build_allele_depth_map, migrate_bam_to_sorted_intervals_grouped};
use phasing::graph_builder::build_phasing_graph as build_graph;
use phasing::structs::HaplotypeConfig;
use phasing::{phase, phasing_input_from_graph};
use haplotype_inspection::identify_misaligned_haps::inspect_haplotypes;
use sdrecall_utils::{Result, SdError};

/// Parameters for one fused FP-control island run, mirroring the per-chunk
/// arguments the Python `realign_filter_per_cov.py` threads through.
#[derive(Clone, Debug)]
pub struct FpControlParams {
    /// Reference genome FASTA (with `.fai`); needed by the allele-depth mpileup
    /// in the graph build.
    pub reference_genome: String,
    /// Per-vertex edge-weight cutoff separating the phasing clique rounds
    /// (Python `edge_weight_cutoff`, default 0.301 in the dumped meta).
    pub edge_weight_cutoff: f32,
    /// Mean read length (drives `HaplotypeConfig`'s score array).
    pub mean_read_length: f32,
    /// MAPQ floor for graph + inspection read filtering (`recall_mq_cutoff`).
    pub mapq_cutoff: u8,
    /// Median-base-quality floor (`basequal_median_cutoff`).
    pub basequal_median_cutoff: u8,
    /// htslib BGZF / collate thread budget.
    pub threads: u8,
    /// Path for the haplotype-comparison meta TSV inspect writes (Python's
    /// `compare_haplotype_meta_tab`). May be empty to skip the dump.
    pub compare_haplotype_meta_tab: String,
}

impl Default for FpControlParams {
    fn default() -> Self {
        Self {
            reference_genome: String::new(),
            edge_weight_cutoff: 0.301,
            mean_read_length: 148.0,
            mapq_cutoff: 10,
            basequal_median_cutoff: 10,
            threads: 4,
            compare_haplotype_meta_tab: String::new(),
        }
    }
}

/// The two disjoint qname sets the FP-control core produces.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FpControlOutput {
    pub correct_qnames: Vec<String>,
    pub mismap_qnames: Vec<String>,
}

impl FpControlOutput {
    /// Sorted for stable output / comparison (the underlying sets are unordered).
    fn from_sets(correct: HashSet<String>, mismap: HashSet<String>) -> Self {
        let mut correct_qnames: Vec<String> = correct.into_iter().collect();
        let mut mismap_qnames: Vec<String> = mismap.into_iter().collect();
        correct_qnames.sort();
        mismap_qnames.sort();
        Self { correct_qnames, mismap_qnames }
    }
}

/// The phasing partition, in the three shapes `inspect_haplotypes` consumes.
/// Kept as a named bundle so the glue between phasing and inspection is one unit
/// rather than three parallel ad-hoc conversions.
struct Partition {
    /// hap_id → [qname]  (Python `hap_qname_info`; values are sets there, lists here).
    hap_qname_info: HashMap<i32, Vec<String>>,
    /// vertex index → hap_id  (Python `qname_hap_info`).
    qname_hap_info: HashMap<i32, i32>,
    /// qname → vertex index  (Python `qname_to_node`).
    qname_to_node: HashMap<String, i32>,
}

/// Derive the three inspection maps from the phasing result + the vertex→qname
/// vector. Replays exactly what `realign_filter_per_cov.py` materialises before
/// calling inspect: `hap_qname_info` (defaultdict(set) → list), `qname_hap_info`
/// keyed by vertex index, and `qname_to_node`.
fn build_partition(vertex_hap: &HashMap<i32, i32>, vertex_qname: &[String]) -> Partition {
    let mut hap_qname_info: HashMap<i32, Vec<String>> = HashMap::new();
    let mut qname_to_node: HashMap<String, i32> = HashMap::new();

    for (&vert, &hid) in vertex_hap {
        if let Some(qname) = vertex_qname.get(vert as usize) {
            hap_qname_info.entry(hid).or_default().push(qname.clone());
            qname_to_node.insert(qname.clone(), vert);
        }
    }
    // Deterministic ordering within each haplotype's qname list.
    for qnames in hap_qname_info.values_mut() {
        qnames.sort();
    }

    Partition {
        hap_qname_info,
        qname_hap_info: vertex_hap.clone(),
        qname_to_node,
    }
}

/// Run the fused Phase-2c FP-control core entirely in Rust.
///
/// `(bam, intrinsic_bam, params) → (correct_qnames, mismap_qnames)`.
///
/// Returns `Ok(None)` when the island is skipped without classification — the
/// two early-outs the Python pipeline also takes:
///   * the phasing graph has ≤ 2 vertices or an empty weight matrix
///     (`build_phasing_graph_rust` returns Python `None` here), and
///   * (handled internally, *not* a `None`) the ≤ 2-haplotype shortcut, where
///     every read is declared correct and nothing is mismapped — matching
///     `realign_filter_per_cov.py:318-320`.
pub fn run_fp_control(
    bam: &str,
    intrinsic_bam: &str,
    params: &FpControlParams,
) -> Result<Option<FpControlOutput>> {
    // ── Stage 1: BAM → read-pair map + allele depths → phasing graph ────────
    // Re-create what `build_phasing_graph_rust` does (it owns the BAM read), but
    // keep the Rust result in-process instead of exporting to Python.
    log::info!("[fp_control] Stage 1: building phasing graph from {bam}");
    let (read_pair_map, header) = migrate_bam_to_sorted_intervals_grouped(
        bam,
        params.mapq_cutoff,
        params.basequal_median_cutoff,
        /* filter_noisy */ true,
        /* use_collate  */ true,
        params.threads,
    )
    .map_err(|e| SdError::Compute(format!("BAM read/pairing failed for {bam}: {e}")))?;

    let allele_depth_map = build_allele_depth_map(
        bam,
        &params.reference_genome,
        params.mapq_cutoff,
        params.basequal_median_cutoff,
    )
    .map_err(|e| SdError::Compute(format!("allele-depth map failed for {bam}: {e}")))?;

    let config = HaplotypeConfig::new(params.mean_read_length);
    let graph = build_graph(&read_pair_map, &allele_depth_map, &header, &config)
        .map_err(|e| SdError::Compute(format!("graph build failed for {bam}: {e}")))?;

    // Early-out #1: same gate as `build_phasing_graph_rust` (None → skip island).
    let weight_matrix = match graph.weight_matrix {
        Some(ref wm) => wm,
        None => {
            log::warn!("[fp_control] no ALT alleles / empty weight matrix for {bam}; skipping");
            return Ok(None);
        }
    };
    if graph.graph.node_count() <= 2 {
        log::warn!(
            "[fp_control] graph has {} vertices (<= 2) for {bam}; skipping",
            graph.graph.node_count()
        );
        return Ok(None);
    }

    // Vertex index → qname (direct correspondence: qname_idx == NodeIndex.index()).
    let n = graph.graph.node_count();
    let mut vertex_qname = vec![String::new(); n];
    for (&qname_idx, rp) in &read_pair_map.readpair_dict {
        if qname_idx < n {
            vertex_qname[qname_idx] = rp.qname.clone();
        }
    }

    // ── Stage 2: phasing (replaces the Python `phasing_realigned_reads` hop) ─
    log::info!("[fp_control] Stage 2: phasing {n} vertices");
    let phasing_input =
        phasing_input_from_graph(&graph, weight_matrix.clone(), params.edge_weight_cutoff);
    let vertex_hap = phase(&phasing_input);
    let partition = build_partition(&vertex_hap, &vertex_qname);

    log::info!(
        "[fp_control] phasing produced {} haplotype clusters over {} vertices",
        partition.hap_qname_info.len(),
        vertex_hap.len()
    );

    // Early-out #2: ≤ 2 haplotype clusters → no choice to make; every read is
    // correct, nothing mismapped (realign_filter_per_cov.py:318-320).
    if partition.hap_qname_info.len() <= 2 {
        let correct: HashSet<String> = partition
            .hap_qname_info
            .values()
            .flat_map(|qs| qs.iter().cloned())
            .collect();
        log::warn!(
            "[fp_control] only {} haplotype clusters for {bam}; declaring all {} reads correct",
            partition.hap_qname_info.len(),
            correct.len()
        );
        return Ok(Some(FpControlOutput::from_sets(correct, HashSet::new())));
    }

    // ── Stage 3: haplotype inspection (consensus / similarity / BILC) ───────
    // inspect_haplotypes re-opens the BAMs itself (single-read fusion deferred).
    log::info!("[fp_control] Stage 3: inspecting haplotypes for {bam}");
    let lowqual: HashSet<String> = graph.lowqual_qnames.iter().cloned().collect();
    let (correct, mismap) = inspect_haplotypes(
        bam,
        intrinsic_bam,
        &partition.hap_qname_info,
        &partition.qname_hap_info,
        &partition.qname_to_node,
        &lowqual,
        &params.compare_haplotype_meta_tab,
        params.mean_read_length as f64,
        params.mapq_cutoff,
        params.basequal_median_cutoff,
    )
    .map_err(|e| SdError::Compute(format!("haplotype inspection failed for {bam}: {e}")))?;

    debug_assert!(
        correct.is_disjoint(&mismap),
        "correct and mismap qname sets overlap"
    );
    log::info!(
        "[fp_control] inspection complete: {} correct, {} mismapped",
        correct.len(),
        mismap.len()
    );

    Ok(Some(FpControlOutput::from_sets(correct, mismap)))
}


#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    /// `build_partition` derives the three inspection maps from a vertex→hap
    /// map and the vertex→qname vector, exactly as the Python pre-inspect glue.
    #[test]
    fn build_partition_derives_three_maps() {
        // vertices 0,1 → hap 5 ; vertex 2 → hap 9
        let mut vertex_hap = HashMap::new();
        vertex_hap.insert(0, 5);
        vertex_hap.insert(1, 5);
        vertex_hap.insert(2, 9);
        let vertex_qname = vec!["qA".to_string(), "qB".to_string(), "qC".to_string()];

        let p = build_partition(&vertex_hap, &vertex_qname);

        assert_eq!(p.qname_hap_info, vertex_hap);
        assert_eq!(p.qname_to_node["qA"], 0);
        assert_eq!(p.qname_to_node["qB"], 1);
        assert_eq!(p.qname_to_node["qC"], 2);
        assert_eq!(p.hap_qname_info[&5], vec!["qA".to_string(), "qB".to_string()]);
        assert_eq!(p.hap_qname_info[&9], vec!["qC".to_string()]);
    }

    /// The ≤2-haplotype shortcut: all reads correct, none mismapped — the same
    /// branch as realign_filter_per_cov.py:318-320. Exercised via the partition
    /// + the from_sets collapse used in `run_fp_control`.
    #[test]
    fn two_haplotype_shortcut_marks_all_correct() {
        let mut vertex_hap = HashMap::new();
        vertex_hap.insert(0, 0);
        vertex_hap.insert(1, 1);
        let vertex_qname = vec!["qA".to_string(), "qB".to_string()];
        let p = build_partition(&vertex_hap, &vertex_qname);
        assert!(p.hap_qname_info.len() <= 2);

        let correct: HashSet<String> = p
            .hap_qname_info
            .values()
            .flat_map(|qs| qs.iter().cloned())
            .collect();
        let out = FpControlOutput::from_sets(correct, HashSet::new());
        assert_eq!(out.correct_qnames, vec!["qA".to_string(), "qB".to_string()]);
        assert!(out.mismap_qnames.is_empty());
    }

    /// `phasing_input_from_graph` flattens `node_read_ids` (drops the `None`
    /// mate), extracts undirected edges, and threads the weight matrix + per-read
    /// vector maps through unchanged. Built on a tiny hand-made graph result.
    #[test]
    fn phasing_input_glue_shapes_match() {
        use phasing::structs::PhasingGraphResult;

        let mut g = PhasingGraphResult::new();
        let a = g.graph.add_node(());
        let b = g.graph.add_node(());
        let c = g.graph.add_node(());
        g.graph.add_edge(a, b, 0.8);
        g.graph.add_edge(b, c, 0.5);

        g.node_read_ids = vec![
            ("qA:65".to_string(), Some("qA:129".to_string())),
            ("qB:65".to_string(), None),
            ("qC:65".to_string(), Some("qC:129".to_string())),
        ];
        g.read_hap_vectors.insert("qA:65".to_string(), vec![1, 1, 1]);
        g.read_error_vectors.insert("qA:65".to_string(), vec![0.01, 0.01, 0.01]);

        let wm = Array2::<f32>::zeros((3, 3));
        let pi = phasing_input_from_graph(&g, wm, 0.301);

        // edges (order follows edge_indices insertion)
        assert_eq!(pi.edges, vec![(0, 1), (1, 2)]);
        // node_read_ids flattened, None mate dropped
        assert_eq!(pi.node_read_ids[0], vec!["qA:65".to_string(), "qA:129".to_string()]);
        assert_eq!(pi.node_read_ids[1], vec!["qB:65".to_string()]);
        // per-read vectors threaded through
        assert_eq!(pi.read_hap["qA:65"], vec![1, 1, 1]);
        assert_eq!(pi.read_err["qA:65"], vec![0.01, 0.01, 0.01]);
        assert_eq!(pi.weight_matrix.dim(), (3, 3));
    }
}
