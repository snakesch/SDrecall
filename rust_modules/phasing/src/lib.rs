//! Standalone BAM phaser for paired short-read data.
//!
//! Given a BAM file and reference FASTA, builds a phasing graph from read-pair
//! overlaps, clusters reads into haplotypes via Greedy-Clique-Expansion, and
//! writes an HP-tagged BAM (`HP:Z:hap{N}` per read).
//!
//! ## High-level API
//!
//! ```no_run
//! use phasing::{phase_bam, PhaserParams};
//! let out = phase_bam("in.bam", "ref.fa", "out.bam", &PhaserParams::default()).unwrap();
//! ```
//!
//! ## Mid-level API (used by fp-control)
//!
//! The individual stages are public modules:
//! - [`bam_reading`] — BAM → ReadPairMap + AlleleDepthMap
//! - [`graph_builder`] — ReadPairMap → PhasingGraphResult
//! - [`phasing`] — weight matrix → vertex→haplotype partition
//! - [`hp_writer`] — partition → HP-tagged BAM

// Graph construction (from build_phasing_graph)
pub mod bam_reading;
pub mod graph_builder;
pub mod haplotype_determination;
pub mod structs;

// Phasing algorithm
pub mod kernels;
pub mod gce;
pub mod phasing;

// HP tag output
pub mod hp_writer;

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

pub use phasing::{phase, qname_partition, PhasingInput, Round};

use sdrecall_utils::{Result, SdError};

/// Parameters for the standalone BAM phaser.
///
/// The reference genome is **not** stored here — it is a required input passed
/// directly to [`phase_bam`], so there is a single source of truth and no way
/// to silently desync the argument from a struct field.
#[derive(Clone, Debug)]
pub struct PhaserParams {
    pub edge_weight_cutoff: f32,
    pub mean_read_length: f32,
    pub mapq_cutoff: u8,
    pub basequal_median_cutoff: u8,
    pub threads: u8,
}

impl Default for PhaserParams {
    fn default() -> Self {
        Self {
            edge_weight_cutoff: 0.301,
            mean_read_length: 148.0,
            mapq_cutoff: 10,
            basequal_median_cutoff: 10,
            threads: 4,
        }
    }
}

/// Output from the standalone BAM phaser.
#[derive(Clone, Debug)]
pub struct PhaserOutput {
    pub output_bam: PathBuf,
    pub n_reads: usize,
    pub n_haplotypes: usize,
}

/// Assemble the phasing-stage input from the graph result.
///
/// Extracts edges, flattens node_read_ids, converts ahash maps to std HashMap.
/// Previously lived in fp-control as glue between the two separate crates.
pub fn phasing_input_from_graph(
    graph: &structs::PhasingGraphResult,
    weight_matrix: ndarray::Array2<f32>,
    edge_weight_cutoff: f32,
) -> PhasingInput {
    let edges: Vec<(i32, i32)> = graph
        .graph
        .edge_indices()
        .filter_map(|e| {
            graph
                .graph
                .edge_endpoints(e)
                .map(|(s, t)| (s.index() as i32, t.index() as i32))
        })
        .collect();

    let node_read_ids: Vec<Vec<String>> = graph
        .node_read_ids
        .iter()
        .map(|(id1, id2)| match id2 {
            Some(id2) => vec![id1.clone(), id2.clone()],
            None => vec![id1.clone()],
        })
        .collect();

    let read_hap: HashMap<String, Vec<i16>> = graph
        .read_hap_vectors
        .iter()
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect();
    let read_err: HashMap<String, Vec<f32>> = graph
        .read_error_vectors
        .iter()
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect();

    PhasingInput {
        weight_matrix,
        edges,
        edge_weight_cutoff,
        node_read_ids,
        read_hap,
        read_err,
    }
}

/// Validate that the BAM contains paired short reads.
///
/// Reads the first ~1000 primary records and checks:
/// - Most reads have the paired flag set
/// - Median sequence length < 1000 bp (short-read threshold)
fn validate_input(bam: &str) -> Result<()> {
    use rust_htslib::bam::{self, Read};

    let mut reader =
        bam::Reader::from_path(bam).map_err(|e| SdError::Htslib(format!("cannot open {bam}: {e}")))?;

    let mut paired_count = 0u32;
    let mut total_count = 0u32;
    let mut lengths: Vec<usize> = Vec::with_capacity(1000);

    for result in reader.records() {
        let rec = result.map_err(|e| SdError::Htslib(format!("BAM read error: {e}")))?;
        if rec.is_secondary() || rec.is_supplementary() {
            continue;
        }
        total_count += 1;
        if rec.is_paired() {
            paired_count += 1;
        }
        lengths.push(rec.seq_len());
        if total_count >= 1000 {
            break;
        }
    }

    if total_count == 0 {
        return Err(SdError::Compute("BAM file contains no primary alignments".into()));
    }

    let paired_frac = paired_count as f64 / total_count as f64;
    if paired_frac < 0.9 {
        return Err(SdError::Compute(format!(
            "BAM does not appear to be paired-end ({:.0}% of first {total_count} reads are paired; need >=90%)",
            paired_frac * 100.0
        )));
    }

    lengths.sort_unstable();
    let median_len = lengths[lengths.len() / 2];
    if median_len >= 1000 {
        return Err(SdError::Compute(format!(
            "BAM appears to contain long reads (median length {median_len} bp >= 1000); \
             this tool is designed for short-read paired-end data"
        )));
    }

    log::info!(
        "[validate_input] OK: {paired_count}/{total_count} paired, median read length {median_len} bp"
    );
    Ok(())
}

/// Build the vertex_idx → qname lookup from a ReadPairMap.
fn build_vertex_qname(
    read_pair_map: &structs::ReadPairMap,
    n: usize,
) -> Vec<String> {
    let mut vertex_qname = vec![String::new(); n];
    for (&qname_idx, rp) in &read_pair_map.readpair_dict {
        if qname_idx < n {
            vertex_qname[qname_idx] = rp.qname.clone();
        }
    }
    vertex_qname
}

/// The phasing partition for one BAM — the shared output of the standalone
/// [`phase_bam`] and the pipeline's `fp_control::run_fp_control`.
///
/// Produced by [`build_and_phase`], which runs the BAM → graph → GCE-partition
/// stages once so neither caller re-implements that wiring.
#[derive(Clone, Debug)]
pub struct PhasedReads {
    /// vertex index → haplotype id (Python `qname_hap_info`).
    pub vertex_hap: HashMap<i32, i32>,
    /// vertex index → qname.
    pub vertex_qname: Vec<String>,
    /// Reads dropped as low quality during BAM reading. The pipeline threads
    /// these into haplotype inspection; the standalone phaser ignores them.
    pub lowqual_qnames: HashSet<String>,
}

/// BAM → phasing graph → GCE partition: the shared core of the standalone phaser
/// ([`phase_bam`]) and the fused FP-control pipeline (`fp_control::run_fp_control`).
///
/// Runs stages 1–2 — read+pair the BAM, build the allele-depth map and phasing
/// graph, then phase the weight matrix. Returns `Ok(None)` for the single
/// early-out both callers special-case: the graph has no weight matrix (no ALT
/// alleles). The downstream gates differ between callers (the phaser writes an
/// unphased BAM; fp-control skips the island / applies the ≤2-vertex and
/// ≤2-haplotype shortcuts), so they stay with the callers.
pub fn build_and_phase(
    bam: &str,
    reference: &str,
    params: &PhaserParams,
) -> Result<Option<PhasedReads>> {
    // Standalone path: no intrinsic BAM, so paralogous-sequence-variant detection
    // is disabled (empty intrinsic AD map → psv_snv_count == 0), matching Python's
    // `intrinsic_ad_dict = {}` default in `graph_build.py`.
    build_and_phase_with_intrinsic(bam, reference, None, params)
}

/// Like [`build_and_phase`], but also threads an optional **intrinsic BAM** whose
/// allele-depth map drives paralogous-sequence-variant (PSV) detection in the
/// edge-weight formula (Python `intrinsic_ad_dict`; see
/// `haplotype_determination::psv_shared_snvs`).
///
/// `intrinsic_bam = None` ⇒ an empty intrinsic AD map ⇒ `psv_snv_count == 0` for
/// every shared SNV (the standalone-phaser behaviour, and Python's default when no
/// intrinsic dict is supplied). The fused FP-control pipeline (which *does* have an
/// intrinsic BAM, mirroring the original `graph_build.py` that built both
/// `nested_ad_dict` and `intrinsic_ad_dict`) should call this with
/// `Some(intrinsic_bam)` so the phasing weights match the Python production
/// pipeline.
///
/// The intrinsic AD map is built by the same `build_allele_depth_map`
/// (`bcftools mpileup … | bcftools query`) used for the main map — i.e. exactly
/// Python's `stat_ad_to_dict` over the intrinsic BAM.
pub fn build_and_phase_with_intrinsic(
    bam: &str,
    reference: &str,
    intrinsic_bam: Option<&str>,
    params: &PhaserParams,
) -> Result<Option<PhasedReads>> {
    let (read_pair_map, header) = bam_reading::migrate_bam_to_sorted_intervals_grouped(
        bam,
        params.mapq_cutoff,
        params.basequal_median_cutoff,
        true,
        true,
        params.threads,
    )
    .map_err(|e| SdError::Compute(format!("BAM read/pairing failed for {bam}: {e}")))?;

    let allele_depth_map = bam_reading::build_allele_depth_map(
        bam,
        reference,
        params.mapq_cutoff,
        params.basequal_median_cutoff,
    )
    .map_err(|e| SdError::Compute(format!("allele-depth map failed for {bam}: {e}")))?;

    // Intrinsic-bam allele-depth map (Python `intrinsic_ad_dict`). Built with the
    // same mpileup filters as the main map (Python `stat_ad_to_dict` hardcodes
    // `-q 10 -Q 10`; PhaserParams default to 10/10). Absent intrinsic BAM (or one
    // with no ALT alleles) → empty map → no PSV deductions.
    let intrinsic_ad_map = match intrinsic_bam {
        Some(ib) => bam_reading::build_allele_depth_map(
            ib,
            reference,
            params.mapq_cutoff,
            params.basequal_median_cutoff,
        )
        .map_err(|e| SdError::Compute(format!("intrinsic allele-depth map failed for {ib}: {e}")))?,
        None => structs::AlleleDepthMap::new(),
    };

    let config = structs::HaplotypeConfig::new(params.mean_read_length);
    let graph = graph_builder::build_phasing_graph(
        &read_pair_map,
        &allele_depth_map,
        &intrinsic_ad_map,
        &header,
        &config,
    )
    .map_err(|e| SdError::Compute(format!("graph build failed for {bam}: {e}")))?;

    let weight_matrix = match graph.weight_matrix {
        Some(ref wm) => wm.clone(),
        None => return Ok(None),
    };

    let n = graph.graph.node_count();
    let vertex_qname = build_vertex_qname(&read_pair_map, n);
    let phasing_input = phasing_input_from_graph(&graph, weight_matrix, params.edge_weight_cutoff);
    let vertex_hap = phase(&phasing_input);

    Ok(Some(PhasedReads {
        vertex_hap,
        vertex_qname,
        lowqual_qnames: graph.lowqual_qnames,
    }))
}

/// High-level API: BAM → HP-tagged BAM.
///
/// No-intrinsic-BAM convenience wrapper over [`phase_bam_with_intrinsic`]
/// (paralogous-sequence-variant detection disabled → Python `intrinsic_ad_dict = {}`).
pub fn phase_bam(
    bam: &str,
    reference: &str,
    output_bam: &str,
    params: &PhaserParams,
) -> Result<PhaserOutput> {
    phase_bam_with_intrinsic(bam, reference, output_bam, None, params)
}

/// Like [`phase_bam`], but threads an optional **intrinsic BAM** so the phasing
/// edge weights include the paralogous-sequence-variant (PSV) term
/// (Python `intrinsic_ad_dict`; see [`build_and_phase_with_intrinsic`] and
/// `haplotype_determination::psv_shared_snvs`). `intrinsic_bam = None` reproduces
/// the standalone-phaser behaviour (no PSV deductions), so the CLI can opt in via
/// `--intrinsic-bam` for Python-parity differential runs.
pub fn phase_bam_with_intrinsic(
    bam: &str,
    reference: &str,
    output_bam: &str,
    intrinsic_bam: Option<&str>,
    params: &PhaserParams,
) -> Result<PhaserOutput> {
    validate_input(bam)?;

    log::info!("[phase_bam] building + phasing graph from {bam}");
    let phased = match build_and_phase_with_intrinsic(bam, reference, intrinsic_bam, params)? {
        Some(p) => p,
        None => {
            log::warn!("[phase_bam] no ALT alleles / empty weight matrix; writing unphased BAM");
            let n_tagged = hp_writer::write_hp_tagged_bam(
                bam,
                output_bam,
                &HashMap::new(),
                &[],
                params.threads,
            )?;
            return Ok(PhaserOutput {
                output_bam: PathBuf::from(output_bam),
                n_reads: n_tagged,
                n_haplotypes: 0,
            });
        }
    };

    let n_haplotypes = {
        let mut hap_ids: Vec<i32> = phased.vertex_hap.values().copied().collect();
        hap_ids.sort_unstable();
        hap_ids.dedup();
        hap_ids.len()
    };
    log::info!(
        "[phase_bam] phasing produced {n_haplotypes} haplotype clusters over {} vertices",
        phased.vertex_hap.len()
    );

    log::info!("[phase_bam] writing HP-tagged BAM to {output_bam}");
    let n_tagged = hp_writer::write_hp_tagged_bam(
        bam,
        output_bam,
        &phased.vertex_hap,
        &phased.vertex_qname,
        params.threads,
    )?;

    Ok(PhaserOutput {
        output_bam: PathBuf::from(output_bam),
        n_reads: n_tagged,
        n_haplotypes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    #[test]
    fn phasing_input_glue_shapes_match() {
        let mut g = structs::PhasingGraphResult::new();
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
        g.read_hap_vectors
            .insert("qA:65".to_string(), vec![1, 1, 1]);
        g.read_error_vectors
            .insert("qA:65".to_string(), vec![0.01, 0.01, 0.01]);

        let wm = Array2::<f32>::zeros((3, 3));
        let pi = phasing_input_from_graph(&g, wm, 0.301);

        assert_eq!(pi.edges, vec![(0, 1), (1, 2)]);
        assert_eq!(
            pi.node_read_ids[0],
            vec!["qA:65".to_string(), "qA:129".to_string()]
        );
        assert_eq!(pi.node_read_ids[1], vec!["qB:65".to_string()]);
        assert_eq!(pi.read_hap["qA:65"], vec![1, 1, 1]);
        assert_eq!(pi.read_err["qA:65"], vec![0.01, 0.01, 0.01]);
        assert_eq!(pi.weight_matrix.dim(), (3, 3));
    }
}
