//! Differential harness: fused `fp-control` vs the current Rust+Python hybrid,
//! on the real per-island dumps under `t2_reval_dump/` + their island BAMs.
//!
//! For each island it:
//!   1. runs the **fused** path (`run_fp_control`: Rust graph → Rust phasing →
//!      inspect) on the real island BAM, and
//!   2. reconstructs the **hybrid** inspect result by feeding the *dumped Python
//!      partition* (`phasing_*_info.json` + `qname_to_node.json` from the hybrid
//!      run) straight into `inspect_haplotypes` on the same BAM,
//!
//! then asserts the two `(correct, mismap)` sets are equal.
//!
//! The phasing stage is already validated partition-identical to Python by the
//! `phasing` T3 harness (270/270 islands), so step 2 is exactly what the hybrid
//! computed. Equality here therefore proves the fuse reproduces the hybrid's
//! end-to-end classification while removing the Rust→Python→Rust round-trip.
//!
//! Usage:
//!   cargo run -p fp-control --example diff_vs_hybrid -- \
//!     --dump <DUMP_ROOT> --reference <hg38.fasta> [--island <N>]

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use phasing::bam_reading::{
    build_allele_depth_map, migrate_bam_to_sorted_intervals_grouped,
};
use phasing::graph_builder::build_phasing_graph;
use phasing::structs::HaplotypeConfig;
use clap::Parser;
use fp_control::{run_fp_control, FpControlParams};
use haplotype_inspection::identify_misaligned_haps::inspect_haplotypes;
use serde::Deserialize;

#[derive(Parser, Debug)]
struct Args {
    /// Dump root (one sub-dir per island).
    #[arg(long)]
    dump: PathBuf,
    /// Reference genome FASTA (with .fai) matching the island BAM contigs.
    #[arg(long)]
    reference: String,
    /// Restrict to a single island number (e.g. 1 → HG002.pooled.raw.deduped.1).
    #[arg(long)]
    island: Option<u32>,
}

#[derive(Deserialize)]
struct Meta {
    bam: String,
    intrinsic_bam: String,
    edge_weight_cutoff: f32,
    mean_read_length: f32,
    recall_mq_cutoff: u8,
    basequal_median_cutoff: u8,
    n_haplotypes: i64,
}

fn read_json<T: for<'de> Deserialize<'de>>(p: &Path) -> T {
    let s = std::fs::read_to_string(p).unwrap_or_else(|e| panic!("read {}: {e}", p.display()));
    serde_json::from_str(&s).unwrap_or_else(|e| panic!("parse {}: {e}", p.display()))
}

/// Recompute the graph's low-qual qname set the way both the fuse and the hybrid
/// do (build_phasing_graph returns it; the hybrid threads it straight into
/// inspect). The dump doesn't store it, so we reproduce it from the BAM.
fn graph_lowqual(meta: &Meta, reference: &str) -> HashSet<String> {
    let (rpm, header) = migrate_bam_to_sorted_intervals_grouped(
        &meta.bam,
        meta.recall_mq_cutoff,
        meta.basequal_median_cutoff,
        true,
        true,
        4,
    )
    .expect("graph BAM read");
    let adm = build_allele_depth_map(
        &meta.bam,
        reference,
        meta.recall_mq_cutoff,
        meta.basequal_median_cutoff,
    )
    .expect("allele depth");
    // Intrinsic AD map (Python `intrinsic_ad_dict`) — built the same way the fuse
    // does in `build_and_phase_with_intrinsic`, so this reproduces the production
    // graph (PSV-aware edge weights) the hybrid was compared against.
    let intrinsic_adm = build_allele_depth_map(
        &meta.intrinsic_bam,
        reference,
        meta.recall_mq_cutoff,
        meta.basequal_median_cutoff,
    )
    .expect("intrinsic allele depth");
    let cfg = HaplotypeConfig::new(meta.mean_read_length);
    let g = build_phasing_graph(&rpm, &adm, &intrinsic_adm, &header, &cfg).expect("graph build");
    g.lowqual_qnames.iter().cloned().collect()
}

/// Reconstruct the hybrid inspect result from the dumped Python partition.
fn hybrid_inspect(
    dir: &Path,
    meta: &Meta,
    lowqual: &HashSet<String>,
) -> (HashSet<String>, HashSet<String>) {
    // hap_id → [qname]
    let hap_qname_str: HashMap<String, Vec<String>> =
        read_json(&dir.join("phasing_hap_qname_info.json"));
    let hap_qname_info: HashMap<i32, Vec<String>> = hap_qname_str
        .into_iter()
        .map(|(k, v)| (k.parse().unwrap(), v))
        .collect();

    // vertex → hap_id
    let qhap_str: HashMap<String, i32> = read_json(&dir.join("phasing_qname_hap_info.json"));
    let qname_hap_info: HashMap<i32, i32> =
        qhap_str.into_iter().map(|(k, v)| (k.parse().unwrap(), v)).collect();

    // qname → vertex
    let qname_to_node: HashMap<String, i32> = read_json(&dir.join("qname_to_node.json"));

    inspect_haplotypes(
        &meta.bam,
        &meta.intrinsic_bam,
        &hap_qname_info,
        &qname_hap_info,
        &qname_to_node,
        lowqual,
        "",
        meta.mean_read_length as f64,
        meta.recall_mq_cutoff,
        meta.basequal_median_cutoff,
    )
    .expect("hybrid inspect_haplotypes failed")
}

fn main() {
    let args = Args::parse();
    let mut dirs: Vec<PathBuf> = std::fs::read_dir(&args.dump)
        .expect("read dump root")
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.is_dir() && p.join("meta.json").exists())
        .collect();
    dirs.sort();

    let mut total = 0usize;
    let mut matched = 0usize;
    let mut skipped = 0usize;

    for dir in &dirs {
        let name = dir.file_name().unwrap().to_str().unwrap();
        if let Some(isl) = args.island {
            if name != format!("HG002.pooled.raw.deduped.{isl}") {
                continue;
            }
        }
        let meta: Meta = read_json(&dir.join("meta.json"));
        if !Path::new(&meta.bam).exists() || !Path::new(&meta.intrinsic_bam).exists() {
            skipped += 1;
            println!("  SKIP   {name}  (BAM missing)");
            continue;
        }
        // Only islands with > 2 haplotypes exercise the real inspect path; the
        // ≤2 shortcut is covered by the lib unit test.
        if meta.n_haplotypes <= 2 {
            skipped += 1;
            continue;
        }

        total += 1;

        let params = FpControlParams {
            reference_genome: args.reference.clone(),
            edge_weight_cutoff: meta.edge_weight_cutoff,
            mean_read_length: meta.mean_read_length,
            mapq_cutoff: meta.recall_mq_cutoff,
            basequal_median_cutoff: meta.basequal_median_cutoff,
            threads: 4,
            compare_haplotype_meta_tab: String::new(),
        };

        let fused = match run_fp_control(&meta.bam, &meta.intrinsic_bam, &params) {
            Ok(Some(o)) => o,
            Ok(None) => {
                skipped += 1;
                println!("  SKIP   {name}  (fuse returned None)");
                continue;
            }
            Err(e) => {
                println!("  ERROR  {name}: {e}");
                continue;
            }
        };
        let fused_correct: HashSet<String> = fused.correct_qnames.into_iter().collect();
        let fused_mismap: HashSet<String> = fused.mismap_qnames.into_iter().collect();

        let lowqual = graph_lowqual(&meta, &args.reference);
        let (hyb_correct, hyb_mismap) = hybrid_inspect(dir, &meta, &lowqual);

        if fused_correct == hyb_correct && fused_mismap == hyb_mismap {
            matched += 1;
            println!(
                "  MATCH  {name}  ({} correct, {} mismap)",
                fused_correct.len(),
                fused_mismap.len()
            );
        } else {
            println!("  MISMATCH {name}");
            println!(
                "    fused  correct={} mismap={}",
                fused_correct.len(),
                fused_mismap.len()
            );
            println!(
                "    hybrid correct={} mismap={}",
                hyb_correct.len(),
                hyb_mismap.len()
            );
            for q in fused_correct.symmetric_difference(&hyb_correct).take(4) {
                println!("    correct-diff: {q}");
            }
            for q in fused_mismap.symmetric_difference(&hyb_mismap).take(4) {
                println!("    mismap-diff:  {q}");
            }
        }
    }

    println!("\nResult: {matched}/{total} islands match (skipped {skipped})");
    if matched != total {
        std::process::exit(1);
    }
}
