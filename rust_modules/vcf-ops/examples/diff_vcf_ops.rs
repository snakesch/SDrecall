//! Differential harness (T5 tier-2): run the Rust `merge` / `inhouse-common` path
//! on a real query+reference VCF pair and compare to a dumped Python output VCF.
//!
//! Per the project's test-verification requirement, run with `RUST_LOG=debug` and
//! redirect to `/paedyl01/disk1/yangyxt/test_tmp/diff_vcf_ops_<path>.log`.
//!
//! Usage:
//!   cargo run --example diff_vcf_ops -- merge \
//!     --query QUERY.vcf.gz --reference REF.vcf.gz --ref-genome REF.fasta \
//!     --python-out PYTHON_MERGED.vcf.gz --out RUST_MERGED.vcf.gz
//!
//!   cargo run --example diff_vcf_ops -- inhouse \
//!     --query QUERY.vcf.gz --cohort COHORT.vcf.gz --ref-genome REF.fasta \
//!     --python-out PYTHON_IC.vcf.gz --out RUST_IC.vcf.gz
//!
//! The comparison is record-level after a common `bcftools norm` on both: sort by
//! (CHROM,POS,REF,ALT), compare CHROM/POS/REF/ALT/FILTER + the GT/AD/GQ FORMAT
//! fields. Reports the first differing record; exit 0 only on zero differences.

use rust_htslib::bcf::{self, Read};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use vcf_ops::{
    annotate_inhouse_common, merge_with_priority, InhouseParams, MergeParams,
};

fn arg(args: &[String], flag: &str) -> Option<String> {
    args.iter().position(|a| a == flag).and_then(|i| args.get(i + 1).cloned())
}

fn main() {
    let level = std::env::var("RUST_LOG")
        .ok()
        .and_then(|l| l.parse().ok())
        .unwrap_or(log::LevelFilter::Info);
    sdrecall_utils::init_console_logger(level);

    let args: Vec<String> = std::env::args().collect();
    let mode = args.get(1).cloned().unwrap_or_default();
    let query = PathBuf::from(arg(&args, "--query").expect("--query"));
    let ref_genome = PathBuf::from(arg(&args, "--ref-genome").expect("--ref-genome"));
    let python_out = PathBuf::from(arg(&args, "--python-out").expect("--python-out"));
    let out = PathBuf::from(arg(&args, "--out").expect("--out"));
    let tmp_dir = PathBuf::from(arg(&args, "--tmp-dir").unwrap_or_else(|| "/tmp".to_string()));

    match mode.as_str() {
        "merge" => {
            let reference = PathBuf::from(arg(&args, "--reference").expect("--reference"));
            merge_with_priority(MergeParams {
                query_vcf: &query,
                reference_vcf: &reference,
                output_vcf: &out,
                ref_genome: &ref_genome,
                added_filter: Some("MISALIGNED"),
                qv_tag: Some("RAW"),
                rv_tag: Some("CLEAN"),
                modify_gt: false,
                threads: 4,
                tmp_dir: &tmp_dir,
            })
            .expect("rust merge failed");
        }
        "inhouse" => {
            let cohort = PathBuf::from(arg(&args, "--cohort").expect("--cohort"));
            annotate_inhouse_common(InhouseParams {
                query_vcf: &query,
                cohort_vcf: &cohort,
                output_vcf: &out,
                ref_genome: &ref_genome,
                added_filter: "INHOUSE_COMMON",
                inhouse_common_cutoff: 0.01,
                conf_level: 0.999,
                threads: 4,
                tmp_dir: &tmp_dir,
            })
            .expect("rust inhouse-common failed");
        }
        other => panic!("unknown mode {other:?} (expected merge|inhouse)"),
    }

    // Normalize both outputs identically, then compare.
    let norm_rust = tmp_dir.join("diff.rust.norm.vcf.gz");
    let norm_py = tmp_dir.join("diff.py.norm.vcf.gz");
    bcftools_norm(&out, &ref_genome, &norm_rust);
    bcftools_norm(&python_out, &ref_genome, &norm_py);

    let rust = load(&norm_rust);
    let py = load(&norm_py);

    let mut diffs = 0usize;
    let all_keys: std::collections::BTreeSet<_> =
        rust.keys().chain(py.keys()).cloned().collect();
    for k in &all_keys {
        match (rust.get(k), py.get(k)) {
            (Some(r), Some(p)) if r == p => {}
            (a, b) => {
                if diffs < 20 {
                    eprintln!("DIFF @ {k:?}\n  rust: {a:?}\n  py:   {b:?}");
                }
                diffs += 1;
            }
        }
    }

    println!(
        "diff_vcf_ops[{mode}]: rust={} py={} total_keys={} diffs={}",
        rust.len(),
        py.len(),
        all_keys.len(),
        diffs
    );
    if diffs == 0 {
        println!("PASS: record-identical after bcftools norm");
    } else {
        eprintln!("FAIL: {diffs} differing records");
        std::process::exit(1);
    }
}

// leaf: external normalization, same as the crate's norm leaf (parity baseline).
fn bcftools_norm(input: &Path, ref_genome: &Path, out: &Path) {
    let script = format!(
        "set -o pipefail; bcftools norm -m -both -f '{rg}' '{inp}' -Ou | \
         bcftools norm -d exact -Ou - | bcftools sort -Oz -o '{out}' - && bcftools index -f '{out}'",
        rg = ref_genome.display(),
        inp = input.display(),
        out = out.display()
    );
    let st = Command::new("bash").arg("-c").arg(&script).status().expect("spawn bcftools");
    assert!(st.success(), "bcftools norm failed for {}", input.display());
}

/// Load a normalized VCF into a map keyed by (chrom,pos,ref,alt) → comparable
/// (FILTER set, GT, AD, GQ) tuple.
fn load(path: &Path) -> BTreeMap<(String, i64, String, String), String> {
    let mut reader = bcf::Reader::from_path(path).expect("open normalized vcf");
    let mut map = BTreeMap::new();
    for res in reader.records() {
        let rec = res.expect("read");
        let hdr = rec.header();
        let chrom = rec
            .rid()
            .map(|r| String::from_utf8_lossy(hdr.rid2name(r).unwrap()).to_string())
            .unwrap_or_default();
        let alleles = rec.alleles();
        let ref_a = String::from_utf8_lossy(alleles[0]).to_string();
        let alt_a = alleles
            .get(1)
            .map(|a| String::from_utf8_lossy(a).to_string())
            .unwrap_or_default();

        // FILTER set (sorted for order-independence at the comparison level — the
        // crate's ordering parity is unit-tested separately).
        let mut filters: Vec<String> =
            rec.filters().map(|id| String::from_utf8_lossy(&hdr.id_to_name(id)).to_string()).collect();
        filters.sort();

        let gt = rec
            .genotypes()
            .ok()
            .map(|g| format!("{:?}", g.get(0)))
            .unwrap_or_default();
        let ad = rec
            .format(b"AD")
            .integer()
            .ok()
            .and_then(|b| b.first().map(|s| format!("{s:?}")))
            .unwrap_or_default();
        let gq = rec
            .format(b"GQ")
            .integer()
            .ok()
            .and_then(|b| b.first().map(|s| format!("{s:?}")))
            .unwrap_or_default();

        map.insert(
            (chrom, rec.pos(), ref_a, alt_a),
            format!("F={filters:?} GT={gt} AD={ad} GQ={gq}"),
        );
    }
    map
}
