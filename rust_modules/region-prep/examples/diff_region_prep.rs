//! Differential harness: Rust `region-prep` vs Python `prepare_masked_align_region`.
//!
//! Runs [`region_prep::prepare_masked_align_region_per_rg`] on the SAME real-data
//! inputs the Python reference generator used (`test_tmp/gen_region_prep_ref.py`)
//! — the HG002 RG0 `RG0_related_homo_regions.bed` (7-col), the CMRG target BED,
//! and a synthetic `.fasta.fai` (huge contig lengths so slop never clamps, making
//! the run reference-independent). It then asserts, per subgroup, that the Rust
//! and Python fc-target / nfc BEDs have **identical sorted-merged interval sets**
//! (the task's pass criterion) AND identical coverage.
//!
//! ## Usage
//!
//! ```text
//! cargo run --example diff_region_prep -- \
//!   --whole-region-bed  <RG0_related_homo_regions.bed> \
//!   --target-bed        <HG002_hg38_CMRG_smallvar_v1.00.bed> \
//!   --ref-genome        <pyref>/synthetic.fasta \
//!   --pyref-manifest    <pyref>/manifest.tsv \
//!   --rust-out-dir      <scratch>
//! ```
//!
//! The Python side is generated out-of-band by
//! `test_tmp/gen_region_prep_ref.py <pyref> <N>` (so this harness has no Python
//! dependency and runs purely against the dumped reference).

use clap::Parser;
use region_prep::prepare_masked_align_region_per_rg;
use sdrecall_io::read_bed;
use sdrecall_utils::GenomicInterval;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

#[derive(Parser)]
#[command(name = "diff_region_prep")]
struct Cli {
    /// The 7-col whole-region BED (RG0_related_homo_regions.bed).
    #[arg(long)]
    whole_region_bed: PathBuf,
    /// The sample target BED (CMRG smallvar).
    #[arg(long)]
    target_bed: PathBuf,
    /// Reference whose <ref>.fasta.fai the Python reference also used.
    #[arg(long)]
    ref_genome: PathBuf,
    /// Python reference manifest: rg<TAB>sub<TAB>fc_bed<TAB>nfc_bed<TAB>fc_size<TAB>nfc_size.
    #[arg(long)]
    pyref_manifest: PathBuf,
    /// Scratch dir for the Rust fc/nfc outputs.
    #[arg(long)]
    rust_out_dir: PathBuf,
    /// RG label (default RG0).
    #[arg(long, default_value = "RG0")]
    rg_label: String,
}

/// Sorted (chrom, start, end) interval set of a BED file as BED3 — the comparison
/// canonical form. Both the Python and Rust outputs are ALREADY final merged BEDs
/// (pybedtools `.merge()` / `merge_bookended`), so we only normalise sort order +
/// strand before set-equality; we do NOT re-merge (that could mask a real
/// interval-count divergence).
fn canonical(path: &Path) -> Vec<GenomicInterval> {
    let mut ivs: Vec<GenomicInterval> = read_bed(path)
        .expect("read bed")
        .into_iter()
        .map(|iv| GenomicInterval::new(iv.chrom, iv.start, iv.end))
        .collect();
    ivs.sort_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });
    ivs
}

fn coverage(ivs: &[GenomicInterval]) -> i64 {
    ivs.iter().map(GenomicInterval::len).sum()
}

fn main() -> ExitCode {
    sdrecall_utils::init_console_logger(log::LevelFilter::Info);
    let cli = Cli::parse();

    // Parse the Python manifest → the subgroup ids + the per-sub reference paths.
    let manifest = std::fs::read_to_string(&cli.pyref_manifest).expect("read manifest");
    let mut rows: Vec<(String, PathBuf, PathBuf, i64, i64)> = Vec::new();
    let mut subids: Vec<String> = Vec::new();
    let mut py_fc_target: Option<PathBuf> = None;
    for line in manifest.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        // rg, sub, fc_bed, nfc_bed, fc_size, nfc_size
        let sub = c[1].to_string();
        let fc_bed = PathBuf::from(c[2]);
        let nfc_bed = PathBuf::from(c[3]);
        let fc_size: i64 = c[4].parse().unwrap();
        let nfc_size: i64 = c[5].parse().unwrap();
        py_fc_target.get_or_insert_with(|| fc_bed.clone());
        subids.push(sub.clone());
        rows.push((sub, fc_bed, nfc_bed, fc_size, nfc_size));
    }
    log::info!("Loaded {} Python reference subgroups", rows.len());

    std::fs::create_dir_all(&cli.rust_out_dir).expect("mkdir rust out");
    let rust_fc_target = cli.rust_out_dir.join("fc_target.bed");
    let rust_nfc_dir = cli.rust_out_dir.join("nfc");

    // Run the Rust path on identical inputs.
    let records = match prepare_masked_align_region_per_rg(
        &cli.rg_label,
        &subids,
        &cli.target_bed,
        &cli.whole_region_bed,
        &cli.ref_genome,
        &rust_fc_target,
        &rust_nfc_dir,
    ) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Rust per_rg failed: {e}");
            return ExitCode::FAILURE;
        }
    };
    let rec_by_sub: std::collections::HashMap<&str, &region_prep::RgSubgroupRecord> =
        records.iter().map(|r| (r.subgroup_id.as_str(), r)).collect();

    // --- FC target comparison (computed once, shared) ---------------------------
    let py_fc = canonical(py_fc_target.as_ref().unwrap());
    let rust_fc = canonical(&rust_fc_target);
    let mut failures = 0usize;
    if py_fc != rust_fc {
        eprintln!(
            "FC TARGET MISMATCH: py {} ivs (cov {}) vs rust {} ivs (cov {})",
            py_fc.len(),
            coverage(&py_fc),
            rust_fc.len(),
            coverage(&rust_fc)
        );
        // show first divergence
        for (i, (p, r)) in py_fc.iter().zip(rust_fc.iter()).enumerate() {
            if p != r {
                eprintln!("  first diff at {i}: py={p:?} rust={r:?}");
                break;
            }
        }
        failures += 1;
    } else {
        log::info!(
            "FC target MATCH: {} intervals, coverage {}",
            rust_fc.len(),
            coverage(&rust_fc)
        );
    }

    // --- per-subgroup NFC comparison -------------------------------------------
    let mut nfc_ok = 0usize;
    for (sub, _py_fc_bed, py_nfc_bed, _py_fc_size, py_nfc_size) in &rows {
        let rec = match rec_by_sub.get(sub.as_str()) {
            Some(r) => r,
            None => {
                eprintln!("sub {sub}: Rust produced no record");
                failures += 1;
                continue;
            }
        };
        let py = canonical(py_nfc_bed);
        let rust = canonical(&rec.nfc_bed);
        let py_cov = coverage(&py);
        let rust_cov = coverage(&rust);
        if py != rust || py_cov != rust_cov || rust_cov != rec.nfc_bed_size {
            eprintln!(
                "NFC MISMATCH sub {sub}: py {} ivs cov {} (manifest {}) vs rust {} ivs cov {} (record {})",
                py.len(),
                py_cov,
                py_nfc_size,
                rust.len(),
                rust_cov,
                rec.nfc_bed_size
            );
            for (i, (p, r)) in py.iter().zip(rust.iter()).enumerate() {
                if p != r {
                    eprintln!("  first diff at {i}: py={p:?} rust={r:?}");
                    break;
                }
            }
            if py.len() != rust.len() {
                eprintln!("  length differs: py {} vs rust {}", py.len(), rust.len());
            }
            failures += 1;
        } else {
            nfc_ok += 1;
        }
    }

    log::info!(
        "DIFFERENTIAL: {nfc_ok}/{} subgroups NFC-identical; {failures} total failures",
        rows.len()
    );
    if failures == 0 {
        println!("PASS: fc target + {nfc_ok} subgroup nfc BEDs identical (Rust == Python)");
        ExitCode::SUCCESS
    } else {
        println!("FAIL: {failures} mismatches");
        ExitCode::FAILURE
    }
}
