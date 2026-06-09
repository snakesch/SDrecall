/// Cross-validation binary: reads Python-generated haplotype_meta.tsv files,
/// feeds the BILC columns into the Rust solver, and compares against Python output.
///
/// Supports both single-file and batch mode (directory or glob via shell expansion).
///
/// The mismap/correct_map columns in the TSV include post-ILP augmentation
/// (lines 1416–1420 of identify_misaligned_haps.py), so we reconstruct
/// the "ILP-only" drop set by subtracting the post-ILP additions.
///
/// Post-ILP additions (haps flagged mismap AFTER the solver):
///   - extreme_vard == True
///   - hap_max_sim_scores > 15
///   - scatter_hap == True AND hap_var_count <= 3
///
/// Usage:
///   # Single file:
///   cargo run --bin validate_bilc_solver -- <path_to_haplotype_meta.tsv>
///
///   # Batch (multiple files):
///   cargo run --bin validate_bilc_solver -- /dir/*.haplotype_meta.tsv
///
///   # Batch with TSV report:
///   cargo run --bin validate_bilc_solver -- --report /path/to/report.tsv /dir/*.haplotype_meta.tsv
use std::collections::BTreeSet;
use std::env;
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::process;

use haplotype_inspection::bilc_solver::lp_solve_remained_haplotypes_with_obj;
use haplotype_inspection::structs::BilcRecord;
use rustc_hash::FxHashSet;

/// Result from validating a single file.
struct ValidationResult {
    file_id: String,              // e.g. "HG002.pooled.raw.deduped.721"
    n_haps: usize,
    n_regions: usize,
    n_rows: usize,
    rust_obj: f64,                // objective (sum of dropped coefficients)
    python_obj: f64,              // reconstructed Python objective
    rust_n_drop: usize,
    python_n_drop: usize,
    ilp_match: bool,
    rust_only: BTreeSet<i32>,     // haps Rust drops but Python doesn't (ILP level)
    python_only: BTreeSet<i32>,   // haps Python drops but Rust doesn't (ILP level)
    full_match: bool,             // after post-ILP augmentation
    n_post_ilp_adds: usize,       // haps added by post-ILP rules
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: validate_bilc_solver [--report <report.tsv>] <file1.tsv> [file2.tsv ...]");
        process::exit(1);
    }

    // Parse --report flag
    let mut report_path: Option<String> = None;
    let mut file_args: Vec<String> = Vec::new();
    let mut i = 1;
    while i < args.len() {
        if args[i] == "--report" {
            if i + 1 >= args.len() {
                eprintln!("--report requires a path argument");
                process::exit(1);
            }
            report_path = Some(args[i + 1].clone());
            i += 2;
        } else {
            file_args.push(args[i].clone());
            i += 1;
        }
    }

    if file_args.is_empty() {
        eprintln!("No input files specified.");
        process::exit(1);
    }

    let mut results: Vec<ValidationResult> = Vec::new();
    let total = file_args.len();

    for (idx, tsv_path) in file_args.iter().enumerate() {
        eprint!("\r[{}/{}] Processing: {} ...", idx + 1, total, short_name(tsv_path));
        match validate_one(tsv_path) {
            Ok(r) => results.push(r),
            Err(e) => eprintln!("\nERROR on {}: {}", tsv_path, e),
        }
    }
    eprintln!("\r{:width$}", " ", width = 120); // clear progress line

    // --- Summary ---
    let n_total = results.len();
    let n_ilp_exact = results.iter().filter(|r| r.ilp_match).count();
    let n_alt_optima = results.iter().filter(|r| !r.ilp_match && r.rust_n_drop == r.python_n_drop).count();
    let n_diff_count = results.iter().filter(|r| !r.ilp_match && r.rust_n_drop != r.python_n_drop).count();
    let n_full_match = results.iter().filter(|r| r.full_match).count();
    let n_full_mismatch = n_total - n_full_match;
    let n_same_obj = results.iter().filter(|r| !r.ilp_match && (r.rust_obj - r.python_obj).abs() < 1e-6).count();

    println!("\n=== BILC Solver Cross-Validation Report ===");
    println!("Files tested: {}", n_total);
    println!("ILP-exact match: {} / {}", n_ilp_exact, n_total);
    println!("ILP alternative optima (same #drops, same obj): {} / {}", n_same_obj, n_total);
    println!("ILP different #drops: {} / {}", n_diff_count, n_total);
    println!("Full (ILP+augmentation) match: {} / {}", n_full_match, n_total);
    println!("Full mismatches: {}", n_full_mismatch);

    // --- Print per-file detail table for ILP-differing cases ---
    let diffs: Vec<&ValidationResult> = results.iter().filter(|r| !r.ilp_match).collect();
    if !diffs.is_empty() {
        println!("\n=== ILP-differing files (detailed) ===");
        println!(
            "{:<30} {:>6} {:>6} {:>7} {:>7} {:>10} {:>10} {:>5} {:<30} {:<30}",
            "File", "#Haps", "#Rgns", "Py#Drp", "Rs#Drp", "Py_Obj", "Rs_Obj", "ObjEq", "Rust-only drops", "Python-only drops"
        );
        println!("{}", "-".repeat(160));
        for r in &diffs {
            let obj_eq = if (r.rust_obj - r.python_obj).abs() < 1e-6 { "YES" } else { "NO" };
            println!(
                "{:<30} {:>6} {:>6} {:>7} {:>7} {:>10.4} {:>10.4} {:>5} {:<30} {:<30}",
                r.file_id,
                r.n_haps,
                r.n_regions,
                r.python_n_drop,
                r.rust_n_drop,
                r.python_obj,
                r.rust_obj,
                obj_eq,
                format!("{:?}", r.rust_only),
                format!("{:?}", r.python_only),
            );
        }
        println!("{}", "-".repeat(160));
        println!(
            "Same objective value: {}/{} differing files → alternative optima confirmed",
            diffs.iter().filter(|r| (r.rust_obj - r.python_obj).abs() < 1e-6).count(),
            diffs.len()
        );
    }

    // --- Write TSV report if requested ---
    if let Some(ref path) = report_path {
        let mut f = fs::File::create(path).expect("Cannot create report file");
        writeln!(
            f,
            "file_id\tn_haps\tn_regions\tn_rows\tpython_n_drop\trust_n_drop\tpython_obj\trust_obj\tilp_match\tobj_match\tfull_match\trust_only_drops\tpython_only_drops\tn_post_ilp_adds"
        ).unwrap();
        for r in &results {
            let obj_match = (r.rust_obj - r.python_obj).abs() < 1e-6;
            writeln!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.file_id,
                r.n_haps,
                r.n_regions,
                r.n_rows,
                r.python_n_drop,
                r.rust_n_drop,
                r.python_obj,
                r.rust_obj,
                r.ilp_match,
                obj_match,
                r.full_match,
                format!("{:?}", r.rust_only),
                format!("{:?}", r.python_only),
                r.n_post_ilp_adds,
            ).unwrap();
        }
        println!("\nTSV report written to: {}", path);
    }

    if n_full_mismatch > 0 {
        process::exit(1);
    }
}

fn short_name(path: &str) -> &str {
    path.rsplit('/').next().unwrap_or(path)
}

fn validate_one(tsv_path: &str) -> Result<ValidationResult, String> {
    let file = fs::File::open(tsv_path).map_err(|e| format!("Cannot open: {}", e))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let header_line = lines
        .next()
        .ok_or("empty file")?
        .map_err(|e| format!("IO: {}", e))?;
    let headers: Vec<&str> = header_line.split('\t').collect();

    let col = |name: &str| -> Result<usize, String> {
        headers
            .iter()
            .position(|h| *h == name)
            .ok_or_else(|| format!("column '{}' not found", name))
    };
    let i_start = col("start")?;
    let i_end = col("end")?;
    let i_hap_id = col("hap_id")?;
    let i_chrom = col("chrom")?;
    let i_var_count = col("var_count")?;
    let i_varc_rank = col("varc_rank")?;
    let i_coefficient = col("coefficient")?;
    let i_mismap = col("mismap")?;
    let i_extreme_vard = col("extreme_vard")?;
    let i_hap_max_sim = col("hap_max_sim_scores")?;
    let i_scatter = col("scatter_hap")?;
    let i_hap_var_count = col("hap_var_count")?;

    let mut records: Vec<BilcRecord> = Vec::new();
    let mut python_mismap_hids: FxHashSet<i32> = FxHashSet::default();
    let mut all_hids: FxHashSet<i32> = FxHashSet::default();
    let mut post_ilp_adds: FxHashSet<i32> = FxHashSet::default();

    // Track unique (hap_id -> first coefficient) for objective computation
    let mut hap_first_coeff: rustc_hash::FxHashMap<i32, f64> = rustc_hash::FxHashMap::default();

    // Track unique regions
    let mut region_keys: FxHashSet<(String, i32, i32)> = FxHashSet::default();
    let mut row_count = 0;

    for line in lines {
        let line = line.map_err(|e| format!("IO: {}", e))?;
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();

        let hap_id: i32 = fields[i_hap_id].parse().map_err(|_| "bad hap_id")?;
        let start: i32 = fields[i_start].parse().map_err(|_| "bad start")?;
        let end: i32 = fields[i_end].parse().map_err(|_| "bad end")?;
        let var_count: i32 = fields[i_var_count].parse().map_err(|_| "bad var_count")?;
        let varc_rank: i32 = fields[i_varc_rank].parse().map_err(|_| "bad varc_rank")?;
        let coefficient: f64 = fields[i_coefficient].parse().map_err(|_| "bad coefficient")?;
        let chrom = fields[i_chrom].to_string();

        let extreme_vard = fields[i_extreme_vard] == "True";
        let hap_max_sim: f64 = fields[i_hap_max_sim].parse().map_err(|_| "bad hap_max_sim")?;
        let scatter = fields[i_scatter] == "True";
        let hap_var_count: i32 = fields[i_hap_var_count].parse().map_err(|_| "bad hap_var_count")?;

        records.push(BilcRecord {
            chrom: chrom.clone(),
            start,
            end,
            hap_id,
            coefficient,
            var_count,
            varc_rank,
        });

        all_hids.insert(hap_id);
        hap_first_coeff.entry(hap_id).or_insert(coefficient);
        region_keys.insert((chrom, start, end));

        if fields[i_mismap] == "True" {
            python_mismap_hids.insert(hap_id);
        }

        // Post-ILP augmentation rules
        if extreme_vard {
            post_ilp_adds.insert(hap_id);
        }
        if hap_max_sim > 15.0 {
            post_ilp_adds.insert(hap_id);
        }
        if scatter && hap_var_count <= 3 {
            post_ilp_adds.insert(hap_id);
        }
        row_count += 1;
    }

    // Reconstruct Python ILP-only drop set
    let python_ilp_drop: FxHashSet<i32> = python_mismap_hids
        .difference(&post_ilp_adds)
        .copied()
        .collect();

    // Compute Python objective value from ILP-only drops
    let python_obj: f64 = python_ilp_drop
        .iter()
        .map(|hid| hap_first_coeff.get(hid).copied().unwrap_or(0.0))
        .sum();

    // Run Rust solver
    let (rust_select, rust_drop, _status, rust_obj) =
        lp_solve_remained_haplotypes_with_obj(&records);

    let ilp_match = rust_drop == python_ilp_drop;

    let rust_only: BTreeSet<i32> = rust_drop.difference(&python_ilp_drop).copied().collect();
    let python_only: BTreeSet<i32> = python_ilp_drop.difference(&rust_drop).copied().collect();

    // Full augmented comparison
    let mut rust_full_drop = rust_drop.clone();
    rust_full_drop.extend(&post_ilp_adds);
    let full_match = rust_full_drop == python_mismap_hids;

    // Extract short file ID
    let file_id = short_name(tsv_path)
        .trim_end_matches(".haplotype_meta.tsv")
        .to_string();

    Ok(ValidationResult {
        file_id,
        n_haps: all_hids.len(),
        n_regions: region_keys.len(),
        n_rows: row_count,
        rust_obj,
        python_obj,
        rust_n_drop: rust_drop.len(),
        python_n_drop: python_ilp_drop.len(),
        ilp_match,
        rust_only,
        python_only,
        full_match,
        n_post_ilp_adds: post_ilp_adds.len(),
    })
}
