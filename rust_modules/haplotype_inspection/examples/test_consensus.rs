/// Standalone integration test for consensus vector extraction and validation.
///
/// Loads a real BAM via Lapper, queries reads overlapping a user-specified interval,
/// extracts hap/err vectors for each read, validates them against the raw CIGAR,
/// groups reads by HP tag, and tests assemble_consensus per HP group.
///
/// Usage:
///   test_consensus <bam_file> <chrom:start-end>
///   RUST_LOG=debug test_consensus input.bam chr1:1633000-1635000
// Use the library crate for cross-module imports
use haplotype_inspection::bam_lappers::{build_lapper_from_bam, query_overlapping_reads};
use haplotype_inspection::identify_misaligned_haps::{
    assemble_consensus, record_hap_err_vectors_per_region,
};
use haplotype_inspection::pairwise_read_inspection::{
    extract_error_vector, extract_hap_vector, is_deletion, is_snv, read_id, HAP_PAD,
};
use ndarray::Array1;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;
use std::collections::HashMap;
use std::env;
use std::io::Write;

fn main() {
    // Initialize env_logger for standalone logging (reads RUST_LOG env var)
    env_logger::init();

    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: {} <bam_file> <chrom:start-end>", args[0]);
        eprintln!("\nExample:");
        eprintln!("  {} input.bam chr1:1633000-1635000", args[0]);
        eprintln!(
            "  RUST_LOG=debug {} input.bam chr1:1633000-1635000",
            args[0]
        );
        std::process::exit(1);
    }

    let bam_path = &args[1];
    let (chrom, start, end) = parse_interval(&args[2]);

    println!("=== Consensus Vector Integration Test ===");
    println!("BAM file: {bam_path}");
    println!("Region: {chrom}:{start}-{end}");
    println!();

    // ── Step 1: Build Lapper from BAM ──────────────────────────────────────
    println!("Building Lapper from BAM...");
    let build_start = std::time::Instant::now();

    let result = match build_lapper_from_bam(bam_path, 0, 0, false, false) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Error building Lapper: {e}");
            std::process::exit(1);
        }
    };

    let build_time = build_start.elapsed();
    println!(
        "Lapper built in {:.2}s ({} qnames retained)",
        build_time.as_secs_f64(),
        result.stats.qnames_retained
    );
    println!();

    // ── Step 2: Query overlapping reads ────────────────────────────────────
    let reads = query_overlapping_reads(&result.lapper_dict, &result.read_dict, &chrom, start, end);
    println!(
        "Queried {}:{}-{} → {} reads",
        chrom,
        start,
        end,
        reads.len()
    );
    println!();

    if reads.is_empty() {
        println!("No reads found in region. Nothing to test.");
        return;
    }

    // ── Step 3: Per-read vector extraction + CIGAR validation ──────────────
    // Print raw CIGAR, hap vector, error vector, base qual for each read
    println!("=== Per-Read Vector Extraction & CIGAR Validation ===");
    let extract_start = std::time::Instant::now();

    let mut total_passed = 0usize;
    let mut total_failed = 0usize;

    for &read in &reads {
        let hap = extract_hap_vector(read).expect("CIGAR must be in =/X mode (--eqx)");
        let err = extract_error_vector(read).expect("CIGAR must be in =/X mode (--eqx)");

        let (hap_ok, hap_details) = validate_hap_against_cigar(read, &hap);
        let (err_ok, err_details) = validate_err_against_cigar(read, &err);

        let passed = hap_ok && err_ok;
        if passed {
            total_passed += 1;
        } else {
            total_failed += 1;
        }

        let qname = String::from_utf8_lossy(read.qname());
        let hp_tag = get_hp_tag(read);
        let status = if passed { "PASS" } else { "FAIL" };
        println!(
            "  [{}] {} flag={} HP={} ref={}..{}",
            status,
            qname,
            read.flags(),
            hp_tag,
            read.pos(),
            read.cigar().end_pos()
        );

        // Print raw CIGAR string
        println!("    CIGAR: {}", cigar_to_string(read));

        // Print raw base quality string (ASCII Phred+33)
        let qual = read.qual();
        let qual_str: String = qual.iter().map(|&q| (q + 33) as char).collect();
        println!("    QUAL:  {qual_str}");

        // Print raw hap vector (compact: show first 80 values, then "...")
        println!("    HAP[{}]: {}", hap.len(), format_vector_i16(&hap, 80));

        // Print raw error vector (compact: first 40 values, then "...")
        println!("    ERR[{}]: {}", err.len(), format_vector_f32(&err, 40));

        if !hap_ok {
            println!("    hap_fail: {hap_details}");
        }
        if !err_ok {
            println!("    err_fail: {err_details}");
        }
    }

    let extract_time = extract_start.elapsed();
    println!();
    println!(
        "Per-read extraction: {} passed, {} failed ({:.3}s)",
        total_passed,
        total_failed,
        extract_time.as_secs_f64()
    );
    println!();

    // ── Step 4: Test record_hap_err_vectors_per_region ─────────────────────
    println!("=== Batch Extraction: record_hap_err_vectors_per_region ===");
    let batch_start = std::time::Instant::now();

    let mut hap_cache: HashMap<String, Array1<i16>> = HashMap::new();
    let mut err_cache: HashMap<String, Array1<f32>> = HashMap::new();

    let (spans, hap_arrays, err_arrays) =
        record_hap_err_vectors_per_region(&reads, &mut hap_cache, &mut err_cache)
            .expect("CIGAR must be in =/X mode (--eqx)");

    let batch_time = batch_start.elapsed();

    // Verify shapes
    let n_reads = reads.len();
    let span_shape_ok = spans.nrows() == n_reads && spans.ncols() == 2;
    let hap_shape_ok = hap_arrays.nrows() == n_reads;
    let err_shape_ok = err_arrays.nrows() == n_reads;
    let max_len = hap_arrays.ncols();

    println!(
        "  spans shape: ({}, {}) {}",
        spans.nrows(),
        spans.ncols(),
        if span_shape_ok { "OK" } else { "MISMATCH" }
    );
    println!(
        "  hap_arrays shape: ({}, {}) {}",
        hap_arrays.nrows(),
        hap_arrays.ncols(),
        if hap_shape_ok { "OK" } else { "MISMATCH" }
    );
    println!(
        "  err_arrays shape: ({}, {}) {}",
        err_arrays.nrows(),
        err_arrays.ncols(),
        if err_shape_ok { "OK" } else { "MISMATCH" }
    );

    // Verify each row's non-padding length matches individual extraction
    let mut row_len_mismatches = 0usize;
    for (i, &read) in reads.iter().enumerate() {
        let individual_hap = extract_hap_vector(read).expect("CIGAR must be in =/X mode (--eqx)");
        let individual_len = individual_hap.len();

        let row = hap_arrays.row(i);
        let row_nonpad = row.iter().filter(|&&v| v != HAP_PAD).count();

        if row_nonpad != individual_len {
            row_len_mismatches += 1;
            let qname = String::from_utf8_lossy(read.qname());
            println!(
                "  ROW MISMATCH: read {qname} row_nonpad={row_nonpad} individual_len={individual_len}"
            );
        }
    }

    // Verify cache was populated
    let unique_rids: std::collections::HashSet<String> =
        reads.iter().map(|&r| read_id(r)).collect();
    let cache_size_ok = hap_cache.len() == unique_rids.len();
    let err_cache_ok = err_cache.len() == unique_rids.len();

    println!(
        "  row_length_mismatches: {} {}",
        row_len_mismatches,
        if row_len_mismatches == 0 {
            "OK"
        } else {
            "FAIL"
        }
    );
    println!(
        "  hap_cache_size: {} (expected {}) {}",
        hap_cache.len(),
        unique_rids.len(),
        if cache_size_ok { "OK" } else { "MISMATCH" }
    );
    println!(
        "  err_cache_size: {} (expected {}) {}",
        err_cache.len(),
        unique_rids.len(),
        if err_cache_ok { "OK" } else { "MISMATCH" }
    );
    println!("  max_vec_len: {max_len}");
    println!("  batch_time: {:.3}s", batch_time.as_secs_f64());
    println!();

    // ── Step 5: Assemble consensus per HP tag ──────────────────────────────
    println!("=== Assemble Consensus per HP Tag ===");

    // Group read indices by HP tag
    let mut hp_groups: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, &read) in reads.iter().enumerate() {
        let hp = get_hp_tag(read);
        hp_groups.entry(hp).or_default().push(i);
    }

    let mut hp_tags_sorted: Vec<&String> = hp_groups.keys().collect();
    hp_tags_sorted.sort();

    let mut consensus_ok = true;
    for hp_tag in &hp_tags_sorted {
        let indices = &hp_groups[*hp_tag];
        println!("  HP={} ({} reads)", hp_tag, indices.len());

        if indices.len() < 2 {
            println!("    Skipping (need >= 2 reads for consensus)");
            continue;
        }

        // Collect reads for this HP group
        let hp_reads: Vec<&Record> = indices.iter().map(|&i| reads[i]).collect();

        // Build batch arrays for this group
        let mut hp_hap_cache: HashMap<String, Array1<i16>> = HashMap::new();
        let mut hp_err_cache: HashMap<String, Array1<f32>> = HashMap::new();
        let (hp_spans, hp_hap_arrays, hp_err_arrays) =
            record_hap_err_vectors_per_region(&hp_reads, &mut hp_hap_cache, &mut hp_err_cache)
                .expect("CIGAR must be in =/X mode (--eqx)");

        // Run assemble_consensus
        let consensus = assemble_consensus(&hp_hap_arrays, &hp_err_arrays, &hp_spans);

        // Print consensus stats
        let consensus_len = consensus.len();
        let n_match = consensus.iter().filter(|&&v| v == 1).count();
        let n_snv = consensus.iter().filter(|&&v| is_snv(v)).count();
        let n_del = consensus.iter().filter(|&&v| is_deletion(v)).count();
        let n_ins = consensus.iter().filter(|&&v| v > 1).count();

        // Compute global ref span for this HP group
        let min_start = (0..hp_spans.nrows())
            .map(|r| hp_spans[[r, 0]])
            .min()
            .unwrap_or(0);
        let max_end = (0..hp_spans.nrows())
            .map(|r| hp_spans[[r, 1]])
            .max()
            .unwrap_or(0);

        println!(
            "    consensus_len={consensus_len} ref_span={min_start}..{max_end} match={n_match} snv={n_snv} del={n_del} ins_marker={n_ins}"
        );
        println!(
            "    consensus[{}]: {}",
            consensus_len,
            format_vector_i16(&consensus, 120)
        );

        // Basic sanity: consensus length should equal max_end - min_start
        let expected_len = (max_end - min_start) as usize;
        if consensus_len != expected_len {
            println!(
                "    WARNING: consensus_len {consensus_len} != expected {expected_len} (max_end - min_start)"
            );
            consensus_ok = false;
        } else {
            println!("    length check: OK");
        }
    }
    println!();

    // ── Summary ────────────────────────────────────────────────────────────
    let all_ok = total_failed == 0
        && span_shape_ok
        && hap_shape_ok
        && err_shape_ok
        && row_len_mismatches == 0
        && cache_size_ok
        && err_cache_ok
        && consensus_ok;

    println!("=== Summary ===");
    println!("Reads processed: {n_reads}");
    println!("Per-read CIGAR validation: {total_passed} passed, {total_failed} failed");
    println!(
        "Batch shape checks: {}",
        if span_shape_ok && hap_shape_ok && err_shape_ok {
            "OK"
        } else {
            "FAIL"
        }
    );
    println!(
        "Row length checks: {}",
        if row_len_mismatches == 0 {
            "OK"
        } else {
            "FAIL"
        }
    );
    println!(
        "Cache checks: {}",
        if cache_size_ok && err_cache_ok {
            "OK"
        } else {
            "FAIL"
        }
    );
    println!(
        "HP-tag consensus: {} groups tested, {}",
        hp_tags_sorted.len(),
        if consensus_ok { "OK" } else { "SOME WARNINGS" }
    );
    println!(
        "Overall: {}",
        if all_ok {
            "ALL PASSED"
        } else {
            "SOME FAILURES"
        }
    );

    // Flush stdout
    std::io::stdout().flush().unwrap();

    if !all_ok {
        std::process::exit(1);
    }
}

// ─── Helper Functions ────────────────────────────────────────────────────────

/// Get the HP tag value from a BAM record, or "no_HP" if absent.
fn get_hp_tag(record: &Record) -> String {
    use rust_htslib::bam::record::Aux;
    match record.aux(b"HP") {
        Ok(Aux::String(s)) => s.to_string(),
        Ok(Aux::I32(v)) => format!("{v}"),
        Ok(Aux::I16(v)) => format!("{v}"),
        Ok(Aux::I8(v)) => format!("{v}"),
        Ok(Aux::U32(v)) => format!("{v}"),
        Ok(Aux::U16(v)) => format!("{v}"),
        Ok(Aux::U8(v)) => format!("{v}"),
        _ => "no_HP".to_string(),
    }
}

/// Convert a Record's CIGAR to a printable string (e.g., "137=1I10=").
fn cigar_to_string(record: &Record) -> String {
    let cigar = record.cigar();
    let mut s = String::new();
    for c in cigar.iter() {
        match c {
            Cigar::Match(n) => s.push_str(&format!("{n}M")),
            Cigar::Ins(n) => s.push_str(&format!("{n}I")),
            Cigar::Del(n) => s.push_str(&format!("{n}D")),
            Cigar::RefSkip(n) => s.push_str(&format!("{n}N")),
            Cigar::SoftClip(n) => s.push_str(&format!("{n}S")),
            Cigar::HardClip(n) => s.push_str(&format!("{n}H")),
            Cigar::Pad(n) => s.push_str(&format!("{n}P")),
            Cigar::Equal(n) => s.push_str(&format!("{n}=")),
            Cigar::Diff(n) => s.push_str(&format!("{n}X")),
        }
    }
    s
}

/// Format an i16 vector compactly: show up to `limit` values, then "...".
fn format_vector_i16(v: &Array1<i16>, limit: usize) -> String {
    let len = v.len();
    if len <= limit {
        let vals: Vec<String> = v.iter().map(|x| x.to_string()).collect();
        format!("[{}]", vals.join(","))
    } else {
        let vals: Vec<String> = v.iter().take(limit).map(|x| x.to_string()).collect();
        format!("[{},...({} more)]", vals.join(","), len - limit)
    }
}

/// Format an f32 vector compactly: show up to `limit` values (3 decimal places), then "...".
fn format_vector_f32(v: &Array1<f32>, limit: usize) -> String {
    let len = v.len();
    if len <= limit {
        let vals: Vec<String> = v.iter().map(|x| format!("{x:.3}")).collect();
        format!("[{}]", vals.join(","))
    } else {
        let vals: Vec<String> = v.iter().take(limit).map(|x| format!("{x:.3}")).collect();
        format!("[{},...({} more)]", vals.join(","), len - limit)
    }
}

// ─── CIGAR Validation Functions ──────────────────────────────────────────────

/// Walk CIGAR and validate hap_vector consistency.
///
/// For each CIGAR operation, checks that the corresponding hap_vector positions
/// contain the expected golden encoding (deletion = -10). Accounts for insertion
/// markers summed onto the first position of the next ref-consuming operation (`> 1`).
///
/// Returns `(passed, details)`.
fn validate_hap_against_cigar(record: &Record, hap: &Array1<i16>) -> (bool, String) {
    let cigar = record.cigar();

    // Compute expected ref_len from CIGAR
    let mut expected_ref_len: usize = 0;
    for c in cigar.iter() {
        match c {
            Cigar::Match(_) => {
                return (
                    false,
                    "CIGAR contains M (op 0); requires =/X mode".to_string(),
                );
            }
            Cigar::Equal(len) | Cigar::Diff(len) | Cigar::Del(len) | Cigar::RefSkip(len) => {
                expected_ref_len += *len as usize;
            }
            _ => {}
        }
    }

    // Length check
    if hap.len() != expected_ref_len {
        return (
            false,
            format!(
                "length mismatch: hap_len={} expected_ref_len={}",
                hap.len(),
                expected_ref_len
            ),
        );
    }

    // Content check: walk CIGAR, tracking ref_pos and pending insertion state
    let mut ref_pos: usize = 0;
    let mut pending_insertion = false;
    let mut errors: Vec<String> = Vec::new();

    for c in cigar.iter() {
        match c {
            Cigar::Equal(len) | Cigar::RefSkip(len) => {
                let n = *len as usize;
                for i in 0..n {
                    let pos = ref_pos + i;
                    if pos >= hap.len() {
                        break;
                    }
                    let val = hap[pos];
                    if i == 0 && pending_insertion {
                        if val <= 1 && val != 1 {
                            errors.push(format!("pos {pos} expected ins_marker or 1, got {val}"));
                        }
                        pending_insertion = false;
                    } else if val != 1 {
                        errors.push(format!("pos {pos} expected 1 (=), got {val}"));
                    }
                }
                ref_pos += n;
            }
            Cigar::Diff(len) => {
                let n = *len as usize;
                for i in 0..n {
                    let pos = ref_pos + i;
                    if pos >= hap.len() {
                        break;
                    }
                    let val = hap[pos];
                    if i == 0 && pending_insertion {
                        if val <= 1 {
                            errors.push(format!("pos {pos} expected ins_marker (>1), got {val}"));
                        }
                        pending_insertion = false;
                    } else if val != -4 {
                        errors.push(format!("pos {pos} expected -4 (X), got {val}"));
                    }
                }
                ref_pos += n;
            }
            Cigar::Del(len) => {
                let n = *len as usize;
                for i in 0..n {
                    let pos = ref_pos + i;
                    if pos >= hap.len() {
                        break;
                    }
                    let val = hap[pos];
                    if i == 0 && pending_insertion {
                        if val <= 1 {
                            errors.push(format!("pos {pos} expected ins_marker (>1), got {val}"));
                        }
                        pending_insertion = false;
                    } else if val != -10 {
                        errors.push(format!("pos {pos} expected -10 (D), got {val}"));
                    }
                }
                ref_pos += n;
            }
            Cigar::Ins(len) => {
                if ref_pos > 0 {
                    pending_insertion = true;
                }
                let _ = len;
            }
            Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
            _ => {}
        }
    }

    if errors.is_empty() {
        (true, String::new())
    } else {
        let detail = if errors.len() <= 5 {
            errors.join("; ")
        } else {
            format!("{} (and {} more)", errors[..5].join("; "), errors.len() - 5)
        };
        (false, detail)
    }
}

/// Walk CIGAR and validate err_vector consistency.
///
/// Checks:
/// - Length matches expected ref_len
/// - D/N positions → 0.0
/// - =/X positions → (0.0, 1.0] unless overwritten by an adjacent insertion
///
/// Returns `(passed, details)`.
fn validate_err_against_cigar(record: &Record, err: &Array1<f32>) -> (bool, String) {
    let cigar = record.cigar();

    let mut expected_ref_len: usize = 0;
    for c in cigar.iter() {
        match c {
            Cigar::Match(_) => {
                return (
                    false,
                    "CIGAR contains M (op 0); requires =/X mode".to_string(),
                );
            }
            Cigar::Equal(len) | Cigar::Diff(len) | Cigar::Del(len) | Cigar::RefSkip(len) => {
                expected_ref_len += *len as usize;
            }
            _ => {}
        }
    }

    if err.len() != expected_ref_len {
        return (
            false,
            format!(
                "length mismatch: err_len={} expected_ref_len={}",
                err.len(),
                expected_ref_len
            ),
        );
    }

    // Collect positions that insertions overwrite (ref_pos - 1 at time of insertion)
    let mut insertion_overwrite_positions: std::collections::HashSet<usize> =
        std::collections::HashSet::new();
    {
        let mut ref_pos: usize = 0;
        for c in cigar.iter() {
            match c {
                Cigar::Equal(len) | Cigar::Diff(len) => {
                    ref_pos += *len as usize;
                }
                Cigar::Del(len) | Cigar::RefSkip(len) => {
                    ref_pos += *len as usize;
                }
                Cigar::Ins(_) => {
                    if ref_pos > 0 {
                        insertion_overwrite_positions.insert(ref_pos - 1);
                    }
                }
                _ => {}
            }
        }
    }

    let mut ref_pos: usize = 0;
    let mut errors: Vec<String> = Vec::new();

    for c in cigar.iter() {
        match c {
            Cigar::Equal(len) | Cigar::Diff(len) => {
                let n = *len as usize;
                for i in 0..n {
                    let pos = ref_pos + i;
                    if pos >= err.len() {
                        break;
                    }
                    let val = err[pos];
                    if insertion_overwrite_positions.contains(&pos) {
                        if val != 0.0 {
                            errors
                                .push(format!("pos {pos} expected 0.0 (ins overwrite), got {val}"));
                        }
                    } else if val <= 0.0 || val > 1.0 {
                        errors.push(format!("pos {pos} expected (0.0, 1.0] (=/X), got {val}"));
                    }
                }
                ref_pos += n;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                let n = *len as usize;
                for i in 0..n {
                    let pos = ref_pos + i;
                    if pos >= err.len() {
                        break;
                    }
                    if err[pos] != 0.0 {
                        errors.push(format!("pos {} expected 0.0 (D/N), got {}", pos, err[pos]));
                    }
                }
                ref_pos += n;
            }
            _ => {}
        }
    }

    if errors.is_empty() {
        (true, String::new())
    } else {
        let detail = if errors.len() <= 5 {
            errors.join("; ")
        } else {
            format!("{} (and {} more)", errors[..5].join("; "), errors.len() - 5)
        };
        (false, detail)
    }
}

// ─── Argument Parsing ────────────────────────────────────────────────────────

/// Parse an interval string like "chr1:1633000-1635000"
fn parse_interval(s: &str) -> (String, u32, u32) {
    let (chrom, range) = s.split_once(':').unwrap_or_else(|| {
        eprintln!("Invalid interval format: '{s}'. Expected chrom:start-end");
        std::process::exit(1);
    });
    let (start_str, end_str) = range.split_once('-').unwrap_or_else(|| {
        eprintln!("Invalid interval format: '{s}'. Expected chrom:start-end");
        std::process::exit(1);
    });
    let start: u32 = start_str.parse().unwrap_or_else(|_| {
        eprintln!("Invalid start position: '{start_str}'");
        std::process::exit(1);
    });
    let end: u32 = end_str.parse().unwrap_or_else(|_| {
        eprintln!("Invalid end position: '{end_str}'");
        std::process::exit(1);
    });
    (chrom.to_string(), start, end)
}
