//! Haplotype misalignment identification.
//!
//! Ports Python functions from `identify_misaligned_haps.py`:
//! - `record_hap_err_vectors_per_region` → batch vector extraction
//! - `assemble_consensus` → consensus assembly
//! - `judge_misalignment_by_extreme_vardensity` → variant density filter
//! - `stat_refseq_similarity` → reference sequence similarity
//! - `cal_similarity_score` → similarity scoring
//! - `extract_continuous_regions_dict` → region grouping
//! - `group_by_dict_optimized` → haplotype grouping
//! - `record_haplotype_rank` → haplotype ranking
//! - `summarize_enclosing_haps` → enclosing haplotype summary
//! - `identify_misalignment_per_region` → per-region orchestrator
//! - `inspect_haplotypes` → main inspection loop

use ndarray::{Array1, Array2};
use rust_htslib::bam::Record;
use rust_lapper::Lapper;
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::{HashMap, HashSet};
use std::collections::hash_map::Entry;
use log::{debug, info, warn};

use crate::bam_lappers::{BamLapperResult, build_lapper_from_bam, query_overlapping_reads};
use crate::pairwise_read_inspection::{
    read_id, extract_hap_vector, extract_error_vector,
    extract_read_qseqs, ReadQseqData,
    count_var, count_continuous_indel_blocks, count_snv,
    HAP_PAD, is_snv, is_indel, is_deletion, insertion_len, is_real,
};
use crate::structs::{BilcRecord, EnrichedRecord, RegionKey};
use crate::bilc_solver::{lp_solve_remained_haplotypes, BilcStatus};


// ─── Data structures for stat_refseq_similarity ─────────────────────────────

/// Per-overlap-region statistics tuple, matching Python's:
/// `(varcount, alt_snv_count, alt_indel_count, shared_psv,
///   verified_shared_snv_pos_abs, verified_shared_indel_pos_abs, overlap_span_size)`
#[derive(Clone, Debug)]
pub struct RegionVarStats {
    pub varcount: i32,
    pub alt_snv_count: i32,
    pub alt_indel_count: i32,
    pub shared_psv: i32,
    pub verified_shared_snv_pos_abs: Vec<i32>,
    pub verified_shared_indel_pos_abs: Vec<i32>,
    pub overlap_span_size: i64,
}

/// Type alias for the nested dict: hid -> homo_refseq_qname -> Vec<RegionVarStats>
pub type VarcountsAmongRefseqs = HashMap<i32, HashMap<String, Vec<RegionVarStats>>>;

// ============================================================================
// Variant Density Filtering
// Ports from identify_misaligned_haps.py
// ============================================================================

/// Calculate sliding window variant density for each position in the array.
///
/// Ports from identify_misaligned_haps.py lines 347-362
///
/// For each position i, calculates the variant density in the window
/// [i-padding_size, i+padding_size]. Density = variant_count / window_size.
///
/// Returns an array of the same length as input, with density values at each position.
///
/// Implementation note: instead of re-slicing and re-counting each window
/// (`O(n * window)`), this builds two prefix-sum arrays once and answers each
/// window in O(1) (`O(n)` overall, no per-position allocation). The result is
/// identical to `count_var` over each clamped window divided by the fixed
/// `window_size`. See `count_window_var_density_reference` in the tests for the
/// straightforward version this is cross-checked against.
pub fn count_window_var_density(array: &Array1<i16>, padding_size: i32) -> Array1<f32> {
    let n = array.len();

    // Check if there are any variants (anything != 1)
    let has_variants = array.iter().any(|&v| v != 1);
    if !has_variants {
        return Array1::zeros(n);
    }

    let pad = padding_size as usize;
    let window_size = (padding_size * 2 + 1) as f32;

    // Prefix sums computed once over the whole sequence (golden decoders, so this
    // stays identical to count_var over each window):
    //   snv_prefix[k]   = count of SNV positions (pure -4 OR compound mismatch-ins
    //                     ending in digit 6) in array[0..k]  — matches count_snv.
    //   block_prefix[k] = count of indel-block STARTS in array[0..k], where a start
    //                     is an indel position (deletion -10 or insertion >1) whose
    //                     left neighbour is not an indel — mirrors count_continuous_blocks.
    let mut snv_prefix = vec![0i32; n + 1];
    let mut block_prefix = vec![0i32; n + 1];
    let mut prev_indel = false;
    for (k, &v) in array.iter().enumerate() {
        let snv = is_snv(v);
        let indel = is_indel(v);
        let is_block_start = indel && !prev_indel;
        snv_prefix[k + 1] = snv_prefix[k] + snv as i32;
        block_prefix[k + 1] = block_prefix[k] + is_block_start as i32;
        prev_indel = indel;
    }

    let is_indel_at = |k: usize| -> bool { is_indel(array[k]) };

    let mut density_arr = Array1::<f32>::zeros(n);
    for i in 0..n {
        let start = i.saturating_sub(pad);
        let end = (i + pad + 1).min(n);

        // SNVs in [start, end): plain prefix difference.
        let snv_count = snv_prefix[end] - snv_prefix[start];

        // Indel blocks in [start, end): block starts strictly inside the window
        // (positions start+1..end) come from block_prefix; the window's own
        // first position counts as a fresh block start if it is an indel
        // (count_continuous_blocks treats each window in isolation).
        let internal_starts = block_prefix[end] - block_prefix[start + 1];
        let boundary_start = is_indel_at(start) as i32;
        let blocks = internal_starts + boundary_start;

        let var_count = (snv_count + blocks) as f32;
        density_arr[i] = var_count / window_size;
    }

    density_arr
}

/// Extract continuous stretches of True values from a boolean array.
/// Returns a Vec of (start_idx, end_idx) tuples (inclusive).
///
/// Ports from numba_operators.py lines 141-170
pub fn extract_true_stretches(bool_array: &Array1<bool>) -> Vec<(usize, usize)> {
    let n = bool_array.len();
    if n == 0 {
        return Vec::new();
    }

    let mut stretches = Vec::new();
    let mut in_stretch = false;
    let mut start_idx = 0;

    for i in 0..n {
        if bool_array[i] {
            if !in_stretch {
                in_stretch = true;
                start_idx = i;
            }
        } else if in_stretch {
            stretches.push((start_idx, i - 1));
            in_stretch = false;
        }
    }

    // Handle the case where the array ends with a True stretch
    if in_stretch {
        stretches.push((start_idx, n - 1));
    }

    stretches
}

/// Judge whether a haplotype sequence is misaligned based on extreme variant density.
///
/// Ports from identify_misaligned_haps.py lines 369-421
///
/// Returns (is_misaligned, max_density):
/// - is_misaligned: true if the sequence shows signs of misalignment
/// - max_density: the maximum local variant density observed across all window sizes
///
/// Algorithm:
/// 1. Calculate variant density using three different window sizes (padding: 42, 65, 74)
/// 2. Track the maximum density across all windows
/// 3. Check three thresholds:
///    - If density >= 5/85 (padding=42) and has >=1 indel block → misaligned
///    - If density >= 6/131 (padding=65) and has >=1 indel block → misaligned
///    - If density >= 10/148 (padding=74) → misaligned
pub fn judge_misalignment_by_extreme_vardensity(seq: &Array1<i16>) -> (bool, f32) {
    // Calculate variant density for three different window sizes
    let five_vard = count_window_var_density(seq, 42);
    let six_vard = count_window_var_density(seq, 65);
    let read_vard = count_window_var_density(seq, 74);

    // Track the overall max local density across all three window sizes.
    // Densities are var-count / window-size (finite); a stray non-finite value is
    // ignored (with a debug note) rather than panicking in `partial_cmp().unwrap()`.
    let mut max_density = 0.0f32;
    for vard in [&five_vard, &six_vard, &read_vard] {
        let mut dropped_non_finite = false;
        let local_max = vard
            .iter()
            .copied()
            .filter(|x| {
                let ok = x.is_finite();
                if !ok {
                    dropped_non_finite = true;
                }
                ok
            })
            .reduce(f32::max);
        if dropped_non_finite {
            debug!("[judge_misalignment_by_extreme_vardensity] ignored non-finite density value(s)");
        }
        if let Some(d) = local_max {
            if d > max_density {
                max_density = d;
            }
        }
    }

    // Check threshold 1: 5/85 with padding=42
    let threshold_1 = 5.0 / 85.0;
    if five_vard.iter().any(|&d| d >= threshold_1) {
        let select_bool = five_vard.mapv(|d| d >= threshold_1);
        let padding = 42;
        let true_segments = extract_true_stretches(&select_bool);

        let mut max_indel_count = 0;
        for (seg_start, seg_end) in true_segments {
            let start = seg_start.saturating_sub(padding);
            let end = (seg_end + padding + 1).min(seq.len());
            let window_seq = seq.slice(ndarray::s![start..end]);
            let indel_count = count_continuous_indel_blocks(&window_seq.to_owned());
            max_indel_count = max_indel_count.max(indel_count);
        }

        if max_indel_count >= 1 {
            return (true, max_density);
        }
    }

    // Check threshold 2: 6/131 with padding=65
    let threshold_2 = 6.0 / 131.0;
    if six_vard.iter().any(|&d| d >= threshold_2) {
        let select_bool = six_vard.mapv(|d| d >= threshold_2);
        let padding = 65;
        let true_segments = extract_true_stretches(&select_bool);

        let mut max_indel_count = 0;
        for (seg_start, seg_end) in true_segments {
            let start = seg_start.saturating_sub(padding);
            let end = (seg_end + padding + 1).min(seq.len());
            let window_seq = seq.slice(ndarray::s![start..end]);
            let indel_count = count_continuous_indel_blocks(&window_seq.to_owned());
            max_indel_count = max_indel_count.max(indel_count);
        }

        if max_indel_count >= 1 {
            return (true, max_density);
        }
    }

    // Check threshold 3: 10/148 with padding=74
    let threshold_3 = 10.0 / 148.0;
    if read_vard.iter().any(|&d| d >= threshold_3) {
        return (true, max_density);
    }

    (false, max_density)
}

// ============================================================================
// Reference Similarity Helper Functions
// Ports from identify_misaligned_haps.py — used by stat_refseq_similarity
// ============================================================================

/// Merge two sorted i32 arrays into a single sorted array with unique values.
///
/// Ports from identify_misaligned_haps.py lines 244-270
/// Replaces: `np.unique(np.concatenate([arr1, arr2]))`
pub fn merge_unique_sorted(arr1: &[i32], arr2: &[i32]) -> Vec<i32> {
    let total = arr1.len() + arr2.len();
    if total == 0 {
        return Vec::new();
    }

    // Merge into one sorted vec, then deduplicate
    let mut merged = Vec::with_capacity(total);
    merged.extend_from_slice(arr1);
    merged.extend_from_slice(arr2);
    merged.sort_unstable();
    merged.dedup();
    merged
}

/// Compute basic reference-genome similarity metrics for a query haplotype vs
/// a homologous genomic haplotype vector.
///
/// Ports from identify_misaligned_haps.py lines 31-53 (`ref_genome_similarity`)
///
/// Returns `(var_count, alt_snv_count, alt_indel_count)`:
/// - `var_count`: total variant count in the query (consensus) vector
/// - `alt_snv_count`: SNV count in the genomic haplotype vector
/// - `alt_indel_count`: continuous indel block count in the genomic haplotype vector
///
/// If the genomic vector is all-reference (all 1s), returns `(0, 0, 0)`.
pub fn ref_genome_similarity(
    query_read_vector: &Array1<i16>,
    genomic_hap_vector: &Array1<i16>,
) -> (i32, i32, i32) {
    // If the genomic hap vector has no variants (all positions == 1), short-circuit
    if genomic_hap_vector.iter().all(|&v| v == 1) {
        return (0, 0, 0);
    }

    let alt_indel_count = count_continuous_indel_blocks(genomic_hap_vector);
    let alt_snv_count = count_snv(genomic_hap_vector);
    let var_size = count_var(query_read_vector);

    (var_size, alt_snv_count, alt_indel_count)
}

/// Find absolute genomic positions of shared variants between two haplotype
/// vectors aligned over the same overlap span.
///
/// Ports from identify_misaligned_haps.py lines 57-108 (`numba_shared_variant_positions`)
///
/// Haplotype vector encoding (golden — decoded via the shared helpers):
/// - `1` = match/reference
/// - `-4` = SNV (pure); a compound mismatch-under-insertion (ends in digit 6) is also an SNV
/// - `-10` = deletion (spans multiple positions)
/// - `>1` = insertion present (`base + 10*L`); the inserted length is `insertion_len`
/// - `-20` = padding / NA ([`HAP_PAD`])
///
/// Returns `(shared_snv_pos, shared_ins_pos, shared_del_pos)`:
/// - `shared_snv_pos`: absolute positions where both vectors are SNVs (pure or compound)
/// - `shared_ins_pos`: absolute positions where both vectors carry an insertion of the
///   **same length** (length match made explicit, not raw-value equality, because the
///   summed point signal differs between match-ins and mismatch-ins of equal length)
/// - `shared_del_pos`: absolute positions of shared deletion spans (exact same start+end in both vectors)
pub fn numba_shared_variant_positions(
    vec1: &Array1<i16>,
    vec2: &Array1<i16>,
    overlap_start: i32,
) -> (Vec<i32>, Vec<i32>, Vec<i32>) {
    let n = vec1.len();
    debug_assert_eq!(n, vec2.len(), "vec1 and vec2 must have the same length");

    let mut snv_out = Vec::new();
    let mut ins_out = Vec::new();
    let mut del_out = Vec::new();

    for i in 0..n {
        let v1 = vec1[i];
        let v2 = vec2[i];

        // Shared SNV: both positions are SNVs (pure -4 or compound mismatch-ins).
        if is_snv(v1) && is_snv(v2) {
            snv_out.push(overlap_start + i as i32);
        }

        // Shared insertion: both carry an insertion of the same inserted length.
        let l1 = insertion_len(v1);
        if l1 > 0 && l1 == insertion_len(v2) {
            ins_out.push(overlap_start + i as i32);
        }
    }

    // Shared deletions: find identical deletion spans (same start and end).
    // A deletion span is a contiguous stretch of HAP_DEL (-10) values.
    let vec1_del_bool = vec1.mapv(is_deletion);
    let vec2_del_bool = vec2.mapv(is_deletion);
    let vec1_del_spans = extract_true_stretches(&vec1_del_bool);
    let vec2_del_spans = extract_true_stretches(&vec2_del_bool);

    for &(s1, e1) in &vec1_del_spans {
        for &(s2, e2) in &vec2_del_spans {
            // Only match if the deletion spans are identical (same start, same end)
            if s1 == s2 && e1 == e2 {
                for i in s1..=e1 {
                    del_out.push(overlap_start + i as i32);
                }
            }
        }
    }

    (snv_out, ins_out, del_out)
}

/// Map absolute genomic positions to the read's encoded base via a
/// reference-offset → query-index mapping.
///
/// Ports from identify_misaligned_haps.py lines 174-193 (`map_positions_to_bases`)
///
/// Base encoding: A=0, T=1, C=2, G=3, N=4.
/// Returns `-1` for positions not covered by the read.
///
/// # Arguments
/// - `shared_pos_abs`: absolute genomic positions to look up
/// - `read_start`: the read's `reference_start`
/// - `ref_qseq_positions`: array mapping reference offset → query index (-1 if no mapping)
/// - `qseq_encoded`: encoded query sequence (i8)
pub fn map_positions_to_bases(
    shared_pos_abs: &[i32],
    read_start: i32,
    ref_qseq_positions: &[i32],
    qseq_encoded: &[i8],
) -> Vec<i8> {
    let n = shared_pos_abs.len();
    let mut out = vec![-1i8; n];

    for i in 0..n {
        let idx = shared_pos_abs[i] - read_start;
        if idx < 0 || idx as usize >= ref_qseq_positions.len() {
            continue;
        }
        let qidx = ref_qseq_positions[idx as usize];
        if qidx >= 0 && (qidx as usize) < qseq_encoded.len() {
            out[i] = qseq_encoded[qidx as usize];
        }
    }

    out
}

/// In-place update of a base-count tally matrix for one read at shared SNV positions.
///
/// Ports from identify_misaligned_haps.py lines 111-136 (`update_tally_for_read`)
///
/// # Arguments
/// - `shared_pos_abs`: absolute genomic positions (aligned with tally rows)
/// - `read_start`: the read's `reference_start`
/// - `ref_qseq_positions`: mapping from reference offset → query index
/// - `qseq_encoded`: encoded query sequence (A=0, T=1, C=2, G=3, N=4)
/// - `tally`: mutable `(num_pos, 5)` matrix; columns 0..3 = A/T/C/G, column 4 = N
pub fn update_tally_for_read(
    shared_pos_abs: &[i32],
    read_start: i32,
    ref_qseq_positions: &[i32],
    qseq_encoded: &[i8],
    tally: &mut [[i32; 5]],
) {
    let num_pos = shared_pos_abs.len();
    for i in 0..num_pos {
        let pos = shared_pos_abs[i];
        let idx = pos - read_start;
        if idx < 0 || idx as usize >= ref_qseq_positions.len() {
            continue;
        }
        let qidx = ref_qseq_positions[idx as usize];
        if qidx < 0 {
            continue;
        }
        let base = qseq_encoded[qidx as usize];
        // Only tally A(0), T(1), C(2), G(3); exclude N(4) and invalid
        if (0..4).contains(&base) {
            tally[i][base as usize] += 1;
        }
    }
}

/// Verify shared SNV positions by checking that the consensus ALT base
/// (argmax of tally ATCG counts) matches the homologous read's base.
///
/// Ports from identify_misaligned_haps.py lines 208-242 (`verify_shared_snv_positions`)
///
/// # Arguments
/// - `tally`: `(num_pos, 5)` base-count matrix (A/T/C/G/N columns)
/// - `h_base`: homolog encoded bases at each shared position (-1 = not covered)
/// - `shared_snv_pos_abs`: absolute genomic positions of shared SNVs
///
/// Returns the subset of `shared_snv_pos_abs` where the consensus ALT matches
/// the homolog's base.
pub fn verify_shared_snv_positions(
    tally: &[[i32; 5]],
    h_base: &[i8],
    shared_snv_pos_abs: &[i32],
) -> Vec<i32> {
    let n = shared_snv_pos_abs.len();
    let mut out = Vec::new();

    for i in 0..n {
        // Sum A/T/C/G counts (columns 0..3) to get coverage
        let cov = tally[i][0] + tally[i][1] + tally[i][2] + tally[i][3];
        if cov <= 0 {
            continue;
        }

        // argmax over columns 0..3
        let mut best_base: i8 = 0;
        let mut best_count = tally[i][0];
        for b in 1..4 {
            if tally[i][b] > best_count {
                best_count = tally[i][b];
                best_base = b as i8;
            }
        }

        // Check homolog base matches consensus
        if h_base[i] >= 0 && best_base == h_base[i] {
            out.push(shared_snv_pos_abs[i]);
        }
    }

    out
}

// ─── stat_refseq_similarity ─────────────────────────────────────────────────

/// Regex for parsing homo_refseq_qname to extract origin chrom:start-end.
/// Matches patterns like "chr1:12345-67890" possibly followed by ":RG<digits>".
/// Python: `re.search(r"([a-zA-Z0-9]+):(\d+)-(\d+)[:RG0-9]*", homo_refseq_qname)`
fn parse_origin_region(qname: &str) -> Option<(String, i64, i64)> {
    // Simple manual parse to avoid regex dependency.
    // Find the pattern: <alphanum>:<digits>-<digits>
    // We look for the FIRST occurrence.
    let bytes = qname.as_bytes();
    let len = bytes.len();

    let mut i = 0;
    while i < len {
        // Look for ':' preceded by alphanumeric characters
        if bytes[i] == b':' && i > 0 {
            // Scan backwards to find start of chrom (alphanumeric run)
            let mut chrom_start = i;
            while chrom_start > 0 && bytes[chrom_start - 1].is_ascii_alphanumeric() {
                chrom_start -= 1;
            }
            if chrom_start == i {
                i += 1;
                continue;
            }
            // Scan forward for digits after ':'
            let mut j = i + 1;
            let digit_start = j;
            while j < len && bytes[j].is_ascii_digit() {
                j += 1;
            }
            if j == digit_start || j >= len || bytes[j] != b'-' {
                i += 1;
                continue;
            }
            let start_str = &qname[digit_start..j];
            // Skip '-'
            j += 1;
            let digit_start2 = j;
            while j < len && bytes[j].is_ascii_digit() {
                j += 1;
            }
            if j == digit_start2 {
                i += 1;
                continue;
            }
            let end_str = &qname[digit_start2..j];

            if let (Ok(start_val), Ok(end_val)) = (start_str.parse::<i64>(), end_str.parse::<i64>()) {
                let chrom = &qname[chrom_start..i];
                return Some((chrom.to_string(), start_val, end_val));
            }
        }
        i += 1;
    }
    None
}

/// Compare a consensus haplotype with overlapping intrinsic (homologous) reference
/// sequences to measure similarity via shared paralogous sequence variants (PSVs).
///
/// Faithfully ports Python's `stat_refseq_similarity`
/// (identify_misaligned_haps.py:888-1018).
///
/// # Arguments
/// - `intrin_lapper` — Lapper result built from the intrinsic BAM (homologous ref seqs)
/// - `chrom` — chromosome name
/// - `span` — `(start, end)` of the continuous region (0-based, half-open)
/// - `hid` — haplotype ID
/// - `consensus_sequence` — assembled consensus for this haplotype region
/// - `reads` — member reads of this haplotype in this region (for SNV verification tally)
/// - `total_genomic_haps` — cache: homo_refseq_id -> hap_vector (mutable, populated on miss)
/// - `qseq_cache` — cache: read_id -> ReadQseqData (mutable, shared across calls)
/// - `varcounts_among_refseqs` — output accumulator (hid -> homo_qname -> Vec<RegionVarStats>)
///
/// # Returns
/// `Ok(())` after mutating `varcounts_among_refseqs`.
///
/// # Errors
/// Propagates [`crate::pairwise_read_inspection::CigarError`] if a homologous refseq
/// or member read has a non-`=`/`X` CIGAR (DivA — fail fast).
#[allow(clippy::too_many_arguments)]
pub fn stat_refseq_similarity(
    intrin_lapper: &BamLapperResult,
    chrom: &str,
    span: (i64, i64),
    hid: i32,
    consensus_sequence: &Array1<i16>,
    reads: &[&Record],
    total_genomic_haps: &mut HashMap<String, Array1<i16>>,
    qseq_cache: &mut HashMap<String, ReadQseqData>,
    varcounts_among_refseqs: &mut VarcountsAmongRefseqs,
) -> Result<(), Box<dyn std::error::Error>> {
    // Query intrinsic Lapper for overlapping homologous reference sequences
    let homo_refseqs = query_overlapping_reads(
        &intrin_lapper.lapper_dict,
        &intrin_lapper.read_dict,
        chrom,
        span.0 as u32,
        span.1 as u32,
    );

    for homo_refseq in &homo_refseqs {
        let homo_refseq_start = homo_refseq.pos();
        let homo_refseq_end = homo_refseq.cigar().end_pos();
        let homo_refseq_qname = String::from_utf8_lossy(homo_refseq.qname()).to_string();

        let overlap_start = homo_refseq_start.max(span.0);
        let overlap_end = homo_refseq_end.min(span.1);
        let overlap_span_size = overlap_end - overlap_start;

        let hregion_str = format!("{chrom}:{homo_refseq_start}-{homo_refseq_end}");
        // Some reference sequences can be mapped to multiple places;
        // to have a unique ID, we append the current region coordinate string.
        let homo_refseq_id = format!("{}:{}", read_id(homo_refseq), hregion_str);

        // ── Get or compute the hap vector for this homologous refseq ──
        let homo_refseq_hap_vector = if let Some(cached) = total_genomic_haps.get(&homo_refseq_id) {
            cached.clone()
        } else {
            // Python calls extract_read_qseqs here to get query_sequence_encoded,
            // then passes it to get_hapvector_from_cigar. The N-base check in
            // get_hapvector_from_cigar is dead code (compares int array slice to
            // string "N"), so extract_hap_vector (CIGAR-only) is equivalent.
            // Cache the qseq data for later use in map_positions_to_bases.
            let homo_rid = read_id(homo_refseq);
            if let Entry::Vacant(slot) = qseq_cache.entry(homo_rid) {
                slot.insert(extract_read_qseqs(homo_refseq)?);
            }
            let hap_vec = extract_hap_vector(homo_refseq)?;
            total_genomic_haps.insert(homo_refseq_id.clone(), hap_vec.clone());
            hap_vec
        };

        // ── Slice the hap vector to the overlap region ──
        let hap_offset_start = (overlap_start - homo_refseq_start) as usize;
        let hap_offset_end = (overlap_end - homo_refseq_start) as usize;
        if hap_offset_end > homo_refseq_hap_vector.len() || hap_offset_start >= hap_offset_end {
            // Python: except IndexError: continue
            continue;
        }
        let interval_genomic_hap = homo_refseq_hap_vector.slice(ndarray::s![hap_offset_start..hap_offset_end]);

        // ── Slice the consensus sequence to the overlap region ──
        let con_offset_start = (overlap_start - span.0) as usize;
        let con_offset_end = (overlap_end - span.0) as usize;
        if con_offset_end > consensus_sequence.len() || con_offset_start >= con_offset_end {
            continue;
        }
        let interval_con_seq = consensus_sequence.slice(ndarray::s![con_offset_start..con_offset_end]);

        // ── Compute ref_genome_similarity ──
        let interval_con_seq_owned = interval_con_seq.to_owned();
        let interval_genomic_hap_owned = interval_genomic_hap.to_owned();
        let (varcount, alt_snv_count, alt_indel_count) =
            ref_genome_similarity(&interval_con_seq_owned, &interval_genomic_hap_owned);

        // ── If zero alt variants, check if this homologous seq is aligned to its origin ──
        if (alt_snv_count + alt_indel_count) == 0 {
            if let Some((origin_chrom, origin_start, origin_end)) = parse_origin_region(&homo_refseq_qname) {
                if origin_chrom == chrom
                    && origin_start <= overlap_start
                    && origin_end >= overlap_end
                {
                    warn!(
                        "The homologous genomic sequence {homo_refseq_qname} overlapping interval {overlap_start}-{overlap_end} \
                         is aligned to its origin. So ignore this homologous sequence."
                    );
                    continue;
                }
            }
        }

        // ── Compute shared variant positions ──
        // Python returns 4 values: (snv, ins, del, del_count).
        // del_count == del_out.len() in Rust since we don't pre-allocate+truncate.
        let (shared_snv_pos_abs, shared_ins_pos_abs, shared_del_pos_abs) =
            numba_shared_variant_positions(
                &interval_con_seq_owned,
                &interval_genomic_hap_owned,
                overlap_start as i32,
            );
        let shared_psv_del = shared_del_pos_abs.len() as i32;

        // ── Verify shared SNV positions via read-level tally ──
        let mut verified_shared_snv_pos_abs: Vec<i32> = Vec::new();

        if !shared_snv_pos_abs.is_empty() {
            let num_pos = shared_snv_pos_abs.len();
            let mut tally = vec![[0i32; 5]; num_pos];

            // Tally consensus ALT codes at shared positions from member reads
            for r in reads {
                let r_id = read_id(r);
                let r_qseq = if let Some(cached) = qseq_cache.get(&r_id) {
                    cached.clone()
                } else {
                    let data = extract_read_qseqs(r)?;
                    qseq_cache.insert(r_id, data.clone());
                    data
                };

                update_tally_for_read(
                    &shared_snv_pos_abs,
                    r.pos() as i32,
                    &r_qseq.ref_to_query,
                    &r_qseq.qseq_encoded,
                    &mut tally,
                );
            }

            // ALT codes from the homologous read (use cache)
            let homo_rid = read_id(homo_refseq);
            let h_qseq = if let Some(cached) = qseq_cache.get(&homo_rid) {
                cached.clone()
            } else {
                let data = extract_read_qseqs(homo_refseq)?;
                qseq_cache.insert(homo_rid, data.clone());
                data
            };
            let h_base = map_positions_to_bases(
                &shared_snv_pos_abs,
                homo_refseq_start as i32,
                &h_qseq.ref_to_query,
                &h_qseq.qseq_encoded,
            );

            // Verify: consensus ALT must match homolog ALT
            verified_shared_snv_pos_abs =
                verify_shared_snv_positions(&tally, &h_base, &shared_snv_pos_abs);
        }

        // ── Count shared PSV insertions ──
        let mut shared_psv_ins: i32 = 0;
        for &ins_pos in &shared_ins_pos_abs {
            let idx = (ins_pos - overlap_start as i32) as usize;
            if idx < interval_genomic_hap_owned.len() {
                let ins_encode_event = interval_genomic_hap_owned[idx];
                if ins_encode_event > 0 {
                    info!(
                        "The homologous genomic sequence aligned at interval {}:{}-{} shared \
                         an insertion at position {}. The encoded event is {}, the alignment \
                         status on consensus sequence is {}",
                        chrom, overlap_start, overlap_end, ins_pos,
                        ins_encode_event,
                        interval_con_seq_owned[idx]
                    );
                    shared_psv_ins += 1;
                } else if ins_encode_event < 0 {
                    warn!(
                        "The homologous genomic sequence aligned at interval {}:{}-{} shared \
                         an insertion at position {}. The encoded event is {}, the alignment \
                         status on consensus sequence is {}",
                        chrom, overlap_start, overlap_end, ins_pos,
                        ins_encode_event,
                        interval_con_seq_owned[idx]
                    );
                }
            }
        }

        let shared_psv = verified_shared_snv_pos_abs.len() as i32 + shared_psv_ins + shared_psv_del;
        let verified_shared_indel_pos_abs = merge_unique_sorted(
            &shared_ins_pos_abs,
            &shared_del_pos_abs,
        );

        // ── Accumulate into varcounts_among_refseqs ──
        let stats = RegionVarStats {
            varcount,
            alt_snv_count,
            alt_indel_count,
            shared_psv,
            verified_shared_snv_pos_abs,
            verified_shared_indel_pos_abs,
            overlap_span_size,
        };

        debug!(
            "[stat_refseq_similarity] hid={hid} homo_qname={homo_refseq_qname} varcount={varcount} alt_snv={alt_snv_count} alt_indel={alt_indel_count} shared_psv={shared_psv}"
        );

        varcounts_among_refseqs
            .entry(hid)
            .or_default()
            .entry(homo_refseq_qname.clone())
            .or_default()
            .push(stats);
    }

    Ok(())
}

// ─── cal_similarity_score ───────────────────────────────────────────────────

/// Per-haplotype result from [`cal_similarity_score`].
#[derive(Clone, Debug)]
pub struct HapSimilarityScore {
    /// Maximum mixed_psv_metric across all reference sequences.
    pub max_sim_score: f64,
    /// total_shared_psv at the best-scoring reference sequence.
    pub max_psv_count: i32,
    /// Sorted unique PSV genomic positions at the best-scoring reference sequence.
    pub max_psv_positions: Vec<i32>,
}

/// Compute per-haplotype maximum similarity scores from the accumulated
/// `varcounts_among_refseqs`.
///
/// Faithfully ports Python `cal_similarity_score` (identify_misaligned_haps.py:1022-1087).
///
/// # Formula
///
/// For each (hid, homo_refseq_qname) combination, aggregate across all region
/// pairs then compute:
///
/// ```text
/// psv_var_ratio     = total_shared_psv / total_varcount       (if total_varcount > 0)
///                   = min(1, total_shared_psv)                 (otherwise)
/// psv_sharing_ratio = total_shared_psv / (alt_snv + alt_indel) (if denom > 0, else 0)
/// non_psv_density_100bp = hid_max_local_density * 100
/// mixed_psv_metric  = sqrt(psv_var_ratio) * total_shared_psv * sqrt(psv_sharing_ratio)
///                     + sqrt(non_psv_density_100bp)
///                     - (4 - psv_var_ratio)
/// ```
///
/// The refseq yielding the **maximum** `mixed_psv_metric` is kept per haplotype.
pub fn cal_similarity_score(
    varcounts_among_refseqs: &VarcountsAmongRefseqs,
    hid_var_count: &HashMap<i32, i32>,
    hid_max_local_density: &HashMap<i32, f64>,
) -> HashMap<i32, HapSimilarityScore> {
    let mut results: HashMap<i32, HapSimilarityScore> = HashMap::new();

    for (&hid, gdict) in varcounts_among_refseqs.iter() {
        let mut max_psv: f64 = -1.0;
        let mut max_psv_c: i32 = 0;
        let mut max_psv_pos: Vec<i32> = Vec::new();

        for (homo_refseq_qname, pairs) in gdict.iter() {
            // Aggregate across all region pairs for this refseq
            let total_shared_psv: i32 = pairs.iter().map(|t| t.shared_psv).sum();
            let alt_snv_count: i32 = pairs.iter().map(|t| t.alt_snv_count).sum();
            let alt_indel_count: i32 = pairs.iter().map(|t| t.alt_indel_count).sum();
            let total_varcount: i32 = *hid_var_count.get(&hid).unwrap_or(&0);

            // Concat all verified shared SNV and indel positions
            let mut all_snv_pos: Vec<i32> = pairs
                .iter()
                .flat_map(|t| t.verified_shared_snv_pos_abs.iter().copied())
                .collect();
            let mut all_indel_pos: Vec<i32> = pairs
                .iter()
                .flat_map(|t| t.verified_shared_indel_pos_abs.iter().copied())
                .collect();

            // Combine, unique, sort → psv_pos_abs  (matches np.unique + np.sort)
            all_snv_pos.append(&mut all_indel_pos);
            all_snv_pos.sort_unstable();
            all_snv_pos.dedup();
            let psv_pos_abs = all_snv_pos; // now sorted & unique

            // Compute psv_var_ratio
            let total_shared_psv_f = total_shared_psv as f64;
            let total_varcount_f = total_varcount as f64;
            let psv_var_ratio: f64 = if total_varcount > 0 {
                total_shared_psv_f / total_varcount_f
            } else {
                // Python: min(1, total_shared_psv)
                (total_shared_psv as f64).min(1.0)
            };

            // Compute psv_sharing_ratio
            let total_psv_count = alt_snv_count + alt_indel_count;
            let psv_sharing_ratio: f64 = if total_psv_count > 0 {
                total_shared_psv_f / (total_psv_count as f64)
            } else {
                0.0
            };

            // Density term
            let max_density = *hid_max_local_density.get(&hid).unwrap_or(&0.0);
            let non_psv_density_100bp = max_density * 100.0;

            // mixed_psv_metric (matches Python exactly)
            let mixed_psv_metric = psv_var_ratio.sqrt()
                * total_shared_psv_f
                * psv_sharing_ratio.sqrt()
                + non_psv_density_100bp.sqrt()
                - (4.0 - psv_var_ratio);

            debug!(
                "[cal_similarity_score] hid={hid} homo_refseq={homo_refseq_qname} psv_var_ratio={psv_var_ratio:.4} \
                 total_shared_psv={total_shared_psv} alt_snv={alt_snv_count} alt_indel={alt_indel_count} psv_sharing_ratio={psv_sharing_ratio:.4} \
                 non_psv_density_100bp={non_psv_density_100bp:.4} mixed_psv_metric={mixed_psv_metric:.6}",
            );

            if mixed_psv_metric > max_psv {
                max_psv_c = total_shared_psv;
                max_psv = mixed_psv_metric;
                max_psv_pos = psv_pos_abs;
            }
        }

        info!(
            "[cal_similarity_score] hid={hid} max_sim_score={max_psv:.6} max_psv_count={max_psv_c} \
             max_psv_positions={max_psv_pos:?}",
        );

        results.insert(
            hid,
            HapSimilarityScore {
                max_sim_score: max_psv,
                max_psv_count: max_psv_c,
                max_psv_positions: max_psv_pos,
            },
        );
    }

    results
}


// ─── ILP Post-Processing Helpers ──────────────────────────────────────────────

/// Rank unique values in a float array by sorted order (1-based).
///
/// Ports from identify_misaligned_haps.py lines 273-285 (numba)
///
/// Algorithm:
/// 1. Extract unique values and sort them
/// 2. For each element, find its position in sorted unique values
/// 3. Return 1-based rank
///
/// Example: `[3.0, 1.0, 3.0, 2.0]` → `[3, 1, 3, 2]`
pub fn rank_unique_values(arr: &[f32]) -> Vec<i32> {
    // Extract unique values and sort them. `total_cmp` is a total order (NaN sorts
    // last, deterministically) so ranking never panics on a non-finite value.
    let mut unique_values: Vec<f32> = arr.to_vec();
    unique_values.sort_by(|a, b| a.total_cmp(b));
    unique_values.dedup();

    // Rank each element by its position in sorted unique values (1-based)
    let mut ranks = vec![0i32; arr.len()];
    for (i, &val) in arr.iter().enumerate() {
        for (j, &uval) in unique_values.iter().enumerate() {
            if val == uval {
                ranks[i] = (j + 1) as i32;
                break;
            }
        }
    }

    ranks
}

/// Calculate per-region coefficient for haplotype ranking.
///
/// Ports from identify_misaligned_haps.py lines 426-433 (numba)
///
/// Input columns (matching the DataFrame column order):
///   [0] start, [1] end, [2] total_depth, [3] hap_id,
///   [4] hap_depth, [5] var_count, [6] indel_count,
///   [7] psv_count, [8] varc_rank, [9] hap_max_sim_scores
///
/// Formula:
///   `coefficient = psv_metric * sqrt(span/100) * sqrt(1 - hap_depth/total_depth)`
///
/// where `psv_metric` is column 9 (`hap_max_sim_scores`),
/// `span` is `(end - start) / 100`, and `depth_frac` is `1 - hap_depth/total_depth`.
pub fn calculate_coefficient(rows: &[&[f32; 10]]) -> Vec<f32> {
    let mut res = Vec::with_capacity(rows.len());
    for row in rows {
        let start = row[0];
        let end = row[1];
        let total_depth = row[2];
        let hap_depth = row[4];
        let psv_metric = row[9];

        let span = (end - start) / 100.0;
        let depth_frac = 1.0 - (hap_depth / total_depth);

        res.push(psv_metric * span.sqrt() * depth_frac.sqrt());
    }
    res
}


// ─── Batch vector collection per region ───────────────────────────────────────

/// Padded per-region vectors produced by [`record_hap_err_vectors_per_region`]:
/// `(read_spans [N×2], hap_vectors [N×max_len], err_vectors [N×max_len])`.
pub type RegionVectors = (Array2<i32>, Array2<i16>, Array2<f32>);

/// Collect haplotype and error vectors for all reads in a region, with caching.
///
/// Ports Python's `record_hap_err_vectors_per_region`
/// (identify_misaligned_haps.py:742-783).
///
/// For each read, extracts (or retrieves from cache) the haplotype and error
/// vectors, then packs them into 2D padded arrays suitable for `assemble_consensus`.
///
/// # Arguments
/// * `records`   - Slice of references: borrows the outer slice and each Record.
///   Avoids copying any Record objects. The Records are owned by
///   `BamLapperResult.read_dict`.
/// * `hap_cache` - Mutable borrow of the cache. Owned `String` keys are required
///   because `HashMap` entries must outlive any individual function call. `Array1`
///   values are owned so the cache can persist across region calls.
/// * `err_cache` - Same ownership semantics as `hap_cache`.
///
/// # Returns
/// * `read_spans`  - `Array2<i32>` shape `(n_reads, 2)` with `[ref_start, ref_end]` per read
/// * `hap_vectors` - `Array2<i16>` shape `(n_reads, max_len)` padded with [`HAP_PAD`] (-20)
/// * `err_vectors` - `Array2<f32>` shape `(n_reads, max_len)` padded with -10.0
///
/// # Errors
/// Propagates [`crate::pairwise_read_inspection::CigarError`] if any read's CIGAR is
/// not in `=`/`X` mode (DivA — fail fast instead of panicking).
pub fn record_hap_err_vectors_per_region(
    records: &[&Record],
    hap_cache: &mut HashMap<String, Array1<i16>>,
    err_cache: &mut HashMap<String, Array1<f32>>,
) -> Result<RegionVectors, Box<dyn std::error::Error>> {
    let n = records.len();

    if n == 0 {
        info!("[record_hap_err_vectors_per_region] no reads to process");
        return Ok((
            Array2::<i32>::zeros((0, 2)),
            Array2::<i16>::zeros((0, 0)),
            Array2::<f32>::zeros((0, 0)),
        ));
    }

    let cache_size_before = hap_cache.len();

    // First pass: ensure each read's hap/err vectors are in the caches, and
    // collect spans + read_ids + the max hap length for padding. We deliberately
    // do NOT clone the cached vectors here — the second pass copies straight from
    // the cache into the output rows, so each Array1 is touched only by reference.
    let mut spans_vec: Vec<[i32; 2]> = Vec::with_capacity(n);
    let mut rids: Vec<String> = Vec::with_capacity(n);
    let mut max_len = 0usize;

    for &record in records {
        let ref_start = record.pos() as i32;
        let ref_end = record.cigar().end_pos() as i32;
        spans_vec.push([ref_start, ref_end]);

        let rid = read_id(record);
        let hap_hit = hap_cache.contains_key(&rid);

        // Populate the caches on miss, propagating CigarError via `?` (DivA: the
        // M-op rejection early-stops here rather than panicking). extract_hap_vector
        // runs first, so a malformed read never reaches extract_error_vector. The
        // Entry API avoids the double lookup that `contains_key` + `insert` incurs.
        if let Entry::Vacant(e) = hap_cache.entry(rid.clone()) {
            e.insert(extract_hap_vector(record)?);
        }
        let hap_len = hap_cache[&rid].len();

        if let Entry::Vacant(e) = err_cache.entry(rid.clone()) {
            e.insert(extract_error_vector(record)?);
        }
        let err_len = err_cache[&rid].len();

        debug!(
            "[record_hap_err_vectors_per_region] read {} cache_hit={} hap_len={} err_len={}",
            String::from_utf8_lossy(record.qname()), hap_hit, hap_len, err_len
        );

        max_len = max_len.max(hap_len);
        rids.push(rid);
    }

    let cache_hits = n - (hap_cache.len() - cache_size_before);
    info!(
        "[record_hap_err_vectors_per_region] n_reads={} max_len={} cache_hits={} cache_misses={} total_cache_size={}",
        n, max_len, cache_hits, n - cache_hits, hap_cache.len()
    );

    // Build read_spans Array2<i32> (n × 2)
    let mut read_spans = Array2::<i32>::zeros((n, 2));
    for (i, span) in spans_vec.iter().enumerate() {
        read_spans[[i, 0]] = span[0];
        read_spans[[i, 1]] = span[1];
    }

    // Second pass: copy each cached vector directly into its padded row.
    // hap and err vectors share the same length (both span the read's reference
    // footprint), so max_len computed from hap lengths bounds both.
    // hap padding is HAP_PAD (-20): golden deletion is -10, so the old -10 padding
    // would collide and drop deletions. err padding stays -10.0 (it does not collide —
    // assemble_consensus truncates qual by the seq-derived count, never testing it).
    let mut hap_vectors = Array2::<i16>::from_elem((n, max_len), HAP_PAD);
    let mut err_vectors = Array2::<f32>::from_elem((n, max_len), -10.0f32);
    for (i, rid) in rids.iter().enumerate() {
        let hap = &hap_cache[rid];
        let hlen = hap.len();
        hap_vectors.row_mut(i).as_slice_mut().unwrap()[..hlen]
            .copy_from_slice(hap.as_slice().unwrap());

        let err = &err_cache[rid];
        let elen = err.len();
        err_vectors.row_mut(i).as_slice_mut().unwrap()[..elen]
            .copy_from_slice(err.as_slice().unwrap());
    }

    Ok((read_spans, hap_vectors, err_vectors))
}

// ─── Consensus assembly ───────────────────────────────────────────────────────

/// Assemble consensus sequence from multiple reads' haplotype and error vectors.
///
/// This is a faithful port of Python's `assemble_consensus` (identify_misaligned_haps.py:288-341).
///
/// All three array arguments are borrowed read-only (`&Array2`); no copy is made.
/// The caller retains ownership and can reuse them after this call.
///
/// # Arguments
/// * `seq_arrays`  - `Array2<i16>` shape `(n_reads, max_len)`, padded with [`HAP_PAD`] (-20).
///   Any value other than `HAP_PAD` is valid (non-NA) — this keeps the golden deletion
///   `-10`, which the old `>= -8` filter would have dropped as padding.
/// * `qual_arrays` - `Array2<f32>` shape `(n_reads, max_len)` (err padding -10.0).
/// * `read_spans`  - `Array2<i32>` shape `(n_reads, 2)` with `[ref_start, ref_end]` per read.
///
/// # Returns
/// Consensus sequence as `Array1<i16>` spanning from `min(ref_start)` to `max(ref_end)`.
/// Default value at each position is `1` (reference match) until overwritten by a read with
/// error probability <= 0.2 that is better (lower) than the incumbent.
pub fn assemble_consensus(
    seq_arrays: &Array2<i16>,
    qual_arrays: &Array2<f32>,
    read_spans: &Array2<i32>,
) -> Array1<i16> {
    // ── Validate inputs ──────────────────────────────────────────────────
    let n_reads = seq_arrays.nrows();
    assert_eq!(
        n_reads,
        qual_arrays.nrows(),
        "seq_arrays and qual_arrays must have the same number of rows"
    );
    assert_eq!(
        n_reads,
        read_spans.nrows(),
        "seq_arrays and read_spans must have the same number of rows"
    );

    if n_reads == 0 {
        return Array1::<i16>::zeros(0);
    }

    // ── Python line 292-293 ──────────────────────────────────────────────
    // start_pos = read_spans[:, 0].min()
    // end_pos   = read_spans[:, 1].max()
    let start_pos = *read_spans.column(0).iter().min().unwrap();
    let end_pos = *read_spans.column(1).iter().max().unwrap();

    let length = (end_pos - start_pos) as usize;

    // ── Python line 298-299 ──────────────────────────────────────────────
    // consensus_seq  = np.ones(end_pos - start_pos, dtype=np.int16)
    // consensus_qual = np.full(end_pos - start_pos, 0.2, dtype=np.float32)
    let mut consensus_seq: Vec<i16> = vec![1i16; length];
    let mut consensus_qual: Vec<f32> = vec![0.2f32; length];

    // ── Python line 301: for i in prange(len(seq_arrays)) ────────────────
    for i in 0..n_reads {
        let seq_row = seq_arrays.row(i);
        let qual_row = qual_arrays.row(i);
        let start = read_spans[[i, 0]];

        // Count real (non-padding) values; padding (HAP_PAD = -20) is always at the
        // tail, so this equals the length of the real prefix. Golden change: was
        // `>= -8`, which would now wrongly drop the deletion signal (-10).
        let non_na_values: usize = seq_row.iter().filter(|&&v| is_real(v)).count();

        // ── Python lines 317-325: strip padding ─────────────────────────
        let nona_seq: Vec<i32> = seq_row.iter().take(non_na_values).map(|&v| v as i32).collect();
        let nona_qual: Vec<f32> = qual_row.iter().take(non_na_values).copied().collect();

        // ── Python lines 328-329 ────────────────────────────────────────
        let rel_start = (start - start_pos) as usize;

        // ── Python lines 334-339 ────────────────────────────────────────
        for j in 0..non_na_values {
            let pos = rel_start + j;
            if pos < length
                && nona_qual[j] <= consensus_qual[pos] && nona_qual[j] <= 0.2 {
                    consensus_seq[pos] = nona_seq[j] as i16;
                    consensus_qual[pos] = nona_qual[j];
                }
        }
    }

    // ── Python line 341: return consensus_seq ────────────────────────────
    Array1::from_vec(consensus_seq)
}

/// Group reads into non-overlapping continuous genomic regions.
///
/// Faithfully ports Python's `extract_continuous_regions_dict`
/// (identify_misaligned_haps.py:477-513).
///
/// # Algorithm
/// 1. Sort reads by `reference_start` (ascending).
/// 2. Walk through the sorted reads.  If a read overlaps the current region
///    (`read_start <= current_end && read_end >= current_start`), extend the
///    region.  Otherwise close the current region and start a new one.
/// 3. Return a `Vec` of `(span, read_indices)` where `span = (start, end)`
///    is 0-based half-open and `read_indices` indexes into the **input** slice
///    (i.e. before sorting).
///
/// # Arguments
/// - `reads` — slice of BAM records (any order; will be sorted internally)
///
/// # Returns
/// `Vec<((i64, i64), Vec<usize>)>` — each entry is a continuous region span
/// paired with the indices into the original `reads` slice that belong to it.
///
/// The ordering of regions is ascending by start coordinate.
pub fn extract_continuous_regions_dict(reads: &[&Record]) -> Vec<((i64, i64), Vec<usize>)> {
    if reads.is_empty() {
        return Vec::new();
    }

    // Build (original_index, ref_start) pairs and sort by ref_start
    let mut order: Vec<usize> = (0..reads.len()).collect();
    order.sort_by_key(|&i| reads[i].pos());

    let mut regions: Vec<((i64, i64), Vec<usize>)> = Vec::new();
    let mut current_start: i64 = 0;
    let mut current_end: i64 = 0;
    let mut current_indices: Vec<usize> = Vec::new();

    for &orig_idx in &order {
        let read = reads[orig_idx];
        let read_start = read.pos(); // 0-indexed
        let read_end = read.cigar().end_pos(); // 0-indexed, one past the last aligned base

        if current_indices.is_empty() {
            // First read — initialise the region
            current_start = read_start;
            current_end = read_end;
            current_indices.push(orig_idx);
        } else if read_start <= current_end && read_end >= current_start {
            // Read overlaps the current region — extend it
            current_start = current_start.min(read_start);
            current_end = current_end.max(read_end);
            current_indices.push(orig_idx);
        } else {
            // No overlap — close the current region and start a new one
            debug!(
                "[extract_continuous_regions_dict] region ({}, {}) with {} reads",
                current_start, current_end, current_indices.len()
            );
            regions.push(((current_start, current_end), current_indices));
            current_start = read_start;
            current_end = read_end;
            current_indices = vec![orig_idx];
        }
    }

    // Flush the last region
    if !current_indices.is_empty() {
        debug!(
            "[extract_continuous_regions_dict] region ({}, {}) with {} reads",
            current_start, current_end, current_indices.len()
        );
        regions.push(((current_start, current_end), current_indices));
    }

    regions
}

// ─── group_by_dict_optimized ────────────────────────────────────────────────

/// Result type for one haplotype group.
///
/// Corresponds to `[vertex_indices, qnames, read_pair_lists]` in Python.
pub struct HapGroup<'a> {
    pub vertex_indices: Vec<i32>,
    pub qnames: Vec<String>,
    pub read_pair_lists: Vec<Vec<&'a Record>>,
}

/// Group overlapping read-pairs by their haplotype label.
///
/// Faithfully ports Python's `group_by_dict_optimized`
/// (identify_misaligned_haps.py:516-527).
///
/// # Arguments
/// - `vprop`    — vertex-index → haplotype-id mapping (Python: `qname_hap_info`)
/// - `vertices` — `{(vertex_idx, qname) → vec_of_reads}` from the overlap query
///
/// # Returns
/// `HashMap<i32, HapGroup>` — keyed by haplotype id
///
/// # Panics
/// Panics if a vertex_idx in `vertices` is not present in `vprop`
/// (matching Python's implicit KeyError behaviour).
pub fn group_by_dict_optimized<'a>(
    vprop: &HashMap<i32, i32>,
    vertices: &HashMap<(i32, String), Vec<&'a Record>>,
) -> HashMap<i32, HapGroup<'a>> {
    let mut grouped: HashMap<i32, HapGroup> = HashMap::new();

    for ((v_idx, qname), reads) in vertices.iter() {
        let label = vprop[v_idx]; // panics if missing — matches Python KeyError
        let group = grouped.entry(label).or_insert_with(|| HapGroup {
            vertex_indices: Vec::new(),
            qnames: Vec::new(),
            read_pair_lists: Vec::new(),
        });
        group.vertex_indices.push(*v_idx);
        group.qnames.push(qname.clone());
        group.read_pair_lists.push(reads.clone());
    }

    debug!(
        "[group_by_dict_optimized] {} vertices → {} groups",
        vertices.len(),
        grouped.len()
    );

    grouped
}

// ─── record_haplotype_rank ──────────────────────────────────────────────────

/// Per-haplotype data needed by `record_haplotype_rank`.
pub struct HaplotypeClusterInfo<'a> {
    /// Consensus sequence for this haplotype within the overlapping span
    pub consensus: Array1<i16>,
    /// Reads belonging to this haplotype cluster
    pub reads: Vec<&'a Record>,
    /// (start, end) span of the overlapping region (0-based, end exclusive)
    pub span: (i64, i64),
    /// Qnames in this haplotype cluster
    pub qnames: Vec<String>,
}

/// Build a 2-D array ranking haplotypes within a region.
///
/// Faithfully ports Python's `record_haplotype_rank`
/// (identify_misaligned_haps.py:530-571).
///
/// # Columns (8 per row):
/// `[start, end, total_depth, hap_id, hap_depth, var_count, indel_count, psv_count]`
///
/// Uses ndarray `Array2<i32>` matching the numpy column_stack in Python.
///
/// # Arguments
/// - `haplotype_dict` — hap_id → `HaplotypeClusterInfo`
/// - `mean_read_length` — not directly used in the calculation (Python default 150)
/// - `hap_max_psv_pos` — hap_id → sorted PSV positions (from `cal_similarity_score`)
pub fn record_haplotype_rank(
    haplotype_dict: &HashMap<i32, HaplotypeClusterInfo>,
    _mean_read_length: i32,
    hap_max_psv_pos: &HashMap<i32, Vec<i32>>,
) -> Array2<i32> {
    let n = haplotype_dict.len();
    if n == 0 {
        return Array2::<i32>::zeros((0, 8));
    }

    let mut starts = Array1::<i32>::zeros(n);
    let mut ends = Array1::<i32>::zeros(n);
    let mut hap_depths = Array1::<i32>::zeros(n);
    let mut hap_ids = Array1::<i32>::zeros(n);
    let mut var_counts = Array1::<i32>::zeros(n);
    let mut indel_counts = Array1::<i32>::zeros(n);
    let mut psv_counts = Array1::<i32>::zeros(n);

    let mut total_depth: i32 = 0;
    // Iterate in sorted hid order for determinism (Python iterates dict insertion order,
    // but the callers reconstruct from DataFrame so ordering doesn't matter for correctness).
    let mut hids: Vec<i32> = haplotype_dict.keys().copied().collect();
    hids.sort_unstable();

    for (i, &hid) in hids.iter().enumerate() {
        let info = &haplotype_dict[&hid];
        let (region_start, _region_end) = info.span;
        let region_start_i32 = region_start as i32;
        let region_len = info.consensus.len();
        let region_end_excl = region_start + region_len as i64;

        starts[i] = region_start_i32;
        ends[i] = info.span.1 as i32;

        // Compute depth: total ref-aligned bases in [region_start, region_end_excl)
        let mut cov_bases: i64 = 0;
        for &r in &info.reads {
            let s = (r.pos()).max(region_start);
            let e = (r.cigar().end_pos()).min(region_end_excl);
            if e > s {
                cov_bases += e - s;
            }
        }
        let depth = if region_len > 0 {
            (cov_bases / region_len as i64) as i32
        } else {
            0
        };
        hap_depths[i] = depth;
        hap_ids[i] = hid;
        var_counts[i] = count_var(&info.consensus);
        indel_counts[i] = count_continuous_indel_blocks(&info.consensus);

        // PSV positions within [region_start, region_end_excl)
        let empty_psv = Vec::new();
        let psv_pos = hap_max_psv_pos.get(&hid).unwrap_or(&empty_psv);
        let psv_in_region = psv_pos
            .iter()
            .filter(|&&p| p >= region_start_i32 && p < region_end_excl as i32)
            .count() as i32;
        psv_counts[i] = psv_in_region;
        total_depth += depth;
    }

    // Build the 2D result: np.column_stack(...)
    let total_depth_col = Array1::from_elem(n, total_depth);
    let mut result = Array2::<i32>::zeros((n, 8));
    result.column_mut(0).assign(&starts);
    result.column_mut(1).assign(&ends);
    result.column_mut(2).assign(&total_depth_col);
    result.column_mut(3).assign(&hap_ids);
    result.column_mut(4).assign(&hap_depths);
    result.column_mut(5).assign(&var_counts);
    result.column_mut(6).assign(&indel_counts);
    result.column_mut(7).assign(&psv_counts);

    debug!(
        "[record_haplotype_rank] {n} haplotypes, total_depth={total_depth}"
    );

    result
}

// ─── summarize_enclosing_haps ───────────────────────────────────────────────

/// Per-span haplotype info returned by [`summarize_enclosing_haps`].
pub struct RegionHapInfo<'a> {
    pub reads: Vec<&'a Record>,
    pub vert_inds: Vec<i32>,  // unique vertex indices from the reads
    pub hap_id: i32,
    pub qnames: Vec<String>,
}

/// Result of [`summarize_enclosing_haps`].
///
/// Returns `Some((region_haplotype_info, overlapping_span))` or `None` if
/// fewer than 2 usable haplotypes are found.
pub type SummarizeResult<'a> = Option<(HashMap<(i64, i64), RegionHapInfo<'a>>, (i64, i64))>;

/// For each haplotype cluster, find continuous regions that **enclose** (or
/// mostly overlap) the target window `[start, end)`.
///
/// Faithfully ports Python's `summarize_enclosing_haps`
/// (identify_misaligned_haps.py:573-636).
///
/// # Recovery logic
/// If fewer than 3 haplotypes fully enclose `[start, end)`:
/// - Try to recover haplotypes with overlap coefficient >= 0.8
/// - Shrink the window to the intersection of recovered spans
/// - If still <= 1 haplotype, returns `None`
///
/// # Arguments
/// - `hap_subgraphs` — output of `group_by_dict_optimized`
/// - `qname_to_node` — qname → vertex_idx mapping
/// - `region`        — `(chrom, start, end)` target window (chrom used only for logging)
pub fn summarize_enclosing_haps<'a>(
    hap_subgraphs: &HashMap<i32, HapGroup<'a>>,
    qname_to_node: &HashMap<String, i32>,
    region: (&str, i64, i64),
) -> SummarizeResult<'a> {
    let (chrom, start, end) = region;

    // {span → RegionHapInfo}  for spans that fully enclose [start, end)
    let mut region_haplotype_info: HashMap<(i64, i64), RegionHapInfo<'a>> = HashMap::new();

    // Partially-overlapping results for potential recovery
    struct InspectCandidate<'b> {
        hap_id: i32,
        qnames: Vec<String>,
        span: (i64, i64),
        reads: Vec<&'b Record>,
        overlap_coef: f64,
    }
    let mut inspect_results: Vec<InspectCandidate<'a>> = Vec::new();

    for (&hap_id, group) in hap_subgraphs.iter() {
        // Flatten all read_pair_lists into a single vec
        let all_reads: Vec<&Record> = group
            .read_pair_lists
            .iter()
            .flat_map(|rpl| rpl.iter().copied())
            .collect();

        // Get continuous covered regions
        let continuous = extract_continuous_regions_dict(&all_reads);

        for ((sp_start, sp_end), indices) in &continuous {
            let span = (*sp_start, *sp_end);
            let sreads: Vec<&Record> = indices.iter().map(|&i| all_reads[i]).collect();

            if span.0 <= start && span.1 >= end {
                // Fully enclosing
                let vert_inds: Vec<i32> = sreads
                    .iter()
                    .filter_map(|r| {
                        let qn = std::str::from_utf8(r.qname()).ok()?;
                        qname_to_node.get(qn).copied()
                    })
                    .collect::<std::collections::HashSet<i32>>()
                    .into_iter()
                    .collect();

                region_haplotype_info.insert(
                    span,
                    RegionHapInfo {
                        reads: sreads,
                        vert_inds,
                        hap_id,
                        qnames: group.qnames.clone(),
                    },
                );
            } else {
                // Partial overlap — candidate for recovery
                let overlap = span.1.min(end) - span.0.max(start);
                let window_len = end - start;
                let overlap_coef = if window_len > 0 {
                    overlap as f64 / window_len as f64
                } else {
                    0.0
                };
                inspect_results.push(InspectCandidate {
                    hap_id,
                    qnames: group.qnames.clone(),
                    span,
                    reads: sreads,
                    overlap_coef,
                });
            }
        }
    }

    // ── Enough enclosing haplotypes? ────────────────────────────────────
    if region_haplotype_info.len() > 2 {
        return Some((region_haplotype_info, (start, end)));
    }

    // ── Recovery: try to rescue partially-overlapping haplotypes ────────
    info!(
        "[summarize_enclosing_haps] At region {}:{}-{}, only found {} enclosing haps. Trying recovery.",
        chrom, start, end, region_haplotype_info.len()
    );

    // Sort by overlap_coef descending. `total_cmp` is a total order (never panics
    // on a non-finite coef; NaN sorts to the end of the descending order).
    inspect_results.sort_by(|a, b| b.overlap_coef.total_cmp(&a.overlap_coef));

    if inspect_results.is_empty() {
        return None;
    }

    if inspect_results[0].overlap_coef < 0.8 {
        // Even the best candidate has < 80% overlap
        return None;
    }

    // Keep candidates with overlap >= 0.8
    let recover_results: Vec<&InspectCandidate> = inspect_results
        .iter()
        .filter(|t| t.overlap_coef >= 0.8)
        .collect();

    // Shrink window: max of candidate starts, min of candidate ends
    let recover_start = recover_results
        .iter()
        .map(|t| t.span.0)
        .max()
        .unwrap();
    let recover_end = recover_results
        .iter()
        .map(|t| t.span.1)
        .min()
        .unwrap();

    let overlapping_span = (recover_start.max(start), recover_end.min(end));

    info!(
        "[summarize_enclosing_haps] Shrinked to ({}, {}), recovering {} haplotypes",
        overlapping_span.0, overlapping_span.1, recover_results.len()
    );

    for cand in &recover_results {
        let vert_inds: Vec<i32> = cand
            .reads
            .iter()
            .filter_map(|r| {
                let qn = std::str::from_utf8(r.qname()).ok()?;
                qname_to_node.get(qn).copied()
            })
            .collect::<std::collections::HashSet<i32>>()
            .into_iter()
            .collect();

        region_haplotype_info.insert(
            cand.span,
            RegionHapInfo {
                reads: cand.reads.clone(),
                vert_inds,
                hap_id: cand.hap_id,
                qnames: cand.qnames.clone(),
            },
        );
    }

    if region_haplotype_info.len() <= 1 {
        info!(
            "[summarize_enclosing_haps] After recovery, still only {} haps at {}:{}-{}",
            region_haplotype_info.len(), chrom, start, end
        );
        return None;
    }

    Some((region_haplotype_info, overlapping_span))
}

// ─── identify_misalignment_per_region ───────────────────────────────────────

/// Result returned by [`identify_misalignment_per_region`].
///
/// Contains the ranked haplotype 2D array (8 columns) and the chromosome name.
/// The 8 columns are: start, end, total_depth, hap_id, hap_depth, var_count,
/// indel_count, psv_count.
pub struct IdentifyMisalignmentResult {
    pub record_2d_arr: Array2<i32>,
    pub chrom: String,
}

/// Rank haplotypes enclosed by a region by sequence similarity metrics.
///
/// Faithfully ports Python's `identify_misalignment_per_region`
/// (`identify_misaligned_haps.py:637-735`).
///
/// # Algorithm
/// 1. Query overlapping reads from BAM lapper.
/// 2. Filter out low-quality reads and reads not in the phasing graph.
/// 3. Group reads by haplotype label via [`group_by_dict_optimized`].
/// 4. Identify enclosing haplotypes via [`summarize_enclosing_haps`].
/// 5. For each haplotype cluster:
///    a. Compute haplotype/error vectors via [`record_hap_err_vectors_per_region`].
///    b. Assemble consensus via [`assemble_consensus`].
///    c. Slice consensus to the overlapping window.
/// 6. Rank all clusters via [`record_haplotype_rank`] and return an `Array2<i32>`.
///
/// # Arguments
/// - `region` — `(chrom, start, end)` tuple for the genomic window.
/// - `lapper_dict` — per-chromosome Lapper intervals from `build_lapper_from_bam`.
/// - `read_dict` — per-qname-index read records from BAM lapper.
/// - `qname_hap_info` — `{vertex_idx → haplotype_id}` from the phasing graph.
/// - `qname_to_node` — `{qname → vertex_idx}` mapping.
/// - `lowqual_qnames` — set of qnames to ignore.
/// - `hap_max_psv_pos` — `{hap_id → [psv_positions]}` for PSV counting.
/// - `hap_cache` — mutable cache for haplotype vectors (String read_id → Array1<i16>).
/// - `err_cache` — mutable cache for error vectors (String read_id → Array1<f32>).
/// - `mean_read_length` — average read length (used to estimate total depth).
///
/// # Returns
/// `Ok(Some(IdentifyMisalignmentResult))` — ranked haplotype array + chrom,
/// `Ok(None)` if the region is skipped (no haplotypes, < 2 enclosing haps, etc.).
///
/// # Errors
/// Propagates [`crate::pairwise_read_inspection::CigarError`] from vector extraction
/// (DivA — a non-`=`/`X` CIGAR fails fast instead of panicking).
#[allow(clippy::too_many_arguments)]
pub fn identify_misalignment_per_region(
    region: (&str, i64, i64),
    lapper_dict: &HashMap<String, Lapper<u32, u32>>,
    read_dict: &FxHashMap<u32, Vec<Record>>,
    qname_hap_info: &HashMap<i32, i32>,
    qname_to_node: &HashMap<String, i32>,
    lowqual_qnames: &HashSet<String>,
    hap_max_psv_pos: &HashMap<i32, Vec<i32>>,
    hap_cache: &mut HashMap<String, Array1<i16>>,
    err_cache: &mut HashMap<String, Array1<f32>>,
    mean_read_length: f64,
) -> Result<Option<IdentifyMisalignmentResult>, Box<dyn std::error::Error>> {
    let (chrom, start, end) = region;
    let region_str = format!("{chrom}:{start}-{end}");

    // ── Step 1: Query overlapping reads ─────────────────────────────────
    let overlap_reads = query_overlapping_reads(
        lapper_dict,
        read_dict,
        chrom,
        start as u32,
        end as u32,
    );

    // ── Step 2: Filter + build vertices map ─────────────────────────────
    // vertices: (vertex_idx, qname) → Vec<&Record>
    let mut vertices: HashMap<(i32, String), Vec<&Record>> = HashMap::new();
    let mut vert_inds: HashSet<i32> = HashSet::new();

    for read in &overlap_reads {
        let qname = String::from_utf8_lossy(read.qname()).to_string();
        if lowqual_qnames.contains(&qname) || !qname_to_node.contains_key(&qname) {
            continue;
        }
        let vert_idx = qname_to_node[&qname];
        vertices
            .entry((vert_idx, qname))
            .or_default()
            .push(read);
        vert_inds.insert(vert_idx);
    }

    debug!(
        "[identify_misalignment_per_region] region {} — {} overlapping reads, {} vertices after filter",
        region_str,
        overlap_reads.len(),
        vert_inds.len()
    );

    // ── Step 3: Group by haplotype label ────────────────────────────────
    let hap_subgraphs = group_by_dict_optimized(qname_hap_info, &vertices);

    if hap_subgraphs.is_empty() {
        warn!(
            "[identify_misalignment_per_region] No haplotype clusters found for region {region_str}. Skipping."
        );
        return Ok(None);
    }

    // ── Step 4: Summarize enclosing haplotypes ──────────────────────────
    let summarize_result = summarize_enclosing_haps(&hap_subgraphs, qname_to_node, region);
    let (region_haplotype_info, overlapping_span) = match summarize_result {
        Some(v) => v,
        None => return Ok(None),
    };

    debug!(
        "[identify_misalignment_per_region] region {}:{}-{} → {} enclosing haplotypes",
        chrom,
        overlapping_span.0,
        overlapping_span.1,
        region_haplotype_info.len()
    );

    // ── Step 5: For each hap cluster → consensus → slice to overlapping window ──
    let mut final_clusters: HashMap<i32, HaplotypeClusterInfo> = HashMap::new();

    for (&(span_start, span_end), hap_info) in &region_haplotype_info {
        let haplotype_idx = hap_info.hap_id;
        let reads = &hap_info.reads;

        // Verify all vertices belong to the same haplotype
        debug_assert!(
            hap_info.vert_inds.iter().all(|vid| qname_hap_info.get(vid) == Some(&haplotype_idx)),
            "Vertices {:?} are not in the same connected component.",
            hap_info.vert_inds
        );

        if reads.is_empty() {
            warn!(
                "[identify_misalignment_per_region] No reads for haplotype across ({span_start}, {span_end}). Skipping cluster."
            );
            continue;
        }

        // Compute hap/err vectors for these reads
        let (read_spans, hap_vectors, err_vectors) =
            record_hap_err_vectors_per_region(reads.as_slice(), hap_cache, err_cache)?;

        // Assemble consensus
        let consensus_sequence = assemble_consensus(&hap_vectors, &err_vectors, &read_spans);

        // Slice consensus to the overlapping window
        // Python: overlapping_con_seq = consensus_sequence[overlapping_span[0] - span[0]:overlapping_span[1] - span[0] + 1]
        let slice_start = (overlapping_span.0 - span_start) as usize;
        let slice_end = (overlapping_span.1 - span_start + 1) as usize;
        let slice_end = slice_end.min(consensus_sequence.len());
        let overlapping_con_seq = if slice_start < consensus_sequence.len() {
            consensus_sequence.slice(ndarray::s![slice_start..slice_end]).to_owned()
        } else {
            warn!(
                "[identify_misalignment_per_region] Consensus too short for hap {} at span ({},{}), overlapping ({},{})",
                haplotype_idx, span_start, span_end, overlapping_span.0, overlapping_span.1
            );
            continue;
        };

        final_clusters.insert(haplotype_idx, HaplotypeClusterInfo {
            consensus: overlapping_con_seq,
            reads: reads.clone(),
            span: overlapping_span,
            qnames: hap_info.qnames.clone(),
        });
    }

    debug!(
        "[identify_misalignment_per_region] {} final clusters for region {}",
        final_clusters.len(),
        region_str
    );

    if final_clusters.is_empty() {
        warn!(
            "[identify_misalignment_per_region] Only {} haplotype clusters for region {}. Skipping.",
            final_clusters.len(),
            region_str
        );
        return Ok(None);
    }

    // ── Step 6: Rank haplotypes ─────────────────────────────────────────
    let record_2d_arr = record_haplotype_rank(
        &final_clusters,
        mean_read_length as i32,
        hap_max_psv_pos,
    );

    debug!(
        "[identify_misalignment_per_region] Haplotype rank array shape: {:?} for region {}",
        record_2d_arr.shape(),
        region_str
    );

    Ok(Some(IdentifyMisalignmentResult {
        record_2d_arr,
        chrom: chrom.to_string(),
    }))
}

/// Select genomic intervals covered by >= `min_haplotypes` distinct haplotypes.
///
/// Ports Python's `select_regions_with_min_haplotypes_from_hapbeds`
/// (identify_misaligned_haps.py:1089-1152).
///
/// In Python this function shells out to BEDOPS CLI tools (sort-bed, bedops --partition,
/// bedmap) operating on per-haplotype BED files on disk.  In Rust we perform the
/// equivalent "genome partitioning + haplotype counting" entirely in memory using a
/// coordinate sweep, which is both faster and avoids temporary files.
///
/// # Algorithm (mirrors BEDOPS --partition + bedmap)
/// 1. Collect **all** interval endpoints across every haplotype into a sorted,
///    deduplicated breakpoint list.
/// 2. For each consecutive pair of breakpoints `[bp_i, bp_{i+1})` (a partition),
///    count how many distinct haplotypes have an interval that fully contains
///    (or overlaps) this partition.
/// 3. Keep partitions where `count >= min_haplotypes`.
/// 4. Merge adjacent partitions on the same chromosome into contiguous intervals.
///
/// # Arguments
/// - `hid_cov_intervals` — per-haplotype coverage.
///   `HashMap<i32, Vec<(String, i64, i64)>>` where each tuple is
///   `(chrom, start, end)` — 0-based half-open, already sorted & merged per
///   haplotype (the Rust caller merges via `extract_continuous_regions_dict`
///   before calling this function).
/// - `min_haplotypes` — minimum number of distinct haplotypes required (default 2).
///
/// # Returns
/// `Option<Vec<(String, i64, i64)>>` — sorted, merged intervals where ≥ `min_haplotypes`
/// haplotypes overlap, or `None` if fewer haplotypes exist than the threshold or
/// no qualifying regions are found.
pub fn select_regions_with_min_haplotypes(
    hid_cov_intervals: &HashMap<i32, Vec<(String, i64, i64)>>,
    min_haplotypes: usize,
) -> Option<Vec<(String, i64, i64)>> {
    // Early return: not enough haplotypes to reach the threshold
    if hid_cov_intervals.len() < min_haplotypes {
        warn!(
            "[select_regions_with_min_haplotypes] Only {} haplotypes have coverage; need at least {}.",
            hid_cov_intervals.len(),
            min_haplotypes
        );
        return None;
    }

    // ── Step 1: Group intervals by chromosome ──────────────────────────
    // Key: chrom, Value: Vec<(start, end, hap_id)>
    let mut chrom_intervals: HashMap<String, Vec<(i64, i64, i32)>> = HashMap::new();
    for (&hid, intervals) in hid_cov_intervals {
        for (chrom, start, end) in intervals {
            chrom_intervals
                .entry(chrom.clone())
                .or_default()
                .push((*start, *end, hid));
        }
    }

    let mut result: Vec<(String, i64, i64)> = Vec::new();

    // ── Step 2: For each chromosome, sweep and partition ───────────────
    // Sort chromosomes for determinism (matches Python's sorted(hid_cov_beds))
    let mut chroms: Vec<&String> = chrom_intervals.keys().collect();
    chroms.sort();

    for chrom in chroms {
        let intervals = &chrom_intervals[chrom];

        // Collect all breakpoints (start and end of every interval)
        let mut breakpoints: Vec<i64> = Vec::with_capacity(intervals.len() * 2);
        for &(start, end, _) in intervals {
            breakpoints.push(start);
            breakpoints.push(end);
        }
        breakpoints.sort_unstable();
        breakpoints.dedup();

        if breakpoints.len() < 2 {
            continue;
        }

        // For each partition [bp_i, bp_{i+1}), count distinct haplotypes
        for window in breakpoints.windows(2) {
            let p_start = window[0];
            let p_end = window[1];
            if p_start >= p_end {
                continue;
            }

            // Count distinct haplotypes whose interval contains this partition
            let mut hap_ids: HashSet<i32> = HashSet::new();
            for &(start, end, hid) in intervals {
                // An interval [start, end) contains partition [p_start, p_end)
                // iff start <= p_start && end >= p_end
                if start <= p_start && end >= p_end {
                    hap_ids.insert(hid);
                }
            }

            if hap_ids.len() >= min_haplotypes {
                result.push((chrom.clone(), p_start, p_end));
            }
        }
    }

    if result.is_empty() {
        info!(
            "[select_regions_with_min_haplotypes] No regions found with >= {min_haplotypes} haplotypes."
        );
        return None;
    }

    // NOTE: Do NOT merge adjacent partitions here.
    // Python's BEDOPS pipeline (sort-bed) keeps individual partitions separate,
    // and downstream `identify_misalignment_per_region` expects fine-grained
    // intervals. Merging collapses many small partitions into large regions
    // that no single haplotype fully encloses, causing the LP solver path
    // to fail and fall back to heuristics.
    info!(
        "[select_regions_with_min_haplotypes] Found {} partition regions with >= {} haplotypes.",
        result.len(),
        min_haplotypes
    );

    Some(result)
}

/// Main haplotype inspection function.
///
/// Ports Python's `inspect_by_haplotypes` (identify_misaligned_haps.py:1153-1441).
///
/// For each haplotype, assembles consensus sequences per continuous region,
/// judges misalignment by variant density and similarity to intrinsic reference
/// sequences, then uses a Binary Integer Linear Programming model to classify
/// haplotypes as correctly mapped or mismapped.
///
/// # Returns
/// `(correct_qnames, mismap_qnames)` — two disjoint sets of read names.
#[allow(clippy::too_many_arguments)]
pub fn inspect_haplotypes(
    bam_path: &str,
    intrinsic_bam_path: &str,
    hap_qname_info: &HashMap<i32, Vec<String>>,
    qname_hap_info: &HashMap<i32, i32>,         // vertex_idx → hap_id
    qname_to_node: &HashMap<String, i32>,        // qname → vertex_idx
    total_lowqual_qnames: &HashSet<String>,
    compare_haplotype_meta_tab: &str,
    mean_read_length: f64,
    mapq_cutoff: u8,
    basequal_median_cutoff: u8,
) -> Result<(HashSet<String>, HashSet<String>), Box<dyn std::error::Error>> {

    // ══════════════════════════════════════════════════════════════════════
    // Phase 0: Build Lapper structures from BAM files
    // ══════════════════════════════════════════════════════════════════════
    info!("[inspect_haplotypes] Building Lapper from input BAM: {bam_path}");
    let bam_lapper = build_lapper_from_bam(bam_path, mapq_cutoff, basequal_median_cutoff, true, true)?;
    info!("[inspect_haplotypes] Building Lapper from intrinsic BAM: {intrinsic_bam_path}");
    let intrin_lapper = build_lapper_from_bam(intrinsic_bam_path, 0, 0, false, false)?;

    // Option A: Build qname → &Vec<Record> index for direct lookup (replaces
    // Python's qname → node → node_read_ids → read_id_read_dict chain).
    let mut qname_to_records: HashMap<&str, &Vec<Record>> = HashMap::new();
    for (qname, &qname_idx) in &bam_lapper.qname_idx_dict {
        if let Some(records) = bam_lapper.read_dict.get(&qname_idx) {
            qname_to_records.insert(qname.as_str(), records);
        }
    }
    info!("[inspect_haplotypes] Built qname→records index with {} entries", qname_to_records.len());

    // Build qname → chrom index from lapper_dict (each chrom's Lapper contains
    // intervals whose val is qname_idx; use qname_dict to resolve to qname string).
    let mut qname_to_chrom: HashMap<&str, String> = HashMap::new();
    for (chrom, lapper) in &bam_lapper.lapper_dict {
        for iv in lapper.iter() {
            if let Some(qname) = bam_lapper.qname_dict.get(&iv.val) {
                qname_to_chrom.entry(qname.as_str()).or_insert_with(|| chrom.clone());
            }
        }
    }

    // ══════════════════════════════════════════════════════════════════════
    // Phase 1: Haplotype iteration loop (Python lines 1169-1257)
    // ══════════════════════════════════════════════════════════════════════
    let mut hid_extreme_vard: HashMap<i32, bool> = HashMap::new();
    let mut hid_var_count: HashMap<i32, i32> = HashMap::new();
    let mut scatter_hid_dict: HashMap<i32, bool> = HashMap::new();
    let mut hid_max_local_density: HashMap<i32, f64> = HashMap::new();
    let mut varcounts_among_refseqs: VarcountsAmongRefseqs = HashMap::new();
    let mut total_qnames: HashSet<String> = HashSet::new();
    let mut hid_cov_beds: HashMap<i32, Vec<(String, i64, i64)>> = HashMap::new();

    // Caches that persist across all haplotypes and regions
    let mut hap_cache: HashMap<String, Array1<i16>> = HashMap::new();
    let mut err_cache: HashMap<String, Array1<f32>> = HashMap::new();
    let mut total_genomic_haps: HashMap<String, Array1<i16>> = HashMap::new();
    let mut qseq_cache: HashMap<String, ReadQseqData> = HashMap::new();

    info!("[inspect_haplotypes] All the haplotype IDs are: {:?}", hap_qname_info.keys().collect::<Vec<_>>());

    for (&hid, qnames_raw) in hap_qname_info {
        // Filter out low-quality qnames
        let qnames: Vec<String> = qnames_raw.iter()
            .filter(|qn| !total_lowqual_qnames.contains(qn.as_str()))
            .cloned()
            .collect();

        total_qnames.extend(qnames.iter().cloned());

        if qnames.len() < 3 {
            scatter_hid_dict.insert(hid, true);
        }

        debug!("[inspect_haplotypes] haplotype {} contains {} read pairs", hid, qnames.len());

        // Collect all Records for this haplotype's qnames (Option A lookup)
        let reads: Vec<&Record> = qnames.iter()
            .flat_map(|qn| {
                qname_to_records.get(qn.as_str())
                    .into_iter()
                    .flat_map(|recs| recs.iter())
            })
            .collect();

        if reads.is_empty() {
            debug!("[inspect_haplotypes] haplotype {hid} has no reads after lookup, skipping");
            continue;
        }

        // Extract continuous regions
        let conregion_dict = extract_continuous_regions_dict(&reads);

        // Build hid_cov_beds for this haplotype
        // Determine chrom from the first qname that has a chrom mapping
        let hap_chrom: Option<String> = qnames.iter()
            .find_map(|qn| qname_to_chrom.get(qn.as_str()).cloned());

        let mut cov_intervals: Vec<(String, i64, i64)> = Vec::new();
        if let Some(ref chrom) = hap_chrom {
            for ((start, end), _) in &conregion_dict {
                cov_intervals.push((chrom.clone(), *start, *end));
            }
        }
        hid_cov_beds.insert(hid, cov_intervals);

        // Per-region inner loop
        let chrom_str = hap_chrom.unwrap_or_default();
        for ((span_start, span_end), ref read_indices) in &conregion_dict {
            let region_reads: Vec<&Record> = read_indices.iter().map(|&i| reads[i]).collect();
            if region_reads.is_empty() {
                continue;
            }
            let span = (*span_start, *span_end);

            debug!("[inspect_haplotypes] haplotype {} — region {}:{}-{} with {} reads",
                   hid, chrom_str, span.0, span.1, region_reads.len());

            // a) Record haplotype/error vectors
            let (read_spans, hap_vectors, err_vectors) =
                record_hap_err_vectors_per_region(&region_reads, &mut hap_cache, &mut err_cache)?;

            // b) Assemble consensus
            let consensus_sequence = assemble_consensus(&hap_vectors, &err_vectors, &read_spans);

            // c) Judge extreme variant density
            let (extreme_vard, region_max_density) =
                judge_misalignment_by_extreme_vardensity(&consensus_sequence);

            let prev_density = *hid_max_local_density.get(&hid).unwrap_or(&0.0);
            hid_max_local_density.insert(hid, prev_density.max(region_max_density as f64));

            // d) Count variants
            let var_count = count_var(&consensus_sequence);
            *hid_var_count.entry(hid).or_insert(0) += var_count;

            let prev_extreme = *hid_extreme_vard.get(&hid).unwrap_or(&false);
            hid_extreme_vard.insert(hid, extreme_vard || prev_extreme);

            debug!("[inspect_haplotypes] hid {} region {}:{}-{} extreme_vard={} var_count={} total_var_count={}",
                   hid, chrom_str, span.0, span.1, extreme_vard, var_count, hid_var_count[&hid]);

            // e) Stat refseq similarity
            stat_refseq_similarity(
                &intrin_lapper,
                &chrom_str,
                span,
                hid,
                &consensus_sequence,
                &region_reads,
                &mut total_genomic_haps,
                &mut qseq_cache,
                &mut varcounts_among_refseqs,
            )?;
        }
    }

    // Post-loop logging
    info!("[inspect_haplotypes] extreme variant density haplotypes: {hid_extreme_vard:?}");
    info!("[inspect_haplotypes] scatter haplotypes: {scatter_hid_dict:?}");
    let scatter_qnames: HashSet<String> = scatter_hid_dict.iter()
        .filter(|(_, &v)| v)
        .flat_map(|(&hid, _)| hap_qname_info.get(&hid).into_iter().flatten().cloned())
        .collect();
    info!("[inspect_haplotypes] {} scatter qnames", scatter_qnames.len());

    // ══════════════════════════════════════════════════════════════════════
    // Phase 2: Similarity scoring (Python lines 1266-1275)
    // ══════════════════════════════════════════════════════════════════════
    let sim_scores = cal_similarity_score(&varcounts_among_refseqs, &hid_var_count, &hid_max_local_density);
    let mut hap_max_sim_scores: HashMap<i32, f64> = HashMap::new();
    let mut hap_max_psvs: HashMap<i32, i32> = HashMap::new();
    let mut hap_max_psv_pos: HashMap<i32, Vec<i32>> = HashMap::new();
    for (hid, score) in &sim_scores {
        hap_max_sim_scores.insert(*hid, score.max_sim_score);
        hap_max_psvs.insert(*hid, score.max_psv_count);
        hap_max_psv_pos.insert(*hid, score.max_psv_positions.clone());
    }
    // Fill missing haplotypes with empty PSV positions
    for hid in hap_qname_info.keys() {
        hap_max_psv_pos.entry(*hid).or_default();
    }
    drop(varcounts_among_refseqs);
    info!("[inspect_haplotypes] similarity scores: {hap_max_sim_scores:?}");
    info!("[inspect_haplotypes] PSV counts: {hap_max_psvs:?}");

    // ══════════════════════════════════════════════════════════════════════
    // Phase 3: Sweep region selection (Python lines 1277-1292)
    // ══════════════════════════════════════════════════════════════════════
    let sweep_regions = select_regions_with_min_haplotypes(&hid_cov_beds, 2);

    // ══════════════════════════════════════════════════════════════════════
    // Phase 4: Early return if no sweep regions (Python lines 1298-1305)
    // ══════════════════════════════════════════════════════════════════════
    if sweep_regions.is_none() {
        let mismap_hids: HashSet<i32> = hid_extreme_vard.iter()
            .filter(|(_, &v)| v).map(|(&k, _)| k).collect();
        let mut mismap_qnames: HashSet<String> = mismap_hids.iter()
            .flat_map(|hid| hap_qname_info.get(hid).into_iter().flatten().cloned())
            .collect();
        mismap_qnames.extend(scatter_qnames);
        let correct_qnames: HashSet<String> = total_qnames.difference(&mismap_qnames).cloned().collect();
        warn!("[inspect_haplotypes] No sweep regions found. Filtering {} mismap, {} correct by extreme_vard only.",
              mismap_qnames.len(), correct_qnames.len());
        return Ok((correct_qnames, mismap_qnames));
    }
    let sweep_regions = sweep_regions.unwrap();
    info!("[inspect_haplotypes] Found {} sweep regions", sweep_regions.len());

    // ══════════════════════════════════════════════════════════════════════
    // Phase 5: Per-region inspection loop (Python lines 1308-1329)
    // ══════════════════════════════════════════════════════════════════════
    let mut record_results: Vec<IdentifyMisalignmentResult> = Vec::new();

    for (chrom, start, end) in &sweep_regions {
        if let Some(result) = identify_misalignment_per_region(
            (chrom.as_str(), *start, *end),
            &bam_lapper.lapper_dict,
            &bam_lapper.read_dict,
            qname_hap_info,
            qname_to_node,
            total_lowqual_qnames,
            &hap_max_psv_pos,
            &mut hap_cache,
            &mut err_cache,
            mean_read_length,
        )? {
            record_results.push(result);
        } else {
            warn!("[inspect_haplotypes] No valid haplotypes for region {chrom}:{start}-{end}");
        }
    }

    // ══════════════════════════════════════════════════════════════════════
    // Phase 6: BILC pre-processing (Python lines 1331-1378)
    // ══════════════════════════════════════════════════════════════════════
    let mut failed_lp = false;
    let mut remove_hids: HashSet<i32> = HashSet::new();

    let mut total_records: Vec<EnrichedRecord> = Vec::new();

    if !record_results.is_empty() {
        // (a) Flatten record_results into EnrichedRecord rows
        for result in &record_results {
            let arr = &result.record_2d_arr;
            let chrom = &result.chrom;
            for row_idx in 0..arr.nrows() {
                let hap_id = arr[[row_idx, 3]];
                let sim_score = *hap_max_sim_scores.get(&hap_id).unwrap_or(&0.0);
                // Round to 1 decimal place (matching Python's round(1))
                let sim_score_rounded = (sim_score * 10.0).round() / 10.0;

                total_records.push(EnrichedRecord {
                    chrom: chrom.clone(),
                    start: arr[[row_idx, 0]],
                    end: arr[[row_idx, 1]],
                    total_depth: arr[[row_idx, 2]],
                    hap_id,
                    hap_depth: arr[[row_idx, 4]],
                    var_count: arr[[row_idx, 5]],
                    indel_count: arr[[row_idx, 6]],
                    psv_count: arr[[row_idx, 7]],
                    extreme_vard: *hid_extreme_vard.get(&hap_id).unwrap_or(&false),
                    scatter_hap: *scatter_hid_dict.get(&hap_id).unwrap_or(&false),
                    hap_var_count: *hid_var_count.get(&hap_id).unwrap_or(&0),
                    hap_max_sim_scores: sim_score_rounded,
                    hap_max_psvs: *hap_max_psvs.get(&hap_id).unwrap_or(&0),
                    varc_rank: 0,
                    rank: 0,
                    interval_coefficient: 0.0,
                    coefficient: 0.0,
                });
            }
        }

        // (b) Determine remove_hids: haps with scatter OR sim>10 OR psvs>=12 OR extreme_vard
        let candidate_remove: HashSet<i32> = total_records.iter()
            .filter(|r| r.scatter_hap || r.hap_max_sim_scores > 10.0
                       || r.hap_max_psvs >= 12 || r.extreme_vard)
            .map(|r| r.hap_id)
            .collect();

        // Except kept_scatter_hids
        let kept_scatter_hids: HashSet<i32> = total_records.iter()
            .filter(|r| r.hap_var_count >= 1 && r.scatter_hap
                       && r.hap_max_psvs < 12 && r.hap_max_sim_scores <= 10.0
                       && !r.extreme_vard)
            .map(|r| r.hap_id)
            .collect();

        remove_hids = candidate_remove.difference(&kept_scatter_hids).copied().collect();
        info!("[inspect_haplotypes] remove_hids (pre-ILP filter): {remove_hids:?}");

        // (c) Filter out remove_hids and low-depth rows
        total_records.retain(|r| !remove_hids.contains(&r.hap_id) && r.total_depth > 5);

        if total_records.is_empty() {
            warn!("[inspect_haplotypes] No records left after filtering. Setting failed_lp=true.");
            failed_lp = true;
        }
    } else {
        warn!("[inspect_haplotypes] No record_results from region inspection. Setting failed_lp=true.");
        failed_lp = true;
    }

    // ══════════════════════════════════════════════════════════════════════
    // Phase 7: Coefficient calculation (Python lines 1380-1407)
    // ══════════════════════════════════════════════════════════════════════
    if !failed_lp {
        // (a) Per-region: rank sim_scores and compute interval_coefficient
        //     Group records by (chrom, start, end)
        let mut region_groups: HashMap<RegionKey, Vec<usize>> = HashMap::new();
        for (i, r) in total_records.iter().enumerate() {
            let key = RegionKey { chrom: r.chrom.clone(), start: r.start, end: r.end };
            region_groups.entry(key).or_default().push(i);
        }

        for indices in region_groups.values() {
            // rank_unique_values on rounded hap_max_sim_scores
            let sim_arr: Vec<f32> = indices.iter()
                .map(|&i| (total_records[i].hap_max_sim_scores * 10.0).round() as f32 / 10.0)
                .collect();
            let ranks = rank_unique_values(&sim_arr);
            for (j, &i) in indices.iter().enumerate() {
                total_records[i].varc_rank = ranks[j];
            }

            // calculate_coefficient: build 10-col f32 rows
            let rows_data: Vec<[f32; 10]> = indices.iter().map(|&i| {
                let r = &total_records[i];
                [
                    r.start as f32, r.end as f32, r.total_depth as f32,
                    r.hap_id as f32, r.hap_depth as f32, r.var_count as f32,
                    r.indel_count as f32, r.psv_count as f32,
                    r.varc_rank as f32, r.hap_max_sim_scores as f32,
                ]
            }).collect();
            let rows_refs: Vec<&[f32; 10]> = rows_data.iter().collect();
            let coefficients = calculate_coefficient(&rows_refs);
            for (j, &i) in indices.iter().enumerate() {
                total_records[i].interval_coefficient = coefficients[j] as f64;
            }
        }

        // (b) Span-weighted average coefficient per hap_id
        let mut hap_coeff_sum: HashMap<i32, f64> = HashMap::new();
        let mut hap_span_sum: HashMap<i32, f64> = HashMap::new();
        for r in &total_records {
            let span = (r.end - r.start) as f64;
            *hap_coeff_sum.entry(r.hap_id).or_insert(0.0) += r.interval_coefficient;
            *hap_span_sum.entry(r.hap_id).or_insert(0.0) += span;
        }
        for r in &mut total_records {
            let span_sum = *hap_span_sum.get(&r.hap_id).unwrap_or(&1.0);
            let coeff_sum = *hap_coeff_sum.get(&r.hap_id).unwrap_or(&0.0);
            r.coefficient = if span_sum > 0.0 { coeff_sum / span_sum } else { 0.0 };
        }

        // (c) Add hap_max_sim_scores to coefficient (Python line 1398)
        for r in &mut total_records {
            r.coefficient += r.hap_max_sim_scores;
        }

        // (d) Rank coefficient per region group (Python lines 1401-1406)
        for indices in region_groups.values() {
            let coeff_arr: Vec<f32> = indices.iter()
                .map(|&i| ((total_records[i].coefficient * 100.0).round() / 100.0) as f32)
                .collect();
            let ranks = rank_unique_values(&coeff_arr);
            // Python stores the coefficient rank in a SEPARATE "rank" column
            // (lines 1401-1406); the BILC solver keys off varc_rank (the
            // sim-score rank set in step (a)), so keep varc_rank intact here.
            for (j, &i) in indices.iter().enumerate() {
                total_records[i].rank = ranks[j];
            }
        }

        // (e) Deduplicate (Python line 1407: drop_duplicates)
        total_records.sort_by(|a, b| {
            a.chrom.cmp(&b.chrom)
                .then(a.start.cmp(&b.start))
                .then(a.end.cmp(&b.end))
                .then(a.hap_id.cmp(&b.hap_id))
        });
        total_records.dedup_by(|a, b| {
            a.chrom == b.chrom && a.start == b.start && a.end == b.end && a.hap_id == b.hap_id
        });
    }

    // ══════════════════════════════════════════════════════════════════════
    // Phase 8: ILP solver (Python lines 1410-1413)
    // ══════════════════════════════════════════════════════════════════════
    let mut select_hids: FxHashSet<i32> = FxHashSet::default();
    let mut drop_hids: FxHashSet<i32> = FxHashSet::default();

    if !failed_lp {
        let bilc_records: Vec<BilcRecord> = total_records.iter().map(|r| BilcRecord {
            chrom: r.chrom.clone(),
            start: r.start,
            end: r.end,
            hap_id: r.hap_id,
            coefficient: r.coefficient,
            var_count: r.var_count,
            varc_rank: r.varc_rank,
        }).collect();

        let (s, d, status) = lp_solve_remained_haplotypes(&bilc_records);
        select_hids = s;
        drop_hids = d;
        info!("[inspect_haplotypes] ILP status: {}, select={}, drop={}", status, select_hids.len(), drop_hids.len());
        if matches!(status, BilcStatus::Infeasible) {
            warn!("[inspect_haplotypes] ILP infeasible. Setting failed_lp=true.");
            failed_lp = true;
        }
    }

    // ══════════════════════════════════════════════════════════════════════
    // Phase 9: Post-ILP augmentation (Python lines 1415-1432)
    // ══════════════════════════════════════════════════════════════════════
    if !failed_lp {
        let mut mismap_hids: HashSet<i32> = drop_hids.into_iter().collect();

        // Augment with extreme_vard
        mismap_hids.extend(hid_extreme_vard.iter().filter(|(_, &v)| v).map(|(&k, _)| k));
        // Augment with sim_score > 15
        mismap_hids.extend(hap_max_sim_scores.iter().filter(|(_, &v)| v > 15.0).map(|(&k, _)| k));
        // Augment with scatter + low var_count
        mismap_hids.extend(scatter_hid_dict.iter()
            .filter(|(&hid, &v)| v && *hid_var_count.get(&hid).unwrap_or(&0) <= 3)
            .map(|(&k, _)| k));
        // Augment with pre-ILP remove_hids
        mismap_hids.extend(&remove_hids);

        let correct_map_hids: HashSet<i32> = select_hids.iter()
            .filter(|hid| !mismap_hids.contains(hid))
            .copied()
            .collect();

        let mismap_qnames: HashSet<String> = mismap_hids.iter()
            .flat_map(|hid| hap_qname_info.get(hid).into_iter().flatten().cloned())
            .collect();
        let all_qnames: HashSet<String> = hap_qname_info.values()
            .flat_map(|qns| qns.iter().cloned())
            .collect();
        let correct_qnames: HashSet<String> = all_qnames.difference(&mismap_qnames).cloned().collect();

        // Write TSV if path is non-empty
        if !compare_haplotype_meta_tab.is_empty() {
            write_haplotype_meta_tsv(&total_records, &mismap_hids, &correct_map_hids, compare_haplotype_meta_tab);
        }

        info!("[inspect_haplotypes] ILP result: {} mismap hids, {} mismap qnames, {} correct qnames",
              mismap_hids.len(), mismap_qnames.len(), correct_qnames.len());
        return Ok((correct_qnames, mismap_qnames));
    }

    // ══════════════════════════════════════════════════════════════════════
    // Phase 10: Failed LP fallback (Python lines 1434-1441)
    // ══════════════════════════════════════════════════════════════════════
    let mut mismap_hids: HashSet<i32> = hid_extreme_vard.iter()
        .filter(|(_, &v)| v).map(|(&k, _)| k).collect();
    mismap_hids.extend(scatter_hid_dict.iter().filter(|(_, &v)| v).map(|(&k, _)| k));
    mismap_hids.extend(hap_max_sim_scores.iter().filter(|(_, &v)| v > 15.0).map(|(&k, _)| k));

    let mismap_qnames: HashSet<String> = mismap_hids.iter()
        .flat_map(|hid| hap_qname_info.get(hid).into_iter().flatten().cloned())
        .collect();
    let all_qnames: HashSet<String> = hap_qname_info.values()
        .flat_map(|qns| qns.iter().cloned())
        .collect();
    let correct_qnames: HashSet<String> = all_qnames.difference(&mismap_qnames).cloned().collect();

    warn!("[inspect_haplotypes] Failed LP fallback: {} mismap hids, {} mismap qnames, {} correct qnames",
          mismap_hids.len(), mismap_qnames.len(), correct_qnames.len());
    Ok((correct_qnames, mismap_qnames))
}

/// Write the enriched haplotype metadata table to a TSV file.
/// Corresponds to Python's `total_record_df.to_csv(compare_haplotype_meta_tab, ...)`.
fn write_haplotype_meta_tsv(
    records: &[EnrichedRecord],
    mismap_hids: &HashSet<i32>,
    correct_map_hids: &HashSet<i32>,
    path: &str,
) {
    use std::io::Write;
    let file = match std::fs::File::create(path) {
        Ok(f) => f,
        Err(e) => {
            warn!("[write_haplotype_meta_tsv] Failed to create {path}: {e}");
            return;
        }
    };
    let mut w = std::io::BufWriter::new(file);

    // Header
    let _ = writeln!(w, "chrom\tstart\tend\ttotal_depth\thap_id\thap_depth\tvar_count\t\
                         indel_count\tpsv_count\textreme_vard\tscatter_hap\thap_var_count\t\
                         hap_max_sim_scores\thap_max_psvs\tcoefficient\tvarc_rank\trank\t\
                         interval_coefficient\tmismap\tcorrect_map");
    for r in records {
        let _ = writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.chrom, r.start, r.end, r.total_depth, r.hap_id, r.hap_depth,
            r.var_count, r.indel_count, r.psv_count,
            r.extreme_vard, r.scatter_hap, r.hap_var_count,
            r.hap_max_sim_scores, r.hap_max_psvs,
            r.coefficient, r.varc_rank, r.rank, r.interval_coefficient,
            mismap_hids.contains(&r.hap_id),
            correct_map_hids.contains(&r.hap_id));
    }
    info!("[write_haplotype_meta_tsv] Saved {} rows to {}", records.len(), path);
}

#[cfg(test)]
mod tests {
    use super::*;
    use log::debug;
    use ndarray::{Array1, Array2, array};

    // ── Array2 helper functions for assemble_consensus tests ─────────


    // ── Helpers to build Array2 from test data ──────────────────────────────

    /// Build Array2<i16> from slices of equal-length rows.
    fn make_seq_array(data: &[&[i16]]) -> Array2<i16> {
        if data.is_empty() {
            return Array2::<i16>::zeros((0, 0));
        }
        let nrows = data.len();
        let ncols = data[0].len();
        let mut arr = Array2::<i16>::zeros((nrows, ncols));
        for (i, row) in data.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                arr[[i, j]] = val;
            }
        }
        arr
    }

    /// Build Array2<f32> from slices of equal-length rows.
    fn make_qual_array(data: &[&[f32]]) -> Array2<f32> {
        if data.is_empty() {
            return Array2::<f32>::zeros((0, 0));
        }
        let nrows = data.len();
        let ncols = data[0].len();
        let mut arr = Array2::<f32>::zeros((nrows, ncols));
        for (i, row) in data.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                arr[[i, j]] = val;
            }
        }
        arr
    }

    /// Build Array2<i32> (n × 2) from `[start, end]` pairs.
    fn make_spans_array(data: &[[i32; 2]]) -> Array2<i32> {
        let nrows = data.len();
        let mut arr = Array2::<i32>::zeros((nrows, 2));
        for (i, span) in data.iter().enumerate() {
            arr[[i, 0]] = span[0];
            arr[[i, 1]] = span[1];
        }
        arr
    }

    /// Build a padded Array2<i16> from variable-length rows (pad with HAP_PAD = -20).
    fn make_padded_seq_array(data: &[Vec<i16>]) -> Array2<i16> {
        if data.is_empty() {
            return Array2::<i16>::zeros((0, 0));
        }
        let nrows = data.len();
        let ncols = data.iter().map(|r| r.len()).max().unwrap();
        let mut arr = Array2::<i16>::from_elem((nrows, ncols), HAP_PAD);
        for (i, row) in data.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                arr[[i, j]] = val;
            }
        }
        arr
    }

    /// Build a padded Array2<f32> from variable-length rows (pad with -10.0).
    fn make_padded_qual_array(data: &[Vec<f32>]) -> Array2<f32> {
        if data.is_empty() {
            return Array2::<f32>::zeros((0, 0));
        }
        let nrows = data.len();
        let ncols = data.iter().map(|r| r.len()).max().unwrap();
        let mut arr = Array2::<f32>::from_elem((nrows, ncols), -10.0f32);
        for (i, row) in data.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                arr[[i, j]] = val;
            }
        }
        arr
    }

    use rust_htslib::bam::record::{Cigar, CigarString};
    use rust_htslib::bam::Record as BamRecord;


    /// Helper: create a BAM record with given CIGAR, sequence, qualities, pos, and name.
    fn make_record(cigar: CigarString, seq: &[u8], qual: &[u8], pos: i64) -> BamRecord {
        make_named_record(b"testread", cigar, seq, qual, pos)
    }

    fn make_named_record(
        name: &[u8],
        cigar: CigarString,
        seq: &[u8],
        qual: &[u8],
        pos: i64,
    ) -> BamRecord {
        let mut record = BamRecord::new();
        record.set(name, Some(&cigar), seq, qual);
        record.set_pos(pos);
        record.set_tid(0);
        record.set_mapq(60);
        record
    }

    // ── Tests for assemble_consensus ─────────────────────────────────────────

    /// Test basic single-read consensus (no padding needed).
    #[test]
    fn test_single_read_no_padding() {
        // One read spanning positions 100..104 with 4 valid bases
        let seq_arrays = make_seq_array(&[&[1i16, -4, 1, 1]]);
        let qual_arrays = make_qual_array(&[&[0.01f32, 0.05, 0.01, 0.01]]);
        let read_spans = make_spans_array(&[[100i32, 104]]);

        let result = assemble_consensus(&seq_arrays, &qual_arrays, &read_spans);
        assert_eq!(result.to_vec(), vec![1i16, -4, 1, 1]);
    }

    /// Test single read with HAP_PAD (-20) padding at the tail.
    #[test]
    fn test_single_read_with_padding() {
        // 4 valid bases + 2 hap-padding values (HAP_PAD = -20). The qual padding
        // stays -10.0 (it does not collide; truncated by the seq-derived count).
        let seq_arrays = make_seq_array(&[&[1i16, -4, 1, 1, -20, -20]]);
        let qual_arrays = make_qual_array(&[&[0.01f32, 0.05, 0.01, 0.01, -10.0, -10.0]]);
        let read_spans = make_spans_array(&[[100i32, 104]]);

        let result = assemble_consensus(&seq_arrays, &qual_arrays, &read_spans);
        // Only 4 valid values; padding (-20) is stripped
        assert_eq!(result.to_vec(), vec![1i16, -4, 1, 1]);
    }

    /// Test two overlapping reads where the better-quality read wins.
    #[test]
    fn test_two_reads_quality_selection() {
        // Read 0: positions 100..103, 3 valid bases
        //   pos 100: seq=1,  qual=0.10
        //   pos 101: seq=-4, qual=0.10
        //   pos 102: seq=1,  qual=0.10
        //
        // Read 1: positions 101..104, 3 valid bases
        //   pos 101: seq=-4, qual=0.05   <-- better quality, wins at pos 101
        //   pos 102: seq=-10, qual=0.15  <-- worse quality than Read 0
        //   pos 103: seq=1,  qual=0.05
        let seq_arrays = make_seq_array(&[
            &[1i16, -4, 1],
            &[-4i16, -10, 1],
        ]);
        let qual_arrays = make_qual_array(&[
            &[0.10f32, 0.10, 0.10],
            &[0.05f32, 0.15, 0.05],
        ]);
        let read_spans = make_spans_array(&[[100i32, 103], [101i32, 104]]);

        let result = assemble_consensus(&seq_arrays, &qual_arrays, &read_spans);
        // Length = 104 - 100 = 4
        // pos 100 (rel 0): only Read 0 → seq=1
        // pos 101 (rel 1): Read 0 qual=0.10, Read 1 qual=0.05 → Read 1 wins → seq=-4
        // pos 102 (rel 2): Read 0 qual=0.10, Read 1 qual=0.15 → Read 0 stays (0.10 < 0.15) → seq=1
        // pos 103 (rel 3): only Read 1 → seq=1
        assert_eq!(result.to_vec(), vec![1i16, -4, 1, 1]);
    }

    /// Test that reads with quality > 0.2 do NOT update the consensus.
    #[test]
    fn test_quality_above_threshold_ignored() {
        let seq_arrays = make_seq_array(&[&[-4i16, -4, -4]]);
        let qual_arrays = make_qual_array(&[&[0.5f32, 0.3, 0.21]]); // all > 0.2
        let read_spans = make_spans_array(&[[100i32, 103]]);

        let result = assemble_consensus(&seq_arrays, &qual_arrays, &read_spans);
        // All positions stay at default (1) since quality is too poor
        assert_eq!(result.to_vec(), vec![1i16, 1, 1]);
    }

    /// Test boundary: quality exactly 0.2 SHOULD update (<=).
    #[test]
    fn test_quality_at_threshold_accepted() {
        let seq_arrays = make_seq_array(&[&[-4i16, -4]]);
        let qual_arrays = make_qual_array(&[&[0.2f32, 0.2]]); // exactly 0.2
        let read_spans = make_spans_array(&[[100i32, 102]]);

        let result = assemble_consensus(&seq_arrays, &qual_arrays, &read_spans);
        // qual=0.2 <= consensus_qual=0.2 AND qual=0.2 <= 0.2 → updates
        assert_eq!(result.to_vec(), vec![-4i16, -4]);
    }

    /// Test empty input returns empty consensus.
    #[test]
    fn test_empty_input() {
        let seq_arrays = Array2::<i16>::zeros((0, 0));
        let qual_arrays = Array2::<f32>::zeros((0, 0));
        let read_spans = Array2::<i32>::zeros((0, 2));

        let result = assemble_consensus(&seq_arrays, &qual_arrays, &read_spans);
        assert_eq!(result.len(), 0);
    }

    /// Test three reads with staggered positions and mixed qualities.
    #[test]
    fn test_three_reads_staggered() {
        // Read 0: pos 10..13, seq=[1, -4, 1], qual=[0.05, 0.15, 0.05]
        // Read 1: pos 11..14, seq=[1,  1, -10], qual=[0.10, 0.01, 0.10]
        // Read 2: pos 12..15, seq=[-4, 1,  1], qual=[0.02, 0.02, 0.02]
        let seq_arrays = make_seq_array(&[
            &[1i16, -4, 1],
            &[1i16, 1, -10],
            &[-4i16, 1, 1],
        ]);
        let qual_arrays = make_qual_array(&[
            &[0.05f32, 0.15, 0.05],
            &[0.10f32, 0.01, 0.10],
            &[0.02f32, 0.02, 0.02],
        ]);
        let read_spans = make_spans_array(&[[10i32, 13], [11i32, 14], [12i32, 15]]);

        let result = assemble_consensus(&seq_arrays, &qual_arrays, &read_spans);
        // Length = 15 - 10 = 5
        // pos 10 (rel 0): Read 0 only → seq=1, qual=0.05
        // pos 11 (rel 1): Read 0 (seq=-4, q=0.15) then Read 1 (seq=1, q=0.10)
        //   → R0: 0.15 <= 0.2 → update → consensus(-4, 0.15)
        //   → R1: 0.10 <= 0.15 → update → consensus(1, 0.10)
        //   → final: seq=1
        // pos 12 (rel 2): Read 0 (seq=1, q=0.05) then R1 (seq=1, q=0.01) then R2 (seq=-4, q=0.02)
        //   → R0: 0.05 <= 0.2 → update → consensus(1, 0.05)
        //   → R1: 0.01 <= 0.05 → update → consensus(1, 0.01)
        //   → R2: 0.02 <= 0.01 → NO (0.02 > 0.01)
        //   → final: seq=1
        // pos 13 (rel 3): Read 1 (seq=-10, q=0.10) then R2 (seq=1, q=0.02)
        //   → R1: 0.10 <= 0.2 → update → consensus(-10, 0.10)
        //   → R2: 0.02 <= 0.10 → update → consensus(1, 0.02)
        //   → final: seq=1
        // pos 14 (rel 4): Read 2 only → seq=1, qual=0.02
        assert_eq!(result.to_vec(), vec![1i16, 1, 1, 1, 1]);
    }

    /// Test that the first valid read always overwrites the default consensus
    /// when its quality is <= 0.2, even if the sequence value equals the default.
    #[test]
    fn test_default_overwrite_with_same_value() {
        // The consensus default is 1 with qual 0.2.
        // A read with seq=1 and qual=0.1 should still "overwrite" (lowering the
        // consensus quality), so a later read with qual=0.15 cannot overwrite.
        let seq_arrays = make_seq_array(&[
            &[1i16],   // same as default
            &[-4i16],  // wants to overwrite
        ]);
        let qual_arrays = make_qual_array(&[
            &[0.1f32],  // better than default 0.2
            &[0.15f32], // worse than 0.1, cannot overwrite
        ]);
        let read_spans = make_spans_array(&[[50i32, 51], [50i32, 51]]);

        let result = assemble_consensus(&seq_arrays, &qual_arrays, &read_spans);
        // Read 0 sets consensus to (1, 0.1).
        // Read 1: 0.15 <= 0.1 → false → no overwrite.
        assert_eq!(result.to_vec(), vec![1i16]);
    }

    /// Golden collision fix: the deletion signal (-10) is real data; only HAP_PAD
    /// (-20) is stripped. The old `>= -8` filter would have wrongly dropped -10.
    #[test]
    fn test_hap_pad_stripped_deletion_kept() {
        // Valid data [1, -4, -10(deletion)], then padding [-20 = HAP_PAD]
        let seq_arrays = make_seq_array(&[&[1i16, -4, -10, -20]]);
        let qual_arrays = make_qual_array(&[&[0.05f32, 0.05, 0.05, -10.0]]);
        let read_spans = make_spans_array(&[[100i32, 103]]);

        let result = assemble_consensus(&seq_arrays, &qual_arrays, &read_spans);
        // non_na_values: count(v != HAP_PAD) = 3 (1, -4, -10) → takes first 3
        assert_eq!(result.to_vec(), vec![1i16, -4, -10]);
    }

    /// Golden: the only NA sentinel is HAP_PAD (-20); a leading deletion (-10) is data.
    #[test]
    fn test_only_hap_pad_is_padding() {
        let seq_arrays = make_seq_array(&[&[-10i16, 1, -20]]);
        let qual_arrays = make_qual_array(&[&[0.05f32, 0.05, -10.0]]);
        let read_spans = make_spans_array(&[[100i32, 102]]);

        let result = assemble_consensus(&seq_arrays, &qual_arrays, &read_spans);
        // -10 != HAP_PAD → real, so non_na_values = 2. Takes first 2: [-10, 1]
        assert_eq!(result.to_vec(), vec![-10i16, 1]);
    }

    // ── Integration: extract vectors → assemble consensus ──────────────────

    #[test]
    fn test_extract_then_assemble() {
        // Two overlapping reads:
        // Read 1: pos 100, CIGAR 5=, qual [10,10,10,10,10] → prob ~0.1
        // Read 2: pos 102, CIGAR 2=1X2=, qual [20,20,20,20,20] → prob ~0.01 (better)
        let cigar1 = CigarString(vec![Cigar::Equal(5)]);
        let cigar2 = CigarString(vec![
            Cigar::Equal(2),
            Cigar::Diff(1),
            Cigar::Equal(2),
        ]);
        let r1 = make_record(cigar1, b"ACGTG", &[10; 5], 100);
        let r2 = make_record(cigar2, b"GCACC", &[20; 5], 102);

        let hap1 = extract_hap_vector(&r1).unwrap();
        let err1 = extract_error_vector(&r1).unwrap();
        let hap2 = extract_hap_vector(&r2).unwrap();
        let err2 = extract_error_vector(&r2).unwrap();

        assert_eq!(hap1.to_vec(), vec![1, 1, 1, 1, 1]);
        assert_eq!(hap2.to_vec(), vec![1, 1, -4, 1, 1]);

        // Pack into padded Array2 for assemble_consensus
        let seq_arrays = make_padded_seq_array(&[hap1.to_vec(), hap2.to_vec()]);
        let qual_arrays = make_padded_qual_array(&[err1.to_vec(), err2.to_vec()]);
        let read_spans = make_spans_array(&[[100i32, 105], [102i32, 107]]);

        let consensus = assemble_consensus(&seq_arrays, &qual_arrays, &read_spans);
        // Span: 100..107 = 7 positions
        // pos 100 (rel 0): R1 only → 1
        // pos 101 (rel 1): R1 only → 1
        // pos 102 (rel 2): R1(1, p10~0.1) then R2(1, p20~0.01) → R2 wins → 1
        // pos 103 (rel 3): R1(1, p10~0.1) then R2(1, p20~0.01) → R2 wins → 1
        // pos 104 (rel 4): R1(1, p10~0.1) then R2(-4, p20~0.01) → R2 wins → -4
        // pos 105 (rel 5): R2 only → 1
        // pos 106 (rel 6): R2 only → 1
        assert_eq!(consensus.len(), 7);
        assert_eq!(consensus[0], 1);
        assert_eq!(consensus[1], 1);
        assert_eq!(consensus[4], -4); // SNV from read 2 wins
        assert_eq!(consensus[5], 1);
        assert_eq!(consensus[6], 1);
    }

    // ── Test record_hap_err_vectors_per_region ─────────────────────────────

    #[test]
    fn test_record_vectors_per_region_basic() {
        let cigar = CigarString(vec![Cigar::Equal(4)]);
        let r1 = make_named_record(b"read1", cigar.clone(), b"ACGT", &[30; 4], 100);
        let r2 = make_named_record(b"read2", cigar, b"TGCA", &[20; 4], 200);

        let records: Vec<&Record> = vec![&r1, &r2];
        let mut hap_cache: HashMap<String, Array1<i16>> = HashMap::new();
        let mut err_cache: HashMap<String, Array1<f32>> = HashMap::new();

        let (spans, haps, errs) =
            record_hap_err_vectors_per_region(&records, &mut hap_cache, &mut err_cache).unwrap();

        // spans is Array2<i32> (2 × 2)
        assert_eq!(spans.nrows(), 2);
        assert_eq!(spans[[0, 0]], 100);
        assert_eq!(spans[[0, 1]], 104);
        assert_eq!(spans[[1, 0]], 200);
        assert_eq!(spans[[1, 1]], 204);

        // haps is Array2<i16> (2 × 4), all matches → all 1s
        assert_eq!(haps.nrows(), 2);
        assert_eq!(haps.ncols(), 4);
        assert_eq!(haps.row(0).to_vec(), vec![1i16; 4]);
        assert_eq!(haps.row(1).to_vec(), vec![1i16; 4]);

        // errs is Array2<f32> (2 × 4)
        assert_eq!(errs.nrows(), 2);
        assert_eq!(errs.ncols(), 4);

        // Caches should be populated
        assert_eq!(hap_cache.len(), 2);
        assert_eq!(err_cache.len(), 2);
    }

    #[test]
    fn test_record_vectors_caching() {
        let cigar = CigarString(vec![Cigar::Equal(3)]);
        let r1 = make_record(cigar, b"ACG", &[30; 3], 100);

        let mut hap_cache: HashMap<String, Array1<i16>> = HashMap::new();
        let mut err_cache: HashMap<String, Array1<f32>> = HashMap::new();

        // First call populates cache
        let records: Vec<&Record> = vec![&r1];
        let (_, haps1, _) =
            record_hap_err_vectors_per_region(&records, &mut hap_cache, &mut err_cache).unwrap();

        // Second call with same record should use cached values
        let (_, haps2, _) =
            record_hap_err_vectors_per_region(&records, &mut hap_cache, &mut err_cache).unwrap();

        assert_eq!(haps1, haps2);
        // Cache still has exactly 1 entry (not duplicated)
        assert_eq!(hap_cache.len(), 1);
    }

    /// Call at the start of every test to enable `RUST_LOG=debug` output.
    fn init_log() {
        let _ = env_logger::try_init();
    }

    // ========================================================================
    // Tests for Variant Density Filtering Functions
    // ========================================================================

    #[test]
    fn test_extract_true_stretches_empty() {
        let bool_arr = array![false, false, false];
        let stretches = extract_true_stretches(&bool_arr);
        assert_eq!(stretches.len(), 0);
    }

    #[test]
    fn test_extract_true_stretches_single_stretch() {
        let bool_arr = array![false, true, true, true, false];
        let stretches = extract_true_stretches(&bool_arr);
        assert_eq!(stretches, vec![(1, 3)]);
    }

    #[test]
    fn test_extract_true_stretches_multiple() {
        let bool_arr = array![true, true, false, true, false, true, true, true];
        let stretches = extract_true_stretches(&bool_arr);
        assert_eq!(stretches, vec![(0, 1), (3, 3), (5, 7)]);
    }

    #[test]
    fn test_extract_true_stretches_ends_with_true() {
        let bool_arr = array![false, true, true];
        let stretches = extract_true_stretches(&bool_arr);
        assert_eq!(stretches, vec![(1, 2)]);
    }

    #[test]
    fn test_extract_true_stretches_starts_with_true() {
        let bool_arr = array![true, true, false];
        let stretches = extract_true_stretches(&bool_arr);
        assert_eq!(stretches, vec![(0, 1)]);
    }

    #[test]
    fn test_extract_true_stretches_all_true() {
        let bool_arr = array![true, true, true, true];
        let stretches = extract_true_stretches(&bool_arr);
        assert_eq!(stretches, vec![(0, 3)]);
    }

    #[test]
    fn test_count_window_var_density_no_variants() {
        // All matches (value = 1), no variants
        let seq = array![1, 1, 1, 1, 1];
        let density = count_window_var_density(&seq, 1);
        assert_eq!(density.len(), 5);
        assert!(density.iter().all(|&d| d == 0.0));
    }

    #[test]
    fn test_count_window_var_density_basic() {
        // Sequence: match, SNV, match, SNV, match
        // Window size = 1 (padding=0 would be window_size=1, but let's use padding=1 for window_size=3)
        let seq = array![1, -4, 1, -4, 1];
        let density = count_window_var_density(&seq, 1);

        // Window size = 2*1 + 1 = 3
        // Position 0: window [0:2] = [1, -4, 1] → 1 SNV → density = 1/3
        // Position 1: window [0:3] = [1, -4, 1] → 1 SNV → density = 1/3
        // Position 2: window [1:4] = [-4, 1, -4] → 2 SNVs → density = 2/3
        // Position 3: window [2:5] = [1, -4, 1] → 1 SNV → density = 1/3
        // Position 4: window [3:5] = [-4, 1] → 1 SNV → density = 1/3

        assert_eq!(density.len(), 5);
        assert!((density[0] - 1.0/3.0).abs() < 0.01);
        assert!((density[2] - 2.0/3.0).abs() < 0.01);
    }

    #[test]
    fn test_count_window_var_density_with_indels() {
        // Sequence with deletion block (golden deletion = -10)
        let seq = array![1, 1, -10, -10, 1, 1];
        let density = count_window_var_density(&seq, 1);

        // Window size = 3
        // The deletion block counts as 1 continuous indel block
        assert_eq!(density.len(), 6);
        // Middle positions should have higher density
        assert!(density[2] > 0.0);
        assert!(density[3] > 0.0);
    }

    /// Straightforward (slow) version: re-slice and re-count each window.
    /// This is the original implementation, kept here only as a reference to
    /// cross-check the optimized prefix-sum `count_window_var_density` against.
    fn count_window_var_density_reference(array: &Array1<i16>, padding_size: i32) -> Array1<f32> {
        let n = array.len();
        let has_variants = array.iter().any(|&v| v != 1);
        if !has_variants {
            return Array1::zeros(n);
        }
        let mut density_arr = Array1::<f32>::zeros(n);
        let window_size = (padding_size * 2 + 1) as f32;
        for i in 0..n {
            let start = i.saturating_sub(padding_size as usize);
            let end = (i + padding_size as usize + 1).min(n);
            let window = array.slice(ndarray::s![start..end]);
            let var_count = count_var(&window.to_owned()) as f32;
            density_arr[i] = var_count / window_size;
        }
        density_arr
    }

    /// Deterministic xorshift used by the cross-check / timing tests.
    fn xorshift(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    #[test]
    fn test_count_window_var_density_matches_reference() {
        init_log();
        // Mix of matches (1), SNVs (-4), deletions (-10), pure insertions (11/21),
        // and compound mismatch+ins (6/16 — counted as BOTH an SNV and an indel),
        // including adjacent mixed indel types so runs straddle window edges. Golden
        // optimized vs. reference must agree even on the double-counted compounds.
        let values = [1i16, 1, 1, -4, -10, -10, 21, 1, -4, 11, -10, 1, 6, -4, -4, 16];
        let mut seed = 0x1234_5678u64;
        for &len in &[0usize, 1, 2, 5, 13, 50, 137, 300] {
            let arr: Array1<i16> = (0..len)
                .map(|_| values[(xorshift(&mut seed) as usize) % values.len()])
                .collect::<Vec<_>>()
                .into();
            for &pad in &[0i32, 1, 2, 5, 42, 65, 74] {
                let got = count_window_var_density(&arr, pad);
                let want = count_window_var_density_reference(&arr, pad);
                assert_eq!(got.len(), want.len(), "len mismatch len={len} pad={pad}");
                for k in 0..got.len() {
                    assert!(
                        (got[k] - want[k]).abs() < 1e-6,
                        "mismatch len={len} pad={pad} idx={k}: got={} want={}",
                        got[k], want[k]
                    );
                }
            }
        }
    }

    #[test]
    #[ignore = "timing benchmark; run with: cargo test --lib -- --ignored --nocapture bench_count_window_var_density"]
    fn bench_count_window_var_density() {
        use std::time::Instant;
        let len = 2000usize;
        let values = [1i16, 1, 1, 1, 1, -4, -10, 21, 1, 1];
        let mut seed = 0x00AB_CDEFu64;
        let arr: Array1<i16> = (0..len)
            .map(|_| values[(xorshift(&mut seed) as usize) % values.len()])
            .collect::<Vec<_>>()
            .into();

        let iters = 200;
        let pads = [42i32, 65, 74];
        let mut sink = 0.0f32;

        let t0 = Instant::now();
        for _ in 0..iters {
            for &pad in &pads {
                sink += count_window_var_density_reference(&arr, pad).sum();
            }
        }
        let old = t0.elapsed();

        let t1 = Instant::now();
        for _ in 0..iters {
            for &pad in &pads {
                sink += count_window_var_density(&arr, pad).sum();
            }
        }
        let new = t1.elapsed();

        eprintln!(
            "count_window_var_density  len={len} x{iters} iters x{} pads:",
            pads.len()
        );
        eprintln!("  reference (re-count each window): {old:?}");
        eprintln!("  prefix-sum (new):                 {new:?}");
        eprintln!(
            "  speedup: {:.1}x  (sink={sink})",
            old.as_secs_f64() / new.as_secs_f64().max(1e-9)
        );
        assert!(new < old, "expected prefix-sum to be faster (old={old:?}, new={new:?})");
    }

    #[test]
    fn test_judge_misalignment_no_variants() {
        // All matches, should not be misaligned
        let seq = array![1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
        let (is_misaligned, max_density) = judge_misalignment_by_extreme_vardensity(&seq);
        assert!(!is_misaligned);
        assert_eq!(max_density, 0.0);
    }

    #[test]
    fn test_judge_misalignment_low_density() {
        // Few variants, low density, should not be misaligned
        let mut seq = Array1::<i16>::ones(200);
        seq[50] = -4; // One SNV
        seq[100] = -4; // Another SNV

        let (is_misaligned, max_density) = judge_misalignment_by_extreme_vardensity(&seq);
        assert!(!is_misaligned);
        assert!(max_density < 0.05); // Should be very low density
    }

    #[test]
    fn test_judge_misalignment_high_density_with_indels() {
        // Create a sequence with high variant density and indels
        // This should trigger misalignment detection
        let mut seq = Array1::<i16>::ones(200);

        // Add a cluster of variants with indels (positions 90-110)
        for i in 90..95 {
            seq[i] = -4; // SNVs
        }
        for i in 95..100 {
            seq[i] = -10; // Deletion block (golden HAP_DEL)
        }
        for i in 100..105 {
            seq[i] = -4; // More SNVs
        }

        let (is_misaligned, max_density) = judge_misalignment_by_extreme_vardensity(&seq);
        // With this many variants clustered together, should be detected as misaligned
        assert!(is_misaligned);
        assert!(max_density > 0.05);
    }

    #[test]
    fn test_judge_misalignment_threshold_3() {
        // Test the third threshold (10/148 ≈ 0.0676)
        // Create a sequence with high density but no indels
        let mut seq = Array1::<i16>::ones(200);

        // Add many SNVs in a cluster to exceed threshold 3
        for i in 90..110 {
            seq[i] = -4; // 20 SNVs in a row
        }

        let (is_misaligned, max_density) = judge_misalignment_by_extreme_vardensity(&seq);
        // Should be detected by threshold 3 (high density regardless of indels)
        assert!(is_misaligned);
        assert!(max_density > 0.06);
    }

    // ========================================================================
    // Tests for rank_unique_values — Python cross-validated
    // Reference: /paedyl01/disk1/yangyxt/test_tmp/cross_validate_helpers.py
    // ========================================================================

    #[test]
    fn test_rank_unique_values_tc1_realistic_sim_scores() {
        // TC1: Realistic hap_max_sim_scores with duplicates, negatives, zero
        // Python: rank_unique_values([3.2, 1.5, 3.2, -1.0, 7.0, 1.5, 0.0, 7.0, -1.0, 3.2])
        //       → [4, 3, 4, 1, 5, 3, 2, 5, 1, 4]
        let arr = vec![3.2f32, 1.5, 3.2, -1.0, 7.0, 1.5, 0.0, 7.0, -1.0, 3.2];
        let ranks = rank_unique_values(&arr);
        assert_eq!(ranks, vec![4, 3, 4, 1, 5, 3, 2, 5, 1, 4]);
    }

    #[test]
    fn test_rank_unique_values_tc2_many_duplicates() {
        // TC2: Many duplicates — simulates haplotypes with same similarity
        // Python: rank_unique_values([2.0, 2.0, 5.0, 5.0, 5.0, 1.0, 3.0, 3.0])
        //       → [2, 2, 4, 4, 4, 1, 3, 3]
        let arr = vec![2.0f32, 2.0, 5.0, 5.0, 5.0, 1.0, 3.0, 3.0];
        let ranks = rank_unique_values(&arr);
        assert_eq!(ranks, vec![2, 2, 4, 4, 4, 1, 3, 3]);
    }

    #[test]
    fn test_rank_unique_values_tc3_all_unique_descending() {
        // TC3: All unique, descending — worst case for ranking
        // Python: rank_unique_values([10.0, 8.0, 6.0, 4.0, 2.0, 0.0, -2.0])
        //       → [7, 6, 5, 4, 3, 2, 1]
        let arr = vec![10.0f32, 8.0, 6.0, 4.0, 2.0, 0.0, -2.0];
        let ranks = rank_unique_values(&arr);
        assert_eq!(ranks, vec![7, 6, 5, 4, 3, 2, 1]);
    }

    #[test]
    fn test_rank_unique_values_tc4_single_element() {
        // TC4: Single element
        // Python: rank_unique_values([42.5]) → [1]
        let arr = vec![42.5f32];
        let ranks = rank_unique_values(&arr);
        assert_eq!(ranks, vec![1]);
    }

    #[test]
    fn test_rank_unique_values_tc5_rounded_precision() {
        // TC5: Rounded to .1 precision (real pipeline does .round(1) before calling)
        // Python: rank_unique_values([4.15, 4.15, 3.80, -1.00, 7.00, 0.42, 0.42])
        //       → [4, 4, 3, 1, 5, 2, 2]
        let arr = vec![4.15f32, 4.15, 3.80, -1.00, 7.00, 0.42, 0.42];
        let ranks = rank_unique_values(&arr);
        assert_eq!(ranks, vec![4, 4, 3, 1, 5, 2, 2]);
    }

    #[test]
    fn test_rank_unique_values_tc6_negatives_and_zero() {
        // TC6: Negative values and zero
        // Python: rank_unique_values([-3.5, 0.0, -3.5, 2.1, 0.0, -1.0, 2.1])
        //       → [1, 3, 1, 4, 3, 2, 4]
        let arr = vec![-3.5f32, 0.0, -3.5, 2.1, 0.0, -1.0, 2.1];
        let ranks = rank_unique_values(&arr);
        assert_eq!(ranks, vec![1, 3, 1, 4, 3, 2, 4]);
    }

    // ========================================================================
    // Tests for calculate_coefficient — Python cross-validated
    // Reference: /paedyl01/disk1/yangyxt/test_tmp/cross_validate_helpers.py
    // ========================================================================
    // columns: [start, end, total_depth, hap_id, hap_depth,
    //           var_count, indel_count, psv_count, varc_rank, hap_max_sim_scores]

    #[test]
    fn test_calculate_coefficient_tc1_4hap_realistic() {
        // TC1: Realistic 4-haplotype region from HG002 chr1:1633000-1635000
        // Python output:
        //   hap 1: 19.7995929718
        //   hap 2:  1.0844353437
        //   hap 3: -6.3245558739
        //   hap 4:  0.0000000000  (depth == total_depth)
        let r1: [f32; 10] = [1633000.0, 1635000.0, 45.0, 1.0, 12.0, 8.0, 2.0, 6.0, 3.0, 5.17];
        let r2: [f32; 10] = [1633000.0, 1635000.0, 45.0, 2.0, 30.0, 15.0, 5.0, 3.0, 1.0, 0.42];
        let r3: [f32; 10] = [1633000.0, 1635000.0, 45.0, 3.0, 5.0, 3.0, 1.0, 0.0, 2.0, -1.50];
        let r4: [f32; 10] = [1633000.0, 1635000.0, 45.0, 4.0, 45.0, 0.0, 0.0, 0.0, 4.0, 7.00];
        let result = calculate_coefficient(&[&r1, &r2, &r3, &r4]);

        assert!((result[0] - 19.799_593).abs() < 0.001, "hap1: got {}", result[0]);
        assert!((result[1] - 1.084_435_3).abs() < 0.001, "hap2: got {}", result[1]);
        assert!((result[2] - (-6.324_556)).abs() < 0.001, "hap3: got {}", result[2]);
        assert!((result[3] - 0.0).abs() < 0.001, "hap4: got {}", result[3]);
    }

    #[test]
    fn test_calculate_coefficient_tc2_span_effect() {
        // TC2: Same psv_metric=3.0, same depth, different spans
        // Tests sqrt(span) scaling
        // Python output:
        //   span=100:   2.5980761051
        //   span=500:   5.8094754219
        //   span=1000:  8.2158384323
        //   span=10000: 25.9807624817
        let r1: [f32; 10] = [0.0, 100.0, 20.0, 1.0, 5.0, 3.0, 1.0, 2.0, 1.0, 3.0];
        let r2: [f32; 10] = [0.0, 500.0, 20.0, 2.0, 5.0, 3.0, 1.0, 2.0, 1.0, 3.0];
        let r3: [f32; 10] = [0.0, 1000.0, 20.0, 3.0, 5.0, 3.0, 1.0, 2.0, 1.0, 3.0];
        let r4: [f32; 10] = [0.0, 10000.0, 20.0, 4.0, 5.0, 3.0, 1.0, 2.0, 1.0, 3.0];
        let result = calculate_coefficient(&[&r1, &r2, &r3, &r4]);

        assert!((result[0] - 2.598_076).abs() < 0.001, "span100: got {}", result[0]);
        assert!((result[1] - 5.809_475_4).abs() < 0.001, "span500: got {}", result[1]);
        assert!((result[2] - 8.215_838).abs() < 0.001, "span1000: got {}", result[2]);
        assert!((result[3] - 25.980_762).abs() < 0.001, "span10000: got {}", result[3]);
    }

    #[test]
    fn test_calculate_coefficient_tc3_depth_frac_effect() {
        // TC3: Same span=400, psv_metric=4.0, different depth ratios
        // Tests sqrt(depth_frac) scaling
        // Python output:
        //   depth_frac=0.90: 7.5894660950
        //   depth_frac=0.75: 6.9282031059
        //   depth_frac=0.50: 5.6568541527
        //   depth_frac=0.10: 2.5298223495
        let r1: [f32; 10] = [100.0, 500.0, 100.0, 1.0, 10.0, 5.0, 1.0, 3.0, 1.0, 4.0];
        let r2: [f32; 10] = [100.0, 500.0, 100.0, 2.0, 25.0, 5.0, 1.0, 3.0, 1.0, 4.0];
        let r3: [f32; 10] = [100.0, 500.0, 100.0, 3.0, 50.0, 5.0, 1.0, 3.0, 1.0, 4.0];
        let r4: [f32; 10] = [100.0, 500.0, 100.0, 4.0, 90.0, 5.0, 1.0, 3.0, 1.0, 4.0];
        let result = calculate_coefficient(&[&r1, &r2, &r3, &r4]);

        assert!((result[0] - 7.589_466).abs() < 0.001, "frac0.90: got {}", result[0]);
        assert!((result[1] - 6.928_203).abs() < 0.001, "frac0.75: got {}", result[1]);
        assert!((result[2] - 5.656_854).abs() < 0.001, "frac0.50: got {}", result[2]);
        assert!((result[3] - 2.529_822_3).abs() < 0.001, "frac0.10: got {}", result[3]);
    }

    #[test]
    fn test_calculate_coefficient_tc4_mixed_psv_metric() {
        // TC4: Mixed psv_metric values: zero, negative, large positive
        // Python output:
        //   psv_metric=0.00:  0.0000000000
        //   psv_metric=-3.17: -6.9451222420
        //   psv_metric=12.50: 27.3861293793
        let r1: [f32; 10] = [200.0, 800.0, 30.0, 1.0, 6.0, 4.0, 2.0, 2.0, 1.0, 0.0];
        let r2: [f32; 10] = [200.0, 800.0, 30.0, 2.0, 6.0, 4.0, 2.0, 2.0, 2.0, -3.17];
        let r3: [f32; 10] = [200.0, 800.0, 30.0, 3.0, 6.0, 4.0, 2.0, 2.0, 3.0, 12.50];
        let result = calculate_coefficient(&[&r1, &r2, &r3]);

        assert!((result[0] - 0.0).abs() < 0.001, "zero: got {}", result[0]);
        assert!((result[1] - (-6.945_122_2)).abs() < 0.001, "neg: got {}", result[1]);
        assert!((result[2] - 27.386_13).abs() < 0.001, "pos: got {}", result[2]);
    }

    #[test]
    fn test_calculate_coefficient_tc5_single_row() {
        // TC5: Single row edge case
        // Python output: 16.1116104126
        let r1: [f32; 10] = [500000.0, 502000.0, 88.0, 7.0, 22.0, 11.0, 3.0, 8.0, 1.0, 4.16];
        let result = calculate_coefficient(&[&r1]);
        assert!((result[0] - 16.111_61).abs() < 0.01, "single: got {}", result[0]);
    }

    #[test]
    fn test_calculate_coefficient_tc6_multi_region_8rows() {
        // TC6: 8 haplotypes across 3 regions — realistic pipeline table
        // Python output:
        //   row 0 (region1 hap1): 19.7995929718
        //   row 1 (region1 hap2): 13.1635856628
        //   row 2 (region1 hap3): -4.0551748276
        //   row 3 (region1 hap4):  1.7260359526
        //   row 4 (region2 hap1): 13.9530639648
        //   row 5 (region2 hap2):  7.6685075760
        //   row 6 (region3 hap1): 33.2039146423
        //   row 7 (region3 hap5): -16.6495361328
        let rows: Vec<[f32; 10]> = vec![
            [1633000.0, 1635000.0, 45.0, 1.0, 12.0,  8.0, 2.0, 6.0, 4.0,  5.17],
            [1633000.0, 1635000.0, 45.0, 2.0, 18.0, 12.0, 4.0, 3.0, 2.0,  3.80],
            [1633000.0, 1635000.0, 45.0, 3.0,  8.0,  3.0, 1.0, 1.0, 3.0, -1.00],
            [1633000.0, 1635000.0, 45.0, 4.0,  7.0,  5.0, 2.0, 0.0, 1.0,  0.42],
            [1700000.0, 1701500.0, 32.0, 1.0,  8.0,  5.0, 1.0, 4.0, 3.0,  4.16],
            [1700000.0, 1701500.0, 32.0, 2.0, 14.0,  9.0, 3.0, 2.0, 1.0,  2.64],
            [1800000.0, 1803000.0, 60.0, 1.0, 15.0,  7.0, 2.0, 5.0, 2.0,  7.00],
            [1800000.0, 1803000.0, 60.0, 5.0, 25.0, 18.0, 6.0, 1.0, 1.0, -3.98],
        ];
        let row_refs: Vec<&[f32; 10]> = rows.iter().collect();
        let result = calculate_coefficient(&row_refs);

        let expected = [19.7995929718, 13.1635856628, -4.0551748276, 1.7260359526,
            13.9530639648, 7.6685075760, 33.2039146423, -16.6495361328];

        for (i, (&got, &exp)) in result.iter().zip(expected.iter()).enumerate() {
            assert!(
                (got - exp as f32).abs() < 0.01,
                "row {i}: got {got} expected {exp}"
            );
        }
    }

    // ========================================================================
    // Tests for Reference Similarity Helper Functions
    // ========================================================================

    // --- merge_unique_sorted ---
    #[test]
    fn test_merge_unique_sorted_both_empty() {
        init_log();
        let result = merge_unique_sorted(&[], &[]);
        debug!("merge_unique_sorted([], []) = {result:?}");
        assert_eq!(result, Vec::<i32>::new());
    }

    #[test]
    fn test_merge_unique_sorted_one_empty() {
        init_log();
        let r1 = merge_unique_sorted(&[1, 3, 5], &[]);
        let r2 = merge_unique_sorted(&[], &[2, 4]);
        debug!("merge_unique_sorted([1,3,5], []) = {r1:?}");
        debug!("merge_unique_sorted([], [2,4]) = {r2:?}");
        assert_eq!(r1, vec![1, 3, 5]);
        assert_eq!(r2, vec![2, 4]);
    }

    #[test]
    fn test_merge_unique_sorted_no_overlap() {
        init_log();
        let result = merge_unique_sorted(&[1, 3, 5], &[2, 4, 6]);
        debug!("merge_unique_sorted([1,3,5], [2,4,6]) = {result:?}");
        assert_eq!(result, vec![1, 2, 3, 4, 5, 6]);
    }

    #[test]
    fn test_merge_unique_sorted_with_duplicates() {
        init_log();
        let result = merge_unique_sorted(&[1, 2, 3], &[2, 3, 4]);
        debug!("merge_unique_sorted([1,2,3], [2,3,4]) = {result:?}");
        assert_eq!(result, vec![1, 2, 3, 4]);
    }

    #[test]
    fn test_merge_unique_sorted_identical() {
        init_log();
        let result = merge_unique_sorted(&[5, 5, 5], &[5, 5]);
        debug!("merge_unique_sorted([5,5,5], [5,5]) = {result:?}");
        assert_eq!(result, vec![5]);
    }

    // --- ref_genome_similarity ---
    #[test]
    fn test_ref_genome_similarity_all_reference() {
        init_log();
        // Genomic hap vector is all matches (1s) → short-circuit to (0,0,0)
        let query = array![1, -4, 1, -10, -10, 1];
        let genomic = array![1, 1, 1, 1, 1, 1];
        let result = ref_genome_similarity(&query, &genomic);
        debug!("ref_genome_similarity(query={:?}, genomic={:?}) = {:?}", query.as_slice().unwrap(), genomic.as_slice().unwrap(), result);
        assert_eq!(result, (0, 0, 0));
    }

    #[test]
    fn test_ref_genome_similarity_basic() {
        init_log();
        // query: 2 SNVs + 1 del block = 3 variants
        // genomic: 1 SNV + 1 del block → alt_snv=1, alt_indel=1
        let query = array![1, -4, -4, 1, -10, -10, 1];
        let genomic = array![1, -4, 1, 1, -10, -10, 1];
        let (var_count, alt_snv, alt_indel) = ref_genome_similarity(&query, &genomic);
        debug!("ref_genome_similarity(query={:?}, genomic={:?}) = (var_count={}, alt_snv={}, alt_indel={})",
               query.as_slice().unwrap(), genomic.as_slice().unwrap(), var_count, alt_snv, alt_indel);
        assert_eq!(var_count, 3); // 2 SNVs + 1 indel block
        assert_eq!(alt_snv, 1);
        assert_eq!(alt_indel, 1);
    }

    #[test]
    fn test_ref_genome_similarity_only_insertions() {
        init_log();
        // query has insertion markers; genomic has insertion markers (golden: base+10*L)
        let query = array![1, 21, 1, 1];  // 1 insertion block (match + 2bp ins)
        let genomic = array![1, 1, 31, 1]; // 1 insertion block (match + 3bp ins)
        let (var_count, alt_snv, alt_indel) = ref_genome_similarity(&query, &genomic);
        debug!("ref_genome_similarity(query={:?}, genomic={:?}) = (var_count={}, alt_snv={}, alt_indel={})",
               query.as_slice().unwrap(), genomic.as_slice().unwrap(), var_count, alt_snv, alt_indel);
        assert_eq!(var_count, 1);   // 1 indel block in query
        assert_eq!(alt_snv, 0);
        assert_eq!(alt_indel, 1);   // 1 indel block in genomic
    }

    // --- numba_shared_variant_positions ---
    #[test]
    fn test_shared_variant_positions_no_shared() {
        init_log();
        // vec1 has SNV at pos 1, vec2 has SNV at pos 3
        let vec1 = array![1, -4, 1, 1, 1];
        let vec2 = array![1, 1, 1, -4, 1];
        let (snvs, ins, dels) = numba_shared_variant_positions(&vec1, &vec2, 100);
        debug!("shared_variant_positions(vec1={:?}, vec2={:?}, overlap_start=100) = (snvs={:?}, ins={:?}, dels={:?})",
               vec1.as_slice().unwrap(), vec2.as_slice().unwrap(), snvs, ins, dels);
        assert!(snvs.is_empty());
        assert!(ins.is_empty());
        assert!(dels.is_empty());
    }

    #[test]
    fn test_shared_variant_positions_shared_snvs() {
        init_log();
        let vec1 = array![1, -4, 1, -4, 1];
        let vec2 = array![1, -4, -4, -4, 1];
        let (snvs, ins, dels) = numba_shared_variant_positions(&vec1, &vec2, 1000);
        debug!("shared_variant_positions(vec1={:?}, vec2={:?}, overlap_start=1000) = (snvs={:?}, ins={:?}, dels={:?})",
               vec1.as_slice().unwrap(), vec2.as_slice().unwrap(), snvs, ins, dels);
        // Shared SNVs at relative indices 1 and 3 → absolute 1001, 1003
        assert_eq!(snvs, vec![1001, 1003]);
        assert!(ins.is_empty());
        assert!(dels.is_empty());
    }

    #[test]
    fn test_shared_variant_positions_shared_insertions() {
        init_log();
        // Golden: shared insertion = same inserted LENGTH (not raw value).
        // 21 = match+2bp, 31 = match+3bp. idx 1 both length 2 (shared); idx 3 length 3 vs 2 (not).
        let vec1 = array![1, 21, 1, 31, 1];
        let vec2 = array![1, 21, 1, 21, 1];
        let (snvs, ins, dels) = numba_shared_variant_positions(&vec1, &vec2, 0);
        debug!("shared_variant_positions(vec1={:?}, vec2={:?}, overlap_start=0) = (snvs={:?}, ins={:?}, dels={:?})",
               vec1.as_slice().unwrap(), vec2.as_slice().unwrap(), snvs, ins, dels);
        assert!(snvs.is_empty());
        assert_eq!(ins, vec![1]); // only idx 1 matches
        assert!(dels.is_empty());
    }

    #[test]
    fn test_shared_variant_positions_shared_deletions() {
        init_log();
        // Two identical deletion spans (golden -10): indices 2..=4
        let vec1 = array![1, 1, -10, -10, -10, 1, 1];
        let vec2 = array![1, 1, -10, -10, -10, 1, 1];
        let (snvs, ins, dels) = numba_shared_variant_positions(&vec1, &vec2, 50);
        debug!("shared_variant_positions(vec1={:?}, vec2={:?}, overlap_start=50) = (snvs={:?}, ins={:?}, dels={:?})",
               vec1.as_slice().unwrap(), vec2.as_slice().unwrap(), snvs, ins, dels);
        assert!(snvs.is_empty());
        assert!(ins.is_empty());
        assert_eq!(dels, vec![52, 53, 54]); // abs positions 50+2, 50+3, 50+4
    }

    #[test]
    fn test_shared_variant_positions_different_deletion_spans() {
        init_log();
        // vec1 has del at 1..=2, vec2 has del at 1..=3 → NOT shared (different span)
        let vec1 = array![1, -10, -10, 1, 1];
        let vec2 = array![1, -10, -10, -10, 1];
        let (_snvs, _ins, dels) = numba_shared_variant_positions(&vec1, &vec2, 0);
        debug!("shared_variant_positions(vec1={:?}, vec2={:?}, overlap_start=0) = dels={:?} (expect empty, different spans)",
               vec1.as_slice().unwrap(), vec2.as_slice().unwrap(), dels);
        assert!(dels.is_empty());
    }

    #[test]
    fn test_shared_variant_positions_mixed() {
        init_log();
        // SNV at idx 0, insertion at idx 2 (both 21 = match+2bp), deletion at idx 4..=5 (-10)
        let vec1 = array![-4, 1, 21, 1, -10, -10, 1];
        let vec2 = array![-4, 1, 21, 1, -10, -10, 1];
        let (snvs, ins, dels) = numba_shared_variant_positions(&vec1, &vec2, 10);
        debug!("shared_variant_positions(vec1={:?}, vec2={:?}, overlap_start=10) = (snvs={:?}, ins={:?}, dels={:?})",
               vec1.as_slice().unwrap(), vec2.as_slice().unwrap(), snvs, ins, dels);
        assert_eq!(snvs, vec![10]);      // abs 10+0
        assert_eq!(ins, vec![12]);       // abs 10+2
        assert_eq!(dels, vec![14, 15]);  // abs 10+4, 10+5
    }

    #[test]
    fn test_shared_variant_positions_compound() {
        init_log();
        // Golden compound base (mismatch + 1bp ins = 6) shared at idx 1 counts as
        // BOTH a shared SNV and a shared insertion (length 1) — the recovered signal.
        let vec1 = array![1, 6, 1];
        let vec2 = array![1, 6, 1];
        let (snvs, ins, dels) = numba_shared_variant_positions(&vec1, &vec2, 100);
        debug!("shared_variant_positions(compound) = (snvs={snvs:?}, ins={ins:?}, dels={dels:?})");
        assert_eq!(snvs, vec![101]);
        assert_eq!(ins, vec![101]);
        assert!(dels.is_empty());
    }

    // --- map_positions_to_bases ---
    #[test]
    fn test_map_positions_to_bases_basic() {
        init_log();
        // Read starts at position 100, 5 ref positions mapping to query indices
        let ref_qseq_positions = vec![0, 1, 2, 3, 4]; // 1:1 mapping
        let qseq_encoded: Vec<i8> = vec![0, 1, 2, 3, 0]; // A, T, C, G, A
        let shared_pos = vec![100, 102, 104];
        let result = map_positions_to_bases(&shared_pos, 100, &ref_qseq_positions, &qseq_encoded);
        debug!("map_positions_to_bases(shared_pos={shared_pos:?}, read_start=100, ref_qseq={ref_qseq_positions:?}, qseq={qseq_encoded:?}) = {result:?}");
        assert_eq!(result, vec![0, 2, 0]); // A, C, A
    }

    #[test]
    fn test_map_positions_to_bases_out_of_range() {
        init_log();
        let ref_qseq_positions = vec![0, 1, 2];
        let qseq_encoded: Vec<i8> = vec![0, 1, 2];
        let shared_pos = vec![99, 100, 105]; // 99 is before read, 105 is after
        let result = map_positions_to_bases(&shared_pos, 100, &ref_qseq_positions, &qseq_encoded);
        debug!("map_positions_to_bases(shared_pos={shared_pos:?}, read_start=100, ref_qseq={ref_qseq_positions:?}, qseq={qseq_encoded:?}) = {result:?}");
        assert_eq!(result, vec![-1, 0, -1]); // only pos 100 maps
    }

    #[test]
    fn test_map_positions_to_bases_with_gaps() {
        init_log();
        // ref_qseq_positions[1] = -1 means deletion at that ref offset
        let ref_qseq_positions = vec![0, -1, 1];
        let qseq_encoded: Vec<i8> = vec![3, 2]; // G, C
        let shared_pos = vec![100, 101, 102];
        let result = map_positions_to_bases(&shared_pos, 100, &ref_qseq_positions, &qseq_encoded);
        debug!("map_positions_to_bases(shared_pos={shared_pos:?}, read_start=100, ref_qseq={ref_qseq_positions:?}, qseq={qseq_encoded:?}) = {result:?}");
        assert_eq!(result, vec![3, -1, 2]); // G, gap, C
    }

    // --- update_tally_for_read ---
    #[test]
    fn test_update_tally_basic() {
        init_log();
        // 3 shared positions, read covers all of them
        let shared_pos = vec![100, 101, 102];
        let ref_qseq_positions = vec![0, 1, 2];
        let qseq_encoded: Vec<i8> = vec![0, 1, 2]; // A, T, C
        let mut tally = vec![[0i32; 5]; 3];

        update_tally_for_read(&shared_pos, 100, &ref_qseq_positions, &qseq_encoded, &mut tally);

        debug!("update_tally(shared_pos={shared_pos:?}, read_start=100, ref_qseq={ref_qseq_positions:?}, qseq={qseq_encoded:?}) => tally={tally:?}");
        assert_eq!(tally[0][0], 1); // A at pos 0
        assert_eq!(tally[1][1], 1); // T at pos 1
        assert_eq!(tally[2][2], 1); // C at pos 2
    }

    #[test]
    fn test_update_tally_multiple_reads() {
        init_log();
        let shared_pos = vec![100, 101];
        let ref_qseq_pos = vec![0, 1];
        let mut tally = vec![[0i32; 5]; 2];

        // Read 1: A, T
        update_tally_for_read(&shared_pos, 100, &ref_qseq_pos, &[0, 1], &mut tally);
        debug!("after read1 (A,T): tally={tally:?}");
        // Read 2: A, C
        update_tally_for_read(&shared_pos, 100, &ref_qseq_pos, &[0, 2], &mut tally);
        debug!("after read2 (A,C): tally={tally:?}");
        // Read 3: G, T
        update_tally_for_read(&shared_pos, 100, &ref_qseq_pos, &[3, 1], &mut tally);
        debug!("after read3 (G,T): tally={tally:?}");

        assert_eq!(tally[0][0], 2); // 2 reads have A at pos 0
        assert_eq!(tally[0][3], 1); // 1 read has G at pos 0
        assert_eq!(tally[1][1], 2); // 2 reads have T at pos 1
        assert_eq!(tally[1][2], 1); // 1 read has C at pos 1
    }

    #[test]
    fn test_update_tally_skips_n_base() {
        init_log();
        let shared_pos = vec![100];
        let ref_qseq_pos = vec![0];
        let qseq_encoded: Vec<i8> = vec![4]; // N
        let mut tally = vec![[0i32; 5]; 1];

        update_tally_for_read(&shared_pos, 100, &ref_qseq_pos, &qseq_encoded, &mut tally);
        debug!("update_tally with N base: qseq={qseq_encoded:?}, tally={tally:?}");
        // N (4) should be excluded — no column updated
        assert_eq!(tally[0], [0, 0, 0, 0, 0]);
    }

    #[test]
    fn test_update_tally_out_of_range() {
        init_log();
        // shared position is outside the read
        let shared_pos = vec![200];
        let ref_qseq_pos = vec![0, 1]; // read covers only 2 ref positions
        let qseq_encoded: Vec<i8> = vec![0, 1];
        let mut tally = vec![[0i32; 5]; 1];

        update_tally_for_read(&shared_pos, 100, &ref_qseq_pos, &qseq_encoded, &mut tally);
        debug!("update_tally out-of-range: shared_pos={:?}, read_start=100, ref_len={}, tally={:?}",
               shared_pos, ref_qseq_pos.len(), tally);
        assert_eq!(tally[0], [0, 0, 0, 0, 0]); // nothing updated
    }

    // --- verify_shared_snv_positions ---
    #[test]
    fn test_verify_shared_snv_all_match() {
        init_log();
        // Tally: consensus = A(10), T(0), C(0), G(0) → argmax = A(0)
        // Homolog base = A(0) → match
        let tally = vec![[10, 0, 0, 0, 0]];
        let h_base: Vec<i8> = vec![0]; // A
        let shared_pos = vec![100];
        let result = verify_shared_snv_positions(&tally, &h_base, &shared_pos);
        debug!("verify_shared_snv(tally={tally:?}, h_base={h_base:?}, shared_pos={shared_pos:?}) = {result:?}");
        assert_eq!(result, vec![100]);
    }

    #[test]
    fn test_verify_shared_snv_no_match() {
        init_log();
        // Consensus = T (argmax of [0, 10, 0, 0, 0])
        // Homolog = G (3)
        let tally = vec![[0, 10, 0, 0, 0]];
        let h_base: Vec<i8> = vec![3]; // G
        let shared_pos = vec![100];
        let result = verify_shared_snv_positions(&tally, &h_base, &shared_pos);
        debug!("verify_shared_snv(tally={tally:?}, h_base={h_base:?}, shared_pos={shared_pos:?}) = {result:?} (expect empty, T!=G)");
        assert!(result.is_empty());
    }

    #[test]
    fn test_verify_shared_snv_zero_coverage() {
        init_log();
        // Zero coverage → skip
        let tally = vec![[0, 0, 0, 0, 0]];
        let h_base: Vec<i8> = vec![0];
        let shared_pos = vec![100];
        let result = verify_shared_snv_positions(&tally, &h_base, &shared_pos);
        debug!("verify_shared_snv zero-coverage: tally={tally:?}, result={result:?} (expect empty)");
        assert!(result.is_empty());
    }

    #[test]
    fn test_verify_shared_snv_homolog_not_covered() {
        init_log();
        // Homolog base = -1 (not covered) → skip
        let tally = vec![[10, 0, 0, 0, 0]];
        let h_base: Vec<i8> = vec![-1];
        let shared_pos = vec![100];
        let result = verify_shared_snv_positions(&tally, &h_base, &shared_pos);
        debug!("verify_shared_snv homolog_not_covered: h_base={h_base:?}, result={result:?} (expect empty)");
        assert!(result.is_empty());
    }

    #[test]
    fn test_verify_shared_snv_multiple_positions() {
        init_log();
        // 3 positions: match, no-match, match
        let tally = vec![
            [5, 1, 0, 0, 0],  // argmax=A(0)
            [0, 0, 8, 2, 0],  // argmax=C(2)
            [1, 0, 0, 7, 0],  // argmax=G(3)
        ];
        let h_base: Vec<i8> = vec![0, 3, 3]; // A, G, G
        let shared_pos = vec![100, 101, 102];
        let result = verify_shared_snv_positions(&tally, &h_base, &shared_pos);
        debug!("verify_shared_snv_multi: tally={tally:?}, h_base={h_base:?}, shared_pos={shared_pos:?} => {result:?}");
        // pos 100: cons A == hom A → verified
        // pos 101: cons C != hom G → not
        // pos 102: cons G == hom G → verified
        assert_eq!(result, vec![100, 102]);
    }

    // ─── parse_origin_region tests ────────────────────────────────────────

    #[test]
    fn test_parse_origin_region_basic() {
        init_log();
        let qname = "chr1:12345-67890";
        let result = parse_origin_region(qname);
        debug!("parse_origin_region({qname:?}) => {result:?}");
        assert_eq!(result, Some(("chr1".to_string(), 12345, 67890)));
    }

    #[test]
    fn test_parse_origin_region_with_rg_suffix() {
        init_log();
        // Python regex: r"([a-zA-Z0-9]+):(\d+)-(\d+)[:RG0-9]*"
        // So "chr1:100-200:RG001" should match "chr1", 100, 200
        let qname = "chr1:100-200:RG001";
        let result = parse_origin_region(qname);
        debug!("parse_origin_region({qname:?}) => {result:?}");
        assert_eq!(result, Some(("chr1".to_string(), 100, 200)));
    }

    #[test]
    fn test_parse_origin_region_no_match() {
        init_log();
        let qname = "READNAME_NO_COORDS";
        let result = parse_origin_region(qname);
        debug!("parse_origin_region({qname:?}) => {result:?}");
        assert_eq!(result, None);
    }

    #[test]
    fn test_parse_origin_region_no_digits_after_colon() {
        init_log();
        let qname = "chr1:abc-200";
        let result = parse_origin_region(qname);
        debug!("parse_origin_region({qname:?}) => {result:?}");
        assert_eq!(result, None);
    }

    #[test]
    fn test_parse_origin_region_no_dash() {
        init_log();
        let qname = "chr1:12345_67890";
        let result = parse_origin_region(qname);
        debug!("parse_origin_region({qname:?}) => {result:?}");
        assert_eq!(result, None);
    }

    #[test]
    fn test_parse_origin_region_chr_with_digits() {
        init_log();
        // Chromosome names like "chr22" should work
        let qname = "chr22:5000000-6000000";
        let result = parse_origin_region(qname);
        debug!("parse_origin_region({qname:?}) => {result:?}");
        assert_eq!(result, Some(("chr22".to_string(), 5000000, 6000000)));
    }

    #[test]
    fn test_parse_origin_region_embedded_in_longer_qname() {
        init_log();
        // qname might have prefix/suffix around the chrom:start-end
        let qname = "sample1_chr1:100-200_hap1";
        let result = parse_origin_region(qname);
        debug!("parse_origin_region({qname:?}) => {result:?}");
        // The first alphanumeric run before ':' that matches the pattern
        // "sample1" is separated by '_', so the chrom_start scan backwards
        // stops at '_'. alphanumeric run would be "chr1" (scanning back from ':')
        // Wait — actually "sample1_chr1:" — scanning back from ':', we get
        // "1" then "r" then "h" then "c" then stop at '_'. So chrom = "chr1".
        assert_eq!(result, Some(("chr1".to_string(), 100, 200)));
    }

    #[test]
    fn test_parse_origin_region_multiple_colons_first_match() {
        init_log();
        // "chr1:100-200:chr2:300-400" should match the first occurrence
        let qname = "chr1:100-200:chr2:300-400";
        let result = parse_origin_region(qname);
        debug!("parse_origin_region({qname:?}) => {result:?}");
        assert_eq!(result, Some(("chr1".to_string(), 100, 200)));
    }

    #[test]
    fn test_parse_origin_region_large_coordinates() {
        init_log();
        let qname = "chrX:100000000-200000000";
        let result = parse_origin_region(qname);
        debug!("parse_origin_region({qname:?}) => {result:?}");
        assert_eq!(result, Some(("chrX".to_string(), 100000000, 200000000)));
    }

    #[test]
    fn test_parse_origin_region_single_digit_coords() {
        init_log();
        let qname = "chr1:1-2";
        let result = parse_origin_region(qname);
        debug!("parse_origin_region({qname:?}) => {result:?}");
        assert_eq!(result, Some(("chr1".to_string(), 1, 2)));
    }

    #[test]
    fn test_parse_origin_region_empty_string() {
        init_log();
        let result = parse_origin_region("");
        debug!("parse_origin_region('') => {result:?}");
        assert_eq!(result, None);
    }

    #[test]
    fn test_parse_origin_region_colon_at_start() {
        init_log();
        // ':' at position 0 means no chrom before it
        let result = parse_origin_region(":100-200");
        debug!("parse_origin_region(':100-200') => {result:?}");
        assert_eq!(result, None);
    }

    #[test]
    fn test_parse_origin_region_missing_end_digits() {
        init_log();
        let result = parse_origin_region("chr1:100-");
        debug!("parse_origin_region('chr1:100-') => {result:?}");
        assert_eq!(result, None);
    }

    // ═══════════════════════════════════════════════════════════════════════
    //  cal_similarity_score tests
    // ═══════════════════════════════════════════════════════════════════════

    /// Helper — build a RegionVarStats with explicit fields.
    fn make_stats(
        varcount: i32,
        alt_snv: i32,
        alt_indel: i32,
        shared_psv: i32,
        snv_pos: Vec<i32>,
        indel_pos: Vec<i32>,
        span: i64,
    ) -> RegionVarStats {
        RegionVarStats {
            varcount,
            alt_snv_count: alt_snv,
            alt_indel_count: alt_indel,
            shared_psv,
            verified_shared_snv_pos_abs: snv_pos,
            verified_shared_indel_pos_abs: indel_pos,
            overlap_span_size: span,
        }
    }

    #[test]
    fn test_cal_similarity_score_empty_input() {
        init_log();
        let varcounts: VarcountsAmongRefseqs = HashMap::new();
        let hid_var_count: HashMap<i32, i32> = HashMap::new();
        let hid_density: HashMap<i32, f64> = HashMap::new();

        let result = cal_similarity_score(&varcounts, &hid_var_count, &hid_density);
        debug!("empty input => {result:?}");
        assert!(result.is_empty());
    }

    #[test]
    fn test_cal_similarity_score_single_hap_single_refseq() {
        init_log();
        // hid=1, one refseq "refA", one region:
        //   shared_psv=10, alt_snv=10, alt_indel=0, total_varcount=10, density=0
        //   This gives ratio=1.0, sharing=1.0, metric=10-3=7.0
        let stats = make_stats(10, 10, 0, 10, vec![100, 200, 300, 400, 500, 600, 700, 800, 900, 1000], vec![], 500);
        let mut varcounts: VarcountsAmongRefseqs = HashMap::new();
        varcounts.entry(1).or_default().insert("refA".to_string(), vec![stats]);

        let mut hid_var_count: HashMap<i32, i32> = HashMap::new();
        hid_var_count.insert(1, 10);
        let hid_density: HashMap<i32, f64> = HashMap::new(); // density=0

        let result = cal_similarity_score(&varcounts, &hid_var_count, &hid_density);
        let r = result.get(&1).unwrap();

        // Hand calculation:
        //   psv_var_ratio = 10/10 = 1.0
        //   psv_sharing_ratio = 10/10 = 1.0
        //   metric = sqrt(1.0)*10*sqrt(1.0) + 0 - (4 - 1.0) = 10 - 3 = 7.0
        let expected = 7.0_f64;
        debug!("single_hap_single_refseq: score={:.10} expected={:.10}", r.max_sim_score, expected);
        assert!((r.max_sim_score - expected).abs() < 1e-10);
        assert_eq!(r.max_psv_count, 10);
        assert_eq!(r.max_psv_positions, vec![100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]);
    }

    #[test]
    fn test_cal_similarity_score_zero_total_varcount_positive_psv() {
        init_log();
        // total_varcount=0 but shared_psv=8 → psv_var_ratio = min(1, 8) = 1.0
        // Use high shared_psv so metric > -1
        let stats = make_stats(0, 8, 2, 8, vec![10, 20, 30, 40], vec![50, 60, 70, 80], 200);
        let mut varcounts: VarcountsAmongRefseqs = HashMap::new();
        varcounts.entry(1).or_default().insert("refA".to_string(), vec![stats]);

        let mut hid_var_count: HashMap<i32, i32> = HashMap::new();
        hid_var_count.insert(1, 0);
        let hid_density: HashMap<i32, f64> = HashMap::new();

        let result = cal_similarity_score(&varcounts, &hid_var_count, &hid_density);
        let r = result.get(&1).unwrap();

        // psv_var_ratio = min(1.0, 8.0) = 1.0
        // psv_sharing_ratio = 8/(8+2) = 0.8
        // metric = sqrt(1.0)*8*sqrt(0.8) + 0 - (4 - 1.0) = 8*0.89443 - 3.0 = 7.1554 - 3 = 4.1554
        let expected = 1.0_f64.sqrt() * 8.0 * (0.8_f64).sqrt() - 3.0;
        debug!("zero_varcount_positive_psv: score={:.10} expected={:.10}", r.max_sim_score, expected);
        assert!((r.max_sim_score - expected).abs() < 1e-10);
        assert_eq!(r.max_psv_count, 8);
    }

    #[test]
    fn test_cal_similarity_score_zero_total_varcount_zero_psv() {
        init_log();
        // total_varcount=0, shared_psv=0 → psv_var_ratio = min(1, 0) = 0.0
        // metric = -4.0 which is < -1, so the floor (-1) wins
        let stats = make_stats(0, 0, 0, 0, vec![], vec![], 100);
        let mut varcounts: VarcountsAmongRefseqs = HashMap::new();
        varcounts.entry(1).or_default().insert("refA".to_string(), vec![stats]);

        let mut hid_var_count: HashMap<i32, i32> = HashMap::new();
        hid_var_count.insert(1, 0);
        let hid_density: HashMap<i32, f64> = HashMap::new();

        let result = cal_similarity_score(&varcounts, &hid_var_count, &hid_density);
        let r = result.get(&1).unwrap();

        // metric = sqrt(0)*0*sqrt(0) + 0 - 4 = -4, but -4 < -1 (initial), so floor wins
        debug!("zero_varcount_zero_psv: score={:.10}", r.max_sim_score);
        assert!((r.max_sim_score - (-1.0)).abs() < 1e-10, "floor should be -1");
        assert_eq!(r.max_psv_count, 0);
        assert!(r.max_psv_positions.is_empty());
    }

    #[test]
    fn test_cal_similarity_score_zero_total_psv_count() {
        init_log();
        // alt_snv=0, alt_indel=0 → psv_sharing_ratio = 0.0
        // metric = 0 - 3.6 = -3.6 < -1, so floor wins
        // This specifically tests the psv_sharing_ratio=0 branch
        let stats = make_stats(5, 0, 0, 2, vec![50], vec![60], 100);
        let mut varcounts: VarcountsAmongRefseqs = HashMap::new();
        varcounts.entry(1).or_default().insert("refA".to_string(), vec![stats]);

        let mut hid_var_count: HashMap<i32, i32> = HashMap::new();
        hid_var_count.insert(1, 5);
        let hid_density: HashMap<i32, f64> = HashMap::new();

        let result = cal_similarity_score(&varcounts, &hid_var_count, &hid_density);
        let r = result.get(&1).unwrap();

        // metric = sqrt(0.4)*2*sqrt(0.0) + 0 - 3.6 = -3.6 < -1 → floor
        debug!("zero_total_psv_count: score={:.10}", r.max_sim_score);
        assert!((r.max_sim_score - (-1.0)).abs() < 1e-10, "floor should be -1");
    }

    #[test]
    fn test_cal_similarity_score_multiple_refseqs_picks_max() {
        init_log();
        // hid=1, two refseqs: "refA" (lower metric) vs "refB" (higher metric)
        let stats_a = make_stats(10, 3, 1, 2, vec![100], vec![], 200);
        let stats_b = make_stats(10, 8, 4, 8, vec![10, 20, 30, 40], vec![50, 60, 70, 80], 300);

        let mut varcounts: VarcountsAmongRefseqs = HashMap::new();
        let inner = varcounts.entry(1).or_default();
        inner.insert("refA".to_string(), vec![stats_a]);
        inner.insert("refB".to_string(), vec![stats_b]);

        let mut hid_var_count: HashMap<i32, i32> = HashMap::new();
        hid_var_count.insert(1, 10);
        let hid_density: HashMap<i32, f64> = HashMap::new();

        let result = cal_similarity_score(&varcounts, &hid_var_count, &hid_density);
        let r = result.get(&1).unwrap();

        // refA: ratio=2/10=0.2, sharing=2/(3+1)=0.5
        //   metric_a = sqrt(0.2)*2*sqrt(0.5) + 0 - 3.8 = 0.4472*2*0.7071 - 3.8 = 0.6324 - 3.8 = -3.1675...
        let metric_a = (0.2_f64).sqrt() * 2.0 * (0.5_f64).sqrt() - (4.0 - 0.2);

        // refB: ratio=8/10=0.8, sharing=8/(8+4)=0.6667
        //   metric_b = sqrt(0.8)*8*sqrt(0.6667) + 0 - 3.2 = 0.8944*8*0.8165 - 3.2 = 5.8423 - 3.2 = 2.6423...
        let metric_b = (0.8_f64).sqrt() * 8.0 * (8.0_f64 / 12.0).sqrt() - (4.0 - 0.8);

        debug!(
            "multiple_refseqs: metric_a={:.10} metric_b={:.10} result={:.10}",
            metric_a, metric_b, r.max_sim_score
        );
        assert!(metric_b > metric_a);
        assert!((r.max_sim_score - metric_b).abs() < 1e-10);
        assert_eq!(r.max_psv_count, 8);
        assert_eq!(r.max_psv_positions, vec![10, 20, 30, 40, 50, 60, 70, 80]);
    }

    #[test]
    fn test_cal_similarity_score_multiple_regions_aggregation() {
        init_log();
        // Two regions (pairs) under one refseq: stats are summed, positions concatenated+deduped
        // Use high values so metric > -1
        let s1 = make_stats(5, 5, 1, 5, vec![100, 200, 300, 400, 500], vec![], 150);
        let s2 = make_stats(5, 5, 1, 5, vec![600, 700, 800, 900, 1000], vec![], 150);

        let mut varcounts: VarcountsAmongRefseqs = HashMap::new();
        varcounts.entry(1).or_default().insert("refA".to_string(), vec![s1, s2]);

        let mut hid_var_count: HashMap<i32, i32> = HashMap::new();
        hid_var_count.insert(1, 12);
        let hid_density: HashMap<i32, f64> = HashMap::new();

        let result = cal_similarity_score(&varcounts, &hid_var_count, &hid_density);
        let r = result.get(&1).unwrap();

        // Aggregated: shared_psv=5+5=10, alt_snv=5+5=10, alt_indel=1+1=2, total_varcount=12
        // positions: [100..1000] (all unique, 10 values)
        // ratio=10/12, sharing=10/(10+2)=10/12
        // metric = sqrt(10/12)*10*sqrt(10/12) + 0 - (4 - 10/12)
        //        = (10/12)*10 - (4 - 10/12) = 100/12 - 38/12 = 62/12 ≈ 5.1667
        let ratio: f64 = 10.0 / 12.0;
        let expected = ratio.sqrt() * 10.0 * ratio.sqrt() - (4.0 - ratio);
        debug!("multiple_regions: score={:.10} expected={:.10}", r.max_sim_score, expected);
        assert!((r.max_sim_score - expected).abs() < 1e-10);
        assert_eq!(r.max_psv_count, 10);
        assert_eq!(r.max_psv_positions, vec![100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]);
    }

    #[test]
    fn test_cal_similarity_score_duplicate_psv_deduped() {
        init_log();
        // Two regions with overlapping positions → deduped.
        // High values so metric > -1 to verify positions at the winning refseq.
        let s1 = make_stats(5, 5, 0, 5, vec![100, 200, 300, 400, 500], vec![], 100);
        let s2 = make_stats(5, 5, 0, 5, vec![200, 400, 600, 700, 800], vec![300], 100);

        let mut varcounts: VarcountsAmongRefseqs = HashMap::new();
        varcounts.entry(1).or_default().insert("refA".to_string(), vec![s1, s2]);

        let mut hid_var_count: HashMap<i32, i32> = HashMap::new();
        hid_var_count.insert(1, 12);
        let hid_density: HashMap<i32, f64> = HashMap::new();

        let result = cal_similarity_score(&varcounts, &hid_var_count, &hid_density);
        let r = result.get(&1).unwrap();

        // SNV pos: [100,200,300,400,500] ∪ [200,400,600,700,800] = 8 unique
        // Indel pos: [] ∪ [300] = [300]
        // Combined unique sorted: [100, 200, 300, 400, 500, 600, 700, 800]
        assert_eq!(r.max_psv_positions, vec![100, 200, 300, 400, 500, 600, 700, 800]);
        debug!("dedup positions: {:?}", r.max_psv_positions);
    }

    #[test]
    fn test_cal_similarity_score_multiple_haplotypes_independent() {
        init_log();
        // hid=1: high PSV → metric > -1
        let stats1 = make_stats(10, 10, 0, 10, vec![100, 200, 300, 400, 500, 600, 700, 800, 900, 1000], vec![], 200);
        // hid=2: also high PSV → metric > -1
        let stats2 = make_stats(20, 15, 5, 12, vec![1100, 1200, 1300, 1400, 1500, 1600], vec![1700, 1800, 1900, 2000, 2100, 2200], 300);

        let mut varcounts: VarcountsAmongRefseqs = HashMap::new();
        varcounts.entry(1).or_default().insert("refA".to_string(), vec![stats1]);
        varcounts.entry(2).or_default().insert("refB".to_string(), vec![stats2]);

        let mut hid_var_count: HashMap<i32, i32> = HashMap::new();
        hid_var_count.insert(1, 10);
        hid_var_count.insert(2, 20);
        let hid_density: HashMap<i32, f64> = HashMap::new();

        let result = cal_similarity_score(&varcounts, &hid_var_count, &hid_density);
        assert_eq!(result.len(), 2);

        // hid=1: ratio=10/10=1.0, sharing=10/10=1.0
        //   metric = 1*10*1 - 3 = 7.0
        let expected_1 = 7.0_f64;
        // hid=2: ratio=12/20=0.6, sharing=12/(15+5)=0.6
        //   metric = sqrt(0.6)*12*sqrt(0.6) - (4-0.6) = 0.6*12 - 3.4 = 7.2 - 3.4 = 3.8
        let expected_2 = (0.6_f64).sqrt() * 12.0 * (0.6_f64).sqrt() - (4.0 - 0.6);

        let r1 = result.get(&1).unwrap();
        let r2 = result.get(&2).unwrap();
        debug!("hid1={:.10} expected={:.10}", r1.max_sim_score, expected_1);
        debug!("hid2={:.10} expected={:.10}", r2.max_sim_score, expected_2);
        assert!((r1.max_sim_score - expected_1).abs() < 1e-10);
        assert!((r2.max_sim_score - expected_2).abs() < 1e-10);
        assert_eq!(r1.max_psv_count, 10);
        assert_eq!(r2.max_psv_count, 12);
    }

    #[test]
    fn test_cal_similarity_score_with_nonzero_density() {
        init_log();
        let stats = make_stats(10, 6, 2, 4, vec![100, 200], vec![300], 200);
        let mut varcounts: VarcountsAmongRefseqs = HashMap::new();
        varcounts.entry(1).or_default().insert("refA".to_string(), vec![stats]);

        let mut hid_var_count: HashMap<i32, i32> = HashMap::new();
        hid_var_count.insert(1, 10);
        let mut hid_density: HashMap<i32, f64> = HashMap::new();
        hid_density.insert(1, 0.05); // max_density = 0.05, so non_psv_density_100bp = 5.0

        let result = cal_similarity_score(&varcounts, &hid_var_count, &hid_density);
        let r = result.get(&1).unwrap();

        // ratio=4/10=0.4, sharing=4/8=0.5, density_100bp=5.0
        // metric = sqrt(0.4)*4*sqrt(0.5) + sqrt(5.0) - (4 - 0.4)
        let expected = (0.4_f64).sqrt() * 4.0 * (0.5_f64).sqrt()
            + (5.0_f64).sqrt()
            - (4.0 - 0.4);
        debug!("with_density: score={:.10} expected={:.10}", r.max_sim_score, expected);
        assert!((r.max_sim_score - expected).abs() < 1e-10);

        // Also verify the density term made a positive difference vs zero density
        let expected_no_density = (0.4_f64).sqrt() * 4.0 * (0.5_f64).sqrt() - (4.0 - 0.4);
        assert!(r.max_sim_score > expected_no_density);
    }

    #[test]
    fn test_cal_similarity_score_negative_metric_floor() {
        init_log();
        // Very few shared PSVs with high total_varcount → metric = -3.98 < -1 → floor wins
        let stats = make_stats(100, 50, 50, 1, vec![42], vec![], 500);
        let mut varcounts: VarcountsAmongRefseqs = HashMap::new();
        varcounts.entry(1).or_default().insert("refA".to_string(), vec![stats]);

        let mut hid_var_count: HashMap<i32, i32> = HashMap::new();
        hid_var_count.insert(1, 100);
        let hid_density: HashMap<i32, f64> = HashMap::new();

        let result = cal_similarity_score(&varcounts, &hid_var_count, &hid_density);
        let r = result.get(&1).unwrap();

        // metric = sqrt(0.01)*1*sqrt(0.01) - 3.99 = 0.01 - 3.99 = -3.98
        // -3.98 < -1 (initial max_psv), so floor of -1.0 wins
        debug!("negative_metric_floor: score={:.10}", r.max_sim_score);
        assert!((r.max_sim_score - (-1.0)).abs() < 1e-10, "floor should be -1");
        // Since no update happened, psv_count remains 0 and positions empty
        assert_eq!(r.max_psv_count, 0);
        assert!(r.max_psv_positions.is_empty());
    }

    /// Helper: create a minimal BAM record at given pos with a simple `n=` CIGAR.
    /// The record spans `[pos, pos + length)`.
    fn make_read(pos: i64, length: u32) -> BamRecord {
        let cigar = CigarString(vec![Cigar::Equal(length)]);
        let seq: Vec<u8> = vec![b'A'; length as usize];
        let qual: Vec<u8> = vec![30; length as usize];
        let mut record = BamRecord::new();
        record.set(b"read", Some(&cigar), &seq, &qual);
        record.set_pos(pos);
        record.set_tid(0);
        record.set_mapq(60);
        record
    }

    /// Helper: create a named read at given pos with a simple `n=` CIGAR.
    fn make_named_read(name: &[u8], pos: i64, length: u32) -> BamRecord {
        let cigar = CigarString(vec![Cigar::Equal(length)]);
        let seq: Vec<u8> = vec![b'A'; length as usize];
        let qual: Vec<u8> = vec![30; length as usize];
        let mut record = BamRecord::new();
        record.set(name, Some(&cigar), &seq, &qual);
        record.set_pos(pos);
        record.set_tid(0);
        record.set_mapq(60);
        record
    }

    #[test]
    fn test_continuous_regions_empty() {
        init_log();
        let reads: Vec<&Record> = vec![];
        let result = extract_continuous_regions_dict(&reads);
        debug!("empty reads => {result:?}");
        assert!(result.is_empty());
    }

    #[test]
    fn test_continuous_regions_single_read() {
        init_log();
        let r = make_read(100, 50); // [100, 150)
        let reads: Vec<&Record> = vec![&r];
        let result = extract_continuous_regions_dict(&reads);
        debug!("single read [100,150) => {result:?}");
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0, (100, 150));
        assert_eq!(result[0].1, vec![0]);
    }

    #[test]
    fn test_continuous_regions_two_overlapping() {
        init_log();
        let r1 = make_read(100, 50); // [100, 150)
        let r2 = make_read(120, 50); // [120, 170)
        let reads: Vec<&Record> = vec![&r1, &r2];
        let result = extract_continuous_regions_dict(&reads);
        debug!("two overlapping [100,150)+[120,170) => {result:?}");
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0, (100, 170));
        // Both reads in the single region; indices are into the original slice
        assert_eq!(result[0].1.len(), 2);
        assert!(result[0].1.contains(&0));
        assert!(result[0].1.contains(&1));
    }

    #[test]
    fn test_continuous_regions_two_disjoint() {
        init_log();
        let r1 = make_read(100, 50);  // [100, 150)
        let r2 = make_read(200, 50);  // [200, 250)
        let reads: Vec<&Record> = vec![&r1, &r2];
        let result = extract_continuous_regions_dict(&reads);
        debug!("two disjoint [100,150)+[200,250) => {result:?}");
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].0, (100, 150));
        assert_eq!(result[0].1, vec![0]);
        assert_eq!(result[1].0, (200, 250));
        assert_eq!(result[1].1, vec![1]);
    }

    #[test]
    fn test_continuous_regions_touching_boundaries() {
        init_log();
        // r1 ends at 150, r2 starts at 150 → read_start <= current_end (150 <= 150) is TRUE
        // so they MERGE into one region (matching Python behavior)
        let r1 = make_read(100, 50);  // [100, 150)
        let r2 = make_read(150, 50);  // [150, 200)
        let reads: Vec<&Record> = vec![&r1, &r2];
        let result = extract_continuous_regions_dict(&reads);
        debug!("touching [100,150)+[150,200) => {result:?}");
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0, (100, 200));
        assert_eq!(result[0].1.len(), 2);
    }

    #[test]
    fn test_continuous_regions_gap_of_one() {
        init_log();
        // r1 ends at 150, r2 starts at 151 → read_start <= current_end (151 <= 150) is FALSE
        // so they are SEPARATE regions
        let r1 = make_read(100, 50);  // [100, 150)
        let r2 = make_read(151, 50);  // [151, 201)
        let reads: Vec<&Record> = vec![&r1, &r2];
        let result = extract_continuous_regions_dict(&reads);
        debug!("gap of 1 [100,150)+[151,201) => {result:?}");
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].0, (100, 150));
        assert_eq!(result[1].0, (151, 201));
    }

    #[test]
    fn test_continuous_regions_unsorted_input() {
        init_log();
        // Reads given out of order — function should sort internally
        let r1 = make_read(200, 50);  // [200, 250)
        let r2 = make_read(100, 50);  // [100, 150)
        let r3 = make_read(120, 50);  // [120, 170)
        let reads: Vec<&Record> = vec![&r1, &r2, &r3];
        let result = extract_continuous_regions_dict(&reads);
        debug!("unsorted input [200,250)+[100,150)+[120,170) => {result:?}");
        // r2 [100,150) and r3 [120,170) overlap → region [100,170)
        // r1 [200,250) is separate → region [200,250)
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].0, (100, 170));
        assert_eq!(result[0].1.len(), 2);
        // Original indices: r1=0, r2=1, r3=2
        // Region [100,170) should contain original indices 1 (r2) and 2 (r3)
        assert!(result[0].1.contains(&1));
        assert!(result[0].1.contains(&2));
        // Region [200,250) should contain original index 0 (r1)
        assert_eq!(result[1].0, (200, 250));
        assert_eq!(result[1].1, vec![0]);
    }

    #[test]
    fn test_continuous_regions_contained_read() {
        init_log();
        // One read fully contains another
        let r1 = make_read(100, 100);  // [100, 200)
        let r2 = make_read(130, 20);   // [130, 150) — fully inside r1
        let reads: Vec<&Record> = vec![&r1, &r2];
        let result = extract_continuous_regions_dict(&reads);
        debug!("contained [100,200)+[130,150) => {result:?}");
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0, (100, 200));
        assert_eq!(result[0].1.len(), 2);
    }

    #[test]
    fn test_continuous_regions_chain_merge() {
        init_log();
        // Three reads forming a chain: each overlaps the next
        let r1 = make_read(100, 50);  // [100, 150)
        let r2 = make_read(140, 50);  // [140, 190)
        let r3 = make_read(180, 50);  // [180, 230)
        let reads: Vec<&Record> = vec![&r1, &r2, &r3];
        let result = extract_continuous_regions_dict(&reads);
        debug!("chain merge [100,150)+[140,190)+[180,230) => {result:?}");
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0, (100, 230));
        assert_eq!(result[0].1.len(), 3);
    }

    #[test]
    fn test_continuous_regions_three_separate() {
        init_log();
        let r1 = make_read(100, 10);   // [100, 110)
        let r2 = make_read(200, 10);   // [200, 210)
        let r3 = make_read(300, 10);   // [300, 310)
        let reads: Vec<&Record> = vec![&r1, &r2, &r3];
        let result = extract_continuous_regions_dict(&reads);
        debug!("three separate regions => {result:?}");
        assert_eq!(result.len(), 3);
        assert_eq!(result[0].0, (100, 110));
        assert_eq!(result[1].0, (200, 210));
        assert_eq!(result[2].0, (300, 310));
    }

    #[test]
    fn test_continuous_regions_many_reads_one_region() {
        init_log();
        // 10 reads all overlapping in one big pile
        let records: Vec<BamRecord> = (0..10)
            .map(|i| make_read(100 + i * 5, 50))
            .collect();
        let reads: Vec<&Record> = records.iter().map(|r| r as &Record).collect();
        let result = extract_continuous_regions_dict(&reads);
        debug!("10 overlapping reads => {result:?}");
        assert_eq!(result.len(), 1);
        // First read starts at 100, last read ends at 100 + 9*5 + 50 = 195
        assert_eq!(result[0].0, (100, 195));
        assert_eq!(result[0].1.len(), 10);
    }

    #[test]
    fn test_continuous_regions_preserves_original_indices() {
        init_log();
        // Verify that returned indices point to the correct reads in the original slice
        let r_a = make_named_read(b"readA", 300, 50);  // idx 0
        let r_b = make_named_read(b"readB", 100, 50);  // idx 1
        let r_c = make_named_read(b"readC", 110, 50);  // idx 2
        let r_d = make_named_read(b"readD", 310, 50);  // idx 3
        let reads: Vec<&Record> = vec![&r_a, &r_b, &r_c, &r_d];
        let result = extract_continuous_regions_dict(&reads);

        debug!("preserve indices => {result:?}");
        // Region 1: [100, 160) — readB(1) + readC(2)
        // Region 2: [300, 360) — readA(0) + readD(3)
        assert_eq!(result.len(), 2);

        let (span1, idxs1) = &result[0];
        assert_eq!(*span1, (100, 160));
        for &idx in idxs1 {
            let name = String::from_utf8_lossy(reads[idx].qname());
            debug!("  region1 idx={idx} name={name}");
            assert!(name == "readB" || name == "readC");
        }

        let (span2, idxs2) = &result[1];
        assert_eq!(*span2, (300, 360));
        for &idx in idxs2 {
            let name = String::from_utf8_lossy(reads[idx].qname());
            debug!("  region2 idx={idx} name={name}");
            assert!(name == "readA" || name == "readD");
        }
    }

    #[test]
    fn test_inspection_placeholder() {
        // Placeholder test - will be expanded during implementation
        assert!(true);
    }

    // ═══════════════════════════════════════════════════════════════════════
    //  group_by_dict_optimized tests
    // ═══════════════════════════════════════════════════════════════════════

    #[test]
    fn test_group_by_empty() {
        init_log();
        let vprop: HashMap<i32, i32> = HashMap::new();
        let vertices: HashMap<(i32, String), Vec<&Record>> = HashMap::new();
        let result = group_by_dict_optimized(&vprop, &vertices);
        debug!("group_by empty => {} groups", result.len());
        assert!(result.is_empty());
    }

    #[test]
    fn test_group_by_single_entry() {
        init_log();
        let r = make_named_read(b"qA", 100, 50);
        let mut vprop: HashMap<i32, i32> = HashMap::new();
        vprop.insert(10, 1); // vertex 10 → hap 1
        let mut vertices: HashMap<(i32, String), Vec<&Record>> = HashMap::new();
        vertices.insert((10, "qA".to_string()), vec![&r]);

        let result = group_by_dict_optimized(&vprop, &vertices);
        debug!("group_by single => {} groups", result.len());
        assert_eq!(result.len(), 1);
        let g = result.get(&1).unwrap();
        assert_eq!(g.vertex_indices, vec![10]);
        assert_eq!(g.qnames, vec!["qA"]);
        assert_eq!(g.read_pair_lists.len(), 1);
        assert_eq!(g.read_pair_lists[0].len(), 1);
    }

    #[test]
    fn test_group_by_two_same_label() {
        init_log();
        let r1 = make_named_read(b"qA", 100, 50);
        let r2 = make_named_read(b"qB", 200, 50);
        let mut vprop: HashMap<i32, i32> = HashMap::new();
        vprop.insert(10, 1);
        vprop.insert(20, 1);
        let mut vertices: HashMap<(i32, String), Vec<&Record>> = HashMap::new();
        vertices.insert((10, "qA".to_string()), vec![&r1]);
        vertices.insert((20, "qB".to_string()), vec![&r2]);

        let result = group_by_dict_optimized(&vprop, &vertices);
        assert_eq!(result.len(), 1);
        let g = result.get(&1).unwrap();
        assert_eq!(g.vertex_indices.len(), 2);
        assert_eq!(g.qnames.len(), 2);
        debug!("two_same_label: verts={:?} qnames={:?}", g.vertex_indices, g.qnames);
    }

    #[test]
    fn test_group_by_two_different_labels() {
        init_log();
        let r1 = make_named_read(b"qA", 100, 50);
        let r2 = make_named_read(b"qB", 200, 50);
        let mut vprop: HashMap<i32, i32> = HashMap::new();
        vprop.insert(10, 1);
        vprop.insert(20, 2);
        let mut vertices: HashMap<(i32, String), Vec<&Record>> = HashMap::new();
        vertices.insert((10, "qA".to_string()), vec![&r1]);
        vertices.insert((20, "qB".to_string()), vec![&r2]);

        let result = group_by_dict_optimized(&vprop, &vertices);
        assert_eq!(result.len(), 2);
        assert!(result.contains_key(&1));
        assert!(result.contains_key(&2));
        debug!("two_different_labels: keys={:?}", result.keys().collect::<Vec<_>>());
    }

    #[test]
    fn test_group_by_multiple_reads_per_entry() {
        init_log();
        let r1 = make_named_read(b"qA", 100, 50);
        let r2 = make_named_read(b"qA", 200, 50); // mate pair
        let mut vprop: HashMap<i32, i32> = HashMap::new();
        vprop.insert(10, 5);
        let mut vertices: HashMap<(i32, String), Vec<&Record>> = HashMap::new();
        vertices.insert((10, "qA".to_string()), vec![&r1, &r2]);

        let result = group_by_dict_optimized(&vprop, &vertices);
        let g = result.get(&5).unwrap();
        assert_eq!(g.read_pair_lists[0].len(), 2);
        debug!("multiple_reads: {} reads in pair list", g.read_pair_lists[0].len());
    }

    // ═══════════════════════════════════════════════════════════════════════
    //  record_haplotype_rank tests
    // ═══════════════════════════════════════════════════════════════════════

    #[test]
    fn test_rank_empty() {
        init_log();
        let hap_dict: HashMap<i32, HaplotypeClusterInfo> = HashMap::new();
        let psv_pos: HashMap<i32, Vec<i32>> = HashMap::new();
        let result = record_haplotype_rank(&hap_dict, 150, &psv_pos);
        debug!("rank_empty => shape {:?}", result.shape());
        assert_eq!(result.shape(), &[0, 8]);
    }

    #[test]
    fn test_rank_single_haplotype() {
        init_log();
        // hid=1: consensus [1,1,-4,1] at span [100,104), two reads covering it
        let r1 = make_named_read(b"qA", 100, 4); // [100,104)
        let r2 = make_named_read(b"qB", 100, 4); // [100,104)
        let consensus = Array1::from(vec![1_i16, 1, -4, 1]);

        let mut hap_dict: HashMap<i32, HaplotypeClusterInfo> = HashMap::new();
        hap_dict.insert(1, HaplotypeClusterInfo {
            consensus: consensus.clone(),
            reads: vec![&r1, &r2],
            span: (100, 104),
            qnames: vec!["qA".to_string(), "qB".to_string()],
        });

        // PSV positions: [101, 102] — 101 is in [100,104), 102 is in [100,104)
        let mut psv_pos: HashMap<i32, Vec<i32>> = HashMap::new();
        psv_pos.insert(1, vec![101, 102]);

        let result = record_haplotype_rank(&hap_dict, 150, &psv_pos);
        debug!("rank_single: {result:?}");
        assert_eq!(result.shape(), &[1, 8]);

        // Col 0: start=100, Col 1: end=104
        assert_eq!(result[[0, 0]], 100);
        assert_eq!(result[[0, 1]], 104);
        // Col 3: hap_id=1
        assert_eq!(result[[0, 3]], 1);
        // Col 4: depth = (4+4)/4 = 2 (two reads each covering 4 bases, region_len=4)
        assert_eq!(result[[0, 4]], 2);
        // Col 5: var_count = count_var([1,1,-4,1]) = 1 (the -4 is an SNV)
        assert_eq!(result[[0, 5]], 1);
        // Col 6: indel_count = 0 (no deletions)
        assert_eq!(result[[0, 6]], 0);
        // Col 7: psv_count = 2 (both 101,102 are in [100,104))
        assert_eq!(result[[0, 7]], 2);
        // Col 2: total_depth = sum of all depths = 2
        assert_eq!(result[[0, 2]], 2);
    }

    #[test]
    fn test_rank_two_haplotypes_total_depth() {
        init_log();
        let r1 = make_named_read(b"q1", 100, 10); // [100,110)
        let r2 = make_named_read(b"q2", 200, 10); // [200,210)
        let con1 = Array1::from(vec![1_i16; 10]);
        let con2 = Array1::from(vec![1_i16, -4, 1, 1, 1, 1, 1, -4, 1, 1]); // 2 SNVs

        let mut hap_dict: HashMap<i32, HaplotypeClusterInfo> = HashMap::new();
        hap_dict.insert(1, HaplotypeClusterInfo {
            consensus: con1,
            reads: vec![&r1],
            span: (100, 110),
            qnames: vec!["q1".to_string()],
        });
        hap_dict.insert(2, HaplotypeClusterInfo {
            consensus: con2,
            reads: vec![&r2],
            span: (200, 210),
            qnames: vec!["q2".to_string()],
        });

        let psv_pos: HashMap<i32, Vec<i32>> = HashMap::new(); // empty

        let result = record_haplotype_rank(&hap_dict, 150, &psv_pos);
        debug!("rank_two: {result:?}");
        assert_eq!(result.shape(), &[2, 8]);

        // Sorted by hid: hid=1 is row 0, hid=2 is row 1
        assert_eq!(result[[0, 3]], 1);
        assert_eq!(result[[1, 3]], 2);
        // hid 1: depth=1, hid 2: depth=1, total_depth=2
        assert_eq!(result[[0, 4]], 1);
        assert_eq!(result[[1, 4]], 1);
        assert_eq!(result[[0, 2]], 2);
        assert_eq!(result[[1, 2]], 2);
        // hid 1: var_count=0, hid 2: var_count=2
        assert_eq!(result[[0, 5]], 0);
        assert_eq!(result[[1, 5]], 2);
        // PSV counts are 0 for both (empty psv_pos)
        assert_eq!(result[[0, 7]], 0);
        assert_eq!(result[[1, 7]], 0);
    }

    #[test]
    fn test_rank_psv_intersection() {
        init_log();
        // Test that only PSV positions within the span are counted
        let r1 = make_named_read(b"qX", 100, 20); // [100,120)
        let con1 = Array1::from(vec![1_i16; 20]);

        let mut hap_dict: HashMap<i32, HaplotypeClusterInfo> = HashMap::new();
        hap_dict.insert(1, HaplotypeClusterInfo {
            consensus: con1,
            reads: vec![&r1],
            span: (100, 120),
            qnames: vec!["qX".to_string()],
        });

        // PSV positions: 90 (before), 100 (in), 110 (in), 119 (in), 120 (at end, exclusive, so out), 130 (after)
        let mut psv_pos: HashMap<i32, Vec<i32>> = HashMap::new();
        psv_pos.insert(1, vec![90, 100, 110, 119, 120, 130]);

        let result = record_haplotype_rank(&hap_dict, 150, &psv_pos);
        debug!("rank_psv_intersection: psv_count={}", result[[0, 7]]);
        // In range [100, 120): positions 100, 110, 119 → 3
        assert_eq!(result[[0, 7]], 3);
    }

    #[test]
    fn test_rank_depth_partial_overlap() {
        init_log();
        // Read only partially overlaps the span
        let r1 = make_named_read(b"qP", 95, 20);  // [95, 115)
        // Span is [100, 110), consensus length = 10
        let con = Array1::from(vec![1_i16; 10]);

        let mut hap_dict: HashMap<i32, HaplotypeClusterInfo> = HashMap::new();
        hap_dict.insert(1, HaplotypeClusterInfo {
            consensus: con,
            reads: vec![&r1],
            span: (100, 110),
            qnames: vec!["qP".to_string()],
        });
        let psv_pos: HashMap<i32, Vec<i32>> = HashMap::new();

        let result = record_haplotype_rank(&hap_dict, 150, &psv_pos);
        // Region: [100, 110), region_len=10
        // Read [95, 115): clamped to [100, 110) → cov_bases=10
        // depth = 10/10 = 1
        debug!("rank_partial_overlap: depth={}", result[[0, 4]]);
        assert_eq!(result[[0, 4]], 1);
    }

    // ═══════════════════════════════════════════════════════════════════════
    //  summarize_enclosing_haps tests
    // ═══════════════════════════════════════════════════════════════════════

    /// Helper: build a HapGroup from a set of named reads
    fn make_hap_group<'a>(
        vertex_indices: Vec<i32>,
        qnames: Vec<String>,
        read_vecs: Vec<Vec<&'a Record>>,
    ) -> HapGroup<'a> {
        HapGroup {
            vertex_indices,
            qnames,
            read_pair_lists: read_vecs,
        }
    }

    #[test]
    fn test_summarize_three_enclosing_haps() {
        init_log();
        // Window [200, 300). Three haplotypes each with reads fully enclosing.
        let r1 = make_named_read(b"q1", 100, 250); // [100, 350) encloses [200,300)
        let r2 = make_named_read(b"q2", 150, 200); // [150, 350) encloses [200,300)
        let r3 = make_named_read(b"q3", 180, 150); // [180, 330) encloses [200,300)

        let mut qname_to_node: HashMap<String, i32> = HashMap::new();
        qname_to_node.insert("q1".to_string(), 10);
        qname_to_node.insert("q2".to_string(), 20);
        qname_to_node.insert("q3".to_string(), 30);

        let mut hap_subgraphs: HashMap<i32, HapGroup> = HashMap::new();
        hap_subgraphs.insert(1, make_hap_group(vec![10], vec!["q1".to_string()], vec![vec![&r1]]));
        hap_subgraphs.insert(2, make_hap_group(vec![20], vec!["q2".to_string()], vec![vec![&r2]]));
        hap_subgraphs.insert(3, make_hap_group(vec![30], vec!["q3".to_string()], vec![vec![&r3]]));

        let result = summarize_enclosing_haps(&hap_subgraphs, &qname_to_node, ("chr1", 200, 300));
        debug!("three_enclosing: {:?}", result.as_ref().map(|(m, s)| (m.len(), s)));
        assert!(result.is_some());
        let (info, span) = result.unwrap();
        assert_eq!(span, (200, 300)); // original window
        assert_eq!(info.len(), 3);
    }

    #[test]
    fn test_summarize_none_enclosing_no_recovery() {
        init_log();
        // Window [200, 300). All reads are too far away.
        let r1 = make_named_read(b"q1", 0, 50); // [0, 50) — no overlap

        let mut qname_to_node: HashMap<String, i32> = HashMap::new();
        qname_to_node.insert("q1".to_string(), 10);

        let mut hap_subgraphs: HashMap<i32, HapGroup> = HashMap::new();
        hap_subgraphs.insert(1, make_hap_group(vec![10], vec!["q1".to_string()], vec![vec![&r1]]));

        let result = summarize_enclosing_haps(&hap_subgraphs, &qname_to_node, ("chr1", 200, 300));
        debug!("none_enclosing: is_some={}", result.is_some());
        assert!(result.is_none());
    }

    #[test]
    fn test_summarize_recovery_with_shrinking() {
        init_log();
        // Window [200, 300). Only 1 hap fully encloses; 2 partial with overlap >= 0.8
        // Hap1: [100, 350) — fully encloses
        let r1 = make_named_read(b"q1", 100, 250);
        // Hap2: [210, 350) — overlap = (300-210)/(300-200) = 90/100 = 0.9
        let r2 = make_named_read(b"q2", 210, 140);
        // Hap3: [190, 290) — overlap = (290-200)/(300-200) = 90/100 = 0.9
        let r3 = make_named_read(b"q3", 190, 100);

        let mut qname_to_node: HashMap<String, i32> = HashMap::new();
        qname_to_node.insert("q1".to_string(), 10);
        qname_to_node.insert("q2".to_string(), 20);
        qname_to_node.insert("q3".to_string(), 30);

        let mut hap_subgraphs: HashMap<i32, HapGroup> = HashMap::new();
        hap_subgraphs.insert(1, make_hap_group(vec![10], vec!["q1".to_string()], vec![vec![&r1]]));
        hap_subgraphs.insert(2, make_hap_group(vec![20], vec!["q2".to_string()], vec![vec![&r2]]));
        hap_subgraphs.insert(3, make_hap_group(vec![30], vec!["q3".to_string()], vec![vec![&r3]]));

        let result = summarize_enclosing_haps(&hap_subgraphs, &qname_to_node, ("chr1", 200, 300));
        debug!("recovery: {:?}", result.as_ref().map(|(m, s)| (m.len(), s)));
        assert!(result.is_some());
        let (info, span) = result.unwrap();
        // Recovery shrinks: max start(210,190)=210, min end(350,290)=290
        // Clamp to window: max(210,200)=210, min(290,300)=290
        // So overlapping_span = (210, 290)
        assert_eq!(span, (210, 290));
        // Should have 3 total haplotypes (1 enclosing + 2 recovered)
        assert!(info.len() >= 2); // at least 2 for non-None return
    }

    #[test]
    fn test_summarize_recovery_fails_low_overlap() {
        init_log();
        // Window [200, 400). Hap1 partially overlaps but coef < 0.8
        // [250, 350) → overlap = (350-250)/(400-200) = 100/200 = 0.5 (< 0.8)
        let r1 = make_named_read(b"q1", 250, 100);
        let r2 = make_named_read(b"q2", 260, 80); // similar range

        let mut qname_to_node: HashMap<String, i32> = HashMap::new();
        qname_to_node.insert("q1".to_string(), 10);
        qname_to_node.insert("q2".to_string(), 20);

        let mut hap_subgraphs: HashMap<i32, HapGroup> = HashMap::new();
        hap_subgraphs.insert(1, make_hap_group(vec![10], vec!["q1".to_string()], vec![vec![&r1]]));
        hap_subgraphs.insert(2, make_hap_group(vec![20], vec!["q2".to_string()], vec![vec![&r2]]));

        let result = summarize_enclosing_haps(&hap_subgraphs, &qname_to_node, ("chr1", 200, 400));
        debug!("recovery_fails: is_some={}", result.is_some());
        assert!(result.is_none());
    }

    #[test]
    fn test_summarize_single_hap_after_recovery_returns_none() {
        init_log();
        // Window [200, 300). Only 1 hap with overlap >= 0.8, no enclosing haps.
        // After recovery still only 1 → return None
        let r1 = make_named_read(b"q1", 195, 100); // [195, 295) overlap=(295-200)/100=0.95

        let mut qname_to_node: HashMap<String, i32> = HashMap::new();
        qname_to_node.insert("q1".to_string(), 10);

        let mut hap_subgraphs: HashMap<i32, HapGroup> = HashMap::new();
        hap_subgraphs.insert(1, make_hap_group(vec![10], vec!["q1".to_string()], vec![vec![&r1]]));

        let result = summarize_enclosing_haps(&hap_subgraphs, &qname_to_node, ("chr1", 200, 300));
        debug!("single_after_recovery: is_some={}", result.is_some());
        assert!(result.is_none()); // <= 1 haplotype after recovery
    }

    // ── Tests for select_regions_with_min_haplotypes ────────────────────

    /// Helper: build hid_cov_intervals from a list of (hid, chrom, start, end)
    fn make_hid_cov(data: &[(i32, &str, i64, i64)]) -> HashMap<i32, Vec<(String, i64, i64)>> {
        let mut m: HashMap<i32, Vec<(String, i64, i64)>> = HashMap::new();
        for &(hid, chrom, start, end) in data {
            m.entry(hid).or_default().push((chrom.to_string(), start, end));
        }
        m
    }

    #[test]
    fn test_select_regions_too_few_haplotypes() {
        init_log();
        // Only 1 haplotype, need 2 → None
        let hid_cov = make_hid_cov(&[(1, "chr1", 100, 500)]);
        let result = select_regions_with_min_haplotypes(&hid_cov, 2);
        debug!("too_few_haplotypes: {result:?}");
        assert!(result.is_none());
    }

    #[test]
    fn test_select_regions_no_overlap() {
        init_log();
        // Two haplotypes, non-overlapping → None
        let hid_cov = make_hid_cov(&[
            (1, "chr1", 100, 200),
            (2, "chr1", 300, 400),
        ]);
        let result = select_regions_with_min_haplotypes(&hid_cov, 2);
        debug!("no_overlap: {result:?}");
        assert!(result.is_none());
    }

    #[test]
    fn test_select_regions_simple_overlap() {
        init_log();
        // Two haplotypes overlap in [200, 300)
        let hid_cov = make_hid_cov(&[
            (1, "chr1", 100, 300),
            (2, "chr1", 200, 400),
        ]);
        let result = select_regions_with_min_haplotypes(&hid_cov, 2);
        debug!("simple_overlap: {result:?}");
        let regions = result.unwrap();
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0], ("chr1".to_string(), 200, 300));
    }

    #[test]
    fn test_select_regions_three_haplotypes_pairwise() {
        init_log();
        // 3 haps, with different overlaps
        // hap1: [100, 400)
        // hap2: [200, 500)
        // hap3: [300, 600)
        // Partitions are NOT merged (mirrors BEDOPS --partition): each kept
        // breakpoint interval is returned separately so downstream per-region
        // inspection sees fine-grained spans.
        //   breakpoints: 100, 200, 300, 400, 500, 600
        //   [100,200): hap1 only → skip
        //   [200,300): hap1+hap2 → keep (2)
        //   [300,400): hap1+hap2+hap3 → keep (3)
        //   [400,500): hap2+hap3 → keep (2)
        //   [500,600): hap3 only → skip
        //   Result: three separate partitions [200,300), [300,400), [400,500)
        let hid_cov = make_hid_cov(&[
            (1, "chr1", 100, 400),
            (2, "chr1", 200, 500),
            (3, "chr1", 300, 600),
        ]);
        let result = select_regions_with_min_haplotypes(&hid_cov, 2);
        debug!("three_haplotypes_pairwise: {result:?}");
        let regions = result.unwrap();
        assert_eq!(regions.len(), 3);
        assert_eq!(regions[0], ("chr1".to_string(), 200, 300));
        assert_eq!(regions[1], ("chr1".to_string(), 300, 400));
        assert_eq!(regions[2], ("chr1".to_string(), 400, 500));
    }

    #[test]
    fn test_select_regions_min_3_haplotypes() {
        init_log();
        // Same 3 haps, but require min 3
        // Only [300,400) has all 3 → single region
        let hid_cov = make_hid_cov(&[
            (1, "chr1", 100, 400),
            (2, "chr1", 200, 500),
            (3, "chr1", 300, 600),
        ]);
        let result = select_regions_with_min_haplotypes(&hid_cov, 3);
        debug!("min_3_haplotypes: {result:?}");
        let regions = result.unwrap();
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0], ("chr1".to_string(), 300, 400));
    }

    #[test]
    fn test_select_regions_multiple_chroms() {
        init_log();
        // 2 haps on chr1 overlap, 2 haps on chr2 overlap
        let hid_cov = make_hid_cov(&[
            (1, "chr1", 100, 300),
            (2, "chr1", 200, 400),
            (1, "chr2", 500, 700),
            (2, "chr2", 600, 800),
        ]);
        let result = select_regions_with_min_haplotypes(&hid_cov, 2);
        debug!("multiple_chroms: {result:?}");
        let regions = result.unwrap();
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0], ("chr1".to_string(), 200, 300));
        assert_eq!(regions[1], ("chr2".to_string(), 600, 700));
    }

    #[test]
    fn test_select_regions_disjoint_overlaps() {
        init_log();
        // Two separate overlap zones that should NOT merge
        // hap1: [100, 200), [500, 600)
        // hap2: [150, 250), [550, 650)
        // Overlaps: [150,200) and [550,600) — two disjoint regions
        let hid_cov = make_hid_cov(&[
            (1, "chr1", 100, 200),
            (1, "chr1", 500, 600),
            (2, "chr1", 150, 250),
            (2, "chr1", 550, 650),
        ]);
        let result = select_regions_with_min_haplotypes(&hid_cov, 2);
        debug!("disjoint_overlaps: {result:?}");
        let regions = result.unwrap();
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0], ("chr1".to_string(), 150, 200));
        assert_eq!(regions[1], ("chr1".to_string(), 550, 600));
    }

    #[test]
    fn test_select_regions_identical_intervals() {
        init_log();
        // All 3 haps cover exact same region
        let hid_cov = make_hid_cov(&[
            (1, "chr1", 100, 500),
            (2, "chr1", 100, 500),
            (3, "chr1", 100, 500),
        ]);
        let result = select_regions_with_min_haplotypes(&hid_cov, 2);
        debug!("identical_intervals: {result:?}");
        let regions = result.unwrap();
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0], ("chr1".to_string(), 100, 500));
    }

    #[test]
    fn test_select_regions_empty_input() {
        init_log();
        let hid_cov: HashMap<i32, Vec<(String, i64, i64)>> = HashMap::new();
        let result = select_regions_with_min_haplotypes(&hid_cov, 2);
        debug!("empty_input: {result:?}");
        assert!(result.is_none());
    }

    #[test]
    fn test_select_regions_contained_interval() {
        init_log();
        // hap1: [100, 600), hap2: [200, 400) (contained within hap1)
        // Overlap: [200, 400)
        let hid_cov = make_hid_cov(&[
            (1, "chr1", 100, 600),
            (2, "chr1", 200, 400),
        ]);
        let result = select_regions_with_min_haplotypes(&hid_cov, 2);
        debug!("contained_interval: {result:?}");
        let regions = result.unwrap();
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0], ("chr1".to_string(), 200, 400));
    }
}
