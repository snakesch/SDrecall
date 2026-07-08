use ahash::AHashMap;
/// Haplotype determination logic
///
/// This module contains the core algorithm for determining whether two read pairs
/// come from the same haplotype by analyzing their sequences and CIGAR operations.
use log::{debug, error, warn};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::Record;

use crate::structs::{
    AlleleDepthMap, HaplotypeConfig, ReadErrorVector, ReadHaplotypeVector, Variant,
};

pub const HAP_DEL: i16 = -10;
pub const INDEL_UNIT: i16 = 10;

#[inline]
fn is_snv_value(v: i16) -> bool {
    v == -4 || (v > 1 && v % INDEL_UNIT == 6)
}

#[inline]
fn is_indel_value(v: i16) -> bool {
    v == HAP_DEL || v > 1
}

/// Result of haplotype compatibility analysis
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HaplotypeResult {
    Same,      // Same haplotype
    Different, // Different haplotypes
    Unknown,   // Cannot determine (insufficient data)
}

/// Extract query sequence and reference positions from a BAM record
///
/// Uses BamRecordExtensions::reference_positions_full() for efficient coordinate mapping
///
/// Returns (query_sequence, ref_positions) where:
/// - query_sequence: Vec<u8> of base-encoded sequence  
/// - ref_positions: Vec<i64> mapping query positions to reference positions
pub fn extract_query_seq(
    record: &Record,
) -> Result<(Vec<u8>, Vec<i64>), Box<dyn std::error::Error>> {
    let query_seq = record.seq().as_bytes();

    // Use BamRecordExtensions for efficient reference position mapping
    let ref_positions: Vec<i64> = record
        .reference_positions_full()
        .map(|pos_opt| pos_opt.map(|p| p as i64).unwrap_or(-1))
        .collect();

    // Convert to our format: A=0, T=1, C=2, G=3, N=4
    // Handle both uppercase and lowercase for robustness (soft-masking, etc.)
    let encoded_seq: Vec<u8> = query_seq
        .iter()
        .map(|&base| {
            match base {
                b'A' | b'a' => 0,
                b'T' | b't' => 1,
                b'C' | b'c' => 2,
                b'G' | b'g' => 3,
                _ => 4, // N or any other base (including ambiguous bases)
            }
        })
        .collect();

    Ok((encoded_seq, ref_positions))
}

/// Slice query sequence to a specific genomic interval using reference positions
///
/// Returns the subsequence (encoded A=0,T=1,C=2,G=3,N=4) that aligns to the given genomic interval
/// [interval_start, interval_end). Insertion bases (ref_positions == -1) are included if they are
/// flanked by aligned positions (nearest prev and next non -1) whose reference coordinates both fall
/// within the interval, indicating an insertion anchored inside the interval.
pub fn slice_seq_to_interval(
    query_seq: &[u8],
    ref_positions: &[i64],
    interval_start: i64,
    interval_end: i64,
) -> Vec<u8> {
    let mut result = Vec::new();
    let n = ref_positions.len();

    for i in 0..n {
        let rp = ref_positions[i];
        if rp >= interval_start && rp < interval_end {
            if i < query_seq.len() {
                result.push(query_seq[i]);
            }
            continue;
        }

        if rp == -1 {
            // Find nearest previous aligned reference position
            let mut prev_ref: Option<i64> = None;
            let mut j = i;
            while j > 0 {
                j -= 1;
                let rpj = ref_positions[j];
                if rpj != -1 {
                    prev_ref = Some(rpj);
                    break;
                }
            }
            // Find nearest next aligned reference position
            let mut next_ref: Option<i64> = None;
            let mut k = i + 1;
            while k < n {
                let rpk = ref_positions[k];
                if rpk != -1 {
                    next_ref = Some(rpk);
                    break;
                }
                k += 1;
            }

            // Include this insertion base only if both flanking aligned positions exist
            // and both lie within the interval
            if let (Some(pr), Some(nr)) = (prev_ref, next_ref) {
                let prev_in = pr >= interval_start && pr < interval_end;
                let next_in = nr >= interval_start && nr < interval_end;
                if prev_in && next_in {
                    if i < query_seq.len() {
                        result.push(query_seq[i]);
                    }
                }
            }
        }
    }

    result
}

/// Compare two sequences for exact match, ignoring N bases
///
/// Returns true if sequences are identical (treating N as wildcard)
pub fn compare_sequences(seq1: &[u8], seq2: &[u8]) -> bool {
    if seq1.len() != seq2.len() {
        return false;
    }

    for (_i, (&base1, &base2)) in seq1.iter().zip(seq2.iter()).enumerate() {
        if base1 != base2 {
            // Allow N (4) to match anything
            if base1 != 4 && base2 != 4 {
                return false;
            }
        }
    }

    true
}

/// Extract haplotype vector from CIGAR operations.
///
/// Returns one entry per reference-consuming base, encoding the read's variants
/// relative to the reference:
/// - `1`: match (`=`, or a lenient `M` — see note)
/// - `-4`: SNV (mismatch, `X`)
/// - `-10`: deletion (`D`)
/// - `base + 10*len` (positive): insertion of `len` bases, summed onto the
///   next reference-consuming base (`11` = match+1I, `6` = SNV+1I)
///
/// # Insertion-marker convention (T2 fix, 2026-06-11)
/// An insertion is recorded by deferring a `pending_ins = len*10` marker to the
/// **next** reference-consuming op, where it is summed onto the base value. This
/// matches the golden haplotype-inspection encoding and preserves SNV+insertion
/// compound bases.
///
/// # Note on `M`
/// Unlike `haplotype_inspection` (which panics on `M` per the golden DivA
/// decision), this crate keeps `M` lenient (treated as a match) for now; under
/// `minimap2 --eqx` no `M` appears, so the two agree on real data. The `M`-policy
/// unification is deferred to the shared `sdrecall-utils` encoder (T0/T4).
pub fn extract_hap_vector(record: &Record) -> Vec<i16> {
    use rust_htslib::bam::record::Cigar;
    let cigar = record.cigar();

    // Pre-size from reference-consuming ops (=, X/M, D, N).
    let ref_len: usize = cigar
        .iter()
        .map(|c| match c {
            Cigar::Match(len)
            | Cigar::Equal(len)
            | Cigar::Diff(len)
            | Cigar::Del(len)
            | Cigar::RefSkip(len) => *len as usize,
            _ => 0,
        })
        .sum();

    let mut hap_vector: Vec<i16> = Vec::with_capacity(ref_len);
    let mut pending_ins: i16 = 0; // 0 = no pending insertion; drained on next ref-consuming op

    for op in cigar.iter() {
        match op {
            // Match (M, lenient) / Equal (=) / RefSkip (N) → match-to-reference (1).
            // N consumes reference like a gap of matches (mirrors haplotype_inspection).
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::RefSkip(len) => {
                let n = *len as usize;
                if pending_ins == 0 {
                    hap_vector.extend(std::iter::repeat_n(1i16, n));
                } else {
                    hap_vector.push(1 + pending_ins);
                    hap_vector.extend(std::iter::repeat_n(1i16, n.saturating_sub(1)));
                    pending_ins = 0;
                }
            }
            // Diff (X) → SNV (-4). A compound insertion+mismatch sums the
            // insertion unit onto the SNV base, producing values ending in 6.
            Cigar::Diff(len) => {
                let n = *len as usize;
                if pending_ins == 0 {
                    hap_vector.extend(std::iter::repeat_n(-4i16, n));
                } else {
                    hap_vector.push(-4 + pending_ins);
                    hap_vector.extend(std::iter::repeat_n(-4i16, n.saturating_sub(1)));
                    pending_ins = 0;
                }
            }
            // Ins (I) → defer the marker to the next ref-consuming op; drop at read start.
            Cigar::Ins(len) => {
                if !hap_vector.is_empty() {
                    pending_ins = (*len as i16) * INDEL_UNIT;
                }
            }
            // Del (D) → deletion (-10).
            Cigar::Del(len) => {
                let n = *len as usize;
                if pending_ins == 0 {
                    hap_vector.extend(std::iter::repeat_n(HAP_DEL, n));
                } else {
                    hap_vector.push(pending_ins);
                    hap_vector.extend(std::iter::repeat_n(HAP_DEL, n.saturating_sub(1)));
                    pending_ins = 0;
                }
            }
            // SoftClip / HardClip / Pad and anything else → no reference contribution.
            _ => {}
        }
    }

    hap_vector
}

/// Count continuous blocks of True values in a boolean array
///
/// Equivalent to Python's count_continuous_blocks function
/// Counts the number of separate continuous regions of True values
fn count_continuous_blocks(bool_array: &[bool]) -> usize {
    if bool_array.is_empty() {
        return 0;
    }

    let mut block_count = 0;
    let mut in_block = false;

    for &is_variant in bool_array {
        if is_variant && !in_block {
            // Starting a new block
            block_count += 1;
            in_block = true;
        } else if !is_variant {
            // Ending a block (if we were in one)
            in_block = false;
        }
    }

    block_count
}

/// Count SNV blocks in a haplotype vector
///
/// Equivalent to Python's count_snv function
fn count_snv_blocks(hap_vector: &[i16]) -> usize {
    // Count SNVs per site, including compound SNV+insertion values.
    hap_vector.iter().filter(|&&val| is_snv_value(val)).count()
}

/// Count indel blocks in a haplotype vector  
///
/// Equivalent to Python's count_continuous_indel_blocks function
/// Counts continuous blocks of deletions (-10) or insertions (>1)
fn count_indel_blocks(hap_vector: &[i16]) -> usize {
    let indel_positions: Vec<bool> = hap_vector.iter().map(|&val| is_indel_value(val)).collect();
    count_continuous_blocks(&indel_positions)
}

/// Count variants in a haplotype vector slice
///
/// Returns (snv_blocks, indel_blocks, total_variant_blocks)
/// Note: This counts continuous blocks, not individual bases
pub fn count_variants(hap_vector: &[i16]) -> (usize, usize, usize) {
    let snv_blocks = count_snv_blocks(hap_vector);
    let indel_blocks = count_indel_blocks(hap_vector);
    let total_blocks = snv_blocks + indel_blocks;

    (snv_blocks, indel_blocks, total_blocks)
}

/// Slice haplotype vector to genomic interval
///
/// Returns the portion of hap_vector corresponding to the genomic interval
pub fn slice_hap_vector(
    hap_vector: &[i16],
    read_start: i64,
    interval_start: i64,
    interval_end: i64,
) -> Vec<i16> {
    let start_offset = (interval_start - read_start) as usize;
    let end_offset = (interval_end - read_start) as usize;

    if start_offset >= hap_vector.len() {
        return Vec::new();
    }

    let actual_end = end_offset.min(hap_vector.len());
    hap_vector[start_offset..actual_end].to_vec()
}

/// Extract error vector from CIGAR operations and base qualities
///
/// **RUST IMPLEMENTATION OF PYTHON'S `get_errorvector_from_cigar`**
///
/// This function creates an error vector representing the probability of sequencing errors
/// at each reference position within the read's alignment. It processes CIGAR operations
/// and base quality scores to compute error probabilities.
///
/// # Key Features:
/// - **CIGAR Processing**: Handles all standard CIGAR operations (M, I, D, S, H, N, =, X)
/// - **Base Quality Conversion**: Converts phred scores to error probabilities using `10^(-qual/10)`
/// - **Indel Handling**: Assigns zero error probability to insertion/deletion positions
/// - **Missing Quality Fallback**: Uses Q20 (1% error) if base qualities are unavailable
///
/// # CIGAR Operation Mapping:
///
/// Operation | Code | Handling
/// ----------|------|----------
/// Match/Equal| 7/8 | Use base quality → error probability  
/// Mismatch  |  8   | Use base quality → error probability
/// Insertion |  1   | Mark previous position, then → 0.0 error
/// Deletion  |  2   | Mark deleted positions → 0.0 error  
/// Soft Clip |  4   | Skip (consume query only)
/// Hard Clip |  5   | Skip (no consumption)
/// RefSkip   |  3   | Mark skipped positions → 0.0 error
///
///
/// # Error Probability Formula:
///
/// error_prob = 10^(-phred_score/10)
///
/// Examples:
/// Q10 → 10% error rate
/// Q20 → 1% error rate  
/// Q30 → 0.1% error rate
///
///
/// # Arguments
/// * `record` - BAM record containing CIGAR string and base qualities
///
/// # Returns  
/// Vector of error probabilities (0.0 to 1.0) for each reference position
///
/// # Python Equivalent
/// ```python
/// def get_errorvector_from_cigar(read, cigar_tuples):
///     # ... (see fp_control/pairwise_read_inspection.py line 265)
/// ```
pub fn extract_error_vector(record: &Record) -> Vec<f32> {
    let ref_length = (record.reference_end() - record.reference_start()) as usize;
    let mut error_vector = vec![0.0; ref_length];

    // Get base qualities (phred scaled)
    let base_qualities = if record.qual().is_empty() {
        // If no qualities, assume reasonable default quality (Q20 = 1% error)
        vec![20u8; record.seq_len()]
    } else {
        record.qual().to_vec()
    };

    let mut query_consume = 0;
    let mut ref_consume = 0;

    // Process each CIGAR operation
    for op in record.cigar().iter() {
        use rust_htslib::bam::record::Cigar;
        match op {
            Cigar::RefSkip(length) => {
                // operation = 3: Skip for reference (mark with placeholder 99)
                for i in 0..*length as usize {
                    if ref_consume + i < error_vector.len() {
                        error_vector[ref_consume + i] = 99.0;
                    }
                }
                ref_consume += *length as usize;
            }
            Cigar::SoftClip(length) => {
                // operation = 4: Soft clipping (consume query only)
                query_consume += *length as usize;
            }
            Cigar::Match(length) | Cigar::Equal(length) => {
                // operation = 7: Match/Equal
                for i in 0..*length as usize {
                    if ref_consume + i < error_vector.len()
                        && query_consume + i < base_qualities.len()
                    {
                        error_vector[ref_consume + i] = base_qualities[query_consume + i] as f32;
                    }
                }
                query_consume += *length as usize;
                ref_consume += *length as usize;
            }
            Cigar::Diff(length) => {
                // operation = 8: Mismatch
                for i in 0..*length as usize {
                    if ref_consume + i < error_vector.len()
                        && query_consume + i < base_qualities.len()
                    {
                        error_vector[ref_consume + i] = base_qualities[query_consume + i] as f32;
                    }
                }
                query_consume += *length as usize;
                ref_consume += *length as usize;
            }
            Cigar::Ins(length) => {
                // operation = 1: Insertion (mark previous reference position with placeholder)
                if ref_consume > 0 && ref_consume - 1 < error_vector.len() {
                    error_vector[ref_consume - 1] = 99.0;
                }
                query_consume += *length as usize;
            }
            Cigar::Del(length) => {
                // operation = 2: Deletion (mark reference positions with placeholder)
                for i in 0..*length as usize {
                    if ref_consume + i < error_vector.len() {
                        error_vector[ref_consume + i] = 99.0;
                    }
                }
                ref_consume += *length as usize;
            }
            Cigar::HardClip(_) => {
                // Hard clipping: no consumption of query or reference
            }
            _ => {
                // Handle any other operations conservatively
                let len = op.len() as usize;
                for i in 0..len {
                    if ref_consume + i < error_vector.len()
                        && query_consume + i < base_qualities.len()
                    {
                        error_vector[ref_consume + i] = base_qualities[query_consume + i] as f32;
                    }
                }
                query_consume += len;
                ref_consume += len;
            }
        }
    }

    // Convert phred scores to error probabilities
    // 99 placeholder → 0.0 (no error for indels)
    // phred score → 10^(-score/10)
    for prob in error_vector.iter_mut() {
        *prob = if *prob == 99.0 {
            0.0 // No error for indels/skips
        } else {
            10.0_f32.powf(-(*prob) / 10.0) // Convert phred to probability
        };
    }

    error_vector
}

/// Get or compute error vector for a read (with caching)
///
/// This function implements caching for error vectors to avoid redundant computation.
/// Error vectors are expensive to compute (CIGAR parsing + quality score conversion),
/// so caching provides significant performance benefits for reads involved in multiple overlaps.
///
/// # Current Usage:
/// - **Caching**: Store computed error vectors to avoid recomputation
/// - **Quality Assessment**: Could be used for dynamic quality filtering  
/// - **Sequencing Error Detection**: Used in `is_sequencing_error()` function
///
/// # Arguments
/// * `record` - BAM record to extract error vector from
/// * `read_error_vectors` - Mutable cache of computed error vectors
///
/// # Returns
/// Vector of error probabilities for this read (cached or newly computed)
///
/// Checks cache first, computes and stores if not found
pub fn get_error_vector(
    record: &Record,
    read_error_vectors: &mut AHashMap<String, ReadErrorVector>,
) -> Result<Vec<f32>, Box<dyn std::error::Error>> {
    let read_id = get_read_id(record)?;

    // Check cache first
    if let Some(cached_vector) = read_error_vectors.get(&read_id) {
        return Ok(cached_vector.clone());
    }

    // Not in cache - compute and store
    let error_vector = extract_error_vector(record);
    read_error_vectors.insert(read_id, error_vector.clone());

    Ok(error_vector)
}

/// Generate unique read ID for caching purposes
///
/// Equivalent to Python's get_read_id function: f"{read.query_name}:{read.flag}"
pub fn get_read_id(record: &Record) -> Result<String, Box<dyn std::error::Error>> {
    let qname = std::str::from_utf8(record.qname())?;
    Ok(format!("{}:{}", qname, record.flags()))
}

/// Get or compute haplotype vector for a read (with caching)
///
/// Checks cache first, computes and stores if not found
pub fn get_hap_vector(
    record: &Record,
    read_hap_vectors: &mut AHashMap<String, ReadHaplotypeVector>,
) -> Result<Vec<i16>, Box<dyn std::error::Error>> {
    let read_id = get_read_id(record)?;

    // Check cache first
    if let Some(cached_vector) = read_hap_vectors.get(&read_id) {
        return Ok(cached_vector.clone());
    }

    // Not in cache - compute and store
    let hap_vector = extract_hap_vector(record);
    read_hap_vectors.insert(read_id, hap_vector.clone());

    Ok(hap_vector)
}

/// Check if a base at a specific position is likely a sequencing error
///
/// **OPTIMIZED VERSION**: Uses BamRecordExtensions::read_pos for O(1) coordinate lookup
/// instead of O(n) linear search through vectors
///
/// Equivalent to Python's seq_err_det_stacked_bases function
/// Uses pre-computed error vector and allele depth to determine sequencing errors
fn is_sequencing_error(
    record: &Record,
    error_vector: &[f32],
    genomic_pos: i64,
    allele_depth_map: &AlleleDepthMap,
    chrom: &str,
) -> bool {
    let qname = std::str::from_utf8(record.qname()).unwrap_or("unknown");
    debug!(
        "[is_sequencing_error] is_sequencing_error: Checking read {} at position {}",
        qname, genomic_pos
    );

    // **COORDINATE LOOKUP**: Use CIGAR directly for position mapping
    // record.cigar() already returns a CigarStringView
    let cigar_view = record.cigar();
    match cigar_view.read_pos(genomic_pos as u32, true, true) {
        Ok(Some(query_pos)) => {
            let qi = query_pos as usize;
            let query_seq = record.seq().as_bytes();

            if qi >= query_seq.len() {
                warn!(
                    "[is_sequencing_error] Query index {} out of bounds (seq_len={})",
                    qi,
                    query_seq.len()
                );
                return false;
            }

            // Use pre-computed error probability (indexed by reference position)
            let read_offset = (genomic_pos - record.reference_start()) as usize;
            if read_offset >= error_vector.len() {
                warn!(
                    "[is_sequencing_error] Position {} out of bounds in error vector (len={})",
                    genomic_pos,
                    error_vector.len()
                );
                return false;
            }

            let error_prob = error_vector[read_offset];
            debug!(
                "[is_sequencing_error] Error probability at position {}: {:.6}",
                genomic_pos, error_prob
            );

            // Convert error probability back to approximate phred score for threshold check
            // error_prob = 10^(-qual/10) → qual = -10 * log10(error_prob)
            let approx_qual = if error_prob > 0.0 {
                (-10.0 * error_prob.log10()) as u8
            } else {
                40 // High quality if error_prob is 0
            };

            debug!(
                "[is_sequencing_error] Approximate quality score: Q{} (from error_prob={:.6})",
                approx_qual, error_prob
            );

            // If base quality >= 20, not a sequencing error
            if approx_qual >= 20 {
                debug!(
                    "[is_sequencing_error] High quality (Q{} >= 20) -> NOT a sequencing error",
                    approx_qual
                );
                return false;
            }

            debug!(
                "[is_sequencing_error] Low quality (Q{} < 20) -> checking allele depth...",
                approx_qual
            );

            // Get allele depth information for this position
            // Use the chromosome name passed from the parent function (already resolved correctly)
            debug!(
                "[is_sequencing_error] Using chromosome: '{}' for position {}",
                chrom, genomic_pos
            );
            if let Some(pos_data) = allele_depth_map.get(&chrom, genomic_pos as u32) {
                let target_base = query_seq[qi];
                // query_seq is raw ASCII (record.seq().as_bytes()), so encode it to the
                // A=0,T=1,C=2,G=3,N=4 allele index — via the SAME base_to_index that
                // build_allele_depth_map used to fill the array — before indexing it.
                // The old `target_base as usize` passed the ASCII code (65..=84, always >=5),
                // so get_allele_depth always returned 0 → af always 0.0 → the allele-frequency
                // test was a no-op. Python seq_err_det_stacked_bases keys the dict by the
                // encoded base, so this restores parity.
                let allele_idx = crate::bam_reading::base_to_index(target_base as char);
                let ad = AlleleDepthMap::get_allele_depth(pos_data, allele_idx);
                let dp = AlleleDepthMap::total_depth(pos_data);

                debug!("[is_sequencing_error] Base at position {}: '{}' -> allele_idx={}, AD={}, DP={}",
                       genomic_pos, target_base as char, allele_idx, ad, dp);

                if dp == 0 {
                    debug!("[is_sequencing_error] Zero depth -> NOT a sequencing error");
                    return false;
                }

                let af = ad as f32 / dp as f32;
                debug!(
                    "[is_sequencing_error] Allele frequency: {:.4} ({}/{})",
                    af, ad, dp
                );

                // Python criteria: (af <= 0.02 or (ad == 1 and dp >= 10)) and base_qual < 13
                let af_criteria = af <= 0.02 || (ad == 1 && dp >= 10);
                let qual_criteria = approx_qual < 13;

                debug!("[is_sequencing_error] Criteria check: AF_criteria={} (af={:.4} <= 0.02 OR (ad={} == 1 AND dp={} >= 10))", 
                       af_criteria, af, ad, dp);
                debug!(
                    "[is_sequencing_error] Criteria check: QUAL_criteria={} (Q{} < 13)",
                    qual_criteria, approx_qual
                );

                let is_error = af_criteria && qual_criteria;
                debug!("[is_sequencing_error] Final decision: {} (AF_criteria={} AND QUAL_criteria={})", 
                       if is_error { "SEQUENCING ERROR" } else { "REAL VARIANT" }, 
                       af_criteria, qual_criteria);

                return is_error;
            } else {
                // No pileup entry at this position → no evidence of a low-AF artifact.
                // Python seq_err_det_stacked_bases hits dp == 0 here and returns False
                // (NOT a sequencing error → the mismatch is treated as a real variant).
                // The old `return true` diverged, over-tolerating mismatches and merging
                // reads Python keeps on separate haplotypes.
                debug!("[is_sequencing_error] No allele depth data for {}:{} -> treat as REAL variant (Python dp==0 -> False)", chrom, genomic_pos);
                return false;
            }
        }
        Ok(None) => {
            // Position is outside the read alignment (before start or after end)
            warn!("[is_sequencing_error] Position {} outside read alignment -> NOT a sequencing error", genomic_pos);
        }
        Err(e) => {
            // Error in read_pos function call
            error!("[is_sequencing_error] read_pos failed for position {}: {} -> NOT a sequencing error", genomic_pos, e);
        }
    }

    false
}

/// **ALGORITHM DESIGN: Two-Stage Mismatch Analysis Sidesteps bcftools Indel Issues**
///
/// This implementation uses a clever two-stage approach that naturally avoids
/// the complexities of bcftools indel representation:
///
/// **STAGE 1: Indel Detection & Immediate Rejection**
/// - Uses haplotype vectors to detect indel mismatches
/// - Indel mismatches → immediate rejection (different haplotypes)
/// - No allele depth queries needed for indels
///
/// **STAGE 2: SNV Sequencing Error Analysis**
/// - Only reached when NO indels are detected
/// - Uses allele depth map for SNV-specific sequencing error detection
/// - bcftools represents SNVs correctly at exact genomic positions
///
/// **Why This Works:**
/// - Indel complexity (anchor positions, complex AD) is bypassed entirely
/// - Allele depth map only used for SNVs where bcftools representation is straightforward
/// - Clean separation of concerns: hap vectors for structural analysis, AD for error rates
///
/// **Result:** The bcftools indel representation issues discussed above are
/// largely irrelevant to this specific algorithm design!

/// Check if mismatches between two reads can be tolerated as sequencing errors
///
/// **UPDATED APPROACH**: Uses haplotype vectors to identify mismatch types
/// but still requires query sequences for sequencing error detection via allele depths
///
/// # Arguments
/// * `hap_vec1` - Sliced haplotype vector for read 1 (interval only)
/// * `hap_vec2` - Sliced haplotype vector for read 2 (interval only)
/// * `seq1` - Query sequence for read 1 (full read, needed for allele lookup)
/// * `ref_pos1` - Reference positions for read 1 (full read)
/// * `seq2` - Query sequence for read 2 (full read, needed for allele lookup)  
/// * `ref_pos2` - Reference positions for read 2 (full read)
/// * `record1` - BAM record for read 1
/// * `error_vec1` - Pre-computed error vector for read 1
/// * `record2` - BAM record for read 2  
/// * `error_vec2` - Pre-computed error vector for read 2
/// * `mismatch_positions` - Genomic positions where haplotype vectors differ
/// * `allele_depth_map` - Allele depth information for sequencing error detection
///
/// # Returns
/// (is_tolerable, tolerated_count) where is_tolerable indicates if all mismatches
/// can be explained as sequencing errors, and tolerated_count is the number tolerated
fn tolerate_mismatches_from_hap_vectors(
    hap_vec1: &[i16],
    hap_vec2: &[i16],
    _seq1: &[u8],
    _ref_pos1: &[i64],
    _seq2: &[u8],
    _ref_pos2: &[i64],
    record1: &Record,
    error_vec1: &[f32],
    record2: &Record,
    error_vec2: &[f32],
    mismatch_positions: &[i64],
    interval_start: i64,
    allele_depth_map: &AlleleDepthMap,
    chrom: &str,
) -> (bool, usize) {
    let qname1 = std::str::from_utf8(record1.qname()).unwrap_or("unknown");
    let qname2 = std::str::from_utf8(record2.qname()).unwrap_or("unknown");
    // Append sam flag of the read to the qname to make a unique read id
    let read1_id = format!("{}_{}", qname1, record1.flags());
    let read2_id = format!("{}_{}", qname2, record2.flags());

    debug!(
        "[tolerate_mismatches_from_hap_vectors] Analyzing {} mismatches between reads {} and {}",
        mismatch_positions.len(),
        read1_id,
        read2_id
    );
    debug!(
        "[tolerate_mismatches_from_hap_vectors] Mismatch positions: {:?}",
        mismatch_positions
    );

    let mut tolerable_count = 0;

    for (idx, &genomic_pos) in mismatch_positions.iter().enumerate() {
        let hap_index = (genomic_pos - interval_start) as usize;

        // Get haplotype values at this position
        let hap1 = if hap_index < hap_vec1.len() {
            hap_vec1[hap_index]
        } else {
            1
        };
        let hap2 = if hap_index < hap_vec2.len() {
            hap_vec2[hap_index]
        } else {
            1
        };

        debug!(
            "[tolerate_mismatches_from_hap_vectors] Mismatch #{} at position {} (hap1={}, hap2={})",
            idx + 1,
            genomic_pos,
            hap1,
            hap2
        );

        // Check if this is an indel mismatch (not tolerable)
        let is_indel1 = is_indel_value(hap1);
        let is_indel2 = is_indel_value(hap2);
        let is_snv1 = is_snv_value(hap1);
        let is_snv2 = is_snv_value(hap2);

        if is_indel1 || is_indel2 {
            debug!("[tolerate_mismatches_from_hap_vectors] Position {} has indel mismatch (hap1={}, hap2={}) - NOT TOLERABLE", 
                   genomic_pos, hap1, hap2);
            debug!("[tolerate_mismatches_from_hap_vectors] tolerate_mismatches: FAILED at indel position {} - returning (false, 0)", genomic_pos);
            return (false, 0);
        }

        // For SNV mismatches, check if they can be explained by sequencing errors
        let error1 = is_sequencing_error(record1, error_vec1, genomic_pos, allele_depth_map, chrom);
        let error2 = is_sequencing_error(record2, error_vec2, genomic_pos, allele_depth_map, chrom);

        debug!("[tolerate_mismatches_from_hap_vectors] Position {}: read1({}) error={}, read2({}) error={}", 
               genomic_pos, read1_id, error1, read2_id, error2);

        // If either read has a sequencing error at this position, we can tolerate it
        if (error1 && is_snv1) || (error2 && is_snv2) {
            tolerable_count += 1;
            debug!("[tolerate_mismatches_from_hap_vectors] Position {} TOLERABLE (sequencing error detected)", genomic_pos);
        } else {
            // This mismatch cannot be explained by sequencing error
            debug!("[tolerate_mismatches_from_hap_vectors] Position {} NOT TOLERABLE (likely real variant)", genomic_pos);
            debug!("[tolerate_mismatches_from_hap_vectors] tolerate_mismatches: FAILED at position {} - returning (false, 0)", genomic_pos);
            return (false, 0);
        }
    }

    // All mismatches can be explained by sequencing errors
    debug!("[tolerate_mismatches_from_hap_vectors] tolerate_mismatches: SUCCESS - all {} mismatches tolerable (sequencing errors)", 
           mismatch_positions.len());
    debug!(
        "[tolerate_mismatches_from_hap_vectors] Total tolerated mismatches: {}",
        tolerable_count
    );
    (true, tolerable_count)
}

/// Find genomic positions where two haplotype vectors differ
///
/// **CORRECT APPROACH**: Uses haplotype vectors instead of query sequences
/// Haplotype vectors have 1-to-1 correspondence with reference genome positions
/// and properly represent alignment status (match, SNV, indel) at each position.
///
/// # Arguments
/// * `hap_vec1` - Sliced haplotype vector for read 1 (interval only)
/// * `hap_vec2` - Sliced haplotype vector for read 2 (interval only)  
/// * `interval_start` - Start of genomic interval
/// * `interval_end` - End of genomic interval
///
/// # Returns
/// Vector of genomic positions where reads have different haplotype values
fn find_mismatch_positions_from_hap_vectors(
    hap_vec1: &[i16],
    hap_vec2: &[i16],
    interval_start: i64,
    interval_end: i64,
) -> Vec<i64> {
    let mut mismatches = Vec::new();
    let interval_length = (interval_end - interval_start) as usize;

    // Ensure both vectors have the same length and cover the interval
    let min_len = hap_vec1.len().min(hap_vec2.len()).min(interval_length);

    // Compare haplotype values position by position
    for i in 0..min_len {
        let hap1 = hap_vec1[i];
        let hap2 = hap_vec2[i];

        // If haplotype values differ, it's an SNV
        if hap1 != hap2 {
            let genomic_pos = interval_start + i as i64;
            mismatches.push(genomic_pos);
        }
    }

    mismatches
}

/// Check if mismatches occur at indel positions using haplotype vectors (which are not tolerable)
///
/// **CORRECT APPROACH**: Uses haplotype vector values to detect indel mismatches
///
/// Haplotype vector values:
/// - 1: Match to reference
/// - -4: SNV (mismatch)
/// - -10: Deletion
/// - Positive values: insertion compound (`base + 10*len`)
///
/// # Arguments
/// * `hap_vec1` - Sliced haplotype vector for read 1 (interval only)
/// * `hap_vec2` - Sliced haplotype vector for read 2 (interval only)
/// * `mismatch_positions` - Genomic positions where reads differ
/// * `interval_start` - Start of genomic interval
///
/// # Returns
/// True if any mismatch involves indels (not tolerable), false if all SNV mismatches
fn has_indel_mismatches_from_hap_vectors(
    hap_vec1: &[i16],
    hap_vec2: &[i16],
    mismatch_positions: &[i64],
    interval_start: i64,
) -> bool {
    for &genomic_pos in mismatch_positions {
        let index = (genomic_pos - interval_start) as usize;

        if index >= hap_vec1.len() || index >= hap_vec2.len() {
            continue;
        }

        let hap1 = hap_vec1[index];
        let hap2 = hap_vec2[index];

        // Check if either value indicates an indel
        let is_indel1 = is_indel_value(hap1);
        let is_indel2 = is_indel_value(hap2);

        if is_indel1 || is_indel2 {
            debug!("[has_indel_mismatches_from_hap_vectors] Found indel mismatch at position {}: hap1={}, hap2={}", 
                   genomic_pos, hap1, hap2);
            return true;
        }
    }
    false
}

/// Determine if two reads come from the same haplotype
///
/// Main function that orchestrates the haplotype comparison workflow
///
/// **OPTIMIZED ALGORITHM FLOW:**
/// ```
/// 1. Extract query sequences from BAM records
/// 2. Slice sequences to genomic interval  
/// 3. Compute haplotype vectors (needed for weight calculation in both paths)
/// 4. Compare sequences for exact match
///    ├─ IDENTICAL → Calculate weight using haplotype vectors
///    └─ DIFFERENT → Mismatch analysis
///       ├─ Compute error vectors (only when needed)
///       ├─ Check sequencing error tolerance
///       └─ If tolerable → Calculate weight using pre-computed haplotype vectors
/// ```
///
/// **KEY OPTIMIZATIONS:**
/// - Haplotype vectors computed once upfront (needed for weight in both paths)
/// - Error vectors computed only for mismatch analysis (when sequences differ)
/// - Efficient caching for both vector types to avoid redundant computation
pub fn determine_same_haplotype(
    read1: &Record,
    read2: &Record,
    start: i64,
    end: i64,
    chrom: &str,
    allele_depth_map: &AlleleDepthMap,
    intrinsic_ad_map: &AlleleDepthMap,
    config: &HaplotypeConfig,
    read_hap_vectors: &mut AHashMap<String, ReadHaplotypeVector>,
    read_error_vectors: &mut AHashMap<String, ReadErrorVector>,
    read_ref_pos_dict: &mut AHashMap<String, (i64, i64)>,
) -> Result<(HaplotypeResult, Option<f32>), Box<dyn std::error::Error>> {
    // `allele_depth_map`  : MAIN-bam AD (Python `nested_ad_dict`) → sequencing-error
    //                       tolerance check (`tolerate_mismatches_*`).
    // `intrinsic_ad_map`  : INTRINSIC-bam AD (Python `intrinsic_ad_dict`) → PSV
    //                       detection inside the edge-weight formula
    //                       (`psv_shared_snvs`). May be empty (standalone phaser /
    //                       Python's `intrinsic_ad_dict = {}`) → psv_snv_count == 0.
    // Extract qnames from records
    let qname1 = std::str::from_utf8(read1.qname())?;
    let qname2 = std::str::from_utf8(read2.qname())?;
    let read1_id = format!("{}_{}", qname1, read1.flags());
    let read2_id = format!("{}_{}", qname2, read2.flags());

    debug!(
        "[determine_same_haplotype] Comparing reads {} and {} at {}:{}-{}",
        read1_id, read2_id, chrom, start, end
    );

    // Step 1: Extract query sequences and reference positions
    let (seq1, ref_pos1) = extract_query_seq(read1)?;
    let (seq2, ref_pos2) = extract_query_seq(read2)?;

    // Step 2: Slice sequences to the genomic interval
    let interval_seq1 = slice_seq_to_interval(&seq1, &ref_pos1, start, end);
    let interval_seq2 = slice_seq_to_interval(&seq2, &ref_pos2, start, end);

    // CRITICAL: Check for empty sequences BEFORE comparison
    // Empty sequences can occur when the overlap interval falls in soft-clipped regions
    // Two empty sequences would incorrectly be treated as "identical" by compare_sequences
    if interval_seq1.is_empty() || interval_seq2.is_empty() {
        warn!("[determine_same_haplotype] Empty sliced sequence(s) detected: read1_len={}, read2_len={} for interval {}:{}-{}. \
               This may indicate the overlap interval falls in soft-clipped regions. Returning UNKNOWN.",
               interval_seq1.len(), interval_seq2.len(), chrom, start, end);
        return Ok((HaplotypeResult::Unknown, None));
    }

    // Step 3: Get haplotype vectors (needed for weight calculation in both paths)
    // OPTIMIZATION: Compute once, use in both identical and tolerable mismatch cases
    let hap_vec1 = get_hap_vector(read1, read_hap_vectors)?;
    let hap_vec2 = get_hap_vector(read2, read_hap_vectors)?;

    // Slice both haplotype vectors to the overlap interval ONCE. The edge-weight
    // formula (overlap_span / indel_num / shared SNVs) and the mismatch analysis
    // in BOTH branches operate on these interval-restricted vectors, so computing
    // them here avoids the previous triple re-slicing.
    let interval_hap1 = slice_hap_vector(&hap_vec1, read1.pos(), start, end);
    let interval_hap2 = slice_hap_vector(&hap_vec2, read2.pos(), start, end);

    // Step 4: Compare sequences for exact match
    if compare_sequences(&interval_seq1, &interval_seq2) {
        // SEQUENCES ARE IDENTICAL: Same haplotype (likely)
        //
        // PERFORMANCE OPTIMIZATION:
        // No error vector computation needed when sequences match exactly.
        // Error vectors are only used for sequencing error detection in mismatches.
        // This saves expensive CIGAR parsing + quality score conversion.

        debug!(
            "[determine_same_haplotype] Sequences IDENTICAL between reads {} and {}",
            read1_id, read2_id
        );
        debug!(
            "[determine_same_haplotype] Interval: {}:{}-{} (length={})",
            chrom,
            start,
            end,
            end - start
        );
        debug!("[determine_same_haplotype] Within this interval, the query seq for {} is {:?}, the query seq for {} is {:?}", read1_id, interval_seq1, read2_id, interval_seq2);

        // Sequences are identical - likely same haplotype
        //
        // Edge weight — faithful port of Python `determine_same_haplotype`
        // (fp_control/pairwise_read_inspection.py lines 689-697):
        //   weight = overlap_span + sum(score_arr[:shared_snv_count])
        //          + mean_read_length * (shared_snv_count - psv_snv_count) * 0.75
        //          + mean_read_length * 3 * indel_num
        // where overlap_span / indel_num are taken over the EQUAL positions of the
        // two interval hap-vectors and shared_snv_count is the base-agnostic count
        // of co-located SNVs (see `compute_edge_weight_base`).
        let weight = compute_edge_weight_base(
            &interval_hap1,
            &interval_hap2,
            start,
            &seq1,
            &ref_pos1,
            &seq2,
            &ref_pos2,
            chrom,
            intrinsic_ad_map,
            config,
        );

        // Normalize weight (Python applies `weight/(mean_read_length*10)` in graph_build).
        let normalized_weight = weight / (config.mean_read_length * 10.0);

        debug!(
            "[determine_same_haplotype] Final weight: {:.6} (normalized from {:.2})",
            normalized_weight, weight
        );

        // No need for error vectors when sequences are identical
        // Error vectors are only used for sequencing error detection in mismatches

        // Store read positions
        read_ref_pos_dict.insert(qname1.to_string(), (start, end));
        read_ref_pos_dict.insert(qname2.to_string(), (start, end));

        debug!(
            "[determine_same_haplotype] SAME HAPLOTYPE (identical) - weight: {:.6}",
            normalized_weight
        );

        Ok((HaplotypeResult::Same, Some(normalized_weight)))
    } else {
        // Sequences differ - need detailed analysis

        // Check if sequences have different lengths (likely indel differences)
        if interval_seq1.len() != interval_seq2.len() {
            debug!("[determine_same_haplotype] Sequence length mismatch: read1_len={}, read2_len={} -> DIFFERENT haplotypes", 
                   interval_seq1.len(), interval_seq2.len());
            return Ok((HaplotypeResult::Different, None));
        }

        // Note: Empty sequence check is now done before compare_sequences() call
        // This branch is only reached when both sequences are non-empty but differ
        // (interval_hap1 / interval_hap2 were sliced once near the top).

        // Find mismatch positions using haplotype vectors (CORRECT APPROACH)
        let mut mismatch_positions =
            find_mismatch_positions_from_hap_vectors(&interval_hap1, &interval_hap2, start, end);

        // Check if mismatches occur at indel positions - these are not tolerable
        if has_indel_mismatches_from_hap_vectors(
            &interval_hap1,
            &interval_hap2,
            &mismatch_positions,
            start,
        ) {
            debug!("[determine_same_haplotype] Found indel mismatches between reads, marking as different haplotypes");
            return Ok((HaplotypeResult::Different, None));
        }

        // Further check the shared SNV positions to see if they share the same ALT allele.
        // `matching_shared_snv_pos` is no longer the weight's shared-SNV source (the
        // weight now uses Python's base-agnostic count via `psv_shared_snvs`); we keep
        // this call only for `discrepant_shared_snv_pos` (in-trans variant detection).
        let (_matching_shared_snv_pos, discrepant_shared_snv_pos) = stat_shared_snv_matches(
            &interval_hap1,
            &interval_hap2,
            start,
            &seq1,
            &ref_pos1,
            &seq2,
            &ref_pos2,
            read1,
            read2,
        )?;

        // If no mismatches found (shouldn't happen since sequences differ), treat as unknown
        if mismatch_positions.is_empty() && discrepant_shared_snv_pos.is_empty() {
            warn!("[determine_same_haplotype] No mismatches found despite sequence differences -> Different Haplotypes for conservative estimation, interval_seq1={:?}, interval_seq2={:?}, interval_hap1={:?}, interval_hap2={:?}. The different hap genomic positions are: {:?}. The shared mismatch positions with different ALT alleles are: {:?}", interval_seq1, interval_seq2, interval_hap1, interval_hap2, mismatch_positions, discrepant_shared_snv_pos);
            return Ok((HaplotypeResult::Different, None));
        }

        if !discrepant_shared_snv_pos.is_empty() {
            debug!("[determine_same_haplotype] Found discrepant shared SNV positions: {:?}, the interval is {}:{}-{}, the two reads compared are {} and {}. Within this overlap interval, their hap vectors are {:?} and {:?}, their query seq are {:?} and {:?}", discrepant_shared_snv_pos, chrom, start, end, read1_id, read2_id, interval_hap1, interval_hap2, interval_seq1, interval_seq2);
        }

        // Merge discrepant_shared_snv_pos into mismatch_positions for downstream analysis
        debug!("[determine_same_haplotype] The different hap genomic positions are: {:?}. The shared mismatch positions with different ALT alleles are: {:?}", mismatch_positions, discrepant_shared_snv_pos);
        mismatch_positions.extend(discrepant_shared_snv_pos);

        if mismatch_positions.len() >= 3 {
            debug!("[determine_same_haplotype] Found {} mismatches, conservatively treat them as from different haplotypes. The interval is {}:{}-{}, the two reads compared are {} and {}. Within this overlap interval, their hap vectors are {:?} and {:?}, their query seq are {:?} and {:?}", mismatch_positions.len(), chrom, start, end, read1_id, read2_id, interval_hap1, interval_hap2, interval_seq1, interval_seq2);
            return Ok((HaplotypeResult::Different, None));
        }

        // SEQUENCES DIFFER: Mismatch analysis required
        //
        // COMPUTE ERROR VECTORS (only when needed for mismatch analysis):
        // Error vectors are expensive to compute (CIGAR parsing + quality conversion)
        // but essential for distinguishing real variants from sequencing errors.
        //
        // Each error vector contains error probabilities for every reference position:
        // - Low error prob (< 0.03) -> likely real variant
        // - High error prob (> 0.03) -> likely sequencing error
        let error_vec1 = get_error_vector(read1, read_error_vectors)?;
        let error_vec2 = get_error_vector(read2, read_error_vectors)?;

        debug!(
            "[determine_same_haplotype] Sequence mismatch detected between reads {} and {}",
            read1_id, read2_id
        );
        debug!(
            "[determine_same_haplotype] Interval: {}:{}-{} (length={})",
            chrom,
            start,
            end,
            end - start
        );
        debug!(
            "[determine_same_haplotype] Sequences: read1_len={}, read2_len={}",
            interval_seq1.len(),
            interval_seq2.len()
        );

        // Check if all mismatches can be explained by sequencing errors
        // Using haplotype vectors to properly identify mismatch types (CORRECT APPROACH)
        let (tolerable, tolerated_count) = tolerate_mismatches_from_hap_vectors(
            &interval_hap1,
            &interval_hap2,
            &seq1,
            &ref_pos1,
            &seq2,
            &ref_pos2,
            read1,
            &error_vec1,
            read2,
            &error_vec2,
            &mismatch_positions,
            start,
            allele_depth_map,
            chrom,
        );

        debug!("[determine_same_haplotype] Mismatch tolerance result: tolerable={}, tolerated_count={}", tolerable, tolerated_count);

        if tolerable {
            // All mismatches are sequencing errors - treat as same haplotype with penalty

            debug!("[determine_same_haplotype] All mismatches tolerable as sequencing errors - calculating weight with penalty");

            // Edge weight: identical Python formula as the IDENTICAL path
            // (lines 689-697), then the tolerated-mismatch penalty
            // `weight = weight - tolerated_count * 20` (Python line 738), clamped at 0.
            let base_weight = compute_edge_weight_base(
                &interval_hap1,
                &interval_hap2,
                start,
                &seq1,
                &ref_pos1,
                &seq2,
                &ref_pos2,
                chrom,
                intrinsic_ad_map,
                config,
            );

            let weight = (base_weight - tolerated_count as f32 * 20.0).max(0.0);

            debug!(
                "[determine_same_haplotype] Penalty applied: base={:.2} - {} x 20.0 -> {:.2}",
                base_weight, tolerated_count, weight
            );

            // Normalize weight
            let normalized_weight = weight / (config.mean_read_length * 10.0);

            debug!(
                "[determine_same_haplotype] Final weight: {:.6} (normalized from {:.2})",
                normalized_weight, weight
            );

            // Update data structures - error vectors already computed above

            read_ref_pos_dict.insert(qname1.to_string(), (start, end));
            read_ref_pos_dict.insert(qname2.to_string(), (start, end));

            debug!(
                "[determine_same_haplotype] SAME HAPLOTYPE (with penalty) - weight: {:.6}",
                normalized_weight
            );

            Ok((HaplotypeResult::Same, Some(normalized_weight)))
        } else {
            // Some mismatches cannot be explained by sequencing errors
            debug!("[determine_same_haplotype] Some mismatches are real variants - marking as DIFFERENT haplotypes");
            debug!(
                "[determine_same_haplotype] DIFFERENT HAPLOTYPES - {} untolerable mismatches",
                mismatch_positions.len()
            );
            Ok((HaplotypeResult::Different, None))
        }
    }
}

/// Raw (un-normalized, pre-penalty) edge weight between two reads in their overlap
/// interval — a faithful port of the weight block of Python `determine_same_haplotype`
/// (fp_control/pairwise_read_inspection.py lines 689-697):
///
/// ```text
/// identical_idx  = positions where interval_hap1 == interval_hap2   (elementwise)
/// identical_part = interval_hap1[identical_idx]
/// overlap_span   = identical_part.size                  # count of EQUAL positions
/// indel_num      = count_continuous_indel_blocks(identical_part)
/// (psv, shared)  = psv_shared_snvs(...)
/// weight = overlap_span + sum(score_arr[:shared])
///        + mean_read_length * (shared - psv) * 0.75
///        + mean_read_length * 3 * indel_num
/// ```
///
/// The three subtleties this encodes (vs. the previous Rust approximation):
/// - the base is `overlap_span` = the number of **EQUAL** positions, NOT the full
///   `end - start` span;
/// - `indel_num` is counted over the **EQUAL** positions (`identical_part`) only;
/// - `shared` is the **base-agnostic** count of co-located SNVs (both hap-vectors
///   `-4`), not `min(snv1, snv2)` nor the alt-base-matching subset.
#[allow(clippy::too_many_arguments)]
fn compute_edge_weight_base(
    interval_hap1: &[i16],
    interval_hap2: &[i16],
    interval_start: i64,
    seq1: &[u8],
    ref_pos1: &[i64],
    seq2: &[u8],
    ref_pos2: &[i64],
    chrom: &str,
    intrinsic_ad_map: &AlleleDepthMap,
    config: &HaplotypeConfig,
) -> f32 {
    let min_len = interval_hap1.len().min(interval_hap2.len());

    // identical_part: the hap values at positions where the two vectors are EQUAL.
    let mut identical_part: Vec<i16> = Vec::with_capacity(min_len);
    for i in 0..min_len {
        if interval_hap1[i] == interval_hap2[i] {
            identical_part.push(interval_hap1[i]);
        }
    }
    let overlap_span = identical_part.len();
    let indel_num = count_indel_blocks(&identical_part);

    let (psv_snv_count, shared_snv_count) = psv_shared_snvs(
        interval_hap1,
        interval_hap2,
        interval_start,
        seq1,
        ref_pos1,
        seq2,
        ref_pos2,
        chrom,
        intrinsic_ad_map,
    );

    let mut weight = overlap_span as f32;
    // numba_sum(score_arr[:shared_snv_count]) — Python slices past the end safely,
    // so cap at the array length.
    let take = shared_snv_count.min(config.score_array.len());
    for s in config.score_array.iter().take(take) {
        weight += *s;
    }
    // Non-PSV shared SNVs (novel variants both reads carry, unsupported by the
    // intrinsic-bam AD) boost the weight; PSVs (paralog-expected) do not.
    weight += config.mean_read_length * (shared_snv_count as f32 - psv_snv_count as f32) * 0.75;
    weight += config.mean_read_length * 3.0 * indel_num as f32;

    debug!(
        "[compute_edge_weight_base] overlap_span={}, shared_snv={}, psv_snv={}, indel_num={} -> raw weight={:.2}",
        overlap_span, shared_snv_count, psv_snv_count, indel_num, weight
    );

    weight
}

/// Faithful port of Python `psv_shared_snvs`
/// (fp_control/pairwise_read_inspection.py lines 475-559).
///
/// Returns `(psv_snv_count, shared_snv_count)`:
/// - `shared_snv_count`: number of positions where BOTH interval hap-vectors are
///   `-4` (a shared SNV site), **base-agnostic** (Python `numba_find_shared_snvs`).
/// - `psv_snv_count`: of those, the count that are *paralogous sequence variants* —
///   both reads carry the SAME non-N alt base AND that alt allele is supported
///   (allele depth > 0) at this position in the **INTRINSIC-bam** allele-depth map.
///
/// `intrinsic_ad_map` mirrors Python's `intrinsic_ad_dict` (built by `stat_ad_to_dict`
/// over the intrinsic BAM). When empty (standalone phaser, or Python's
/// `intrinsic_ad_dict = {}`), `psv_snv_count` is always 0 and every shared SNV gets
/// the full `mean_read_length * 0.75` boost.
///
/// Representational note: Python's inner dict only stores ALT (and DP) keys, never
/// the REF base, whereas the Rust `PositionAlleleDepth` array also carries the REF
/// slot. This never diverges because `alt1` here is the read's base at a `-4` (SNV)
/// site — i.e. a mismatch to the reference — so it is never the REF allele; we only
/// ever query a non-REF slot.
#[allow(clippy::too_many_arguments)]
fn psv_shared_snvs(
    interval_hap1: &[i16],
    interval_hap2: &[i16],
    interval_start: i64,
    seq1: &[u8],
    ref_pos1: &[i64],
    seq2: &[u8],
    ref_pos2: &[i64],
    chrom: &str,
    intrinsic_ad_map: &AlleleDepthMap,
) -> (usize, usize) {
    let min_len = interval_hap1.len().min(interval_hap2.len());

    let mut shared_snv_count = 0usize;
    let mut psv_snv_count = 0usize;

    for i in 0..min_len {
        // Shared SNV site iff both hap-vectors carry an SNV signal
        // (including SNV+insertion compound values).
        if !is_snv_value(interval_hap1[i]) || !is_snv_value(interval_hap2[i]) {
            continue;
        }
        shared_snv_count += 1;

        let genomic_pos = interval_start + i as i64;
        if genomic_pos < 0 {
            continue;
        }

        // Read base at this position in each read (encoded A=0,T=1,C=2,G=3,N=4).
        let (alt1, alt2) = match (
            get_base_at_position(seq1, ref_pos1, genomic_pos),
            get_base_at_position(seq2, ref_pos2, genomic_pos),
        ) {
            (Some(a), Some(b)) => (a, b),
            // Position falls in a deletion / is not covered in one read → skip,
            // mirroring Python skipping when get_interval_seq returns empty.
            _ => continue,
        };

        // Both reads must carry the same, non-N alt base (Python: alt1 == alt2 != N).
        if alt1 != alt2 || alt1 == 4 {
            continue;
        }

        // PSV iff this alt allele is supported (AD > 0) in the intrinsic-bam map.
        if let Some(pos_data) = intrinsic_ad_map.get(chrom, genomic_pos as u32) {
            if AlleleDepthMap::get_allele_depth(pos_data, alt1 as usize) > 0 {
                psv_snv_count += 1;
            }
        }
    }

    (psv_snv_count, shared_snv_count)
}

/// Extract positions where two overlapping reads share the same SNV and require that the
/// alternate bases are exactly the same in both reads. No allele depth data is used here.
///
/// This is a Rust-specific in-trans-variant detector — distinct from the actual
/// `psv_shared_snvs` port above. Its `discrepant_shared_snv_positions` output feeds
/// the mismatch analysis (shared SNV sites with *different* alt bases = in-trans),
/// while the edge weight's shared-SNV count comes from `psv_shared_snvs` instead.
///
/// Returns (matching_shared_snv_positions, discrepant_shared_snv_positions) where:
/// - matching_shared_snv_positions: genomic positions where both reads have SNVs with identical alt bases
/// - discrepant_shared_snv_positions: genomic positions where both reads have SNVs but with different alt bases
pub fn stat_shared_snv_matches(
    interval_hap1: &[i16],
    interval_hap2: &[i16],
    interval_start: i64,
    seq1: &[u8],
    ref_pos1: &[i64],
    seq2: &[u8],
    ref_pos2: &[i64],
    record1: &Record,
    record2: &Record,
) -> Result<(Vec<i64>, Vec<i64>), Box<dyn std::error::Error>> {
    let qname1 = std::str::from_utf8(record1.qname())?;
    let qname2 = std::str::from_utf8(record2.qname())?;

    debug!(
        "[stat_shared_snv_matches] Analyzing shared SNVs between reads {} and {}",
        qname1, qname2
    );

    // Find indices where both vectors have SNVs, including compound SNV+insertion values.
    let mut all_shared_snv_positions = Vec::new();
    let min_len = interval_hap1.len().min(interval_hap2.len());

    for i in 0..min_len {
        if is_snv_value(interval_hap1[i]) && is_snv_value(interval_hap2[i]) {
            let genomic_pos = interval_start + i as i64;
            all_shared_snv_positions.push(genomic_pos);
        }
    }

    debug!(
        "[stat_shared_snv_matches] Found {} positions with shared SNVs",
        all_shared_snv_positions.len()
    );

    if all_shared_snv_positions.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }

    let mut matching_shared_snv_positions = Vec::new();
    let mut discrepant_shared_snv_positions = Vec::new();

    // Process each shared SNV position and check if alt bases are identical
    for genomic_pos in &all_shared_snv_positions {
        // Extract bases at this position using reference position mapping
        let alt_base1_opt = get_base_at_position(seq1, ref_pos1, *genomic_pos);
        let alt_base2_opt = get_base_at_position(seq2, ref_pos2, *genomic_pos);

        let (alt_base1, alt_base2) = match (alt_base1_opt, alt_base2_opt) {
            (Some(b1), Some(b2)) => (b1, b2),
            _ => {
                debug!(
                    "[stat_shared_snv_matches] Cannot extract bases at position {} - skipping",
                    genomic_pos
                );
                continue;
            }
        };

        // Check if both reads have the same alt base (and not N=4)
        if alt_base1 != 4 && alt_base1 == alt_base2 {
            matching_shared_snv_positions.push(*genomic_pos);
            debug!(
                "[stat_shared_snv_matches] Both reads have same alt base {} at position {}",
                alt_base1, genomic_pos
            );
        } else if alt_base1 != 4 && alt_base2 != 4 && alt_base1 != alt_base2 {
            discrepant_shared_snv_positions.push(*genomic_pos);
            debug!(
                "[stat_shared_snv_matches] Alt bases differ at position {}: {} vs {}",
                genomic_pos, alt_base1, alt_base2
            );
        } else {
            debug!(
                "[stat_shared_snv_matches] One or both alt bases are N at position {}: {} vs {}",
                genomic_pos, alt_base1, alt_base2
            );
        }
    }

    debug!("[stat_shared_snv_matches] Final result: {} matching positions, {} discrepant positions out of {} total shared SNV positions", 
           matching_shared_snv_positions.len(), discrepant_shared_snv_positions.len(), all_shared_snv_positions.len());
    Ok((
        matching_shared_snv_positions,
        discrepant_shared_snv_positions,
    ))
}

/// Helper function to extract base at a specific genomic position from a read
///
/// Returns the encoded base (0=A, 1=T, 2=C, 3=G, 4=N) at the given genomic position
/// Returns None if the position is not covered by this read
fn get_base_at_position(query_seq: &[u8], ref_positions: &[i64], genomic_pos: i64) -> Option<u8> {
    // Find query index for this genomic position
    for (query_idx, &ref_pos) in ref_positions.iter().enumerate() {
        if ref_pos == genomic_pos {
            if query_idx < query_seq.len() {
                return Some(query_seq[query_idx]);
            }
        }
    }
    None
}

/// Check if variants at the same position are compatible
pub fn are_variants_compatible(vars1: &[&Variant], vars2: &[&Variant]) -> bool {
    // For now, simple check: SNVs must have the same alt base
    for v1 in vars1 {
        for v2 in vars2 {
            match (v1, v2) {
                (Variant::Snv { alt_base: b1, .. }, Variant::Snv { alt_base: b2, .. }) => {
                    if b1 != b2 {
                        return false;
                    }
                }
                // Different variant types at same position = incompatible
                _ => return false,
            }
        }
    }
    true
}

#[cfg(test)]
mod hap_vector_tests {
    use super::extract_hap_vector;
    use rust_htslib::bam::record::{Cigar, CigarString};
    use rust_htslib::bam::Record;

    /// Build a record with the given CIGAR (seq/qual sized to the query length).
    fn rec(cigar: Vec<Cigar>) -> Record {
        let qlen: usize = cigar
            .iter()
            .map(|c| match c {
                Cigar::Match(l)
                | Cigar::Ins(l)
                | Cigar::SoftClip(l)
                | Cigar::Equal(l)
                | Cigar::Diff(l) => *l as usize,
                _ => 0,
            })
            .sum();
        let seq = vec![b'A'; qlen];
        let qual = vec![30u8; qlen];
        let cs = CigarString(cigar);
        let mut r = Record::new();
        r.set(b"r", Some(&cs), &seq, &qual);
        r.set_pos(100);
        r.set_tid(0);
        r.set_mapq(60);
        r
    }

    #[test]
    fn insertion_marker_is_deferred_to_next_ref_position() {
        // 3=1I3=: marker at index 3 (the position AFTER the 3 matches), not index 2.
        // This is the T2 fix — the old code put it at index 2 (the position before).
        let v = extract_hap_vector(&rec(vec![Cigar::Equal(3), Cigar::Ins(1), Cigar::Equal(3)]));
        assert_eq!(v, vec![1, 1, 1, 11, 1, 1]);
    }

    #[test]
    fn insertion_at_read_start_is_dropped() {
        // 1I5=: a leading insertion has no prior ref position → dropped (Python `if index > 0`).
        let v = extract_hap_vector(&rec(vec![Cigar::Ins(1), Cigar::Equal(5)]));
        assert_eq!(v, vec![1, 1, 1, 1, 1]);
    }

    #[test]
    fn pending_insertion_drains_into_refskip() {
        // 3=1I2N3=: the deferred marker drains onto the first N base. The old `_=>`
        // catch-all pushed N as matches WITHOUT draining a pending insertion
        // (the latent second off-by-one) — this guards the fix.
        let v = extract_hap_vector(&rec(vec![
            Cigar::Equal(3),
            Cigar::Ins(1),
            Cigar::RefSkip(2),
            Cigar::Equal(3),
        ]));
        assert_eq!(v, vec![1, 1, 1, 11, 1, 1, 1, 1]);
    }

    #[test]
    fn compound_insertion_then_mismatch_preserves_snv_and_insertion() {
        // 3=1I1X2=: the compound base is -4 + 10 = 6, preserving both signals.
        let v = extract_hap_vector(&rec(vec![
            Cigar::Equal(3),
            Cigar::Ins(1),
            Cigar::Diff(1),
            Cigar::Equal(2),
        ]));
        assert_eq!(v, vec![1, 1, 1, 6, 1, 1]);
    }

    #[test]
    fn plain_match_snv_deletion() {
        assert_eq!(
            extract_hap_vector(&rec(vec![Cigar::Equal(5)])),
            vec![1, 1, 1, 1, 1]
        );
        assert_eq!(
            extract_hap_vector(&rec(vec![Cigar::Equal(2), Cigar::Diff(1), Cigar::Equal(2)])),
            vec![1, 1, -4, 1, 1]
        );
        assert_eq!(
            extract_hap_vector(&rec(vec![Cigar::Equal(3), Cigar::Del(2), Cigar::Equal(3)])),
            vec![1, 1, 1, -10, -10, 1, 1, 1]
        );
    }
}

/// Tests for the Python-faithful edge-weight formula (T3 fixes #6 + #4):
///   weight = overlap_span                                  # EQUAL-position count
///          + sum(score_arr[:shared_snv_count])             # shared_snv base-agnostic
///          + mean_read_length * (shared_snv - psv) * 0.75  # intrinsic-AD PSV term
///          + mean_read_length * 3 * indel_num              # indels over EQUAL positions
/// covering `compute_edge_weight_base` / `psv_shared_snvs` directly and both
/// branches of `determine_same_haplotype` end-to-end.
#[cfg(test)]
mod weight_tests {
    use super::*;
    use crate::structs::{AlleleDepthMap, HaplotypeConfig};
    use ahash::AHashMap;
    use rust_htslib::bam::record::{Cigar, CigarString};
    use rust_htslib::bam::Record;

    const MRL: f32 = 148.0;

    fn cfg() -> HaplotypeConfig {
        HaplotypeConfig::new(MRL)
    }

    /// Encoded read fully aligned (one ref-consuming base per query base) starting
    /// at `start`: builds (encoded_seq, reference_positions_full) the way
    /// `extract_query_seq` would for an all-`=`/`X` CIGAR.
    fn aligned(bases: &[u8], start: i64) -> (Vec<u8>, Vec<i64>) {
        let ref_pos: Vec<i64> = (0..bases.len() as i64).map(|i| start + i).collect();
        (bases.to_vec(), ref_pos)
    }

    /// Intrinsic AD map carrying a single supported allele (`base_idx`, depth>0) at
    /// `(chrom, pos)` — the Rust analogue of one `intrinsic_ad_dict[pos][base]` entry.
    fn intrinsic_with(chrom: &str, pos: u32, base_idx: usize, depth: u32) -> AlleleDepthMap {
        let mut m = AlleleDepthMap::new();
        let mut data = AlleleDepthMap::new_position_data(depth.max(1));
        AlleleDepthMap::set_allele_depth(&mut data, base_idx, depth);
        m.insert(chrom, pos, data);
        m
    }

    // ---- compute_edge_weight_base (the shared core of BOTH paths) -------------

    /// FIX #6 (base = overlap_span): the base is the count of EQUAL positions, NOT
    /// the full `end - start` span. Here 3 of 4 positions are equal (idx1 differs
    /// 1 vs -4), so the base is 3, not 4. No shared SNVs / indels here.
    #[test]
    fn base_is_equal_position_count_not_full_span() {
        let h1 = [1, 1, 1, 1];
        let h2 = [1, -4, 1, 1];
        let (s1, r1) = aligned(&[0, 0, 0, 0], 100);
        let (s2, r2) = aligned(&[0, 0, 0, 0], 100);
        let w = compute_edge_weight_base(
            &h1,
            &h2,
            100,
            &s1,
            &r1,
            &s2,
            &r2,
            "chr1",
            &AlleleDepthMap::new(),
            &cfg(),
        );
        assert!(
            (w - 3.0).abs() < 1e-3,
            "expected overlap_span base 3.0, got {w}"
        );
    }

    /// FIX #6 (indel_num over EQUAL positions): two indel blocks (a `-10` run and a
    /// positive insertion marker) within the identical part → `mrl*3*2`.
    #[test]
    fn indel_blocks_counted_over_identical_part() {
        let h1 = [1, 1, -10, -10, 1, 11, 1];
        let h2 = [1, 1, -10, -10, 1, 11, 1];
        let (s, r) = aligned(&[0; 7], 100);
        let w = compute_edge_weight_base(
            &h1,
            &h2,
            100,
            &s,
            &r,
            &s,
            &r,
            "chr1",
            &AlleleDepthMap::new(),
            &cfg(),
        );
        // overlap_span 7 + 0 SNV + 0 PSV term + mrl*3*2
        let expect = 7.0 + MRL * 3.0 * 2.0;
        assert!((w - expect).abs() < 1e-2, "expected {expect}, got {w}");
    }

    /// FIX #6 (shared_snv base-agnostic) + FIX #4 (PSV term, intrinsic supported):
    /// both positions are co-located SNVs (-4/-4) so shared_snv_count == 2 even
    /// though only one position has matching alt bases. With the intrinsic AD
    /// supporting the matching alt base, psv == 1, so the PSV term is mrl*(2-1)*0.75.
    #[test]
    fn shared_snv_base_agnostic_with_psv_deduction() {
        let h1 = [-4, -4];
        let h2 = [-4, -4];
        // pos100: both T(1) (matching); pos101: T(1) vs C(2) (mismatched alt)
        let (s1, r1) = aligned(&[1, 1], 100);
        let (s2, r2) = aligned(&[1, 2], 100);
        let intr = intrinsic_with("chr1", 100, 1, 5); // T supported at pos100
        let w = compute_edge_weight_base(&h1, &h2, 100, &s1, &r1, &s2, &r2, "chr1", &intr, &cfg());
        // overlap_span 2 + score_arr[0]+score_arr[1] + mrl*(2-1)*0.75 + 0
        let expect = 2.0 + (MRL + 2.0 * MRL) + MRL * 1.0 * 0.75;
        assert!((w - expect).abs() < 1e-2, "expected {expect}, got {w}");
    }

    /// FIX #4 (empty intrinsic ⇒ psv == 0): the same shared SNVs with NO intrinsic
    /// support give the full `mrl*(2-0)*0.75` boost (standalone-phaser behaviour /
    /// Python `intrinsic_ad_dict = {}`).
    #[test]
    fn empty_intrinsic_gives_zero_psv() {
        let h1 = [-4, -4];
        let h2 = [-4, -4];
        let (s1, r1) = aligned(&[1, 1], 100);
        let (s2, r2) = aligned(&[1, 1], 100);
        let w = compute_edge_weight_base(
            &h1,
            &h2,
            100,
            &s1,
            &r1,
            &s2,
            &r2,
            "chr1",
            &AlleleDepthMap::new(),
            &cfg(),
        );
        let expect = 2.0 + (MRL + 2.0 * MRL) + MRL * 2.0 * 0.75;
        assert!((w - expect).abs() < 1e-2, "expected {expect}, got {w}");
    }

    /// Unequal-length interval hap vectors are compared over `min_len` only (a
    /// guard against a truncated slice). hap2 is shorter, so only its 3 positions
    /// are scored: overlap_span 3, one shared SNV at idx0 (T/T, unsupported), no
    /// indels → 3 + score_arr[0] + mrl*(1-0)*0.75.
    #[test]
    fn unequal_length_hap_vectors_compare_over_min_len() {
        let h1 = [-4, 1, 1, 1, 1];
        let h2 = [-4, 1, 1];
        let (s1, r1) = aligned(&[1, 0, 0, 0, 0], 100);
        let (s2, r2) = aligned(&[1, 0, 0], 100);
        let w = compute_edge_weight_base(
            &h1,
            &h2,
            100,
            &s1,
            &r1,
            &s2,
            &r2,
            "chr1",
            &AlleleDepthMap::new(),
            &cfg(),
        );
        let expect = 3.0 + MRL + MRL * 0.75;
        assert!((w - expect).abs() < 1e-2, "expected {expect}, got {w}");
    }

    // ---- psv_shared_snvs (counts) --------------------------------------------

    /// shared_snv_count is base-agnostic (both -4) = 2; psv_snv_count counts only
    /// the matching-alt position that is intrinsic-supported = 1.
    #[test]
    fn psv_counts_shared_and_supported() {
        let h1 = [-4, -4];
        let h2 = [-4, -4];
        let (s1, r1) = aligned(&[1, 1], 100);
        let (s2, r2) = aligned(&[1, 2], 100);
        let intr = intrinsic_with("chr1", 100, 1, 3);
        let (psv, shared) = psv_shared_snvs(&h1, &h2, 100, &s1, &r1, &s2, &r2, "chr1", &intr);
        assert_eq!((psv, shared), (1, 2));
        // Without intrinsic support, psv collapses to 0 (shared unchanged).
        let (psv0, shared0) = psv_shared_snvs(
            &h1,
            &h2,
            100,
            &s1,
            &r1,
            &s2,
            &r2,
            "chr1",
            &AlleleDepthMap::new(),
        );
        assert_eq!((psv0, shared0), (0, 2));
    }

    /// N alt bases and supported-but-discordant alt bases never become PSVs, but
    /// they are still counted in shared_snv_count (base-agnostic).
    #[test]
    fn psv_skips_n_and_discordant_alt() {
        let h1 = [-4, -4];
        let h2 = [-4, -4];
        // pos100: N(4) vs N(4) → skipped for PSV; pos101: T(1) vs T(1) supported → PSV
        let (s1, r1) = aligned(&[4, 1], 100);
        let (s2, r2) = aligned(&[4, 1], 100);
        let intr = intrinsic_with("chr1", 101, 1, 7);
        let (psv, shared) = psv_shared_snvs(&h1, &h2, 100, &s1, &r1, &s2, &r2, "chr1", &intr);
        assert_eq!((psv, shared), (1, 2));
    }

    // ---- determine_same_haplotype end-to-end (both paths) --------------------

    fn make_record(
        qname: &[u8],
        cigar: Vec<Cigar>,
        seq: &[u8],
        quals: &[u8],
        pos: i64,
        flags: u16,
    ) -> Record {
        let cs = CigarString(cigar);
        let mut r = Record::new();
        r.set(qname, Some(&cs), seq, quals);
        r.set_pos(pos);
        r.set_tid(0);
        r.set_mapq(60);
        r.set_flags(flags);
        r
    }

    fn empty_maps() -> (
        AHashMap<String, Vec<i16>>,
        AHashMap<String, Vec<f32>>,
        AHashMap<String, (i64, i64)>,
    ) {
        (AHashMap::new(), AHashMap::new(), AHashMap::new())
    }

    /// IDENTICAL path: two reads with an identical `5=1X4=` alignment (one shared
    /// SNV, no intrinsic support). The normalized weight equals the Python formula
    /// `(overlap_span 10 + score_arr[0] + mrl*1*0.75) / (mrl*10)`.
    #[test]
    fn identical_path_weight_matches_formula() {
        let cigar = vec![Cigar::Equal(5), Cigar::Diff(1), Cigar::Equal(4)];
        let seq = b"AAAAATAAAA"; // alt 'T' at the X (idx 5)
        let quals = [30u8; 10];
        let r1 = make_record(b"rA", cigar.clone(), seq, &quals, 100, 65);
        let r2 = make_record(b"rB", cigar, seq, &quals, 100, 65);

        let ad = AlleleDepthMap::new();
        let intr = AlleleDepthMap::new();
        let config = cfg();
        let (mut hv, mut ev, mut rp) = empty_maps();

        let (res, w) = determine_same_haplotype(
            &r1, &r2, 100, 110, "chr1", &ad, &intr, &config, &mut hv, &mut ev, &mut rp,
        )
        .expect("determine_same_haplotype failed");

        assert_eq!(res, HaplotypeResult::Same);
        let raw = 10.0 + MRL + MRL * 0.75; // 269.0
        let expect = raw / (MRL * 10.0);
        let got = w.expect("identical path should return a weight");
        assert!((got - expect).abs() < 1e-5, "expected {expect}, got {got}");
    }

    /// TOLERATED path: read1 carries a single low-quality SNV (`15=1X14=`, baseQ 5)
    /// vs read2's plain `30=`; the main-bam AD flags it a sequencing error so the
    /// pair is tolerated. Weight = (overlap_span 29 − 1×20) / (mrl*10).
    #[test]
    fn tolerated_path_weight_applies_penalty() {
        let mut seq1 = vec![b'A'; 30];
        seq1[15] = b'T'; // alt at the X
        let mut quals1 = [30u8; 30];
        quals1[15] = 5; // low base quality → candidate sequencing error
        let r1 = make_record(
            b"rA",
            vec![Cigar::Equal(15), Cigar::Diff(1), Cigar::Equal(14)],
            &seq1,
            &quals1,
            100,
            65,
        );
        let seq2 = vec![b'A'; 30];
        let quals2 = [30u8; 30];
        let r2 = make_record(b"rB", vec![Cigar::Equal(30)], &seq2, &quals2, 100, 129);

        // Main-bam AD: at 0-based pos 115 the alt 'T'(idx1) has AD 1 over DP 100
        // (af 0.01 ≤ 0.02 and ad==1, dp≥10) → is_sequencing_error == true.
        let mut ad = AlleleDepthMap::new();
        let mut data = AlleleDepthMap::new_position_data(100);
        AlleleDepthMap::set_allele_depth(&mut data, 1, 1);
        ad.insert("chr1", 115, data);

        let intr = AlleleDepthMap::new();
        let config = cfg();
        let (mut hv, mut ev, mut rp) = empty_maps();

        let (res, w) = determine_same_haplotype(
            &r1, &r2, 100, 130, "chr1", &ad, &intr, &config, &mut hv, &mut ev, &mut rp,
        )
        .expect("determine_same_haplotype failed");

        assert_eq!(res, HaplotypeResult::Same, "mismatch should be tolerated");
        let raw = (29.0 - 20.0_f32).max(0.0); // base overlap_span 29, minus 1 tolerated*20
        let expect = raw / (MRL * 10.0);
        let got = w.expect("tolerated path should return a weight");
        assert!((got - expect).abs() < 1e-5, "expected {expect}, got {got}");
    }
}
