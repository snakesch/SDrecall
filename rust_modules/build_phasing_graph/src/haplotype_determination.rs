/// Haplotype determination logic
/// 
/// This module contains the core algorithm for determining whether two read pairs
/// come from the same haplotype by analyzing their sequences and CIGAR operations.

use log::{debug, warn, error, info};
use ahash::AHashMap;
use rust_htslib::bam::Record;
use rust_htslib::bam::ext::BamRecordExtensions;


use crate::structs::{
    AlleleDepthMap, HaplotypeConfig, ReadHaplotypeVector, 
    ReadErrorVector, Variant
};

/// Result of haplotype compatibility analysis
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HaplotypeResult {
    Same,       // Same haplotype
    Different,  // Different haplotypes  
    Unknown,    // Cannot determine (insufficient data)
}

/// Extract query sequence and reference positions from a BAM record
/// 
/// Uses BamRecordExtensions::reference_positions_full() for efficient coordinate mapping
/// 
/// Returns (query_sequence, ref_positions) where:
/// - query_sequence: Vec<u8> of base-encoded sequence  
/// - ref_positions: Vec<i64> mapping query positions to reference positions
pub fn extract_query_seq(record: &Record) -> Result<(Vec<u8>, Vec<i64>), Box<dyn std::error::Error>> {
    let query_seq = record.seq().as_bytes();
    
    // Use BamRecordExtensions for efficient reference position mapping
    let ref_positions: Vec<i64> = record.reference_positions_full()
        .map(|pos_opt| pos_opt.map(|p| p as i64).unwrap_or(-1))
        .collect();
    
    // Convert to our format: A=0, T=1, C=2, G=3, N=4
    let encoded_seq: Vec<u8> = query_seq.iter().map(|&base| {
        match base {
            b'A' | b'a' => 0,
            b'T' | b't' => 1, 
            b'C' | b'c' => 2,
            b'G' | b'g' => 3,
            _ => 4, // N or any other base
        }
    }).collect();
    
    Ok((encoded_seq, ref_positions))
}

/// Slice query sequence to a specific genomic interval using reference positions
/// 
/// Returns the subsequence that aligns to the given genomic interval
pub fn slice_seq_to_interval(
    query_seq: &[u8], 
    ref_positions: &[i64],
    interval_start: i64, 
    interval_end: i64
) -> Vec<u8> {
    let mut result = Vec::new();
    
    for (i, &ref_pos) in ref_positions.iter().enumerate() {
        if ref_pos >= interval_start && ref_pos < interval_end {
            if i < query_seq.len() {
                result.push(query_seq[i]);
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

/// Extract haplotype vector from CIGAR operations
/// 
/// Returns a vector representing variants relative to reference:
/// - 1: Match
/// - -4: SNV (mismatch)  
/// - -6: Deletion
/// - positive values: Insertion length
pub fn extract_hap_vector(record: &Record) -> Vec<i16> {
    let cigar = record.cigar();
    let mut hap_vector = Vec::new();
    
    for &op in cigar.iter() {
        use rust_htslib::bam::record::Cigar;
        match op {
            Cigar::Match(len) | Cigar::Equal(len) => {
                // Match to reference
                for _ in 0..len {
                    hap_vector.push(1);
                }
            }
            Cigar::Diff(len) => {
                // Mismatch (SNV)
                for _ in 0..len {
                    hap_vector.push(-4);
                }
            }
            Cigar::Ins(len) => {
                // Insertion - mark previous position with insertion length
                if !hap_vector.is_empty() {
                    let last_idx = hap_vector.len() - 1;
                    hap_vector[last_idx] = len as i16 * 4; // Encode insertion
                }
            }
            Cigar::Del(len) => {
                // Deletion
                for _ in 0..len {
                    hap_vector.push(-6);
                }
            }
            Cigar::SoftClip(_) | Cigar::HardClip(_) => {
                // Skip clipped bases
            }
            _ => {
                // Handle other operations as matches for now
                let len = op.len();
                for _ in 0..len {
                    hap_vector.push(1);
                }
            }
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
    let snv_positions: Vec<bool> = hap_vector.iter().map(|&val| val == -4).collect();
    count_continuous_blocks(&snv_positions)
}

/// Count indel blocks in a haplotype vector  
/// 
/// Equivalent to Python's count_continuous_indel_blocks function
/// Counts continuous blocks of deletions (-6) or insertions (>1)
fn count_indel_blocks(hap_vector: &[i16]) -> usize {
    let indel_positions: Vec<bool> = hap_vector.iter().map(|&val| val == -6 || val > 1).collect();
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
    interval_end: i64
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
                    if ref_consume + i < error_vector.len() && query_consume + i < base_qualities.len() {
                        error_vector[ref_consume + i] = base_qualities[query_consume + i] as f32;
                    }
                }
                query_consume += *length as usize;
                ref_consume += *length as usize;
            }
            Cigar::Diff(length) => {
                // operation = 8: Mismatch
                for i in 0..*length as usize {
                    if ref_consume + i < error_vector.len() && query_consume + i < base_qualities.len() {
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
                    if ref_consume + i < error_vector.len() && query_consume + i < base_qualities.len() {
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
            0.0  // No error for indels/skips
        } else {
            10.0_f32.powf(-(*prob) / 10.0)  // Convert phred to probability
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
fn get_error_vector(
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
fn get_read_id(record: &Record) -> Result<String, Box<dyn std::error::Error>> {
    let qname = std::str::from_utf8(record.qname())?;
    Ok(format!("{}:{}", qname, record.flags()))
}

/// Get or compute haplotype vector for a read (with caching)
/// 
/// Checks cache first, computes and stores if not found
fn get_hap_vector(
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
) -> bool {
    let qname = std::str::from_utf8(record.qname()).unwrap_or("unknown");
    debug!("[ANALYZE] is_sequencing_error: Checking read {} at position {}", qname, genomic_pos);
    
    // **OPTIMIZED COORDINATE LOOKUP**: Use CigarStringView::read_pos 
    // This provides O(1) or O(log n) lookup instead of O(n) linear search
    if let Some(cigar_view) = record.cigar_cached() {
        match cigar_view.read_pos(genomic_pos as u32, true, true) {
        Ok(Some(query_pos)) => {
            let qi = query_pos as usize;
            let query_seq = record.seq().as_bytes();
            
            if qi >= query_seq.len() {
                warn!("[WARNING] Query index {} out of bounds (seq_len={})", qi, query_seq.len());
                return false;
            }
            
            // Use pre-computed error probability (indexed by reference position)
            let read_offset = (genomic_pos - record.reference_start()) as usize;
            if read_offset >= error_vector.len() {
                warn!("[WARNING] Position {} out of bounds in error vector (len={})", genomic_pos, error_vector.len());
                return false;
            }
            
            let error_prob = error_vector[read_offset];
            debug!("[STATS] Error probability at position {}: {:.6}", genomic_pos, error_prob);
            
            // Convert error probability back to approximate phred score for threshold check
            // error_prob = 10^(-qual/10) → qual = -10 * log10(error_prob)
            let approx_qual = if error_prob > 0.0 {
                (-10.0 * error_prob.log10()) as u8
            } else {
                40 // High quality if error_prob is 0
            };
            
            debug!("[QUALITY] Approximate quality score: Q{} (from error_prob={:.6})", approx_qual, error_prob);
            
            // If base quality >= 15, not a sequencing error
            if approx_qual >= 15 {
                debug!("[SUCCESS] High quality (Q{} >= 15) -> NOT a sequencing error", approx_qual);
                return false;
            }
            
            warn!("[WARNING] Low quality (Q{} < 15) -> checking allele depth...", approx_qual);
            
            // Get allele depth information for this position
            // Note: For full functionality, we'd need header access to convert tid to chromosome name
            // For now, use a placeholder approach - this may need header parameter in future
            let tid = record.tid();
            let chrom = format!("chr{}", tid); // Temporary placeholder - proper implementation needs header
            if let Some(pos_data) = allele_depth_map.get(&chrom, genomic_pos as u32) {
                let target_base = query_seq[qi];
                let ad = AlleleDepthMap::get_allele_depth(pos_data, target_base as usize);
                let dp = AlleleDepthMap::total_depth(pos_data);
                
                debug!("[DNA] Base at position {}: encoded={}, which is {}, AD={}, DP={}", 
                       genomic_pos, target_base, target_base as char, ad, dp);
                
                if dp == 0 {
                    warn!("[WARNING] Zero depth -> NOT a sequencing error");
                    return false;
                }
                
                let af = ad as f32 / dp as f32;
                debug!("[ALLELE] Allele frequency: {:.4} ({}/{})", af, ad, dp);
                
                // Python criteria: (af <= 0.02 or (ad == 1 and dp >= 10)) and base_qual < 13
                let af_criteria = af <= 0.02 || (ad == 1 && dp >= 10);
                let qual_criteria = approx_qual < 13;
                
                debug!("[CHECK] Criteria check: AF_criteria={} (af={:.4} <= 0.02 OR (ad={} == 1 AND dp={} >= 10))", 
                       af_criteria, af, ad, dp);
                debug!("[CHECK] Criteria check: QUAL_criteria={} (Q{} < 13)", qual_criteria, approx_qual);
                
                let is_error = af_criteria && qual_criteria;
                debug!("[RESULT] Final decision: {} (AF_criteria={} AND QUAL_criteria={})", 
                       if is_error { "SEQUENCING ERROR" } else { "REAL VARIANT" }, 
                       af_criteria, qual_criteria);
                
                return is_error;
            } else {
                warn!("[WARNING] No allele depth data for {}:{} -> NOT a sequencing error", chrom, genomic_pos);
            }
        }
        Ok(None) => {
            // Position is outside the read alignment (before start or after end)
            warn!("[WARNING] Position {} outside read alignment -> NOT a sequencing error", genomic_pos);
        }
            Err(e) => {
                // Error in read_pos function call
                error!("[ERROR] read_pos failed for position {}: {} -> NOT a sequencing error", genomic_pos, e);
            }
        }
    } else {
        // No cached CIGAR view available
        warn!("[WARNING] No cached CIGAR view available for position {} -> NOT a sequencing error", genomic_pos);
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
    hap_vec1: &[i16], hap_vec2: &[i16],
    _seq1: &[u8], _ref_pos1: &[i64], 
    _seq2: &[u8], _ref_pos2: &[i64],
    record1: &Record, error_vec1: &[f32],
    record2: &Record, error_vec2: &[f32],
    mismatch_positions: &[i64],
    interval_start: i64,
    allele_depth_map: &AlleleDepthMap,
) -> (bool, usize) {
    let qname1 = std::str::from_utf8(record1.qname()).unwrap_or("unknown");
    let qname2 = std::str::from_utf8(record2.qname()).unwrap_or("unknown");
    
    debug!("[TOLERATE] Analyzing {} mismatches between reads {} and {}", 
           mismatch_positions.len(), qname1, qname2);
    debug!("[POSITIONS] Mismatch positions: {:?}", mismatch_positions);
    
    let mut tolerable_count = 0;
    
    for (idx, &genomic_pos) in mismatch_positions.iter().enumerate() {
        let hap_index = (genomic_pos - interval_start) as usize;
        
        // Get haplotype values at this position
        let hap1 = if hap_index < hap_vec1.len() { hap_vec1[hap_index] } else { 1 };
        let hap2 = if hap_index < hap_vec2.len() { hap_vec2[hap_index] } else { 1 };
        
        debug!("[ANALYZE] Mismatch #{} at position {} (hap1={}, hap2={})", 
               idx + 1, genomic_pos, hap1, hap2);
        
        // Check if this is an indel mismatch (not tolerable)
        let is_indel1 = hap1 == -6 || hap1 > 1;  // Deletion or insertion
        let is_indel2 = hap2 == -6 || hap2 > 1;  // Deletion or insertion
        
        if is_indel1 || is_indel2 {
            debug!("[INDEL] Position {} has indel mismatch (hap1={}, hap2={}) - NOT TOLERABLE", 
                   genomic_pos, hap1, hap2);
            debug!("[REJECT] tolerate_mismatches: FAILED at indel position {} - returning (false, 0)", genomic_pos);
            return (false, 0);
        }
        
        // For SNV mismatches, check if they can be explained by sequencing errors
        let error1 = is_sequencing_error(
            record1, error_vec1, genomic_pos, allele_depth_map
        );
        let error2 = is_sequencing_error(
            record2, error_vec2, genomic_pos, allele_depth_map
        );
        
        debug!("[STATS] Position {}: read1({}) error={}, read2({}) error={}", 
               genomic_pos, qname1, error1, qname2, error2);
        
        // If either read has a sequencing error at this position, we can tolerate it
        if error1 || error2 {
            tolerable_count += 1;
            debug!("[SUCCESS] Position {} TOLERABLE (sequencing error detected)", genomic_pos);
        } else {
            // This mismatch cannot be explained by sequencing error
            debug!("[FAIL] Position {} NOT TOLERABLE (likely real variant)", genomic_pos);
            debug!("[REJECT] tolerate_mismatches: FAILED at position {} - returning (false, 0)", genomic_pos);
            return (false, 0);
        }
    }
    
    // All mismatches can be explained by sequencing errors
    debug!("[SUCCESS] tolerate_mismatches: SUCCESS - all {} mismatches tolerable (sequencing errors)", 
           mismatch_positions.len());
    debug!("[STATS] Total tolerated mismatches: {}", tolerable_count);
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
        
        // If haplotype values differ, it's a mismatch
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
/// - -6: Deletion
/// - Positive values: Insertion length
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
        let is_indel1 = hap1 == -6 || hap1 > 1;  // Deletion or insertion
        let is_indel2 = hap2 == -6 || hap2 > 1;  // Deletion or insertion
        
        if is_indel1 || is_indel2 {
            debug!("[INDEL] Found indel mismatch at position {}: hap1={}, hap2={}", 
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
    config: &HaplotypeConfig,
    read_hap_vectors: &mut AHashMap<String, ReadHaplotypeVector>,
    read_error_vectors: &mut AHashMap<String, ReadErrorVector>, 
    read_ref_pos_dict: &mut AHashMap<String, (i64, i64)>,
) -> Result<(HaplotypeResult, Option<f32>), Box<dyn std::error::Error>> {
    // Extract qnames from records
    let qname1 = std::str::from_utf8(read1.qname())?;
    let qname2 = std::str::from_utf8(read2.qname())?;
    
    debug!("[HAPLOTYPE] Comparing reads {} and {} at {}:{}-{}", 
           qname1, qname2, chrom, start, end);
    
    // Step 1: Extract query sequences and reference positions
    let (seq1, ref_pos1) = extract_query_seq(read1)?;
    let (seq2, ref_pos2) = extract_query_seq(read2)?;
    
    // Step 2: Slice sequences to the genomic interval
    let interval_seq1 = slice_seq_to_interval(&seq1, &ref_pos1, start, end);
    let interval_seq2 = slice_seq_to_interval(&seq2, &ref_pos2, start, end);
    
    // Step 3: Get haplotype vectors (needed for weight calculation in both paths)
    // OPTIMIZATION: Compute once, use in both identical and tolerable mismatch cases
    let hap_vec1 = get_hap_vector(read1, read_hap_vectors)?;
    let hap_vec2 = get_hap_vector(read2, read_hap_vectors)?;
    
    // Step 4: Compare sequences for exact match
    if compare_sequences(&interval_seq1, &interval_seq2) {
        // SEQUENCES ARE IDENTICAL: Same haplotype (likely)
        // 
        // PERFORMANCE OPTIMIZATION: 
        // No error vector computation needed when sequences match exactly.
        // Error vectors are only used for sequencing error detection in mismatches.
        // This saves expensive CIGAR parsing + quality score conversion.
        
        debug!("[IDENTICAL] Sequences IDENTICAL between reads {} and {}", qname1, qname2);
        info!("[INFO] Interval: {}:{}-{} (length={})", chrom, start, end, end - start);
        
        // Sequences are identical - likely same haplotype
        
        // Step 5: Slice haplotype vectors to interval (using pre-computed vectors)
        let interval_hap1 = slice_hap_vector(&hap_vec1, read1.pos(), start, end);
        let interval_hap2 = slice_hap_vector(&hap_vec2, read2.pos(), start, end);
        
        // Step 6: Count shared variants for weight calculation
        let (snv1, indel1, _) = count_variants(&interval_hap1);
        let (snv2, indel2, _) = count_variants(&interval_hap2);
        
        // Use minimum counts as shared variants (conservative estimate)
        let shared_snvs = snv1.min(snv2);
        let shared_indels = indel1.min(indel2);
        
        debug!("[STATS] Variant counts: read1(snv={}, indel={}), read2(snv={}, indel={}), shared(snv={}, indel={})", 
               snv1, indel1, snv2, indel2, shared_snvs, shared_indels);
        
        // Calculate weight based on overlap and shared variants
        let overlap_length = end - start;
        let mut weight = overlap_length as f32;
        
        // Add weight for shared variants using score array
        for i in 0..shared_snvs.min(config.score_array.len()) {
            weight += config.score_array[i];
        }
        
        // Add weight for shared indels
        weight += shared_indels as f32 * config.mean_read_length * 3.0;
        
        debug!("[WEIGHT] Weight calculation: base={} + variant_bonus={:.2} + indel_bonus={:.2} = {:.2}", 
               overlap_length,
               (0..shared_snvs.min(config.score_array.len())).map(|i| config.score_array[i]).sum::<f32>(),
               shared_indels as f32 * config.mean_read_length * 3.0,
               weight);
        
        // Normalize weight
        let normalized_weight = weight / (config.mean_read_length * 10.0);
        
        debug!("[RESULT] Final weight: {:.6} (normalized from {:.2})", normalized_weight, weight);
        
        // No need for error vectors when sequences are identical
        // Error vectors are only used for sequencing error detection in mismatches
        
        // Store read positions
        read_ref_pos_dict.insert(qname1.to_string(), (start, end));
        read_ref_pos_dict.insert(qname2.to_string(), (start, end));
        
        debug!("[SUCCESS] SAME HAPLOTYPE (identical) - weight: {:.6}", normalized_weight);
        
        Ok((HaplotypeResult::Same, Some(normalized_weight)))
        
    } else {
        // Sequences differ - need detailed analysis
        
        // Check if sequences have different lengths (likely indel differences)
        if interval_seq1.len() != interval_seq2.len() {
            debug!("[DIFF] Sequence length mismatch: read1_len={}, read2_len={} -> DIFFERENT haplotypes", 
                   interval_seq1.len(), interval_seq2.len());
            return Ok((HaplotypeResult::Different, None));
        }
        
        // Check if either sequence is empty
        if interval_seq1.is_empty() || interval_seq2.is_empty() {
            warn!("[WARNING] Empty sequence(s): read1_len={}, read2_len={} -> UNKNOWN", 
                   interval_seq1.len(), interval_seq2.len());
            return Ok((HaplotypeResult::Unknown, None));
        }
        
        // Slice haplotype vectors to the interval for mismatch analysis
        let interval_hap1 = slice_hap_vector(&hap_vec1, read1.pos(), start, end);
        let interval_hap2 = slice_hap_vector(&hap_vec2, read2.pos(), start, end);
        
        // Find mismatch positions using haplotype vectors (CORRECT APPROACH)
        let mismatch_positions = find_mismatch_positions_from_hap_vectors(
            &interval_hap1, &interval_hap2, start, end
        );
        
        // If no mismatches found (shouldn't happen since sequences differ), treat as unknown
        if mismatch_positions.is_empty() {
            warn!("[WARNING] No mismatches found despite sequence differences -> UNKNOWN");
            return Ok((HaplotypeResult::Unknown, None));
        }
        
        // Check if mismatches occur at indel positions - these are not tolerable
        if has_indel_mismatches_from_hap_vectors(&interval_hap1, &interval_hap2, &mismatch_positions, start) {
            debug!("[INDEL] Found indel mismatches between reads, marking as different haplotypes");
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
        
        debug!("[MISMATCH] Sequence mismatch detected between reads {} and {}", qname1, qname2);
        info!("[INFO] Interval: {}:{}-{} (length={})", chrom, start, end, end - start);
        debug!("[CHECK] Sequences: read1_len={}, read2_len={}", interval_seq1.len(), interval_seq2.len());
        
        // Check if all mismatches can be explained by sequencing errors
        // Using haplotype vectors to properly identify mismatch types (CORRECT APPROACH)
        let (tolerable, tolerated_count) = tolerate_mismatches_from_hap_vectors(
            &interval_hap1, &interval_hap2,
            &seq1, &ref_pos1, 
            &seq2, &ref_pos2,
            read1, &error_vec1,
            read2, &error_vec2,
            &mismatch_positions,
            start,
            allele_depth_map,
        );
        
        debug!("[RESULT] Mismatch tolerance result: tolerable={}, tolerated_count={}", tolerable, tolerated_count);
        
        if tolerable {
            // All mismatches are sequencing errors - treat as same haplotype with penalty
            
            debug!("[SUCCESS] All mismatches tolerable as sequencing errors - calculating weight with penalty");
            
            // Slice haplotype vectors for weight calculation (using pre-computed vectors)
            let interval_hap1 = slice_hap_vector(&hap_vec1, read1.pos(), start, end);
            let interval_hap2 = slice_hap_vector(&hap_vec2, read2.pos(), start, end);
            
            // Count shared variants
            let (snv1, indel1, _) = count_variants(&interval_hap1);
            let (snv2, indel2, _) = count_variants(&interval_hap2);
            let shared_snvs = snv1.min(snv2);
            let shared_indels = indel1.min(indel2);
            
            debug!("[STATS] Variant counts: read1(snv={}, indel={}), read2(snv={}, indel={}), shared(snv={}, indel={})", 
                   snv1, indel1, snv2, indel2, shared_snvs, shared_indels);
            
            // Calculate weight with penalty for tolerated mismatches
            let overlap_length = end - start;
            let mut weight = overlap_length as f32;
            
            // Add weight for shared variants
            for i in 0..shared_snvs.min(config.score_array.len()) {
                weight += config.score_array[i];
            }
            weight += shared_indels as f32 * config.mean_read_length * 3.0;
            
            debug!("[WEIGHT] Weight before penalty: {:.2} (overlap={}, variant_bonus={:.2}, indel_bonus={:.2})", 
                   weight, overlap_length, 
                   (0..shared_snvs.min(config.score_array.len())).map(|i| config.score_array[i]).sum::<f32>(),
                   shared_indels as f32 * config.mean_read_length * 3.0);
            
            // Apply penalty for tolerated sequencing errors (like Python: weight - tolerated_count * 20)
            weight -= tolerated_count as f32 * 20.0;
            
            debug!("[PENALTY] Penalty applied: -{} x 20.0 = -{:.1}", tolerated_count, tolerated_count as f32 * 20.0);
            
            // Normalize weight
            let normalized_weight = weight / (config.mean_read_length * 10.0);
            
            debug!("[RESULT] Final weight: {:.6} (normalized from {:.2})", normalized_weight, weight);
            
            // Update data structures - error vectors already computed above
            
            read_ref_pos_dict.insert(qname1.to_string(), (start, end));
            read_ref_pos_dict.insert(qname2.to_string(), (start, end));
            
            debug!("[SUCCESS] SAME HAPLOTYPE (with penalty) - weight: {:.6}", normalized_weight);
            
            Ok((HaplotypeResult::Same, Some(normalized_weight)))
            
        } else {
            // Some mismatches cannot be explained by sequencing errors
            debug!("[FAIL] Some mismatches are real variants - marking as DIFFERENT haplotypes");
            debug!("[REJECT] DIFFERENT HAPLOTYPES - {} untolerable mismatches", mismatch_positions.len());
            Ok((HaplotypeResult::Different, None))
        }
    }
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