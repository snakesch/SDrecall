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
    // Handle both uppercase and lowercase for robustness (soft-masking, etc.)
    let encoded_seq: Vec<u8> = query_seq.iter().map(|&base| {
        match base {
            b'A' | b'a' => 0,
            b'T' | b't' => 1, 
            b'C' | b'c' => 2,
            b'G' | b'g' => 3,
            _ => 4, // N or any other base (including ambiguous bases)
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
    chrom: &str,
) -> bool {
    let qname = std::str::from_utf8(record.qname()).unwrap_or("unknown");
    debug!("[is_sequencing_error] is_sequencing_error: Checking read {} at position {}", qname, genomic_pos);
    
    // **COORDINATE LOOKUP**: Use CIGAR directly for position mapping
    // record.cigar() already returns a CigarStringView
    let cigar_view = record.cigar();
    match cigar_view.read_pos(genomic_pos as u32, true, true) {
        Ok(Some(query_pos)) => {
            let qi = query_pos as usize;
            let query_seq = record.seq().as_bytes();
            
            if qi >= query_seq.len() {
                warn!("[is_sequencing_error] Query index {} out of bounds (seq_len={})", qi, query_seq.len());
                return false;
            }
            
            // Use pre-computed error probability (indexed by reference position)
            let read_offset = (genomic_pos - record.reference_start()) as usize;
            if read_offset >= error_vector.len() {
                warn!("[is_sequencing_error] Position {} out of bounds in error vector (len={})", genomic_pos, error_vector.len());
                return false;
            }
            
            let error_prob = error_vector[read_offset];
            debug!("[is_sequencing_error] Error probability at position {}: {:.6}", genomic_pos, error_prob);
            
            // Convert error probability back to approximate phred score for threshold check
            // error_prob = 10^(-qual/10) → qual = -10 * log10(error_prob)
            let approx_qual = if error_prob > 0.0 {
                (-10.0 * error_prob.log10()) as u8
            } else {
                40 // High quality if error_prob is 0
            };
            
            debug!("[is_sequencing_error] Approximate quality score: Q{} (from error_prob={:.6})", approx_qual, error_prob);
            
            // If base quality >= 20, not a sequencing error
            if approx_qual >= 20 {
                debug!("[is_sequencing_error] High quality (Q{} >= 20) -> NOT a sequencing error", approx_qual);
                return false;
            }
            
            debug!("[is_sequencing_error] Low quality (Q{} < 20) -> checking allele depth...", approx_qual);
            
            // Get allele depth information for this position
            // Use the chromosome name passed from the parent function (already resolved correctly)
            debug!("[is_sequencing_error] Using chromosome: '{}' for position {}", chrom, genomic_pos);
            if let Some(pos_data) = allele_depth_map.get(&chrom, genomic_pos as u32) {
                let target_base = query_seq[qi];
                let ad = AlleleDepthMap::get_allele_depth(pos_data, target_base as usize);
                let dp = AlleleDepthMap::total_depth(pos_data);
                
                debug!("[is_sequencing_error] Base at position {}: encoded={}, which is {}, AD={}, DP={}", 
                       genomic_pos, target_base, target_base as char, ad, dp);
                
                if dp == 0 {
                    debug!("[is_sequencing_error] Zero depth -> NOT a sequencing error");
                    return false;
                }
                
                let af = ad as f32 / dp as f32;
                debug!("[is_sequencing_error] Allele frequency: {:.4} ({}/{})", af, ad, dp);
                
                // Python criteria: (af <= 0.02 or (ad == 1 and dp >= 10)) and base_qual < 13
                let af_criteria = af <= 0.02 || (ad == 1 && dp >= 10);
                let qual_criteria = approx_qual < 13;
                
                debug!("[is_sequencing_error] Criteria check: AF_criteria={} (af={:.4} <= 0.02 OR (ad={} == 1 AND dp={} >= 10))", 
                       af_criteria, af, ad, dp);
                debug!("[is_sequencing_error] Criteria check: QUAL_criteria={} (Q{} < 13)", qual_criteria, approx_qual);
                
                let is_error = af_criteria && qual_criteria;
                debug!("[is_sequencing_error] Final decision: {} (AF_criteria={} AND QUAL_criteria={})", 
                       if is_error { "SEQUENCING ERROR" } else { "REAL VARIANT" }, 
                       af_criteria, qual_criteria);
                
                return is_error;
            } else {
                debug!("[is_sequencing_error] No allele depth data for {}:{}, indicating 0 depth of the allele, CAN BE a sequencing error", chrom, genomic_pos);
                return true;
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
    hap_vec1: &[i16], hap_vec2: &[i16],
    _seq1: &[u8], _ref_pos1: &[i64], 
    _seq2: &[u8], _ref_pos2: &[i64],
    record1: &Record, error_vec1: &[f32],
    record2: &Record, error_vec2: &[f32],
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
    
    debug!("[tolerate_mismatches_from_hap_vectors] Analyzing {} mismatches between reads {} and {}", 
           mismatch_positions.len(), read1_id, read2_id);
    debug!("[tolerate_mismatches_from_hap_vectors] Mismatch positions: {:?}", mismatch_positions);
    
    let mut tolerable_count = 0;
    
    for (idx, &genomic_pos) in mismatch_positions.iter().enumerate() {
        let hap_index = (genomic_pos - interval_start) as usize;
        
        // Get haplotype values at this position
        let hap1 = if hap_index < hap_vec1.len() { hap_vec1[hap_index] } else { 1 };
        let hap2 = if hap_index < hap_vec2.len() { hap_vec2[hap_index] } else { 1 };
        
        debug!("[tolerate_mismatches_from_hap_vectors] Mismatch #{} at position {} (hap1={}, hap2={})", 
               idx + 1, genomic_pos, hap1, hap2);
        
        // Check if this is an indel mismatch (not tolerable)
        let is_indel1 = hap1 == -6 || hap1 > 1;  // Deletion or insertion
        let is_indel2 = hap2 == -6 || hap2 > 1;  // Deletion or insertion
        let is_snv1 = hap1 == -4;
        let is_snv2 = hap2 == -4;
        
        if is_indel1 || is_indel2 {
            debug!("[tolerate_mismatches_from_hap_vectors] Position {} has indel mismatch (hap1={}, hap2={}) - NOT TOLERABLE", 
                   genomic_pos, hap1, hap2);
            debug!("[tolerate_mismatches_from_hap_vectors] tolerate_mismatches: FAILED at indel position {} - returning (false, 0)", genomic_pos);
            return (false, 0);
        }
        
        // For SNV mismatches, check if they can be explained by sequencing errors
        let error1 = is_sequencing_error(
            record1, error_vec1, genomic_pos, allele_depth_map, chrom
        );
        let error2 = is_sequencing_error(
            record2, error_vec2, genomic_pos, allele_depth_map, chrom
        );
        
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
    debug!("[tolerate_mismatches_from_hap_vectors] Total tolerated mismatches: {}", tolerable_count);
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
    config: &HaplotypeConfig,
    read_hap_vectors: &mut AHashMap<String, ReadHaplotypeVector>,
    read_error_vectors: &mut AHashMap<String, ReadErrorVector>, 
    read_ref_pos_dict: &mut AHashMap<String, (i64, i64)>,
) -> Result<(HaplotypeResult, Option<f32>), Box<dyn std::error::Error>> {
    // Extract qnames from records
    let qname1 = std::str::from_utf8(read1.qname())?;
    let qname2 = std::str::from_utf8(read2.qname())?;
    let read1_id = format!("{}_{}", qname1, read1.flags());
    let read2_id = format!("{}_{}", qname2, read2.flags());
    
    debug!("[determine_same_haplotype] Comparing reads {} and {} at {}:{}-{}", 
           read1_id, read2_id, chrom, start, end);
    
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
        
        debug!("[determine_same_haplotype] Sequences IDENTICAL between reads {} and {}", read1_id, read2_id);
        debug!("[determine_same_haplotype] Interval: {}:{}-{} (length={})", chrom, start, end, end - start);
        
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
        
        debug!("[determine_same_haplotype] Variant counts: read1(snv={}, indel={}), read2(snv={}, indel={}), shared(snv={}, indel={})", 
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
        
        debug!("[determine_same_haplotype] Weight calculation: base={} + variant_bonus={:.2} + indel_bonus={:.2} = {:.2}", 
               overlap_length,
               (0..shared_snvs.min(config.score_array.len())).map(|i| config.score_array[i]).sum::<f32>(),
               shared_indels as f32 * config.mean_read_length * 3.0,
               weight);
        
        // Normalize weight
        let normalized_weight = weight / (config.mean_read_length * 10.0);
        
        debug!("[determine_same_haplotype] Final weight: {:.6} (normalized from {:.2})", normalized_weight, weight);
        
        // No need for error vectors when sequences are identical
        // Error vectors are only used for sequencing error detection in mismatches
        
        // Store read positions
        read_ref_pos_dict.insert(qname1.to_string(), (start, end));
        read_ref_pos_dict.insert(qname2.to_string(), (start, end));
        
        debug!("[determine_same_haplotype] SAME HAPLOTYPE (identical) - weight: {:.6}", normalized_weight);
        
        Ok((HaplotypeResult::Same, Some(normalized_weight)))
        
    } else {
        // Sequences differ - need detailed analysis
        
        // Check if sequences have different lengths (likely indel differences)
        if interval_seq1.len() != interval_seq2.len() {
            debug!("[determine_same_haplotype] Sequence length mismatch: read1_len={}, read2_len={} -> DIFFERENT haplotypes", 
                   interval_seq1.len(), interval_seq2.len());
            return Ok((HaplotypeResult::Different, None));
        }
        
        // Check if either sequence is empty
        if interval_seq1.is_empty() || interval_seq2.is_empty() {
            warn!("[determine_same_haplotype] Empty sequence(s): read1_len={}, read2_len={} -> UNKNOWN", 
                   interval_seq1.len(), interval_seq2.len());
            return Ok((HaplotypeResult::Unknown, None));
        }
        
        // Slice haplotype vectors to the interval for mismatch analysis
        let interval_hap1 = slice_hap_vector(&hap_vec1, read1.pos(), start, end);
        let interval_hap2 = slice_hap_vector(&hap_vec2, read2.pos(), start, end);
        
        // Find mismatch positions using haplotype vectors (CORRECT APPROACH)
        let mut mismatch_positions = find_mismatch_positions_from_hap_vectors(
            &interval_hap1, &interval_hap2, start, end
        );
        
        // Check if mismatches occur at indel positions - these are not tolerable
        if has_indel_mismatches_from_hap_vectors(&interval_hap1, &interval_hap2, &mismatch_positions, start) {
            debug!("[determine_same_haplotype] Found indel mismatches between reads, marking as different haplotypes");
            return Ok((HaplotypeResult::Different, None));
        }

        // Further check the shared SNV positions to see if they share the same ALT allele
        let (matching_shared_snv_pos, discrepant_shared_snv_pos) = stat_shared_snv_matches(
            &interval_hap1, &interval_hap2, start, &seq1, &ref_pos1, &seq2, &ref_pos2, read1, read2
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
        
        debug!("[determine_same_haplotype] Sequence mismatch detected between reads {} and {}", read1_id, read2_id);
        debug!("[determine_same_haplotype] Interval: {}:{}-{} (length={})", chrom, start, end, end - start);
        debug!("[determine_same_haplotype] Sequences: read1_len={}, read2_len={}", interval_seq1.len(), interval_seq2.len());
        
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
            chrom,
        );
        
        debug!("[determine_same_haplotype] Mismatch tolerance result: tolerable={}, tolerated_count={}", tolerable, tolerated_count);
        
        if tolerable {
            // All mismatches are sequencing errors - treat as same haplotype with penalty
            
            debug!("[determine_same_haplotype] All mismatches tolerable as sequencing errors - calculating weight with penalty");
            
            // Slice haplotype vectors for weight calculation (using pre-computed vectors)
            let interval_hap1 = slice_hap_vector(&hap_vec1, read1.pos(), start, end);
            let interval_hap2 = slice_hap_vector(&hap_vec2, read2.pos(), start, end);
            
            // Here we really need to identify the shared SNVs and indels
            // Not by pick the smaller SNV/Indel count in either haplotype
            let (snv1, indel1, _) = count_variants(&interval_hap1);
            let (snv2, indel2, _) = count_variants(&interval_hap2);
            let shared_snvs = matching_shared_snv_pos.len();
            let shared_indels = indel1.min(indel2);
            
            debug!("[determine_same_haplotype] Variant counts: read1(snv={}, indel={}), read2(snv={}, indel={}), shared(snv={}, indel={})", 
                   snv1, indel1, snv2, indel2, shared_snvs, shared_indels);
            
            // Calculate weight with penalty for tolerated mismatches
            let overlap_length = end - start;
            let mut weight = overlap_length as f32;
            
            // Add weight for shared variants
            for i in 0..shared_snvs.min(config.score_array.len()) {
                weight += config.score_array[i];
            }
            weight += shared_indels as f32 * config.mean_read_length * 3.0;
            
            debug!("[determine_same_haplotype] Weight before penalty: {:.2} (overlap={}, variant_bonus={:.2}, indel_bonus={:.2})", 
                   weight, overlap_length, 
                   (0..shared_snvs.min(config.score_array.len())).map(|i| config.score_array[i]).sum::<f32>(),
                   shared_indels as f32 * config.mean_read_length * 3.0);
            
            // Apply penalty for tolerated sequencing errors (like Python: weight - tolerated_count * 20)
            weight -= tolerated_count as f32 * 20.0;
            // Pick max between 0 and weight
            weight = weight.max(0.0);
            
            debug!("[determine_same_haplotype] Penalty applied: -{} x 20.0 = -{:.1}", tolerated_count, tolerated_count as f32 * 20.0);
            
            // Normalize weight
            let normalized_weight = weight / (config.mean_read_length * 10.0);
            
            debug!("[determine_same_haplotype] Final weight: {:.6} (normalized from {:.2})", normalized_weight, weight);
            
            // Update data structures - error vectors already computed above
            
            read_ref_pos_dict.insert(qname1.to_string(), (start, end));
            read_ref_pos_dict.insert(qname2.to_string(), (start, end));
            
            debug!("[determine_same_haplotype] SAME HAPLOTYPE (with penalty) - weight: {:.6}", normalized_weight);
            
            Ok((HaplotypeResult::Same, Some(normalized_weight)))
            
        } else {
            // Some mismatches cannot be explained by sequencing errors
            debug!("[determine_same_haplotype] Some mismatches are real variants - marking as DIFFERENT haplotypes");
            debug!("[determine_same_haplotype] DIFFERENT HAPLOTYPES - {} untolerable mismatches", mismatch_positions.len());
            Ok((HaplotypeResult::Different, None))
        }
    }
}


/// **RUST IMPLEMENTATION OF PYTHON'S `stat_shared_snv_matches()` FUNCTION**
/// 
/// Extract positions where two overlapping reads share the same SNV and require that the
/// alternate bases are exactly the same in both reads. No allele depth data is used here.
/// This is a critical function for proper edge weight calculation that was missing in the simplified Rust version.
/// 
/// Equivalent to Python's psv_shared_snvs() function in fp_control/pairwise_read_inspection.py lines 471-555
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
    
    debug!("[stat_shared_snv_matches] Analyzing shared SNVs between reads {} and {}", qname1, qname2);
    
    // Find indices where both vectors have SNVs (-4)
    let mut all_shared_snv_positions = Vec::new();
    let min_len = interval_hap1.len().min(interval_hap2.len());
    
    for i in 0..min_len {
        if interval_hap1[i] == -4 && interval_hap2[i] == -4 {
            let genomic_pos = interval_start + i as i64;
            all_shared_snv_positions.push(genomic_pos);
        }
    }
    
    debug!("[stat_shared_snv_matches] Found {} positions with shared SNVs", all_shared_snv_positions.len());
    
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
                debug!("[stat_shared_snv_matches] Cannot extract bases at position {} - skipping", genomic_pos);
                continue;
            }
        };
        
        // Check if both reads have the same alt base (and not N=4)
        if alt_base1 != 4 && alt_base1 == alt_base2 {
            matching_shared_snv_positions.push(*genomic_pos);
            debug!("[stat_shared_snv_matches] Both reads have same alt base {} at position {}", alt_base1, genomic_pos);
        } else if alt_base1 != 4 && alt_base2 != 4 && alt_base1 != alt_base2 {
            discrepant_shared_snv_positions.push(*genomic_pos);
            debug!("[stat_shared_snv_matches] Alt bases differ at position {}: {} vs {}", 
                   genomic_pos, alt_base1, alt_base2);
        } else {
            debug!("[stat_shared_snv_matches] One or both alt bases are N at position {}: {} vs {}", 
                   genomic_pos, alt_base1, alt_base2);
        }
    }
    
    debug!("[stat_shared_snv_matches] Final result: {} matching positions, {} discrepant positions out of {} total shared SNV positions", 
           matching_shared_snv_positions.len(), discrepant_shared_snv_positions.len(), all_shared_snv_positions.len());
    Ok((matching_shared_snv_positions, discrepant_shared_snv_positions))
}

/// Helper function to extract base at a specific genomic position from a read
/// 
/// Returns the encoded base (0=A, 1=T, 2=C, 3=G, 4=N) at the given genomic position
/// Returns None if the position is not covered by this read
fn get_base_at_position(
    query_seq: &[u8],
    ref_positions: &[i64],
    genomic_pos: i64,
) -> Option<u8> {
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