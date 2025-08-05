/// Haplotype determination logic
/// 
/// This module contains the core algorithm for determining whether two read pairs
/// come from the same haplotype by analyzing their sequences and CIGAR operations.

use std::collections::{HashSet, HashMap};
use log::{info, debug, warn};
use ahash::AHashMap;
use rust_htslib::bam::Record;
use rust_htslib::bam::Read;

use crate::structs::{
    ReadPair, AlleleDepthMap, HaplotypeConfig, ReadHaplotypeVector, 
    ReadErrorVector, Variant, VariantCompatibility
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
/// Returns (query_sequence, ref_positions) where:
/// - query_sequence: Vec<u8> of base-encoded sequence  
/// - ref_positions: Vec<i64> mapping query positions to reference positions
pub fn extract_query_seq(record: &Record) -> Result<(Vec<u8>, Vec<i64>), Box<dyn std::error::Error>> {
    let query_seq = record.seq().as_bytes();
    let ref_positions = record.reference_positions_full();
    
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
    
    // Convert reference positions, handling None as -1
    let ref_pos_vec: Vec<i64> = ref_positions.iter().map(|&pos_opt| {
        pos_opt.map(|p| p as i64).unwrap_or(-1)
    }).collect();
    
    Ok((encoded_seq, ref_pos_vec))
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
    
    for (i, (&base1, &base2)) in seq1.iter().zip(seq2.iter()).enumerate() {
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
        use rust_htslib::bam::record::CigarOp;
        match op {
            CigarOp::Match(len) | CigarOp::Equal(len) => {
                // Match to reference
                for _ in 0..len {
                    hap_vector.push(1);
                }
            }
            CigarOp::Diff(len) => {
                // Mismatch (SNV)
                for _ in 0..len {
                    hap_vector.push(-4);
                }
            }
            CigarOp::Ins(len) => {
                // Insertion - mark previous position with insertion length
                if !hap_vector.is_empty() {
                    let last_idx = hap_vector.len() - 1;
                    hap_vector[last_idx] = len as i16 * 4; // Encode insertion
                }
            }
            CigarOp::Del(len) => {
                // Deletion
                for _ in 0..len {
                    hap_vector.push(-6);
                }
            }
            CigarOp::SoftClip(_) | CigarOp::HardClip(_) => {
                // Skip clipped bases
            }
            _ => {
                // Handle other operations as matches for now
                if let Some(len) = op.len() {
                    for _ in 0..len {
                        hap_vector.push(1);
                    }
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

/// Determine if two reads come from the same haplotype
/// 
/// Main function that orchestrates the haplotype comparison workflow
pub fn determine_same_haplotype(
    read1: &Record,
    read2: &Record,
    start: i64,
    end: i64,
    chrom: &str,
    allele_depth_map: &AlleleDepthMap,
    intrinsic_allele_depth_map: Option<&AlleleDepthMap>,
    config: &HaplotypeConfig,
    read_hap_vectors: &mut AHashMap<String, ReadHaplotypeVector>,
    read_error_vectors: &mut AHashMap<String, ReadErrorVector>, 
    read_ref_pos_dict: &mut AHashMap<String, (i64, i64)>,
) -> Result<(HaplotypeResult, Option<f32>), Box<dyn std::error::Error>> {
    // Extract qnames from records
    let qname1 = std::str::from_utf8(read1.qname())?;
    let qname2 = std::str::from_utf8(read2.qname())?;
    
    // Step 1: Extract query sequences and reference positions
    let (seq1, ref_pos1) = extract_query_seq(read1)?;
    let (seq2, ref_pos2) = extract_query_seq(read2)?;
    
    // Step 2: Slice sequences to the genomic interval
    let interval_seq1 = slice_seq_to_interval(&seq1, &ref_pos1, start, end);
    let interval_seq2 = slice_seq_to_interval(&seq2, &ref_pos2, start, end);
    
    // Step 3: Compare sequences for exact match
    if compare_sequences(&interval_seq1, &interval_seq2) {
        // Sequences are identical - likely same haplotype
        
        // Step 4: Extract haplotype vectors for weight calculation
        let hap_vec1 = extract_hap_vector(read1);
        let hap_vec2 = extract_hap_vector(read2);
        
        // Step 5: Slice haplotype vectors to interval
        let interval_hap1 = slice_hap_vector(&hap_vec1, read1.pos(), start, end);
        let interval_hap2 = slice_hap_vector(&hap_vec2, read2.pos(), start, end);
        
        // Step 6: Count shared variants for weight calculation
        let (snv1, indel1, _) = count_variants(&interval_hap1);
        let (snv2, indel2, _) = count_variants(&interval_hap2);
        
        // Use minimum counts as shared variants (conservative estimate)
        let shared_snvs = snv1.min(snv2);
        let shared_indels = indel1.min(indel2);
        
        // Calculate weight based on overlap and shared variants
        let overlap_length = end - start;
        let mut weight = overlap_length as f32;
        
        // Add weight for shared variants using score array
        for i in 0..shared_snvs.min(config.score_array.len()) {
            weight += config.score_array[i];
        }
        
        // Add weight for shared indels
        weight += shared_indels as f32 * config.mean_read_length * 3.0;
        
        // Normalize weight
        let normalized_weight = weight / (config.mean_read_length * 10.0);
        
        // Update data structures
        let overlap_size = overlap_length as usize;
        let read1_id = format!("{}:{}", qname1, read1.flags());
        let read2_id = format!("{}:{}", qname2, read2.flags());
        
        read_hap_vectors.insert(read1_id, vec![1; overlap_size]);
        read_error_vectors.insert(read2_id, vec![0.01; overlap_size]);
        
        // Store read positions
        read_ref_pos_dict.insert(qname1.to_string(), (start, end));
        read_ref_pos_dict.insert(qname2.to_string(), (start, end));
        
        Ok((HaplotypeResult::Same, Some(normalized_weight)))
        
    } else {
        // Sequences differ - different haplotypes or insufficient data
        if interval_seq1.is_empty() || interval_seq2.is_empty() {
            Ok((HaplotypeResult::Unknown, None))
        } else {
            Ok((HaplotypeResult::Different, None))
        }
    }
}

/// Compare variants between two specific reads
pub fn compare_variants(
    variants1: &[Variant],
    variants2: &[Variant],
    allele_depth_map: &AlleleDepthMap,
    intrinsic_ad_map: Option<&AlleleDepthMap>,
) -> Result<VariantCompatibility, Box<dyn std::error::Error>> {
    // Group variants by position
    let mut pos_to_var1: HashMap<i64, Vec<&Variant>> = HashMap::new();
    let mut pos_to_var2: HashMap<i64, Vec<&Variant>> = HashMap::new();
    
    for var in variants1 {
        pos_to_var1.entry(var.position()).or_default().push(var);
    }
    
    for var in variants2 {
        pos_to_var2.entry(var.position()).or_default().push(var);
    }
    
    let mut shared_count = 0;
    let mut incompatible = false;
    
    // Check each position where either read has a variant
    let all_positions: HashSet<i64> = pos_to_var1.keys()
        .chain(pos_to_var2.keys())
        .cloned()
        .collect();
    
    for pos in all_positions {
        let vars1 = pos_to_var1.get(&pos);
        let vars2 = pos_to_var2.get(&pos);
        
        match (vars1, vars2) {
            (Some(v1), Some(v2)) => {
                // Both have variants at this position - check compatibility
                if are_variants_compatible(v1, v2) {
                    shared_count += 1;
                } else {
                    incompatible = true;
                    break;
                }
            }
            (Some(_), None) | (None, Some(_)) => {
                // Only one has a variant - could be heterozygous
                // This is generally okay unless we have strong evidence otherwise
            }
            _ => unreachable!(),
        }
    }
    
    if incompatible {
        Ok(VariantCompatibility::Incompatible)
    } else if shared_count > 0 {
        Ok(VariantCompatibility::Compatible(shared_count))
    } else {
        Ok(VariantCompatibility::NoData)
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