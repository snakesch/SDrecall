/// Haplotype determination logic
/// 
/// This module contains the core algorithm for determining whether two read pairs
/// come from the same haplotype by analyzing their variants and overlap regions.

use std::collections::{HashSet, HashMap};
use log::{info, debug, warn};
use ahash::AHashMap;
use rust_htslib::bam::Record;

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

/// Determine if two read pairs come from the same haplotype
/// 
/// This is the main logic for comparing two read pairs in a specific genomic region
/// to determine if they likely originated from the same haplotype. This is a critical
/// function that implements the core algorithm for phasing.
/// 
/// # Arguments
/// 
/// * `read_pair1` - First read pair
/// * `read_pair2` - Second read pair  
/// * `start` - Start position of analysis region
/// * `end` - End position of analysis region
/// * `chrom` - Chromosome name
/// * `allele_depth_map` - Variant depth information
/// * `intrinsic_allele_depth_map` - Optional intrinsic variant data
/// * `config` - Analysis configuration
/// * `read_hap_vectors` - Mutable reference to haplotype vectors
/// * `read_error_vectors` - Mutable reference to error vectors
/// * `read_ref_pos_dict` - Mutable reference to position dictionary
/// * `low_qual_qnames` - Set of low quality read names (will be updated)
///
/// # Returns
///
/// * `(HaplotypeResult, Option<f32>)` - Result of compatibility analysis and optional weight
pub fn determine_same_haplotype(
    read_pair1: &ReadPair,
    read_pair2: &ReadPair,
    start: i64,
    end: i64,
    chrom: &str,
    allele_depth_map: &AlleleDepthMap,
    intrinsic_allele_depth_map: Option<&AlleleDepthMap>,
    config: &HaplotypeConfig,
    read_hap_vectors: &mut AHashMap<String, ReadHaplotypeVector>,
    read_error_vectors: &mut AHashMap<String, ReadErrorVector>, 
    read_ref_pos_dict: &mut AHashMap<String, (i64, i64)>,
    low_qual_qnames: &mut HashSet<String>,
) -> Result<(HaplotypeResult, Option<f32>), Box<dyn std::error::Error>> {
    use rust_htslib::bam::Read;
    use std::cmp::{min, max};
    
    // Get read IDs for logging and storage
    let read1_id = format!("{}:{}", 
        std::str::from_utf8(read_pair1.read1.qname())?, 
        read_pair1.read1.flags()
    );
    let read2_id = read_pair1.read2.as_ref().map(|r| 
        format!("{}:{}", std::str::from_utf8(r.qname()).unwrap_or(""), r.flags())
    );
    
    let other_read1_id = format!("{}:{}", 
        std::str::from_utf8(read_pair2.read1.qname())?, 
        read_pair2.read1.flags()
    );
    let other_read2_id = read_pair2.read2.as_ref().map(|r| 
        format!("{}:{}", std::str::from_utf8(r.qname()).unwrap_or(""), r.flags())
    );
    
    // Process all reads in both pairs
    let reads1 = if let Some(ref r2) = read_pair1.read2 {
        vec![&read_pair1.read1, r2]
    } else {
        vec![&read_pair1.read1]
    };
    
    let reads2 = if let Some(ref r2) = read_pair2.read2 {
        vec![&read_pair2.read1, r2]
    } else {
        vec![&read_pair2.read1]
    };
    
    // Track variants found in the overlap region
    let mut pair1_variants = Vec::new();
    let mut pair2_variants = Vec::new();
    
    // Analyze reads from first pair
    for read in &reads1 {
        let read_start = read.pos();
        let read_end = read.cigar().end_pos();
        
        // Skip if read doesn't overlap with analysis region
        if read_end <= start || read_start >= end {
            continue;
        }
        
        // Extract variants in the overlap region
        let variants = extract_variants_in_region(read, max(start, read_start), min(end, read_end))?;
        
        // Check read quality
        if is_noisy_read(read, &variants) {
            low_qual_qnames.insert(read_pair1.qname.clone());
            return Ok((HaplotypeResult::Unknown, None));
        }
        
        pair1_variants.extend(variants);
    }
    
    // Analyze reads from second pair
    for read in &reads2 {
        let read_start = read.pos();
        let read_end = read.cigar().end_pos();
        
        // Skip if read doesn't overlap with analysis region
        if read_end <= start || read_start >= end {
            continue;
        }
        
        // Extract variants in the overlap region
        let variants = extract_variants_in_region(read, max(start, read_start), min(end, read_end))?;
        
        // Check read quality
        if is_noisy_read(read, &variants) {
            low_qual_qnames.insert(read_pair2.qname.clone());
            return Ok((HaplotypeResult::Unknown, None));
        }
        
        pair2_variants.extend(variants);
    }
    
    // Compare variants between the two pairs
    let compatibility = compare_variants(&pair1_variants, &pair2_variants, allele_depth_map, intrinsic_allele_depth_map)?;
    
    // Calculate weight based on overlap and shared variants
    let overlap_length = end - start;
    let base_weight = overlap_length as f32;
    
    // Adjust weight based on shared variants
    let variant_weight = match compatibility {
        VariantCompatibility::Compatible(shared_count) => {
            // Add weight for each shared variant
            let mut weight = base_weight;
            for i in 0..shared_count.min(config.score_array.len()) {
                weight += config.score_array[i];
            }
            weight
        }
        VariantCompatibility::Incompatible => {
            // Incompatible variants indicate different haplotypes
            return Ok((HaplotypeResult::Different, None));
        }
        VariantCompatibility::NoData => {
            // Not enough data to make a determination
            return Ok((HaplotypeResult::Unknown, None));
        }
    };
    
    // Normalize weight
    let normalized_weight = variant_weight / (config.mean_read_length * 10.0);
    
    // Update data structures
    // Store haplotype vectors (simplified for now)
    read_hap_vectors.insert(read1_id.clone(), vec![1; overlap_length as usize]);
    if let Some(id) = read2_id {
        read_hap_vectors.insert(id, vec![1; overlap_length as usize]);
    }
    
    // Store error vectors (placeholder)
    read_error_vectors.insert(read1_id.clone(), vec![0.01; overlap_length as usize]);
    
    // Store read positions
    read_ref_pos_dict.insert(read_pair1.qname.clone(), (start, end));
    read_ref_pos_dict.insert(read_pair2.qname.clone(), (start, end));
    
    Ok((HaplotypeResult::Same, Some(normalized_weight)))
}

/// Extract variants from a read in a specific region
pub fn extract_variants_in_region(
    read: &Record, 
    region_start: i64, 
    region_end: i64
) -> Result<Vec<Variant>, Box<dyn std::error::Error>> {
    let mut variants = Vec::new();
    let cigar = read.cigar();
    let seq = read.seq();
    let qual = read.qual();
    
    let mut ref_pos = read.pos();
    let mut read_pos = 0;
    
    // Parse CIGAR to find variants
    for &op in cigar.iter() {
        use rust_htslib::bam::record::CigarOp;
        match op {
            CigarOp::Match(len) | CigarOp::Equal(len) => {
                ref_pos += len as i64;
                read_pos += len as u32;
            }
            CigarOp::Diff(len) => {
                // Mismatch - potential SNV
                for i in 0..len {
                    let pos = ref_pos + i as i64;
                    if pos >= region_start && pos < region_end {
                        let base = seq[(read_pos + i) as usize];
                        let quality = qual[(read_pos + i) as usize];
                        if base != b'N' && quality >= 15 {  // Quality threshold
                            variants.push(Variant::Snv { 
                                position: pos, 
                                alt_base: base,
                                quality,
                            });
                        }
                    }
                }
                ref_pos += len as i64;
                read_pos += len;
            }
            CigarOp::Ins(len) => {
                // Insertion
                if ref_pos >= region_start && ref_pos < region_end {
                    variants.push(Variant::Insertion { 
                        position: ref_pos,
                        length: len,
                    });
                }
                read_pos += len;
            }
            CigarOp::Del(len) => {
                // Deletion
                if ref_pos >= region_start && ref_pos < region_end {
                    variants.push(Variant::Deletion { 
                        position: ref_pos,
                        length: len,
                    });
                }
                ref_pos += len as i64;
            }
            CigarOp::SoftClip(len) => {
                read_pos += len;
            }
            _ => {} // Skip other operations
        }
    }
    
    Ok(variants)
}

/// Check if a read has too many errors to be reliable
pub fn is_noisy_read(read: &Record, variants: &[Variant]) -> bool {
    // Count high-quality mismatches
    let high_qual_mismatches = variants.iter()
        .filter(|v| matches!(v, Variant::Snv { quality, .. } if *quality > 20))
        .count();
    
    // If more than 3 high-quality SNVs, consider it noisy
    high_qual_mismatches >= 3
}

/// Compare variants between two read pairs
pub fn compare_variants(
    variants1: &[Variant],
    variants2: &[Variant],
    allele_depth_map: &AlleleDepthMap,
    intrinsic_ad_map: Option<&AlleleDepthMap>,
) -> Result<VariantCompatibility, Box<dyn std::error::Error>> {
    use std::collections::HashMap;
    
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