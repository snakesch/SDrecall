use std::collections::HashSet;
use log::{info, debug, warn};
use petgraph::Graph;
use petgraph::graph::NodeIndex;
use crate::structs::{
    ReadPairMap, AlleleDepthMap, PhasingGraphResult, HaplotypeConfig, 
    FastIntervals, OverlapInterval, ReadHaplotypeVector, ReadErrorVector
};
use rustc_hash::FxHashMap;
use ahash::AHashMap;

/// Build a phasing graph from BAM data for efficient haplotype identification.
///
/// This function builds a graph where each vertex represents a read pair and each edge
/// represents the confidence that two read pairs originated from the same haplotype.
/// The graph is used for subsequent haplotype identification through clique finding.
///
/// # Arguments
///
/// * `read_pair_map` - The ReadPairMap containing processed BAM data with interval trees
/// * `allele_depth_map` - Allele depth information from variant calling
/// * `intrinsic_allele_depth_map` - Optional intrinsic allele depth data 
/// * `config` - Configuration parameters for haplotype determination
/// * `low_qual_qnames` - Set of query names identified as low quality
///
/// # Returns
///
/// * `Result<PhasingGraphResult, Box<dyn std::error::Error>>` - The constructed phasing graph
///   and associated data structures, or an error if construction fails
///
/// # Notes
///
/// - The function uses petgraph library for efficient graph operations (equivalent to graph-tool)
/// - Edge weights are calculated based on two factors:
///   1. The total overlapping span between two read pairs
///   2. The number of variants (SNVs and small indels) shared by two read pairs
/// - Node indices in petgraph are automatically assigned and may differ from qname_idx
/// - The resulting graph is used for subsequent haplotype identification through clique finding
///
/// # Panics
///
/// The function may panic if there are inconsistencies in the input data structures
pub fn build_phasing_graph(
    read_pair_map: &ReadPairMap,
    allele_depth_map: &AlleleDepthMap,
    intrinsic_allele_depth_map: Option<&AlleleDepthMap>,
    config: &HaplotypeConfig,
    low_qual_qnames: &HashSet<String>,
) -> Result<PhasingGraphResult, Box<dyn std::error::Error>> {
    
    let total_read_pairs = read_pair_map.readpair_dict.len();
    info!("Building phasing graph with {} read pairs, mean read length: {}", 
          total_read_pairs, config.mean_read_length);
    
    let mut result = PhasingGraphResult::new();
    
    // Track which node pairs we've already checked to avoid duplicates
    // Using a more efficient representation than Python's qname_check_dict
    let mut checked_pairs: HashSet<(NodeIndex, NodeIndex)> = HashSet::new();
    
    // Create weight matrix with initial capacity
    // In Rust, we use HashMap instead of numpy array for sparse representation
    result.weight_matrix.reserve(total_read_pairs * total_read_pairs / 4); // Rough estimate
    
    // Iterate through all read pairs to build the graph
    // This is equivalent to Python's: for qname_idx, paired_reads in ncls_read_dict.items()
    for (&qname_idx, read_pair) in &read_pair_map.readpair_dict {
        let qname = &read_pair.qname;
        
        // Add node to graph if it doesn't exist
        // In petgraph, NodeIndex is automatically assigned (unlike Python where we control it)
        let node_idx = if let Some(&existing_node) = result.qname_to_node.get(qname) {
            existing_node
        } else {
            let new_node = result.graph.add_node(qname.clone());
            result.qname_to_node.insert(qname.clone(), new_node);
            new_node
        };
        
        // Get chromosome and position info for overlap queries
        let chrom = get_read_chromosome(read_pair)?;
        let (start, end) = get_read_pair_span(read_pair);
        
        // Find overlapping read pairs using the interval tree
        // This replaces Python's overlap_qname_idx_iterator
        let overlapping_qname_indices = read_pair_map
            .find_overlapping_qname_indices(&chrom, start, end)?;
        
        // Process each overlapping read pair
        for other_qname_idx in overlapping_qname_indices {
            // Skip self-comparison
            if qname_idx == other_qname_idx {
                continue;
            }
            
            // Get the other read pair
            let other_read_pair = match read_pair_map.readpair_dict.get(&other_qname_idx) {
                Some(pair) => pair,
                None => {
                    warn!("Could not find read pair for qname_idx: {}", other_qname_idx);
                    continue;
                }
            };
            
            let other_qname = &other_read_pair.qname;
            
            // Add other node to graph if it doesn't exist
            let other_node_idx = if let Some(&existing_node) = result.qname_to_node.get(other_qname) {
                existing_node
            } else {
                let new_node = result.graph.add_node(other_qname.clone());
                result.qname_to_node.insert(other_qname.clone(), new_node);
                new_node
            };
            
            // Check if we've already processed this pair
            // Use ordered pair to avoid checking (A,B) and (B,A) separately
            let pair_key = if node_idx < other_node_idx {
                (node_idx, other_node_idx)
            } else {
                (other_node_idx, node_idx)
            };
            
            if checked_pairs.contains(&pair_key) {
                continue;
            }
            checked_pairs.insert(pair_key);
            
            // Skip low quality read pairs (equivalent to Python's total_lowqual_qnames check)
            if low_qual_qnames.contains(qname) || low_qual_qnames.contains(other_qname) {
                // Mark as negative weight in matrix (indicates low quality)
                result.weight_matrix.insert((node_idx.index(), other_node_idx.index()), -1.0);
                result.weight_matrix.insert((other_node_idx.index(), node_idx.index()), -1.0);
                debug!("Skipping low quality read pair comparison: {} vs {}", qname, other_qname);
                continue;
            }
            
            // Find overlap intervals between the two read pairs
            // This replaces Python's get_overlap_intervals function
            let overlap_intervals = find_overlap_intervals(read_pair, other_read_pair)?;
            
            if overlap_intervals.is_empty() {
                debug!("No overlap intervals found between {} and {}", qname, other_qname);
                continue;
            }
            
            debug!("Found {} overlap intervals between {} and {}", 
                   overlap_intervals.len(), qname, other_qname);
            
            // Process overlapping regions to determine haplotype compatibility
            let edge_weight = process_overlap_regions(
                &overlap_intervals,
                read_pair,
                other_read_pair,
                allele_depth_map,
                intrinsic_allele_depth_map,
                config,
                &mut result.read_hap_vectors,
                &mut result.read_error_vectors,
                &mut result.read_ref_pos_dict,
            )?;
            
            // Add edge to graph based on compatibility
            match edge_weight {
                Some(weight) if weight > config.edge_weight_cutoff => {
                    // Reads are compatible - add edge to graph
                    let edge = result.graph.add_edge(node_idx, other_node_idx, weight);
                    result.weight_matrix.insert((node_idx.index(), other_node_idx.index()), weight);
                    result.weight_matrix.insert((other_node_idx.index(), node_idx.index()), weight);
                    debug!("Added edge between {} and {} with weight {}", qname, other_qname, weight);
                }
                Some(_) => {
                    // Weight too low - don't add edge but record the weight
                    debug!("Edge weight too low between {} and {}", qname, other_qname);
                }
                None => {
                    // Reads are incompatible - mark with negative weight
                    result.weight_matrix.insert((node_idx.index(), other_node_idx.index()), -1.0);
                    result.weight_matrix.insert((other_node_idx.index(), node_idx.index()), -1.0);
                    debug!("Reads {} and {} are incompatible", qname, other_qname);
                }
            }
        }
    }
    
    info!("Finished building phasing graph: {} vertices, {} edges", 
          result.vertex_count(), result.edge_count());
    
    Ok(result)
}

/// Extract chromosome name from a read pair
fn get_read_chromosome(read_pair: &crate::structs::ReadPair) -> Result<String, Box<dyn std::error::Error>> {
    // Use rust-htslib to get chromosome name
    // This is equivalent to Python's paired_reads[0].reference_name
    if let Some(tid) = read_pair.read1.tid() {
        // In a real implementation, you'd need access to the BAM header to convert tid to chromosome name
        // For now, we'll use a placeholder - this needs to be passed from the calling function
        Ok(format!("chr{}", tid)) // Placeholder - needs proper header lookup
    } else {
        Err("Read has no valid chromosome assignment".into())
    }
}

/// Calculate the genomic span covered by a read pair
fn get_read_pair_span(read_pair: &crate::structs::ReadPair) -> (i64, i64) {
    let read1_start = read_pair.read1.reference_start();
    let read1_end = read_pair.read1.reference_end();
    
    let (read2_start, read2_end) = if let Some(ref read2) = read_pair.read2 {
        (read2.reference_start(), read2.reference_end())
    } else {
        (read1_start, read1_end) // If no read2, use read1 coordinates
    };
    
    // Return the span covering both reads
    // Equivalent to Python's: start = min(r.reference_start for r in paired_reads)
    //                         end = max(r.reference_end for r in paired_reads)
    (read1_start.min(read2_start), read1_end.max(read2_end))
}

/// Find overlap intervals between two read pairs
/// This replaces Python's get_overlap_intervals function
fn find_overlap_intervals(
    read_pair1: &crate::structs::ReadPair,
    read_pair2: &crate::structs::ReadPair,
) -> Result<Vec<OverlapInterval>, Box<dyn std::error::Error>> {
    let mut overlaps = Vec::new();
    
    // Check all possible combinations of reads between the two pairs
    let reads1 = vec![&read_pair1.read1];
    let reads1 = if let Some(ref read2) = read_pair1.read2 {
        vec![&read_pair1.read1, read2]
    } else {
        vec![&read_pair1.read1]
    };
    
    let reads2 = if let Some(ref read2) = read_pair2.read2 {
        vec![&read_pair2.read1, read2]
    } else {
        vec![&read_pair2.read1]
    };
    
    for read1 in &reads1 {
        for read2 in &reads2 {
            // Check if reads overlap
            let start1 = read1.reference_start();
            let end1 = read1.reference_end();
            let start2 = read2.reference_start();
            let end2 = read2.reference_end();
            
            // Calculate overlap
            let overlap_start = start1.max(start2);
            let overlap_end = end1.min(end2);
            
            if overlap_start < overlap_end {
                overlaps.push(OverlapInterval {
                    start: overlap_start,
                    end: overlap_end,
                    read1_ref: format!("{}:{}-{}", read1.tid().unwrap_or(0), start1, end1),
                    read2_ref: format!("{}:{}-{}", read2.tid().unwrap_or(0), start2, end2),
                });
            }
        }
    }
    
    Ok(overlaps)
}

/// Process overlap regions to determine haplotype compatibility
/// This is the main logic equivalent to Python's nested loop processing
fn process_overlap_regions(
    overlap_intervals: &[OverlapInterval],
    read_pair1: &crate::structs::ReadPair,
    read_pair2: &crate::structs::ReadPair,
    allele_depth_map: &AlleleDepthMap,
    intrinsic_allele_depth_map: Option<&AlleleDepthMap>,
    config: &HaplotypeConfig,
    read_hap_vectors: &mut AHashMap<String, ReadHaplotypeVector>,
    read_error_vectors: &mut AHashMap<String, ReadErrorVector>,
    read_ref_pos_dict: &mut AHashMap<String, (i64, i64)>,
) -> Result<Option<f32>, Box<dyn std::error::Error>> {
    
    // Track haplotype compatibility results
    // Equivalent to Python's qname_bools array
    let mut compatibility_results = Vec::new();
    let mut total_weight = 0.0f32;
    
    // Track inspected overlaps to avoid double-counting
    // Equivalent to Python's inspected_overlaps = FastIntervals(max_size = len(overlap_intervals))
    let mut inspected_overlaps = FastIntervals::new(overlap_intervals.len());
    
    for overlap in overlap_intervals {
        // Find uncovered regions within this overlap
        // This replaces Python's find_uncovered_regions_numba call
        let uncovered_regions = find_uncovered_regions(&inspected_overlaps, overlap.start, overlap.end);
        
        for (uncovered_start, uncovered_end) in uncovered_regions {
            // Determine if reads come from same haplotype in this region
            // This is the PLACEHOLDER for the complex determine_same_haplotype function
            let haplotype_result = determine_same_haplotype_placeholder(
                read_pair1,
                read_pair2,
                uncovered_start,
                uncovered_end,
                allele_depth_map,
                intrinsic_allele_depth_map,
                config,
                read_hap_vectors,
                read_error_vectors,
                read_ref_pos_dict,
            )?;
            
            match haplotype_result {
                HaplotypeCompatibility::Same(weight) => {
                    compatibility_results.push(1);
                    total_weight += weight;
                    debug!("Reads compatible in region {}:{}-{} with weight {}", 
                           overlap.start, uncovered_start, uncovered_end, weight);
                }
                HaplotypeCompatibility::Different => {
                    compatibility_results.push(-1);
                    debug!("Reads incompatible in region {}:{}-{}", 
                           overlap.start, uncovered_start, uncovered_end);
                }
                HaplotypeCompatibility::Unknown => {
                    compatibility_results.push(0);
                    debug!("Unknown compatibility in region {}:{}-{}", 
                           overlap.start, uncovered_start, uncovered_end);
                }
            }
        }
        
        // Mark this overlap as inspected
        inspected_overlaps.add(overlap.start, overlap.end);
    }
    
    // Determine final compatibility based on all results
    // Equivalent to Python's any_false_numba(qname_bools) check
    if compatibility_results.iter().any(|&result| result == -1) {
        // Any incompatible region means overall incompatibility
        Ok(None)
    } else if compatibility_results.iter().all(|&result| result >= 0) {
        // All regions compatible or unknown
        let normalized_weight = if total_weight > 0.0 {
            total_weight / (config.mean_read_length * 10.0)
        } else {
            1e-4  // Minimum weight for compatible pairs
        };
        Ok(Some(normalized_weight))
    } else {
        // Mixed results - treat as incompatible for safety
        Ok(None)
    }
}

/// Result of haplotype compatibility analysis
#[derive(Debug, Clone)]
enum HaplotypeCompatibility {
    Same(f32),    // Same haplotype with confidence weight
    Different,    // Different haplotypes
    Unknown,      // Cannot determine (insufficient data)
}

/// PLACEHOLDER: Determine if two read pairs come from the same haplotype
/// 
/// This is a complex function that needs to be implemented based on the original Python logic.
/// It involves variant calling, quality assessment, and statistical analysis.
/// 
/// # Arguments
/// 
/// * `read_pair1` - First read pair
/// * `read_pair2` - Second read pair  
/// * `start` - Start position of analysis region
/// * `end` - End position of analysis region
/// * `allele_depth_map` - Variant depth information
/// * `intrinsic_allele_depth_map` - Optional intrinsic variant data
/// * `config` - Analysis configuration
/// * `read_hap_vectors` - Mutable reference to haplotype vectors
/// * `read_error_vectors` - Mutable reference to error vectors
/// * `read_ref_pos_dict` - Mutable reference to position dictionary
///
/// # Returns
///
/// * `HaplotypeCompatibility` - Result of compatibility analysis
fn determine_same_haplotype_placeholder(
    _read_pair1: &crate::structs::ReadPair,
    _read_pair2: &crate::structs::ReadPair,
    _start: i64,
    _end: i64,
    _allele_depth_map: &AlleleDepthMap,
    _intrinsic_allele_depth_map: Option<&AlleleDepthMap>,
    _config: &HaplotypeConfig,
    _read_hap_vectors: &mut AHashMap<String, ReadHaplotypeVector>,
    _read_error_vectors: &mut AHashMap<String, ReadErrorVector>, 
    _read_ref_pos_dict: &mut AHashMap<String, (i64, i64)>,
) -> Result<HaplotypeCompatibility, Box<dyn std::error::Error>> {
    
    // PLACEHOLDER IMPLEMENTATION
    // TODO: Implement the full logic from Python's determine_same_haplotype function
    // This should include:
    // 1. Variant extraction from reads in the specified region
    // 2. Quality assessment and filtering
    // 3. Allele depth comparison with reference data
    // 4. Statistical analysis to determine compatibility
    // 5. Weight calculation based on evidence strength
    
    warn!("Using placeholder haplotype determination - implement full logic!");
    
    // For now, return a dummy result based on simple overlap
    let overlap_length = _end - _start;
    if overlap_length > 50 {
        Ok(HaplotypeCompatibility::Same(overlap_length as f32 / 100.0))
    } else {
        Ok(HaplotypeCompatibility::Unknown)
    }
}

/// Find uncovered regions within an overlap interval
/// This replaces Python's find_uncovered_regions_numba function
fn find_uncovered_regions(
    inspected_overlaps: &FastIntervals,
    query_start: i64,
    query_end: i64,
) -> Vec<(i64, i64)> {
    let (starts, ends) = inspected_overlaps.get_intervals();
    
    if starts.is_empty() {
        // No previous overlaps - entire region is uncovered
        return vec![(query_start, query_end)];
    }
    
    let mut uncovered = Vec::new();
    let mut current_pos = query_start;
    
    // Sort intervals by start position for processing
    let mut intervals: Vec<(i64, i64)> = starts.iter().zip(ends.iter()).map(|(&s, &e)| (s, e)).collect();
    intervals.sort_by_key(|&(start, _)| start);
    
    for (start, end) in intervals {
        // Skip intervals that don't overlap with our query region
        if end <= query_start || start >= query_end {
            continue;
        }
        
        // Clip interval to query region
        let clipped_start = start.max(query_start);
        let clipped_end = end.min(query_end);
        
        // Add uncovered region before this interval
        if current_pos < clipped_start {
            uncovered.push((current_pos, clipped_start));
        }
        
        // Update current position
        current_pos = clipped_end.max(current_pos);
    }
    
    // Add final uncovered region
    if current_pos < query_end {
        uncovered.push((current_pos, query_end));
    }
    
    uncovered
}
