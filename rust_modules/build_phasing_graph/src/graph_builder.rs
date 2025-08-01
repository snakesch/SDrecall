use std::collections::HashSet;
use log::{info, debug, warn};
use petgraph::Graph;
use petgraph::graph::NodeIndex;
use crate::structs::{
    ReadPairMap, AlleleDepthMap, PhasingGraphResult, HaplotypeConfig, 
    FastIntervals, OverlapInterval, ReadHaplotypeVector, ReadErrorVector
};
use crate::haplotype_determination::{determine_same_haplotype, HaplotypeResult};
use rustc_hash::FxHashMap;
use ahash::AHashMap;
use rust_htslib::bam::Record;

/// Build a phasing graph from BAM data for efficient haplotype identification.
///
/// This function constructs a graph where each vertex represents a read pair and each edge
/// represents the confidence that two read pairs originated from the same haplotype.
/// The graph is used for subsequent haplotype identification through clique finding.
///
/// # Arguments
///
/// * `read_pair_map` - The ReadPairMap containing processed BAM data with interval trees
/// * `allele_depth_map` - Allele depth information from variant calling
/// * `intrinsic_allele_depth_map` - Optional intrinsic allele depth data 
/// * `config` - Configuration parameters for haplotype determination
/// * `low_qual_qnames` - Mutable set of query names identified as low quality (will be updated)
///
/// # Returns
///
/// * `Result<PhasingGraphResult, Box<dyn std::error::Error>>` - The constructed phasing graph
///   and associated data structures, or an error if construction fails
///
/// # Notes
///
/// - The function uses petgraph library for efficient graph operations (equivalent to Python's graph-tool)
/// - Edge weights are calculated based on two factors:
///   1. The total overlapping span between two read pairs
///   2. The number of variants (SNVs and small indels) shared by two read pairs
/// - Node indices in petgraph are automatically assigned and may differ from qname_idx
/// - We maintain a mapping between qname_idx and NodeIndex for cross-referencing
/// - The resulting graph is used for subsequent haplotype identification through clique finding (Greedy-Clique-Expansion algorithm)
///
/// # Differences from Python
///
/// - Python uses graph-tool which allows custom vertex indices, but petgraph auto-assigns NodeIndex
/// - Python uses numpy arrays for weight matrix, we use FxHashMap for sparse representation
/// - Python's qname_check_dict becomes checked_pairs HashSet with NodeIndex tuples
pub fn build_phasing_graph(
    read_pair_map: &ReadPairMap,
    allele_depth_map: &AlleleDepthMap,
    intrinsic_allele_depth_map: Option<&AlleleDepthMap>,
    config: &HaplotypeConfig,
    low_qual_qnames: &mut HashSet<String>,
) -> Result<PhasingGraphResult, Box<dyn std::error::Error>> {
    
    let total_read_pairs = read_pair_map.readpair_dict.len();
    info!("There are totally {} pair of reads, mean read length is {:.1}. with adequate mapping or base quality which can be used to build the graph", 
          total_read_pairs, config.mean_read_length);
    
    // Check if no ALT alleles found (equivalent to Python's early return)
    if allele_depth_map.is_empty() {
        warn!("No ALT allele found in this BAM file. Skip this entire script");
        return Ok(PhasingGraphResult::new());
    }
    
    let mut result = PhasingGraphResult::new();
    
    // Create a mapping from qname_idx to NodeIndex for cross-referencing
    // This is necessary because petgraph doesn't let us control node indices
    let mut qname_idx_to_node: FxHashMap<usize, NodeIndex> = FxHashMap::default();
    
    // Create a set to track which node pairs we've already checked
    // Equivalent to Python's qname_check_dict but using NodeIndex pairs
    let mut checked_pairs: HashSet<(NodeIndex, NodeIndex)> = HashSet::new();
    
    // Initialize weight matrix with estimated capacity
    // Using sparse representation (HashMap) instead of dense numpy array
    let estimated_edges = total_read_pairs * total_read_pairs / 4; // Rough estimate
    result.weight_matrix.reserve(estimated_edges);
    
    // Use NCLS read dict to iterate through all the reads to build a graph
    // One good thing is that every key: value stores a pair of read objects
    for (&qname_idx, read_pair) in &read_pair_map.readpair_dict {
        let qname = &read_pair.qname;
        
        // Verify all reads in the pair have the same query name (assertion from Python)
        debug_assert!(
            read_pair.read2.as_ref().map_or(true, |r2| 
                std::str::from_utf8(read_pair.read1.qname()).ok() == 
                std::str::from_utf8(r2.qname()).ok()
            ),
            "The query names of the reads in the paired_reads are not the same"
        );
        
        // Check if the query name already exists in the graph
        let node_idx = if let Some(&existing_node) = result.qname_to_node.get(qname) {
            existing_node
        } else {
            // Add a new node to the graph
            let new_node = result.graph.add_node(qname.clone());
            result.qname_to_node.insert(qname.clone(), new_node);
            qname_idx_to_node.insert(qname_idx, new_node);
            debug!("Added a new node {:?} for qname {} to the graph", new_node, qname);
            new_node
        };
        
        // Get chromosome and position info for overlap queries
        let chrom = get_read_chromosome(read_pair, &read_pair_map)?;
        let (start, end) = get_read_pair_span(read_pair);
        
        // Find overlapping read pairs using the interval tree
        // This replaces Python's overlap_qname_idx_iterator function
        let overlapping_qname_indices = read_pair_map
            .find_overlapping_qname_indices(&chrom, start, end)?;
        
        // Iterate through the overlapping reads
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
            
            // Check if the other query name already exists in the graph
            let other_node_idx = if let Some(&existing_node) = result.qname_to_node.get(other_qname) {
                existing_node
            } else {
                // Add a new node to the graph
                let new_node = result.graph.add_node(other_qname.clone());
                result.qname_to_node.insert(other_qname.clone(), new_node);
                qname_idx_to_node.insert(other_qname_idx, new_node);
                debug!("Added a new node {:?} for qname {} to the graph", new_node, other_qname);
                new_node
            };
            
            // Check if we've already inspected this edge
            // Use ordered pair to avoid checking (A,B) and (B,A) separately
            let inspect_res = check_edge(node_idx, other_node_idx, &checked_pairs);
            if inspect_res {
                continue;
            }
            
            // Mark this pair as checked
            checked_pairs.insert((node_idx.min(other_node_idx), node_idx.max(other_node_idx)));
            
            // Check if either read is low quality
            if low_qual_qnames.contains(qname) {
                // Mark as negative weight in matrix (indicates low quality)
                result.weight_matrix.insert((node_idx.index(), other_node_idx.index()), -1.0);
                result.weight_matrix.insert((other_node_idx.index(), node_idx.index()), -1.0);
                debug!("The read pair {} is low quality, skip inspecting this pair of fragments", qname);
                continue;
            }
            
            if low_qual_qnames.contains(other_qname) {
                // Mark as negative weight in matrix (indicates low quality)
                result.weight_matrix.insert((node_idx.index(), other_node_idx.index()), -1.0);
                result.weight_matrix.insert((other_node_idx.index(), node_idx.index()), -1.0);
                debug!("The read pair {} is low quality, skip inspecting this pair of fragments", other_qname);
                continue;
            }
            
            // Inspect the overlap between the two pairs of reads
            let overlap_intervals = get_overlap_intervals(read_pair, other_read_pair)?;
            debug!("Found {} overlap intervals between {} and {}", 
                   overlap_intervals.len(), qname, other_qname);
            
            if overlap_intervals.is_empty() {
                continue;
            }
            
            // Process overlapping regions to determine haplotype compatibility
            // Using more descriptive variable names than Python's qname_bools
            let mut compatibility_results: Vec<i32> = Vec::with_capacity(4);
            let mut pair_weight: Option<f32> = None;
            
            // Track inspected overlaps to avoid redundant processing
            // Clear all the contents inside the inspected_overlaps
            let mut inspected_overlaps = FastIntervals::new(overlap_intervals.len());
            
            for overlap in &overlap_intervals {
                // Find uncovered regions within this overlap
                // This replaces Python's find_uncovered_regions_numba call
                let uncovered_regions = find_uncovered_regions_numba(
                    &inspected_overlaps,
                    overlap.start,
                    overlap.end
                );
                
                debug!("Found {} uncovered regions for overlap {:?}", 
                       uncovered_regions.len(), overlap);
                
                for uncovered_region in uncovered_regions {
                    let (uncovered_start, uncovered_end) = (uncovered_region[0], uncovered_region[1]);
                    
                    // Check if either read became low quality during processing
                    if low_qual_qnames.contains(other_qname) {
                        debug!("The read pair {} is low quality, skip this pair of reads", other_qname);
                        break;
                    }
                    if low_qual_qnames.contains(qname) {
                        debug!("The read pair {} is low quality, skip this pair of reads", qname);
                        break;
                    }
                    
                    // Determine if reads come from same haplotype in this region
                    // This calls the placeholder function that will be implemented later
                    let (haplotype_result, read_weight) = determine_same_haplotype(
                        read_pair,
                        other_read_pair,
                        uncovered_start,
                        uncovered_end,
                        &chrom,
                        allele_depth_map,
                        intrinsic_allele_depth_map,
                        config,
                        &mut result.read_hap_vectors,
                        &mut result.read_error_vectors,
                        &mut result.read_ref_pos_dict,
                        low_qual_qnames,
                    )?;
                    
                    // Process the haplotype comparison result
                    match haplotype_result {
                        HaplotypeResult::Same => {
                            compatibility_results.push(1);
                            if let Some(weight) = read_weight {
                                info!("Weight is {}. Found the reads {} and {} are in the same haplotype, overlap region ({}-{})", 
                                      weight, qname, other_qname, overlap.start, overlap.end);
                            }
                        }
                        HaplotypeResult::Different => {
                            compatibility_results.push(-1);
                            info!("Found the reads {} and {} are in different haplotypes, overlap region ({}-{})", 
                                  qname, other_qname, overlap.start, overlap.end);
                        }
                        HaplotypeResult::Unknown => {
                            compatibility_results.push(0);
                        }
                    }
                    
                    // Update pair weight
                    if let Some(weight) = read_weight {
                        let weight = weight.max(0.0); // Ensure non-negative
                        let norm_weight = weight / (config.mean_read_length * 10.0);
                        
                        match pair_weight {
                            None => pair_weight = Some(norm_weight),
                            Some(ref mut pw) => *pw += norm_weight,
                        }
                    }
                }
                
                // Mark this overlap as inspected
                inspected_overlaps.add(overlap.start, overlap.end);
            }
            
            // Check if reads became low quality during processing
            if low_qual_qnames.contains(other_qname) || low_qual_qnames.contains(qname) {
                result.weight_matrix.insert((node_idx.index(), other_node_idx.index()), -1.0);
                result.weight_matrix.insert((other_node_idx.index(), node_idx.index()), -1.0);
                continue;
            }
            
            // Determine final compatibility based on all results
            // Equivalent to Python's any_false_numba(qname_bools) check
            let final_weight = if any_false(&compatibility_results) {
                // Any incompatible region means overall incompatibility
                info!("Qname_bools are {:?}, Found two pairs {} and {} are in different haplotypes", 
                      compatibility_results, qname, other_qname);
                result.weight_matrix.insert((node_idx.index(), other_node_idx.index()), -1.0);
                result.weight_matrix.insert((other_node_idx.index(), node_idx.index()), -1.0);
                None
            } else {
                // All regions compatible or unknown
                let weight = pair_weight.unwrap_or(0.0).max(1e-4); // Minimum weight for compatible pairs
                debug!("Between {} and {}, the pair weight is {}", qname, other_qname, weight);
                
                result.weight_matrix.insert((node_idx.index(), other_node_idx.index()), weight);
                result.weight_matrix.insert((other_node_idx.index(), node_idx.index()), weight);
                
                // Add edge to graph
                let edge = result.graph.add_edge(node_idx, other_node_idx, weight);
                info!("Added edge {:?} between {} and {} with weight {}", edge, qname, other_qname, weight);
                
                Some(weight)
            };
        }
    }
    
    info!("Now we finished building up the edges in the graph. There are currently {} vertices and {} edges in the graph", 
          result.vertex_count(), result.edge_count());
    
    Ok(result)
}

/// Check if an edge between two nodes has already been inspected
/// 
/// Equivalent to Python's check_edge function but adapted for Rust's HashSet
/// Returns true if the edge was already checked, false otherwise
fn check_edge(u: NodeIndex, v: NodeIndex, checked_pairs: &HashSet<(NodeIndex, NodeIndex)>) -> bool {
    // Check both orderings since we store ordered pairs
    let ordered = if u < v { (u, v) } else { (v, u) };
    checked_pairs.contains(&ordered)
}

/// Extract chromosome name from a read pair
/// 
/// This function needs access to the BAM header to convert tid to chromosome name
/// In Python this is simple: paired_reads[0].reference_name
/// In Rust, we need to handle the conversion properly
fn get_read_chromosome(
    read_pair: &crate::structs::ReadPair,
    read_pair_map: &ReadPairMap,
) -> Result<String, Box<dyn std::error::Error>> {
    // For now, we'll extract the chromosome from the interval trees
    // This assumes the read was already properly indexed in the interval trees
    
    // Find which chromosome this read pair belongs to by checking interval trees
    for (chrom, _) in &read_pair_map.interval_trees {
        if read_pair_map.find_overlapping_qname_indices(
            chrom, 
            read_pair.read1.reference_start(), 
            read_pair.read1.reference_end()
        )?.contains(&read_pair.qname_idx) {
            return Ok(chrom.clone());
        }
    }
    
    Err("Read pair not found in any chromosome interval tree".into())
}

/// Calculate the genomic span covered by a read pair
/// 
/// Returns the minimum start position and maximum end position across all reads in the pair
/// Equivalent to Python's:
///   start = min(r.reference_start for r in paired_reads)
///   end = max(r.reference_end for r in paired_reads)
fn get_read_pair_span(read_pair: &crate::structs::ReadPair) -> (i64, i64) {
    let read1_start = read_pair.read1.reference_start();
    let read1_end = read_pair.read1.reference_end();
    
    let (read2_start, read2_end) = if let Some(ref read2) = read_pair.read2 {
        (read2.reference_start(), read2.reference_end())
    } else {
        (read1_start, read1_end) // If no read2, use read1 coordinates
    };
    
    // Return the span covering both reads
    (read1_start.min(read2_start), read1_end.max(read2_end))
}

/// Extract the overlapping intervals between two pairs of reads
/// 
/// This function finds all genomic regions where reads from two different pairs overlap
/// Returns a vector of OverlapInterval structs containing the overlap coordinates
/// 
/// Equivalent to Python's get_overlap_intervals function
fn get_overlap_intervals(
    read_pair1: &crate::structs::ReadPair,
    read_pair2: &crate::structs::ReadPair,
) -> Result<Vec<OverlapInterval>, Box<dyn std::error::Error>> {
    let mut overlaps = Vec::new();
    
    // Collect all reads from both pairs
    let reads1: Vec<&Record> = if let Some(ref read2) = read_pair1.read2 {
        vec![&read_pair1.read1, read2]
    } else {
        vec![&read_pair1.read1]
    };
    
    let reads2: Vec<&Record> = if let Some(ref read2) = read_pair2.read2 {
        vec![&read_pair2.read1, read2]
    } else {
        vec![&read_pair2.read1]
    };
    
    // Check all possible combinations of reads between the two pairs
    for (i, read1) in reads1.iter().enumerate() {
        for (j, read2) in reads2.iter().enumerate() {
            let start1 = read1.reference_start();
            let end1 = read1.reference_end();
            let start2 = read2.reference_start();
            let end2 = read2.reference_end();
            
            // Calculate overlap
            let overlap_start = start1.max(start2);
            let overlap_end = end1.min(end2);
            
            if overlap_start < overlap_end {
                // Create read identifiers for debugging
                let read1_id = format!("{}_{}", 
                    std::str::from_utf8(read1.qname()).unwrap_or("unknown"), 
                    if i == 0 { "R1" } else { "R2" }
                );
                let read2_id = format!("{}_{}", 
                    std::str::from_utf8(read2.qname()).unwrap_or("unknown"),
                    if j == 0 { "R1" } else { "R2" }
                );
                
                overlaps.push(OverlapInterval {
                    start: overlap_start,
                    end: overlap_end,
                    read1_ref: read1_id,
                    read2_ref: read2_id,
                });
                
                debug!("Found overlap interval {}-{} between {} and {}", 
                       overlap_start, overlap_end, 
                       overlaps.last().unwrap().read1_ref,
                       overlaps.last().unwrap().read2_ref);
            }
        }
    }
    
    Ok(overlaps)
}

/// Find uncovered regions within an overlap interval using binary search
/// 
/// This is a Rust implementation of Python's find_uncovered_regions_numba function
/// Uses efficient algorithms to find regions not yet inspected within a query interval
/// 
/// # Arguments
/// * `inspected_overlaps` - FastIntervals tracking already processed regions
/// * `query_start` - Start of the query interval
/// * `query_end` - End of the query interval
/// 
/// # Returns
/// Vec of [start, end] arrays representing uncovered regions
fn find_uncovered_regions_numba(
    inspected_overlaps: &FastIntervals,
    query_start: i64,
    query_end: i64,
) -> Vec<[i64; 2]> {
    let (existing_starts, existing_ends) = inspected_overlaps.get_intervals();
    
    if existing_starts.is_empty() {
        return vec![[query_start, query_end]];
    }
    
    // Find overlapping intervals using binary search logic
    let overlaps: Vec<bool> = existing_starts.iter().zip(existing_ends.iter())
        .map(|(&start, &end)| start <= query_end && end >= query_start)
        .collect();
    
    if !overlaps.iter().any(|&x| x) {
        return vec![[query_start, query_end]];
    }
    
    // Pre-allocate result vector
    let mut result = Vec::new();
    let mut current_start = query_start;
    
    for (i, &overlaps_i) in overlaps.iter().enumerate() {
        if overlaps_i {
            if existing_starts[i] > current_start {
                result.push([current_start, existing_starts[i]]);
            }
            current_start = current_start.max(existing_ends[i]);
        }
    }
    
    if current_start < query_end {
        result.push([current_start, query_end]);
    }
    
    result
}

/// Check if any element in the array is false (-1 in this case)
/// 
/// Equivalent to Python's any_false_numba function
/// Returns true if any element equals -1, false otherwise
fn any_false(arr: &[i32]) -> bool {
    arr.iter().any(|&x| x == -1)
}
