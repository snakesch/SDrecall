use std::collections::HashSet;
use log::{info, debug, warn};

use petgraph::graph::NodeIndex;
use half::f16;
use crate::structs::{
    ReadPairMap, AlleleDepthMap, PhasingGraphResult, HaplotypeConfig, 
    FastIntervals, OverlapInterval
};
use crate::haplotype_determination::{determine_same_haplotype, HaplotypeResult};
use rust_htslib::bam::Record;
use rust_htslib::bam::ext::BamRecordExtensions;

/*
================================================================================
                            ALGORITHM WORKFLOW DIAGRAM
================================================================================

┌─────────────────────────────────────────────────────────────────────────────┐
│                         BUILD_PHASING_GRAPH ALGORITHM                      │
└─────────────────────────────────────────────────────────────────────────────┘

Step 1: Initialize and Validate Input
┌─────────────────────────────────────┐
│ • Check total_read_pairs count      │
│ • Validate allele_depth_map         │
│ • Create PhasingGraphResult         │
│ • Log initial statistics            │
└─────────────────────────────────────┘

Step 2: Pre-allocate Graph Structure
┌─────────────────────────────────────┐
│ For each qname_idx (0..num_nodes):  │
│                                     │
│ ┌─ Node Creation ──────────────┐   │
│ │                              │   │
│ │ • Add node to petgraph       │   │
│ │ • Store qname_idx as weight  │   │
│ │ • Update qname_to_node map   │   │
│ └──────────────────────────────┘   │
└─────────────────────────────────────┘

Step 3: Initialize Weight Matrix and Tracking
┌─────────────────────────────────────────────────────────────┐
│                                                             │
│ • Initialize num_nodes × num_nodes weight matrix           │
│ • Create checked_pairs HashSet for edge deduplication     │
│ • Set up interval tracking structures                     │
│                                                             │
└─────────────────────────────────────────────────────────────┘

Step 4: Main Processing Loop - Read Pair Iteration
┌─────────────────────────────────────────────────────────────┐
│ For each (qname_idx, read_pair) in readpair_dict:          │
│                                                             │
│ ┌─ Spatial Query ──────────────────────────────────────┐   │
│ │ • Extract chromosome and genomic span               │   │
│ │ • Query interval tree for overlapping reads        │   │
│ │ • Get list of overlapping_qname_indices             │   │
│ └─────────────────────────────────────────────────────┘   │
│                                                             │
│ ┌─ Overlap Processing ─────────────────────────────────┐   │
│ │ For each other_qname_idx in overlapping reads:      │   │
│ │ • Skip self-comparison                              │   │
│ │ • Check edge deduplication (checked_pairs)         │   │
│ │ • Proceed to pairwise analysis                      │   │
│ └─────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘

Step 5: Pairwise Read Compatibility Analysis
┌─────────────────────────────────────────────────────────────┐
│                                                             │
│ Input: Two read pairs (current vs overlapping)             │
│ Note: All reads pre-filtered for basic quality in BAM stage│
│                                                             │
│ ┌─ Overlap Interval Detection ─────────────────────────┐   │
│ │ • Find all read-vs-read overlaps (R1-R1, R1-R2...)  │   │
│ │ • Create OverlapInterval structures                 │   │
│ │ • Skip if no overlaps found                         │   │
│ └─────────────────────────────────────────────────────┘   │
│                                                             │
│ ┌─ Region-by-Region Haplotype Analysis ────────────────┐   │
│ │ For each overlap_interval:                          │   │
│ │ • Find uncovered regions (avoid double-counting)   │   │
│ │ • Dynamic quality check during analysis            │   │
│ │ • Call determine_same_haplotype() for each region  │   │
│ │ • Collect compatibility results [1, -1, 0]         │   │
│ │ • Accumulate pair weights                          │   │
│ └─────────────────────────────────────────────────────┘   │
│                                                             │
│ Output: compatibility_results[], pair_weight               │
└─────────────────────────────────────────────────────────────┘

Step 6: Final Edge Weight Determination
┌─────────────────────────────────────┐
│ Analyze compatibility_results:      │
│                                     │
│ ┌─ Any False (-1)? ──────────┐   │
│ │                              │   │
│ │ YES: Set weight = -1.0       │   │
│ │      (Different haplotype)   │   │
│ │                              │   │
│ │ NO:  Calculate final weight  │   │
│ │      Add edge to graph       │   │
│ │      Store in weight matrix  │   │
│ └──────────────────────────────┘   │
└─────────────────────────────────────┘

Step 7: Graph Finalization
┌─────────────────────────────────────┐
│ • Log final graph statistics        │
│ • Return completed PhasingGraphResult│
│ • Graph ready for clique finding    │
└─────────────────────────────────────┘

Data Flow Timeline:
BAM reads ──▶ [Spatial Index] ──▶ [Overlap Detection] ──▶ [Haplotype Analysis] ──▶ Phasing Graph
    │              │                    │                        │                    │
    ▼              ▼                    ▼                        ▼                    ▼
[ReadPairs]   [Interval Trees]    [Overlap Regions]        [Compatibility]      [Weighted Edges]



Complexity: O(n × k × r) where n=read_pairs, k=average_overlaps_per_read, r=regions_per_overlap
Graph Structure: Vertices=ReadPairs, Edges=HaplotypeCompatibility, Weights=OverlapConfidence
*/

// Rust idioms used in this code:
// • ? operator: Error propagation - automatically returns errors up the call stack
// • ok_or(): Converts Option<T> to Result<T, E> - None becomes Err with provided message  
// • * operator: Dereferencing - accesses the value a reference points to
// • as_ref(): Converts &Option<T> to Option<&T> for comparison without taking ownership
// • unwrap_or(): Returns the contained value or a provided default if Option is None

/*
================================================================================
                           PERFORMANCE OPTIMIZATIONS
================================================================================

KEY OPTIMIZATION: Direct qname_idx ↔ NodeIndex Correspondence
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Instead of maintaining expensive HashMap mappings, we leverage petgraph's 
sequential NodeIndex assignment guarantee:

  qname_idx = 0 → NodeIndex(0) → Weight Matrix[0][0]
  qname_idx = 1 → NodeIndex(1) → Weight Matrix[1][1]  
  qname_idx = n → NodeIndex(n) → Weight Matrix[n][n]

Benefits:
• ELIMINATED: node_indices: Vec<NodeIndex>        (saves 8*n bytes)
• ELIMINATED: qname_to_node: HashMap<String, NodeIndex>  (saves ~64*n bytes)
• O(1) NodeIndex lookup instead of O(1) HashMap lookup (but no hashing overhead)
• Direct matrix indexing: matrix[qname_idx][other_qname_idx]
• Simplified memory layout and better cache locality

Memory Savings: ~72 bytes per read pair + reduced allocation overhead
*/

/// Create NodeIndex directly from qname_idx using petgraph's sequential assignment
/// 
/// # Safety Note
/// This function relies on the fact that we add nodes sequentially (0..n) 
/// during graph initialization, so NodeIndex(k) corresponds to the kth node added.
/// 
/// # Arguments
/// * `qname_idx` - The continuous qname index (0..num_reads-1)
/// 
/// # Returns
/// NodeIndex that can be used directly with the graph
#[inline]
fn qname_idx_to_node_index(qname_idx: usize) -> NodeIndex {
    NodeIndex::new(qname_idx)
}

pub fn build_phasing_graph(
    read_pair_map: &ReadPairMap,
    allele_depth_map: &AlleleDepthMap,
    header: &rust_htslib::bam::HeaderView,
    config: &HaplotypeConfig) -> Result<PhasingGraphResult, Box<dyn std::error::Error>> {
    
    // ========== STEP 1: Initialize and Validate Input ==========
    let total_read_pairs = read_pair_map.readpair_dict.len();
    info!("[build_phasing_graph] There are totally {} pair of reads, mean read length is {:.1}. with adequate mapping or base quality which can be used to build the graph", 
          total_read_pairs, config.mean_read_length);
    
    if allele_depth_map.is_empty() {
        warn!("[build_phasing_graph] No ALT allele found in this BAM file. Skip this entire script");
        return Ok(PhasingGraphResult::new());
    }
    
    let mut result = PhasingGraphResult::new();
    
    // ========== STEP 2: Pre-allocate Graph Structure ==========
    // With direct correspondence: qname_idx = NodeIndex.index()
    // We only need to ensure we have num_nodes sequential nodes (0..num_nodes-1)
    let num_nodes = read_pair_map.readpair_dict.len();
    
    // Pre-allocate graph nodes ensuring sequential assignment
    // petgraph should assign: first node = NodeIndex(0), second = NodeIndex(1), etc.
    for qname_idx in 0..num_nodes {
        // Add node with unit data - we rely on sequential assignment for correspondence
        let node = result.graph.add_node(()); // Empty unit data - no storage needed
        
        // Verify petgraph's sequential assignment - return error instead of panic
        if qname_idx != node.index() {
            return Err(format!(
                "Graph node assignment failed: expected NodeIndex({}), got NodeIndex({}). \
                This breaks the required 1-to-1 correspondence between qname_idx and NodeIndex.",
                qname_idx, node.index()
            ).into());
        }
        
        debug!("[build_phasing_graph] Verified node correspondence: qname_idx {} = NodeIndex({})", qname_idx, node.index());
    }
    
    // ========== STEP 3: Initialize Weight Matrix and Tracking ==========
    // Initialize weight matrix with known size
    result.initialize_weight_matrix(num_nodes);
    info!("[build_phasing_graph] Initialized {}x{} weight matrix", num_nodes, num_nodes);
    
    // Create a set to track which node pairs we've already checked
    let mut checked_pairs: HashSet<(usize, usize)> = HashSet::new();
    
    // ========== STEP 4: Main Processing Loop - Read Pair Iteration ==========
    // Use read dict to iterate through all the reads to build a graph
    for (&qname_idx, read_pair) in &read_pair_map.readpair_dict {
        let qname = &read_pair.qname;
        
        // Direct NodeIndex construction from qname_idx (leveraging sequential guarantee)
        let node_idx = qname_idx_to_node_index(qname_idx);
        
        // ---------- Step 4a: Spatial Query ----------
        // Get chromosome for overlap queries
        let chrom = get_read_chromosome(read_pair, header)?;
        
        // Query overlapping reads for actual coverage regions (not gaps)
        let overlapping_qname_indices = get_overlapping_qname_indices_for_pair(
            read_pair_map, 
            read_pair, 
            &chrom
        )?;
        
        // ---------- Step 4b: Overlap Processing ----------
        // Iterate through the overlapping reads
        for other_qname_idx in overlapping_qname_indices {
            if qname_idx == other_qname_idx {
                continue;
            }
            
            // Check if we've already inspected this edge (use qname_idx directly)
            let edge_key = (qname_idx.min(other_qname_idx), qname_idx.max(other_qname_idx));
            if checked_pairs.contains(&edge_key) {
                continue;
            }
            checked_pairs.insert(edge_key);
            
            // Get the other read pair
            let other_read_pair = match read_pair_map.readpair_dict.get(&other_qname_idx) {
                Some(pair) => pair, // HashMap get() returns Option<&ReadPair> (not owned)
                None => {
                    warn!("[build_phasing_graph] Could not find read pair for qname_idx: {}", other_qname_idx);
                    continue;
                }
            };
            
            let other_qname = &other_read_pair.qname;
            
            // Direct NodeIndex construction for the other read
            let other_node_idx = qname_idx_to_node_index(other_qname_idx);
            
            // ========== STEP 5: Pairwise Read Compatibility Analysis ==========
            
            // Note: Initial quality filtering already done in bam_reading.rs
            // All reads in readpair_dict are pre-filtered for basic quality issues
            
            // ---------- Step 5a: Overlap Interval Detection ----------
            // Inspect the overlap between the two pairs of reads
            let overlap_intervals = get_overlap_intervals(read_pair, other_read_pair)?;
            debug!("[build_phasing_graph] Found {} overlap intervals between {} and {}", 
                   overlap_intervals.len(), qname, other_qname);
            
            if overlap_intervals.is_empty() {
                continue;
            }
            
            // ---------- Step 5b: Region-by-Region Haplotype Analysis ----------
            // Process overlapping regions to determine haplotype compatibility
            // Using more descriptive variable names than Python's qname_bools
            let mut compatibility_results: Vec<i32> = Vec::with_capacity(4);
            let mut pair_weight: Option<f32> = None;
            
            // Track inspected overlaps to avoid redundant processing
            // Clear all the contents inside the inspected_overlaps
            let mut inspected_overlaps = FastIntervals::new(overlap_intervals.len());
            
            for overlap in &overlap_intervals {
                // Find uncovered regions within this overlap
                // This replaces Python's find_uncovered_regions call
                let uncovered_regions = find_uncovered_regions(
                    &inspected_overlaps,
                    overlap.start,
                    overlap.end
                );
                
                debug!("[build_phasing_graph] Found {} uncovered regions for overlap {:?}", 
                       uncovered_regions.len(), overlap);
                
                for uncovered_region in uncovered_regions {
                    let (uncovered_start, uncovered_end) = (uncovered_region[0], uncovered_region[1]);
                    
                    // Now we have direct access to the specific reads from the overlap
                    let (haplotype_result, read_weight) = determine_same_haplotype(
                        overlap.read1,
                        overlap.read2,
                        uncovered_start,
                        uncovered_end,
                        &chrom,
                        allele_depth_map,
                        config,
                        &mut result.read_hap_vectors,
                        &mut result.read_error_vectors,
                        &mut result.read_ref_pos_dict,
                    )?;
                    
                    // Process the haplotype comparison result
                    match haplotype_result {
                        HaplotypeResult::Same => {
                            compatibility_results.push(1);
                            if let Some(weight) = read_weight {
                                debug!("[build_phasing_graph] Weight is {}. Found the reads {} and {} are in the same haplotype, overlap region ({}:{}-{})", 
                                      weight, qname, other_qname, chrom, overlap.start, overlap.end);
                            }
                        }
                        HaplotypeResult::Different => {
                            compatibility_results.push(-1);
                            debug!("[build_phasing_graph] Found the reads {} and {} are in different haplotypes, overlap region ({}:{}-{})", 
                                  qname, other_qname, chrom, overlap.start, overlap.end);
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
            
            // ========== STEP 6: Final Edge Weight Determination ==========
            // Determine final compatibility based on all results
            // Equivalent to Python's any_false_numba(qname_bools) check
            let _final_weight = if any_false(&compatibility_results) {
                // Any incompatible region means overall incompatibility
                debug!("[build_phasing_graph] Qname_bools are {:?}, Found two pairs {} and {} are in different haplotypes\n", 
                      compatibility_results, qname, other_qname);
                result.set_weight(qname_idx, other_qname_idx, -1.0)?;
                None
            } else {
                // All regions compatible or unknown
                let weight = pair_weight.unwrap_or(0.0).max(1e-4); // Minimum weight for compatible pairs
                debug!("[build_phasing_graph] Between {} and {}, the pair weight is {}\n", qname, other_qname, weight);
                
                // Set weight in matrix using qname_idx directly
                result.set_weight(qname_idx, other_qname_idx, weight)?;
                
                // Add edge to graph with f16 weight
                let weight_f16 = f16::from_f32(weight);
                result.graph.add_edge(node_idx, other_node_idx, weight_f16);
                
                Some(weight)
            };
        }
    }
    
    // ========== STEP 7: Graph Finalization ==========
    // Extract lowqual_qnames from ReadPairMap to include in result
    result.lowqual_qnames = read_pair_map.noisy_qnames.keys().cloned().collect();
    
    info!("[build_phasing_graph] Finished building graph with {} vertices and {} edges", 
          result.vertex_count(), result.edge_count());
    info!("[build_phasing_graph] Identified {} low-quality qnames during processing\n", result.lowqual_qnames.len());
    
    Ok(result)
}



/// Extract chromosome name from a read pair using BAM header
/// 
/// This function directly accesses the chromosome name from the BAM record's tid
/// using rust-htslib's header.tid2name() function for O(1) lookup
/// Equivalent to Python's paired_reads[0].reference_name
fn get_read_chromosome(
    read_pair: &crate::structs::ReadPair,
    header: &rust_htslib::bam::HeaderView,
) -> Result<String, Box<dyn std::error::Error>> {
    let tid = read_pair.read1.tid();
    
    // Handle unmapped reads
    if tid < 0 {
        return Err("Read is unmapped - no reference chromosome".into());
    }
    
    // Convert tid to chromosome name using BAM header
    let chrom_bytes = header.tid2name(tid as u32);
    let chrom = std::str::from_utf8(chrom_bytes)
        .map_err(|e| format!("Invalid chromosome name encoding: {}", e))?;
    
    Ok(chrom.to_string())
}

/// Query overlapping qname indices for a read pair using actual coverage regions
/// 
/// This function performs separate interval tree queries for each read in the pair,
/// avoiding the inclusion of gaps between reads that don't have actual sequence coverage.
/// This is more accurate than the single-span approach used in the Python implementation.
/// 
/// Benefits over single-span approach:
/// - More accurate: Only queries actual read coverage regions
/// - Better performance: Avoids querying irrelevant gap regions  
/// - Biological correctness: Matches actual sequenced bases, not gaps
/// - Scalability: Important for long insert-size paired reads with large gaps
/// 
/// # Arguments
/// * `read_pair_map` - The read pair map containing interval trees
/// * `read_pair` - The read pair to query overlaps for
/// * `chrom` - The chromosome name for the queries
/// 
/// # Returns
/// A vector of qname_indices that overlap with any read in the pair
fn get_overlapping_qname_indices_for_pair(
    read_pair_map: &ReadPairMap,
    read_pair: &crate::structs::ReadPair,
    chrom: &str,
) -> Result<Vec<usize>, Box<dyn std::error::Error>> {
    use std::collections::HashSet;
    
    let mut all_overlapping_qnames = HashSet::new();
    
    // Query overlaps for read1
    let read1_overlaps = read_pair_map.find_overlapping_qname_indices(
        chrom,
        read_pair.read1.reference_start(),
        read_pair.read1.reference_end()
    )?;
    all_overlapping_qnames.extend(read1_overlaps);
    
    // Query overlaps for read2 if it exists
    if let Some(ref read2) = read_pair.read2 {
        let read2_overlaps = read_pair_map.find_overlapping_qname_indices(
            chrom,
            read2.reference_start(),
            read2.reference_end()
        )?;
        all_overlapping_qnames.extend(read2_overlaps);
    }
    
    // Convert HashSet back to Vec for compatibility with existing code
    Ok(all_overlapping_qnames.into_iter().collect())
}

/// Extract the overlapping intervals between two pairs of reads
/// 
/// This function finds all genomic regions where reads from two different pairs overlap
/// Returns a vector of OverlapInterval structs containing the overlap coordinates
/// and direct references to the overlapping reads
/// 
/// Equivalent to Python's get_overlap_intervals function
fn get_overlap_intervals<'a>(
    read_pair1: &'a crate::structs::ReadPair,
    read_pair2: &'a crate::structs::ReadPair,
) -> Result<Vec<OverlapInterval<'a>>, Box<dyn std::error::Error>> {
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
    for read1 in &reads1 {
        for read2 in &reads2 {
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
                    read1,  // Direct reference to the record
                    read2,  // Direct reference to the record
                });
                
                debug!("[get_overlap_intervals] Found overlap interval {}-{} between reads", 
                       overlap_start, overlap_end);
            }
        }
    }
    
    Ok(overlaps)
}

/// Find uncovered regions within an overlap interval using binary search
/// 
/// This is a Rust implementation of Python's find_uncovered_regions function
/// Uses efficient algorithms to find regions not yet inspected within a query interval
/// 
/// # Arguments
/// * `inspected_overlaps` - FastIntervals tracking already processed regions
/// * `query_start` - Start of the query interval
/// * `query_end` - End of the query interval
/// 
/// # Returns
/// Vec of [start, end] arrays representing uncovered regions
fn find_uncovered_regions(
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
