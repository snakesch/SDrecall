use std::collections::HashSet;
use rust_htslib::bam::Record;
use rustc_hash::FxHashMap;
use ahash::AHashMap;
use petgraph::Graph;
use petgraph::Undirected;
use ndarray::Array2;
use half::f16; // Add f16 import

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Interval {
    pub start: i64,
    pub end: i64,
    pub qname_idx: usize,
}

impl Interval {
    pub fn new(start: i64, end: i64, qname_idx: usize) -> Self {
        Self { start, end, qname_idx }
    }
}

pub struct SortedVecIntervals {
    // Use FxHashMap for integer keys (qname_idx)
    interval_map: FxHashMap<usize, [Option<Interval>; 2]>,
    // For querying phase: pre-allocated vector for fast overlap queries
    sorted_intervals: Option<Vec<Interval>>,
    is_finalized: bool,
}

impl SortedVecIntervals {
    pub fn new() -> Self {
        Self {
            interval_map: FxHashMap::default(),
            sorted_intervals: None,
            is_finalized: false,
        }
    }

    /// Add an interval for a qname_idx (must have exactly 2 intervals for paired reads)
    pub fn add_interval(&mut self, start: i64, end: i64, qname_idx: usize) -> Result<(), String> {
        if self.is_finalized {
            return Err("Cannot add intervals after finalization".to_string());
        }
        
        let interval = Interval::new(start, end, qname_idx);
        let intervals = self.interval_map.entry(qname_idx).or_insert([None, None]);
        
        // Find first empty slot
        if intervals[0].is_none() {
            intervals[0] = Some(interval);
            Ok(())
        } else if intervals[1].is_none() {
            intervals[1] = Some(interval);
            Ok(())
        } else {
            Err(format!("Cannot add more than 2 intervals for qname_idx {}", qname_idx))
        }
    }

    /// Remove all intervals for a specific qname_idx (O(1) operation)
    pub fn remove_qname_idx(&mut self, qname_idx: usize) -> Result<(), String> {
        if self.is_finalized {
            return Err("Cannot remove intervals after finalization".to_string());
        }
        
        self.interval_map.remove(&qname_idx);
        Ok(())
    }

    /// Remove incomplete pairs (qname_idx with != 2 intervals)
    pub fn remove_incomplete_pairs(&mut self) -> Result<(), String> {
        if self.is_finalized {
            return Err("Cannot remove intervals after finalization".to_string());
        }
        
        self.interval_map.retain(|_, intervals| {
            intervals[0].is_some() && intervals[1].is_some()
        });
        Ok(())
    }

    /// Finalize the structure for efficient querying
    pub fn finalize(&mut self) -> Result<(), String> {
        if self.is_finalized {
            return Ok(());
        }

        // First remove incomplete pairs
        self.remove_incomplete_pairs()?;

        // Pre-allocate Vec with known capacity (2 intervals per qname_idx)
        let total_intervals = self.interval_map.len() * 2;
        let mut all_intervals = Vec::with_capacity(total_intervals);
        
        for intervals_array in self.interval_map.values() {
            for interval_opt in intervals_array {
                if let Some(interval) = interval_opt {
                    all_intervals.push(*interval);
                }
            }
        }
        
        // Sort by start position, then by end position
        all_intervals.sort_by_key(|interval| (interval.start, interval.end));
        
        self.sorted_intervals = Some(all_intervals);
        self.is_finalized = true;
        
        // Clear the HashMap to save memory
        self.interval_map.clear();
        Ok(())
    }

    /// Find overlapping qname_indices (deduplicated) - only works after finalization
    pub fn find_overlaps(&self, query_start: i64, query_end: i64) -> Result<Vec<usize>, String> {
        if !self.is_finalized {
            return Err("Must call finalize() before querying".to_string());
        }

        let intervals = self.sorted_intervals.as_ref()
            .ok_or("Sorted intervals not initialized")?;
        
        // Use HashSet to ensure uniqueness, then convert to Vec
        let mut unique_qname_indices = std::collections::HashSet::new();
        
        // Binary search for first interval that could overlap
        let start_pos = intervals.iter().position(|interval| interval.end > query_start)
            .unwrap_or(intervals.len());
        
        // Collect overlapping intervals and ensure uniqueness
        for interval in &intervals[start_pos..] {
            if interval.start >= query_end {
                break; // No more overlaps possible
            }
            
            if interval.end > query_start && interval.start < query_end {
                unique_qname_indices.insert(interval.qname_idx);
            }
        }
        
        // Convert HashSet to Vec - guaranteed unique qname_indices
        Ok(unique_qname_indices.into_iter().collect())
    }

    pub fn len(&self) -> usize {
        if self.is_finalized {
            self.sorted_intervals.as_ref().map_or(0, |v| v.len())
        } else {
            self.interval_map.values().map(|v| {
                v.iter().filter(|interval| interval.is_some()).count()
            }).sum()
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

#[derive(Clone, Debug)]
pub struct ReadPair {
    pub read1: Record,
    pub read2: Option<Record>,  // Optional until complete
    pub qname: String,
    pub qname_idx: usize,
}

impl ReadPair {
    pub fn new_incomplete(read1: Record, qname: String, qname_idx: usize) -> Self {
        Self {
            read1,
            read2: None,
            qname,
            qname_idx,
        }
    }
    
    pub fn new_complete(read1: Record, read2: Record, qname: String, qname_idx: usize) -> Self {
        Self {
            read1,
            read2: Some(read2),
            qname,
            qname_idx,
        }
    }
    
    pub fn complete_with_read2(&mut self, read2: Record) {
        self.read2 = Some(read2);
    }
    
    pub fn is_complete(&self) -> bool {
        self.read2.is_some()
    }
    
    pub fn update_read1(&mut self, new_record: Record) {
        self.read1 = new_record;
    }
    
    pub fn update_read2(&mut self, new_record: Record) {
        self.read2 = Some(new_record);
    }
}

pub struct ReadPairMap {
    // Use AHashMap for string keys (chromosome names)
    pub interval_trees: AHashMap<String, SortedVecIntervals>,
    // Use FxHashMap for integer keys (qname_idx)
    pub readpair_dict: FxHashMap<usize, ReadPair>,
    // Use AHashMap for string keys (qname)
    pub qname_to_idx: AHashMap<String, usize>,
    // Use FxHashMap for integer keys (qname_idx)
    pub idx_to_qname: FxHashMap<usize, String>,
    // Use AHashMap for string keys (qname)
    pub noisy_qnames: AHashMap<String, ()>, // Use HashMap as HashSet with () values
}

impl ReadPairMap {
    pub fn new() -> Self {
        Self {
            interval_trees: AHashMap::new(),
            readpair_dict: FxHashMap::default(),
            qname_to_idx: AHashMap::new(),
            idx_to_qname: FxHashMap::default(),
            noisy_qnames: AHashMap::new(),
        }
    }

    /// Find overlapping qname indices using efficient interval search
    pub fn find_overlapping_qname_indices(&self, chrom: &str, start: i64, end: i64) -> Result<Vec<usize>, String> {
        match self.interval_trees.get(chrom) {
            Some(tree) => tree.find_overlaps(start, end),
            None => Ok(Vec::new()),
        }
    }

    /// Get ReadPair by qname index - O(1) HashMap lookup
    pub fn get_readpair(&self, qname_idx: usize) -> Option<&ReadPair> {
        self.readpair_dict.get(&qname_idx)
    }

    /// Iterator for overlapping ReadPairs (combines interval search + HashMap lookup)
    pub fn overlapping_readpairs(&self, chrom: &str, start: i64, end: i64) -> impl Iterator<Item = &ReadPair> {
        self.find_overlapping_qname_indices(chrom, start, end)
            .unwrap_or_default()
            .into_iter()
            .filter_map(move |idx| self.readpair_dict.get(&idx))
    }
}

/// Simple array for allele depths: [A, T, C, G, N, total_depth]
/// Much simpler than a custom struct for just storing 6 u16 values
pub type PositionAlleleDepth = [u16; 6];

/// Efficient nested structure for chromosome -> position -> allele depths
/// Uses specialized HashMaps for better performance
pub struct AlleleDepthMap {
    /// Outer map: chromosome name -> inner map
    /// Using AHashMap for string keys (chromosome names)
    /// Using FxHashMap with u32 for genomic positions (max ~250M for human genome)
    chromosomes: AHashMap<String, FxHashMap<u32, PositionAlleleDepth>>,
}

impl AlleleDepthMap {
    pub fn new() -> Self {
        Self {
            chromosomes: AHashMap::new(),
        }
    }
    
    pub fn insert(&mut self, chrom: &str, pos: u32, data: PositionAlleleDepth) {
        self.chromosomes
            .entry(chrom.to_string())
            .or_insert_with(FxHashMap::default)
            .insert(pos, data);
    }
    
    pub fn get(&self, chrom: &str, pos: u32) -> Option<&PositionAlleleDepth> {
        self.chromosomes.get(chrom)?.get(&pos)
    }
    
    pub fn get_mut(&mut self, chrom: &str, pos: u32) -> Option<&mut PositionAlleleDepth> {
        self.chromosomes.get_mut(chrom)?.get_mut(&pos)
    }
    
    /// Iterator over all positions in a chromosome
    pub fn chromosome_positions(&self, chrom: &str) -> Option<impl Iterator<Item = (&u32, &PositionAlleleDepth)>> {
        self.chromosomes.get(chrom).map(|pos_map| pos_map.iter())
    }
}

/// Helper functions for working with PositionAlleleDepth arrays
impl AlleleDepthMap {
    /// Create a new allele depth array with total depth
    pub fn new_position_data(total_depth: u16) -> PositionAlleleDepth {
        [0, 0, 0, 0, 0, total_depth] // [A, T, C, G, N, total]
    }
    
    /// Set allele depth at given index
    pub fn set_allele_depth(data: &mut PositionAlleleDepth, allele_index: usize, depth: u16) {
        if allele_index < 5 {
            data[allele_index] = depth;
        }
    }
    
    /// Get allele depth at given index
    pub fn get_allele_depth(data: &PositionAlleleDepth, allele_index: usize) -> u16 {
        if allele_index < 5 {
            data[allele_index]
        } else {
            0
        }
    }
    
    /// Get total depth
    pub fn total_depth(data: &PositionAlleleDepth) -> u16 {
        data[5]
    }
    
    /// Check if the map is empty (no chromosomes or no data)
    pub fn is_empty(&self) -> bool {
        self.chromosomes.is_empty() || 
        self.chromosomes.values().all(|pos_map| pos_map.is_empty())
    }
    
    /// Get the number of chromosomes in the map
    pub fn chromosome_count(&self) -> usize {
        self.chromosomes.len()
    }
}

/// Overlap interval between two reads
/// Represents a genomic region where two reads overlap
#[derive(Debug, Clone)]
pub struct OverlapInterval<'a> {
    pub start: i64,
    pub end: i64,
    pub read1: &'a Record,  // Direct reference to first read
    pub read2: &'a Record,  // Direct reference to second read
}

/// Read haplotype vector - stores variant information for a read
/// This is equivalent to the Python read_hap_vectors dictionary
pub type ReadHaplotypeVector = Vec<i16>;  // Could be variant calls, error flags, etc.

/// Read error vector - stores error information for a read  
pub type ReadErrorVector = Vec<f32>;     // Error probabilities or quality scores

/// Fast interval storage for tracking inspected overlaps
/// Equivalent to Python's FastIntervals class
#[derive(Debug)]
pub struct FastIntervals {
    starts: Vec<i64>,
    ends: Vec<i64>,
    max_size: usize,
    current_size: usize,
}

impl FastIntervals {
    pub fn new(max_size: usize) -> Self {
        Self {
            starts: Vec::with_capacity(max_size),
            ends: Vec::with_capacity(max_size),
            max_size,
            current_size: 0,
        }
    }
    
    pub fn add(&mut self, start: i64, end: i64) {
        if self.current_size < self.max_size {
            self.starts.push(start);
            self.ends.push(end);
            self.current_size += 1;
        }
    }
    
    pub fn clear(&mut self) {
        self.starts.clear();
        self.ends.clear();
        self.current_size = 0;
    }
    
    pub fn get_intervals(&self) -> (&[i64], &[i64]) {
        (&self.starts[..self.current_size], &self.ends[..self.current_size])
    }
}

/// Phasing graph structure using petgraph with unit nodes and f16 weights
/// Node data is () since we use direct correspondence: qname_idx = NodeIndex.index()
/// Edge weight is f16 for memory efficiency
/// This ensures direct mapping: qname_idx = NodeIndex = matrix index
pub type PhasingGraph = Graph<(), f16, Undirected>;

/// Complete phasing graph result
pub struct PhasingGraphResult {
    /// The main graph structure with unit nodes (direct qname_idx ↔ NodeIndex correspondence)
    pub graph: PhasingGraph,
    
    /// Weight matrix using half-precision for 50% memory reduction
    pub weight_matrix: Option<Array2<f16>>,
    
    /// Read haplotype vectors
    pub read_hap_vectors: AHashMap<String, ReadHaplotypeVector>,
    
    /// Read error vectors  
    pub read_error_vectors: AHashMap<String, ReadErrorVector>,
    
    /// Read reference position dictionary
    pub read_ref_pos_dict: AHashMap<String, (i64, i64)>,
    
    /// Low quality/noisy qnames that were filtered out during BAM processing
    pub lowqual_qnames: HashSet<String>,
}

impl PhasingGraphResult {
    pub fn new() -> Self {
        Self {
            graph: Graph::new_undirected(),
            weight_matrix: None,
            read_hap_vectors: AHashMap::new(),
            read_error_vectors: AHashMap::new(),
            read_ref_pos_dict: AHashMap::new(),
            lowqual_qnames: HashSet::new(),
        }
    }
    
    /// Initialize the weight matrix with f16 precision
    pub fn initialize_weight_matrix(&mut self, size: usize) {
        if self.weight_matrix.is_none() {
            // Create identity matrix with f16 precision
            let mut matrix = Array2::<f16>::zeros((size, size));
            
            // Set diagonal to 1.0 (f16 can represent this exactly)
            for i in 0..size {
                matrix[[i, i]] = f16::ONE;
            }
            
            self.weight_matrix = Some(matrix);
        }
    }
    
    /// Set weight between two nodes with f16 precision
    pub fn set_weight(&mut self, i: usize, j: usize, weight: f32) -> Result<(), String> {
        match &mut self.weight_matrix {
            Some(ref mut matrix) => {
                // Check bounds
                if i >= matrix.nrows() || j >= matrix.ncols() {
                    return Err(format!("Index out of bounds: ({}, {}) for matrix of size {}x{}", 
                                     i, j, matrix.nrows(), matrix.ncols()));
                }
                
                // Convert f32 to f16 - this is where precision reduction happens
                let weight_f16 = f16::from_f32(weight);
                
                // Set symmetric values
                matrix[[i, j]] = weight_f16;
                matrix[[j, i]] = weight_f16;
                Ok(())
            }
            None => Err("Weight matrix not initialized. Call initialize_weight_matrix() first.".to_string())
        }
    }
    
    /// Get weight as f32 for calculations
    pub fn get_weight(&self, i: usize, j: usize) -> Option<f32> {
        self.weight_matrix.as_ref()?.get([i, j]).map(|&w| w.to_f32())
    }
    
    /// Get reference to weight matrix
    pub fn weight_matrix(&self) -> Option<&Array2<f16>> {
        self.weight_matrix.as_ref()
    }
    
    /// Get number of vertices in the graph
    pub fn vertex_count(&self) -> usize {
        self.graph.node_count()
    }
    
    /// Get number of edges in the graph  
    pub fn edge_count(&self) -> usize {
        self.graph.edge_count()
    }
}

/// Configuration for haplotype determination
/// Groups all the parameters needed for the determine_same_haplotype function
#[derive(Debug, Clone)]
pub struct HaplotypeConfig {
    pub mean_read_length: f32,
    pub edge_weight_cutoff: f32,
    pub score_array: Vec<f32>,  // Equivalent to Python's score_arr
}

impl HaplotypeConfig {
    pub fn new(mean_read_length: f32) -> Self {
        // Create score array equivalent to Python: [mean_read_length + mean_read_length * i for i in range(50)]
        let score_array: Vec<f32> = (0..50)
            .map(|i| mean_read_length + mean_read_length * i as f32)
            .collect();
            
        Self {
            mean_read_length,
            edge_weight_cutoff: 0.201,  // Default heuristic cutoff from Python
            score_array,
        }
    }
}

/// Variant types found in reads
#[derive(Debug, Clone)]
pub enum Variant {
    Snv { position: i64, alt_base: u8, quality: u8 },
    Insertion { position: i64, length: u32 },
    Deletion { position: i64, length: u32 },
}

impl Variant {
    pub fn position(&self) -> i64 {
        match self {
            Variant::Snv { position, .. } |
            Variant::Insertion { position, .. } |
            Variant::Deletion { position, .. } => *position,
        }
    }
}

/// Result of variant compatibility analysis
pub enum VariantCompatibility {
    Compatible(usize),  // Number of shared variants
    Incompatible,       // Conflicting variants found
    NoData,            // Not enough data to determine
}
