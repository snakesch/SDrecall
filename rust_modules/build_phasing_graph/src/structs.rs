use std::collections::HashMap;
use rust_htslib::bam::Record;

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
    // Use fixed-size arrays instead of Vec for exactly 2 intervals per qname_idx
    interval_map: HashMap<usize, [Option<Interval>; 2]>,
    // For querying phase: pre-allocated vector for fast overlap queries
    sorted_intervals: Option<Vec<Interval>>,
    is_finalized: bool,
}

impl SortedVecIntervals {
    pub fn new() -> Self {
        Self {
            interval_map: HashMap::new(),
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
    pub interval_trees: HashMap<String, SortedVecIntervals>,
    pub readpair_dict: HashMap<usize, ReadPair>,
    pub qname_to_idx: HashMap<String, usize>,
    pub idx_to_qname: HashMap<usize, String>,
    pub noisy_qnames: HashSet<String>,
}

impl ReadPairMap {
    pub fn new() -> Self {
        Self {
            interval_trees: HashMap::new(),
            readpair_dict: HashMap::new(),
            qname_to_idx: HashMap::new(),
            idx_to_qname: HashMap::new(),
            noisy_qnames: HashSet::new(),
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