// Potential Rust implementation of pairwise_read_inspection.py, can be used as a replacement for the Python version, not tested yet
// File: fp_control/pairwise_read_inspection.rs
// This module creates a Python extension using Rust for high-performance read comparison

// === IMPORTS SECTION ===
// `use` is like Python's `import`. It brings items into scope
use pyo3::prelude::*;  // PyO3 main functionality - the * means "import everything"
use pyo3::types::{PyDict, PySet};  // Specific Python types we'll work with
use numpy::{PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};  // NumPy array handling
use ndarray::{Array1, Array2, ArrayView1, s};  // Rust's n-dimensional arrays (like NumPy)
use std::collections::{HashMap, HashSet};  // Rust's dictionary and set types
use rayon::prelude::*;  // Parallel iterator traits

// === TYPE ALIASES ===
// `type` creates shortcuts for complex types (like Python's type hints)
type HapVector = Array1<i16>;      // 1D array of 16-bit integers
type ErrorVector = Array1<f32>;    // 1D array of 32-bit floats
type SequenceArray = Array1<i8>;   // 1D array of 8-bit integers
type QualityArray = Array1<i8>;    // 1D array of 8-bit integers

// === CONSTANTS ===
// `const` defines compile-time constants (never change)
const MATCH: i16 = 1;
const SNV: i16 = -4;
const DELETION: i16 = -6;
const AMBIGUOUS_BASE: i8 = 4;  // 'N' base

// === STRUCT DEFINITION ===
// A `struct` is like a Python class that only holds data (no __init__ method)
// #[derive(Clone)] automatically implements the Clone trait (allows .clone())
#[derive(Clone)]
struct ReadData {
    read_id: String,                           // String type (owned, not borrowed)
    reference_start: i32,                      // 32-bit integer
    hap_vector: HapVector,                     // Using our type alias
    error_vector: ErrorVector,
    ref_positions: Array1<i32>,                // Generic Array1 with i32 elements
    qseq_ref_positions: Array1<i32>,
    query_sequence_encoded: SequenceArray,
    query_sequence_qualities: QualityArray,
}

// === MAIN CLASS EXPOSED TO PYTHON ===
// #[pyclass] tells PyO3 to make this struct available to Python
#[pyclass]
pub struct PairwiseInspectionEngine {
    // HashMap is Rust's dictionary, <String, ReadData> means String keys, ReadData values
    read_data_cache: HashMap<String, ReadData>,
    mean_read_length: f32,  // f32 is 32-bit float
    score_array: Array1<f32>,
}

// === PYTHON-CALLABLE METHODS ===
// #[pymethods] marks the following impl block as Python methods
#[pymethods]
impl PairwiseInspectionEngine {
    // #[new] marks this as the Python __init__ method
    #[new]
    fn new(mean_read_length: f32) -> Self {  // Self means PairwiseInspectionEngine
        // Create score array using iterator and collect
        // (0..50) creates a range from 0 to 49
        // .map() transforms each value (like Python's map())
        // .collect() gathers results into the specified type
        let score_array = (0..50)
            .map(|i| mean_read_length + mean_read_length * i as f32)  // `as` casts types
            .collect::<Array1<f32>>();  // ::<Type> specifies what to collect into
        
        // Return a new instance (no need for 'return' keyword)
        PairwiseInspectionEngine {
            read_data_cache: HashMap::new(),  // :: accesses associated functions
            mean_read_length,  // Shorthand for mean_read_length: mean_read_length
            score_array,
        }
    }

    // Main method exposed to Python
    // #[pyo3(signature = ...)] defines the Python function signature
    #[pyo3(signature = (
        read1, read2, overlap_start, overlap_end, 
        read_hap_vectors, read_error_vectors, 
        nested_ad_dict, read_ref_pos_dict,
        total_lowqual_qnames, empty_dict,
        intrinsic_ad_dict, base_dict
    ))]
    fn determine_same_haplotype(
        &mut self,  // &mut self means this method can modify the struct
        py: Python,  // Python GIL token - needed for Python API calls
        read1: PyObject,  // PyObject is a generic Python object
        read2: PyObject,
        overlap_start: i32,
        overlap_end: i32,
        read_hap_vectors: &PyDict,  // & means borrowed reference (doesn't own it)
        read_error_vectors: &PyDict,
        nested_ad_dict: &PyDict,
        read_ref_pos_dict: &PyDict,
        total_lowqual_qnames: &PySet,
        empty_dict: &PyDict,
        intrinsic_ad_dict: &PyDict,
        base_dict: &PyDict,
    ) -> PyResult<(PyObject, Py<PyDict>, Py<PyDict>, Py<PyDict>, Py<PySet>, Option<f32>)> {
        // PyResult is like Python's Union[ReturnType, Exception]
        // The return type is a tuple with 6 elements
        // Option<f32> is like Optional[float] - can be Some(value) or None
        
        // Extract read data - ? operator propagates errors (like Python's raise)
        let read1_data = self.extract_or_cache_read_data(
            py, &read1, read_hap_vectors, read_error_vectors, 
            read_ref_pos_dict, base_dict
        )?;  // ? means "if error, return error immediately"
        
        // Check if reads are noisy - destructuring tuple return
        let (read1_noisy, _) = self.check_noisy_read(&read1_data)?;  // _ ignores value
        let (read2_noisy, _) = self.check_noisy_read(&read2_data)?;
        
        // Clone the PySet (makes a copy)
        let mut updated_lowqual = total_lowqual_qnames.clone();  // mut means mutable
        
        // Conditional logic
        if read1_noisy {
            let qname = read1.getattr(py, "query_name")?;  // Like Python's getattr()
            updated_lowqual.add(qname)?;
        }
        if read2_noisy {
            let qname = read2.getattr(py, "query_name")?;
            updated_lowqual.add(qname)?;
        }

        // Get haplotype vectors for overlap region
        let interval_hap1 = self.slice_vector(
            &read1_data.hap_vector, 
            overlap_start, 
            overlap_end, 
            read1_data.reference_start
        );
        let interval_hap2 = self.slice_vector(
            &read2_data.hap_vector,
            overlap_start,
            overlap_end,
            read2_data.reference_start
        );

        // Find differences between haplotype vectors
        let diff_indices = self.find_diff_indices(&interval_hap1, &interval_hap2);
        
        // Early rejection if too many differences
        if diff_indices.len() >= 3 {
            // Return tuple - Ok() wraps successful result
            return Ok((
                py.False().into(),  // .into() converts types
                read_ref_pos_dict.into(),  // Converts &PyDict to Py<PyDict>
                read_hap_vectors.into(), 
                read_error_vectors.into(),
                updated_lowqual.into(),
                None  // None in Option<f32>
            ));
        }

        // Extract sequences in overlap region
        let (seq1, qr_idx1) = self.get_interval_seq(
            &read1_data,
            overlap_start,
            overlap_end - 1
        );
        let (seq2, qr_idx2) = self.get_interval_seq(
            &read2_data,
            overlap_start,
            overlap_end - 1
        );

        // Compare sequences
        let sequences_match = self.compare_sequences(&seq1, &seq2, AMBIGUOUS_BASE);
        
        // Calculate shared variants and weight
        let (psv_snv_count, shared_snv_count) = self.count_shared_snvs(
            &interval_hap1, &interval_hap2, overlap_start, overlap_end,
            &read1_data, &read2_data, intrinsic_ad_dict
        )?;

        // Calculate weight
        let weight = self.calculate_weight(
            &interval_hap1, &interval_hap2,
            shared_snv_count, psv_snv_count,
            overlap_start, overlap_end
        );

        // Pattern matching on boolean
        if sequences_match {
            return Ok((
                py.True().into(),
                read_ref_pos_dict.into(),
                read_hap_vectors.into(),
                read_error_vectors.into(),
                updated_lowqual.into(),
                Some(weight)  // Some() wraps a value in Option
            ));
        }

        // Check if we can tolerate mismatches
        if seq1.len() != seq2.len() || interval_hap1.len() != interval_hap2.len() {
            return Ok((
                py.False().into(),
                read_ref_pos_dict.into(),
                read_hap_vectors.into(),
                read_error_vectors.into(),
                updated_lowqual.into(),
                None
            ));
        }

        // Find sequence differences and check tolerance
        let seq_diff_indices = self.find_diff_indices_seq(&seq1, &seq2);
        // Method chaining with iterator
        let ref_diff_indices = seq_diff_indices.iter()  // Create iterator
            .map(|&i| qr_idx1[i])  // |&i| is a closure (lambda) that destructures
            .collect::<Vec<_>>();  // Vec<_> lets Rust infer the type

        // Check if any difference is on an indel
        if ref_diff_indices.iter().any(|&idx| idx == -1) {
            return Ok((
                py.False().into(),
                read_ref_pos_dict.into(),
                read_hap_vectors.into(),
                read_error_vectors.into(),
                updated_lowqual.into(),
                None
            ));
        }

        // Check if we can tolerate the mismatches
        let (can_tolerate, tolerated_count) = self.tolerate_mismatches_two_seq(
            &read1_data, &read2_data,
            &ref_diff_indices, nested_ad_dict, empty_dict, py
        )?;

        // Conditional expression
        let final_weight = if can_tolerate {
            Some(weight - tolerated_count as f32 * 20.0)
        } else {
            None
        };

        // Ternary-like expression
        Ok((
            if can_tolerate { py.True() } else { py.False() }.into(),
            read_ref_pos_dict.into(),
            read_hap_vectors.into(),
            read_error_vectors.into(),
            updated_lowqual.into(),
            final_weight
        ))
    }
}

// === PRIVATE IMPLEMENTATION ===
// This impl block is not exposed to Python (no #[pymethods])
impl PairwiseInspectionEngine {
    // Private method - note no 'pub' keyword
    fn extract_or_cache_read_data(
        &mut self,
        py: Python,
        pysam_read: &PyObject,
        hap_vectors: &PyDict,
        error_vectors: &PyDict,
        ref_pos_dict: &PyDict,
        base_dict: &PyDict,
    ) -> PyResult<ReadData> {
        // Convert pysam to Rust struct once
        let rust_read = PysameRead::from_pysam(py, pysam_read)?;
        let read_id = rust_read.get_read_id();
        
        // Check cache
        if let Some(cached) = self.read_data_cache.get(&read_id) {
            return Ok(cached.clone());
        }

        // Extract read attributes
        let reference_start: i32 = read.getattr(py, "reference_start")?.extract(py)?;
        
        // Get or compute haplotype vector
        let hap_vector = if let Ok(cached_hap) = hap_vectors.get_item(&read_id) {
            // Convert from Python numpy array
            let py_array: &PyArray1<i16> = cached_hap.downcast(py)?;
            py_array.to_owned_array()
        } else {
            // Compute haplotype vector from CIGAR
            let cigar_tuples = read.getattr(py, "cigartuples")?;
            let query_seq = read.getattr(py, "query_sequence")?.extract::<String>(py)?;
            self.compute_haplotype_vector(py, cigar_tuples, &query_seq, base_dict)?
        };

        // Get or compute error vector
        let error_vector = if let Ok(cached_err) = error_vectors.get_item(&read_id) {
            let py_array: &PyArray1<f32> = cached_err.downcast(py)?;
            py_array.to_owned_array()
        } else {
            let cigar_tuples = read.getattr(py, "cigartuples")?;
            self.compute_error_vector(py, read, cigar_tuples)?
        };

        // Get reference positions
        let (ref_positions, qseq_ref_positions, query_sequence_encoded, query_sequence_qualities) = 
            if let Ok(cached_data) = ref_pos_dict.get_item(&read_id) {
                // Unpack tuple from cache
                let tuple = cached_data.downcast::<pyo3::types::PyTuple>(py)?;
                (
                    tuple.get_item(0)?.downcast::<PyArray1<i32>>()?.to_owned_array(),
                    tuple.get_item(1)?.downcast::<PyArray1<i32>>()?.to_owned_array(),
                    tuple.get_item(2)?.downcast::<PyArray1<i8>>()?.to_owned_array(),
                    tuple.get_item(3)?.downcast::<PyArray1<i8>>()?.to_owned_array(),
                )
            } else {
                self.extract_read_sequences(py, read, base_dict)?
            };

        let read_data = ReadData {
            read_id: read_id.clone(),
            reference_start: rust_read.reference_start,
            hap_vector: Array1::from(hap_vector),
            error_vector: Array1::from(error_vector),
            ref_positions,
            qseq_ref_positions: Array1::zeros(0),
            query_sequence_encoded: Array1::zeros(0),
            query_sequence_qualities: Array1::zeros(0),
        };

        // Cache the data
        self.read_data_cache.insert(read_id, read_data.clone());
        Ok(read_data)
    }

    /// Generate read ID from pysam read object
    fn get_read_id(&self, py: Python, read: &PyObject) -> PyResult<String> {
        let query_name: String = read.getattr(py, "query_name")?.extract(py)?;
        let flag: u16 = read.getattr(py, "flag")?.extract(py)?;
        Ok(format!("{}:{}", query_name, flag))
    }

    /// Slice a vector based on genomic coordinates
    fn slice_vector(&self, vec: &Array1<i16>, start: i32, end: i32, read_start: i32) -> Array1<i16> {
        let start_idx = (start - read_start) as usize;
        let end_idx = (end - read_start) as usize;
        vec.slice(s![start_idx..end_idx]).to_owned()
    }

    /// Find indices where two arrays differ
    fn find_diff_indices(&self, arr1: &Array1<i16>, arr2: &Array1<i16>) -> Vec<usize> {
        arr1.iter()
            .zip(arr2.iter())
            .enumerate()
            .filter_map(|(i, (a, b))| if a != b { Some(i) } else { None })
            .collect()
    }

    /// Find indices where two sequence arrays differ
    fn find_diff_indices_seq(&self, arr1: &Array1<i8>, arr2: &Array1<i8>) -> Vec<usize> {
        arr1.iter()
            .zip(arr2.iter())
            .enumerate()
            .filter_map(|(i, (a, b))| if a != b { Some(i) } else { None })
            .collect()
    }

    /// Check if a read is noisy based on error probabilities
    fn check_noisy_read(&self, read_data: &ReadData) -> PyResult<(bool, usize)> {
        // Count SNVs with high error probability
        let snv_positions: Vec<usize> = read_data.hap_vector.iter()
            .enumerate()
            .filter_map(|(i, &val)| if val == SNV { Some(i) } else { None })
            .collect();
        
        let high_error_snvs = snv_positions.iter()
            .filter(|&&pos| read_data.error_vector[pos] > 0.03)
            .count();
        
        Ok((high_error_snvs >= 3, high_error_snvs))
    }

    /// Compare two sequences, ignoring ambiguous bases
    fn compare_sequences(&self, seq1: &Array1<i8>, seq2: &Array1<i8>, except_char: i8) -> bool {
        if seq1.len() != seq2.len() {
            return false;
        }
        
        seq1.iter().zip(seq2.iter()).all(|(a, b)| {
            a == b || *a == except_char || *b == except_char
        })
    }

    /// Extract interval sequence from read data
    fn get_interval_seq(
        &self,
        read_data: &ReadData,
        interval_start: i32,
        interval_end: i32,
    ) -> (Array1<i8>, Array1<i32>) {
        let read_start = read_data.reference_start;
        let ref_positions = &read_data.ref_positions;
        let qseq_ref_positions = &read_data.qseq_ref_positions;
        let query_sequence = &read_data.query_sequence_encoded;
        
        // Find start index in query sequence
        let mut start_idx = interval_start - read_start;
        let mut interval_start_qidx = ref_positions[start_idx as usize];
        
        while interval_start_qidx < 0 && start_idx <= interval_end - read_start {
            start_idx += 1;
            interval_start_qidx = ref_positions[start_idx as usize];
        }
        
        if start_idx > interval_end - read_start {
            return (Array1::zeros(0), Array1::zeros(0));
        }
        
        // Find end index
        let mut end_idx = interval_end - read_start;
        let mut interval_end_qidx = ref_positions[end_idx as usize];
        
        while interval_end_qidx < 0 {
            end_idx -= 1;
            interval_end_qidx = ref_positions[end_idx as usize];
        }
        
        // Extract sequence and positions
        let start = interval_start_qidx as usize;
        let end = (interval_end_qidx + 1) as usize;
        
        (
            query_sequence.slice(s![start..end]).to_owned(),
            qseq_ref_positions.slice(s![start..end]).to_owned()
        )
    }

    /// Count shared SNVs between two reads
    fn count_shared_snvs(
        &self,
        hap1: &Array1<i16>,
        hap2: &Array1<i16>,
        overlap_start: i32,
        overlap_end: i32,
        read1_data: &ReadData,
        read2_data: &ReadData,
        intrinsic_ad_dict: &PyDict,
    ) -> PyResult<(i32, i32)> {
        let mut shared_snvs = 0;
        let mut psv_snvs = 0;
        
        // Find positions where both have SNVs
        for (i, (&val1, &val2)) in hap1.iter().zip(hap2.iter()).enumerate() {
            if val1 == SNV && val2 == SNV {
                shared_snvs += 1;
                let genomic_pos = overlap_start + i as i32;
                
                // Check if supported by AD data
                if let Ok(pos_dict) = intrinsic_ad_dict.get_item(genomic_pos) {
                    if let Ok(pos_dict) = pos_dict.downcast::<PyDict>() {
                        // Extract bases at this position
                        let (seq1, _) = self.get_interval_seq(
                            read1_data, genomic_pos, genomic_pos + 1
                        );
                        let (seq2, _) = self.get_interval_seq(
                            read2_data, genomic_pos, genomic_pos + 1
                        );
                        
                        if seq1.len() > 0 && seq2.len() > 0 && seq1[0] == seq2[0] {
                            let alt_base = seq1[0];
                            if let Ok(ad) = pos_dict.get_item(alt_base) {
                                if let Ok(ad_val) = ad.extract::<i16>() {
                                    if ad_val > 0 {
                                        psv_snvs += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        Ok((psv_snvs, shared_snvs))
    }

    /// Calculate weight for read pair similarity
    fn calculate_weight(
        &self,
        hap1: &Array1<i16>,
        hap2: &Array1<i16>,
        shared_snv_count: i32,
        psv_snv_count: i32,
        overlap_start: i32,
        overlap_end: i32,
    ) -> f32 {
        // Find identical positions
        let identical_positions: Vec<usize> = hap1.iter()
            .zip(hap2.iter())
            .enumerate()
            .filter_map(|(i, (a, b))| if a == b { Some(i) } else { None })
            .collect();
        
        let overlap_span = identical_positions.len() as f32;
        
        // Count indel blocks in identical part
        let mut indel_blocks = 0;
        let mut in_indel = false;
        for &pos in &identical_positions {
            let is_indel = hap1[pos] > 1 || hap1[pos] == DELETION;
            if is_indel && !in_indel {
                indel_blocks += 1;
                in_indel = true;
            } else if !is_indel {
                in_indel = false;
            }
        }
        
        // Calculate weight using score array
        let mut weight = overlap_span;
        
        // Add scores for shared SNVs
        for i in 0..shared_snv_count.min(50) {
            weight += self.score_array[i as usize];
        }
        
        // Adjust for PSV vs non-PSV SNVs
        weight += self.mean_read_length * (shared_snv_count - psv_snv_count) as f32 * 0.75;
        
        // Add weight for indels
        weight += self.mean_read_length * 3.0 * indel_blocks as f32;
        
        weight
    }

    /// Check if mismatches can be tolerated based on sequencing artifacts
    fn tolerate_mismatches_two_seq(
        &self,
        read1_data: &ReadData,
        read2_data: &ReadData,
        diff_indices: &[i32],
        nested_ad_dict: &PyDict,
        empty_dict: &PyDict,
        py: Python,
    ) -> PyResult<(bool, usize)> {
        let mut tolerate_count = 0;
        
        for &diff_ind in diff_indices {
            let seq_err1 = self.check_seq_error_at_position(
                read1_data, diff_ind, nested_ad_dict, empty_dict, py
            )?;
            let seq_err2 = self.check_seq_error_at_position(
                read2_data, diff_ind, nested_ad_dict, empty_dict, py
            )?;
            
            if seq_err1 || seq_err2 {
                tolerate_count += 1;
            } else {
                return Ok((false, 0));
            }
        }
        
        Ok((true, tolerate_count))
    }

    /// Check if a position has sequencing error based on AD data
    fn check_seq_error_at_position(
        &self,
        read_data: &ReadData,
        position: i32,
        nested_ad_dict: &PyDict,
        empty_dict: &PyDict,
        py: Python,
    ) -> PyResult<bool> {
        let ridx = (position - read_data.reference_start) as usize;
        let qidx = read_data.ref_positions[ridx];
        
        if qidx < 0 {
            return Ok(false);
        }
        
        let base = read_data.query_sequence_encoded[qidx as usize];
        let qual = read_data.query_sequence_qualities[qidx as usize];
        
        if qual >= 15 {
            return Ok(false);
        }
        
        // Get AD dict for this position
        let pos_dict = if let Ok(dict) = nested_ad_dict.get_item(position) {
            dict.downcast::<PyDict>(py)?
        } else {
            empty_dict
        };
        
        // Get allele depth and total depth
        let ad = if let Ok(ad_val) = pos_dict.get_item(base) {
            ad_val.extract::<i16>(py)?
        } else {
            0
        };
        
        let dp = if let Ok(dp_val) = pos_dict.get_item(5i8) {  // 5 is DP key
            dp_val.extract::<i16>(py)?
        } else {
            return Ok(false);
        };
        
        if dp == 0 {
            return Ok(false);
        }
        
        let af = ad as f32 / dp as f32;
        
        Ok((af <= 0.02 || (ad == 1 && dp >= 10)) && qual < 13)
    }

    /// Compute haplotype vector from CIGAR string
    fn compute_haplotype_vector(
        &self,
        py: Python,
        cigar_tuples: PyObject,
        query_sequence: &str,
        base_dict: &PyDict,
    ) -> PyResult<Array1<i16>> {
        // This is a simplified version - you'd need to implement full CIGAR parsing
        // For now, returning a placeholder
        let cigar_list = cigar_tuples.extract::<Vec<(u8, u32)>>(py)?;
        
        // Calculate reference length
        let ref_length: u32 = cigar_list.iter()
            .filter(|(op, _)| matches!(op, 0 | 7 | 8 | 2 | 3))
            .map(|(_, len)| len)
            .sum();
        
        let mut hap_vector = Array1::zeros(ref_length as usize);
        let mut ref_pos = 0;
        let mut query_pos = 0;
        let mut insertion = 0;
        
        for (operation, length) in cigar_list {
            match operation {
                7 => { // Match
                    for i in 0..*length as usize {
                        if insertion > 0 {
                            hap_vector[ref_pos] = insertion;
                            insertion = 0;
                            if i > 0 {
                                hap_vector[ref_pos + 1..ref_pos + *length as usize].fill(MATCH);
                            }
                        } else {
                            hap_vector[ref_pos..ref_pos + *length as usize].fill(MATCH);
                        }
                    }
                    ref_pos += *length as usize;
                    query_pos += *length as usize;
                },
                8 => { // Mismatch
                    for i in 0..*length as usize {
                        let base = &query_sequence[query_pos..query_pos + 1];
                        if base == "N" {
                            hap_vector[ref_pos + i] = MATCH;
                        } else {
                            hap_vector[ref_pos + i] = SNV;
                        }
                    }
                    ref_pos += *length as usize;
                    query_pos += *length as usize;
                },
                1 => { // Insertion
                    if ref_pos > 0 {
                        insertion = (*length as i16) * 4;
                    }
                    query_pos += *length as usize;
                },
                2 => { // Deletion
                    hap_vector[ref_pos..ref_pos + *length as usize].fill(DELETION);
                    ref_pos += *length as usize;
                },
                3 => { // Skip
                    hap_vector[ref_pos..ref_pos + *length as usize].fill(MATCH);
                    ref_pos += *length as usize;
                },
                4 => { // Soft clip
                    query_pos += *length as usize;
                },
                _ => {}
            }
        }
        
        Ok(hap_vector)
    }

    /// Compute error vector from CIGAR and quality scores
    fn compute_error_vector(
        &self,
        py: Python,
        read: &PyObject,
        cigar_tuples: PyObject,
    ) -> PyResult<Array1<f32>> {
        let reference_start: i32 = read.getattr(py, "reference_start")?.extract(py)?;
        let reference_end: i32 = read.getattr(py, "reference_end")?.extract(py)?;
        let qualities = read.getattr(py, "query_qualities")?
            .extract::<Vec<u8>>(py)?;
        
        let cigar_list = cigar_tuples.extract::<Vec<(u8, u32)>>(py)?;
        let ref_length = (reference_end - reference_start) as usize;
        
        let mut error_vector = Array1::zeros(ref_length);
        let mut ref_pos = 0;
        let mut query_pos = 0;
        
        for (operation, length) in cigar_list {
            match operation {
                7 | 8 => { // Match or mismatch
                    for i in 0..*length as usize {
                        let qual = qualities[query_pos + i] as f32;
                        error_vector[ref_pos + i] = 10f32.powf(-qual / 10.0);
                    }
                    ref_pos += *length as usize;
                    query_pos += *length as usize;
                },
                1 => { // Insertion
                    if ref_pos > 0 {
                        error_vector[ref_pos - 1] = 0.0;
                    }
                    query_pos += *length as usize;
                },
                2 => { // Deletion
                    for i in 0..*length as usize {
                        error_vector[ref_pos + i] = 0.0;
                    }
                    ref_pos += *length as usize;
                },
                3 => { // Skip
                    for i in 0..*length as usize {
                        error_vector[ref_pos + i] = 0.0;
                    }
                    ref_pos += *length as usize;
                },
                4 => { // Soft clip
                    query_pos += *length as usize;
                },
                _ => {}
            }
        }
        
        Ok(error_vector)
    }

    /// Extract read sequences and positions
    fn extract_read_sequences(
        &self,
        py: Python,
        read: &PyObject,
        base_dict: &PyDict,
    ) -> PyResult<(Array1<i32>, Array1<i32>, Array1<i8>, Array1<i8>)> {
        // Get reference positions from pysam
        let ref_positions_method = read.getattr(py, "get_reference_positions")?;
        let kwargs = pyo3::types::PyDict::new(py);
        kwargs.set_item("full_length", true)?;
        let ref_positions_list = ref_positions_method.call((), Some(kwargs))?;
        
        // Convert to array, replacing None with -1
        let qseq_ref_positions: Array1<i32> = ref_positions_list
            .extract::<Vec<Option<i32>>>(py)?
            .into_iter()
            .map(|x| x.unwrap_or(-1))
            .collect();
        
        // Get query sequence
        let query_sequence = read.getattr(py, "query_sequence")?.extract::<String>(py)?;
        
        // Encode sequence
        let query_sequence_encoded: Array1<i8> = query_sequence
            .chars()
            .map(|c| {
                let c_str = c.to_string();
                base_dict.get_item(c_str)
                    .and_then(|x| x.extract::<i8>().ok())
                    .unwrap_or(4)
            })
            .collect();
        
        // Get quality scores
        let qualities = read.getattr(py, "query_qualities")?;
        let query_sequence_qualities: Array1<i8> = if qualities.is_none(py) {
            Array1::zeros(query_sequence.len())
        } else {
            qualities.extract::<Vec<u8>>(py)?
                .into_iter()
                .map(|x| x as i8)
                .collect()
        };
        
        // Prepare reference to query index map
        let reference_start: i32 = read.getattr(py, "reference_start")?.extract(py)?;
        let max_ref_pos = qseq_ref_positions.iter().max().copied().unwrap_or(0);
        let ref_positions = self.prepare_ref_query_idx_map(
            &qseq_ref_positions, 
            reference_start,
            max_ref_pos
        );
        
        Ok((ref_positions, qseq_ref_positions, query_sequence_encoded, query_sequence_qualities))
    }

    /// Prepare reference to query index mapping
    fn prepare_ref_query_idx_map(
        &self,
        qseq_ref_positions: &Array1<i32>,
        read_ref_start: i32,
        max_ref_pos: i32,
    ) -> Array1<i32> {
        let size = (max_ref_pos - read_ref_start + 1) as usize;
        let mut ref_positions = Array1::from_elem(size, -1i32);
        
        for (qi, &ri) in qseq_ref_positions.iter().enumerate() {
            if ri >= 0 {
                ref_positions[(ri - read_ref_start) as usize] = qi as i32;
            }
        }
        
        ref_positions
    }
}

/// Create the Python module
#[pymodule]
fn pairwise_inspection(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PairwiseInspectionEngine>()?;
    Ok(())
}