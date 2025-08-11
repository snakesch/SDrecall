use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use numpy::{PyArray1, PyArray2, ToPyArray};

use log::{LevelFilter, info, debug, Record, Level, Metadata};
use std::io::{self, Write};

use crate::structs::{HaplotypeConfig, PhasingGraphResult};
use crate::graph_builder::build_phasing_graph;

/// Simple custom logger that provides detailed logging without SIMD dependencies
struct SimpleEnhancedLogger;

impl log::Log for SimpleEnhancedLogger {
    fn enabled(&self, _metadata: &Metadata) -> bool {
        true
    }

    fn log(&self, record: &Record) {
        if self.enabled(record.metadata()) {
            // Get current time using standard library (no external dependencies)
            let now = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap_or_default()
                .as_secs();
            
            let module = record.module_path().unwrap_or("unknown_module");
            let file = record.file().unwrap_or("unknown_file");
            let line = record.line().unwrap_or(0);
            
            eprintln!(
                "[{}] [{}] [{}] [{}:{}] {}",
                now,
                record.level(),
                module,
                file,
                line,
                record.args()
            );
        }
    }

    fn flush(&self) {
        io::stderr().flush().ok();
    }
}

/// Configure comprehensive Rust logging with file names, line numbers, and function stacks
/// 
/// This function sets up detailed logging that includes:
/// - File names and line numbers where the log was generated
/// - Function names and call stack context
/// - Proper level filtering based on Python logging level
/// 
/// Uses minimal custom logger to avoid any SIMD compatibility issues
/// 
/// # Arguments
/// * `level_str` - Python logging level as string ("DEBUG", "INFO", "WARNING", "ERROR")
fn configure_rust_logging(level_str: &str) {
    let rust_level = match level_str.to_uppercase().as_str() {
        "DEBUG" => LevelFilter::Debug,
        "INFO" => LevelFilter::Info, 
        "WARNING" | "WARN" => LevelFilter::Warn,
        "ERROR" => LevelFilter::Error,
        "CRITICAL" => LevelFilter::Error,
        _ => LevelFilter::Info, // Default fallback
    };
    
    // Initialize our custom logger (only once)
    static INIT: std::sync::Once = std::sync::Once::new();
    INIT.call_once(|| {
        log::set_boxed_logger(Box::new(SimpleEnhancedLogger))
            .map(|()| log::set_max_level(rust_level))
            .ok();
    });
    
    // Update the max level for this configuration
    log::set_max_level(rust_level);
    
    info!("[configure_rust_logging] Rust logging configured with level: {:?} (from Python: {})", rust_level, level_str);
    debug!("[configure_rust_logging] Enhanced logging with source location is active (minimal SIMD-free implementation)");
}

/// Python-exposed function for building the phasing graph
/// 
/// This function serves as the bridge between Python and Rust, handling all necessary
/// data conversions and returning results in a format compatible with graph-tool.
#[pyfunction]
#[pyo3(signature = (
    bam_file_path,
    reference_genome,
    mean_read_length,
    _edge_weight_cutoff=0.201,
    mapq_filter=10,
    basequal_median_filter=10,
    filter_noisy=true,
    use_collate=true,
    threads=4,
    log_level=None
))]
pub fn build_phasing_graph_rust(
    py: Python<'_>,
    bam_file_path: &str,
    reference_genome: &str,
    mean_read_length: f32,
    _edge_weight_cutoff: f32,
    mapq_filter: u8,
    basequal_median_filter: u8,
    filter_noisy: bool,
    use_collate: bool,
    threads: u8,
    log_level: Option<&str>,
) -> PyResult<PyObject> {
    info!("[build_phasing_graph_rust] Starting Rust graph building with mean_read_length={}", mean_read_length);
    
    // Configure Rust logging level based on Python logger level
    if let Some(level_str) = log_level {
        configure_rust_logging(level_str);
    } else {
        log::set_max_level(LevelFilter::Debug);
    }
    
    info!("[build_phasing_graph_rust] Processing BAM file: {}", bam_file_path);
    
    // Read BAM file and build read pair map directly in Rust
    info!("[build_phasing_graph_rust] Step 1: Reading BAM file and building read pair map");
    let (read_pair_map, header) = crate::bam_reading::migrate_bam_to_sorted_intervals_grouped(
        bam_file_path,
        mapq_filter,
        basequal_median_filter,
        filter_noisy,
        use_collate,
        threads,
    ).map_err(|e| {
        let error_msg = format!("BAM reading failed: {}", e);
        log::error!("[build_phasing_graph_rust] Step 1 error: {}", error_msg);
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(error_msg)
    })?;
    
    info!("[build_phasing_graph_rust] Step 1 completed: Read pair map contains {} pairs", read_pair_map.readpair_dict.len());
    
    // Build allele depth map
    info!("[build_phasing_graph_rust] Step 2: Building allele depth map");
    let allele_depth_map = crate::bam_reading::build_allele_depth_map(
        bam_file_path,
        reference_genome,
        mapq_filter,
        basequal_median_filter,
    ).map_err(|e| {
        let error_msg = format!("Allele depth mapping failed: {}", e);
        log::error!("[build_phasing_graph_rust] Step 2 error: {}", error_msg);
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(error_msg)
    })?;
    
    info!("[build_phasing_graph_rust] Step 2 completed: Allele depth map has {} chromosomes", allele_depth_map.chromosome_count());
    
    // Create configuration
    info!("[build_phasing_graph_rust] Step 3: Creating haplotype configuration");
    let config = HaplotypeConfig::new(mean_read_length);
    
    // Call the Rust graph building function
    info!("[build_phasing_graph_rust] Step 4: Building phasing graph");
    let result = build_phasing_graph(
        &read_pair_map,
        &allele_depth_map,
        &header,
        &config,
    ).map_err(|e| {
        let error_msg = format!("Graph building failed: {}", e);
        log::error!("[build_phasing_graph_rust] Step 4 error: {}", error_msg);
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(error_msg)
    })?;
    
    info!("[build_phasing_graph_rust] Step 4 completed: Graph has {} vertices and {} edges", 
              result.graph.node_count(), result.graph.edge_count());
    
    // Convert the result back to Python format
    info!("[build_phasing_graph_rust] Step 5: Converting result to Python format");
    export_graph_result_to_python(py, result)
}

// Note: BAM reading and allele depth mapping are now handled directly in Rust
// The unimplemented conversion functions have been removed since we process
// BAM files directly rather than converting from Python data structures

/// Export the Rust graph result to Python format compatible with graph-tool
fn export_graph_result_to_python(
    py: Python<'_>,
    result: PhasingGraphResult,
) -> PyResult<PyObject> {
    let dict = PyDict::new_bound(py);
    
    // 1. Export edges with integer node indices
    let num_edges = result.graph.edge_count();
    let mut edges = Vec::with_capacity(num_edges);
    let mut weights = Vec::with_capacity(num_edges);

    for edge in result.graph.edge_indices() {
        let (source, target) = result.graph.edge_endpoints(edge).unwrap();
        // With direct correspondence: qname_idx = NodeIndex.index()
        let source_idx = source.index();
        let target_idx = target.index();
        edges.push([source_idx as u32, target_idx as u32]);
        
        // Convert f16 weight to f32 for Python
        let weight_f16 = result.graph[edge];
        weights.push(weight_f16.to_f32());
    }
    
    // Convert to numpy arrays
    let edges_vec: Vec<Vec<u32>> = edges.iter().map(|&[a, b]| vec![a, b]).collect();
    let edge_array = PyArray2::from_vec2_bound(py, &edges_vec)?;
    let weight_array = PyArray1::from_vec_bound(py, weights);
    
    dict.set_item("edges", edge_array)?;
    dict.set_item("weights", weight_array)?;
    
    // 2. Export node names in the correct order (index i contains qname for node i)
    // With direct correspondence, we need to get qnames from read_pair_map
    let mut node_names = vec![String::new(); result.graph.node_count()];
    // TODO: This would need the read_pair_map passed in to get qnames by index
    // For now, we'll create placeholder names based on index
    for i in 0..result.graph.node_count() {
        node_names[i] = format!("read_pair_{}", i);
    }
    let node_list = PyList::new_bound(py, node_names);
    dict.set_item("vertex_names", node_list)?;
    
    // 3. Export weight matrix as numpy array - much simpler with ndarray!
    match result.weight_matrix() {
        Some(matrix) => {
            // Convert f16 matrix to f32 for Python compatibility
            let f32_matrix: ndarray::Array2<f32> = matrix.mapv(|x| x.to_f32());
            let weight_matrix_array = f32_matrix.to_pyarray_bound(py);
            dict.set_item("weight_matrix", weight_matrix_array)?;
        }
        None => {
            // Create empty f32 matrix for Python
            let empty_matrix = ndarray::Array2::<f32>::zeros((0, 0));
            let weight_matrix_array = empty_matrix.to_pyarray_bound(py);
            dict.set_item("weight_matrix", weight_matrix_array)?;
        }
    }
    
    // 4. Note: qname_to_node mapping removed due to direct correspondence optimization
    // Direct mapping: qname_idx = NodeIndex.index()
    
    // 5. Export read_hap_vectors
    let hap_vectors_dict = PyDict::new_bound(py);
    for (qname, vec) in &result.read_hap_vectors {
        let py_array = PyArray1::from_vec_bound(py, vec.clone());
        hap_vectors_dict.set_item(qname, py_array)?;
    }
    dict.set_item("read_hap_vectors", hap_vectors_dict)?;
    
    // 6. Export read_error_vectors
    let error_vectors_dict = PyDict::new_bound(py);
    for (qname, vec) in &result.read_error_vectors {
        let py_array = PyArray1::from_vec_bound(py, vec.clone());
        error_vectors_dict.set_item(qname, py_array)?;
    }
    dict.set_item("read_error_vectors", error_vectors_dict)?;
    
    // 7. Export read_ref_pos_dict
    let ref_pos_dict = PyDict::new_bound(py);
    for (qname, (start, end)) in &result.read_ref_pos_dict {
        let tuple = PyTuple::new_bound(py, &[start, end]);
        ref_pos_dict.set_item(qname, tuple)?;
    }
    dict.set_item("read_ref_pos_dict", ref_pos_dict)?;
    
    // 8. Export low_qual_qnames (noisy qnames identified during BAM processing)
    let low_qual_vec: Vec<&str> = result.lowqual_qnames.iter().map(|s| s.as_str()).collect();
    let low_qual_list = PyList::new_bound(py, low_qual_vec);
    dict.set_item("low_qual_qnames", low_qual_list)?;
    
    // 9. Add metadata
    dict.set_item("num_vertices", result.graph.node_count())?;
    dict.set_item("num_edges", result.graph.edge_count())?;
    
    Ok(dict.into())
}