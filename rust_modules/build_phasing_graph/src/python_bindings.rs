use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use numpy::{PyArray1, PyArray2, ToPyArray};

use log::{LevelFilter, info, debug, warn, error};
use pyo3_log;

use crate::structs::{HaplotypeConfig, PhasingGraphResult, ReadPairMap};
use crate::graph_builder::build_phasing_graph;


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
        "ERROR" | "CRITICAL" => LevelFilter::Error,
        _ => LevelFilter::Info,
    };

    // Initialize pyo3-log once per process; safe to call multiple times
    let _ = pyo3_log::try_init();
    log::set_max_level(rust_level);
	// Print the log level to stderr

    info!("[configure_rust_logging] Rust logging configured with level: {:?} (from Python: {})", rust_level, level_str);
    debug!("[configure_rust_logging] pyo3-log forwarding active");
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
    
    // Configure Rust logging level based on Python logger level
    if let Some(level_str) = log_level {
        configure_rust_logging(level_str);
    } else {
		println!("Rust logging configured with level: Debug (from Python: None), configure_rust_logging is not called");
        log::set_max_level(LevelFilter::Debug);
    }
    
    info!("[build_phasing_graph_rust] Starting Rust graph building with mean_read_length={}, the log level is {}", mean_read_length, log_level.unwrap_or("unknown"));
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
        error!("[build_phasing_graph_rust] Step 1 error: {}", error_msg);
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
        error!("[build_phasing_graph_rust] Step 2 error: {}", error_msg);
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
        error!("[build_phasing_graph_rust] Step 4 error: {}", error_msg);
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(error_msg)
    })?;

    // If result.weight_matrix is None, return None early
    if result.weight_matrix.is_none() {
        warn!("[build_phasing_graph_rust] Step 4 completed: Graph has {} vertices and {} edges. Skipping this region.", 
              result.graph.node_count(), result.graph.edge_count());
        return Ok(py.None().into_py(py));
    }

    // If result.graph.node_count() <= 2, return None early
    if result.graph.node_count() <= 2 {
        warn!("[build_phasing_graph_rust] Step 4 completed: Graph has {} vertices and {} edges. Skipping this region.", 
              result.graph.node_count(), result.graph.edge_count());
        return Ok(py.None().into_py(py));
    }
    
    info!("[build_phasing_graph_rust] Step 4 completed: Graph has {} vertices and {} edges", 
              result.graph.node_count(), result.graph.edge_count());
    
    // Convert the result back to Python format
    info!("[build_phasing_graph_rust] Step 5: Converting result to Python format");
    export_graph_result_to_python(py, result, &read_pair_map)
}

// Note: BAM reading and allele depth mapping are now handled directly in Rust
// The unimplemented conversion functions have been removed since we process
// BAM files directly rather than converting from Python data structures

/// Export the Rust graph result to Python format compatible with graph-tool
fn export_graph_result_to_python(
    py: Python<'_>,
    result: PhasingGraphResult,
    read_pair_map: &ReadPairMap,
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
        
        // Get f32 weight directly from graph
        let weight = result.graph[edge];
        weights.push(weight);
    }
    
    // Convert to numpy arrays
    let edges_vec: Vec<Vec<u32>> = edges.iter().map(|&[a, b]| vec![a, b]).collect();
    let edge_array = PyArray2::from_vec2_bound(py, &edges_vec)?;
    let weight_array = PyArray1::from_vec_bound(py, weights);
    
    dict.set_item("edges", edge_array)?;
    dict.set_item("weights", weight_array)?;
    
    // 2. Export node names in the correct order (index i contains qname for node i)
    // With direct correspondence: qname_idx = NodeIndex.index()
    let mut node_names = vec![String::new(); result.graph.node_count()];
    
    // Get actual qnames from read_pair_map using direct correspondence
    for (qname_idx, read_pair) in &read_pair_map.readpair_dict {
        if *qname_idx < node_names.len() {
            node_names[*qname_idx] = read_pair.qname.clone();
        }
    }
    
    let node_list = PyList::new_bound(py, node_names);
    dict.set_item("vertex_names", node_list)?;

    // 2b. Export per-node read IDs as list of (id1, id2_or_None) tuples
    let node_ids_py = PyList::empty_bound(py);
    for (id1, id2_opt) in &result.node_read_ids {
        let items = vec![
            id1.clone().into_py(py),
            match id2_opt {
                Some(s) => s.clone().into_py(py),
                None => py.None().into_py(py),
            },
        ];
        let tuple = PyTuple::new_bound(py, &items);
        node_ids_py.append(tuple)?;
    }
    dict.set_item("node_read_ids", node_ids_py)?;
    
    // 3. Export weight matrix as numpy array - direct f32 export
    match result.weight_matrix() {
        Some(matrix) => {
            // Weight matrix is already f32, export directly
            let weight_matrix_array = matrix.to_pyarray_bound(py);
            dict.set_item("weight_matrix", weight_matrix_array)?;
        }
        None => {
            // Return None early instead of Error
            warn!("[export_graph_result_to_python] Weight matrix is None");
            return Ok(py.None().into_py(py));
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