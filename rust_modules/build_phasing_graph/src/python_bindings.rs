use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use numpy::{PyArray1, PyArray2, ToPyArray};
use std::collections::HashMap;
use ahash::AHashMap;
use rustc_hash::FxHashMap;
use log::info;

use crate::structs::{ReadPairMap, AlleleDepthMap, HaplotypeConfig, PhasingGraphResult};
use crate::graph_builder::build_phasing_graph;

/// Python-exposed function for building the phasing graph
/// 
/// This function serves as the bridge between Python and Rust, handling all necessary
/// data conversions and returning results in a format compatible with graph-tool.
#[pyfunction]
#[pyo3(signature = (
    read_pair_map_data,
    allele_depth_data, 
    intrinsic_allele_depth_data,
    mean_read_length,
    edge_weight_cutoff=0.201
))]
pub fn build_phasing_graph_rust(
    py: Python<'_>,
    read_pair_map_data: &PyDict,
    allele_depth_data: &PyDict,
    intrinsic_allele_depth_data: Option<&PyDict>,
    mean_read_length: f32,
    edge_weight_cutoff: f32,
) -> PyResult<PyObject> {
    info!("Starting Rust graph building with mean_read_length={}", mean_read_length);
    
    // Convert Python data to Rust structures
    let read_pair_map = convert_read_pair_map(py, read_pair_map_data)?;
    let allele_depth_map = convert_allele_depth_map(allele_depth_data)?;
    let intrinsic_ad_map = intrinsic_allele_depth_data
        .map(|d| convert_allele_depth_map(d))
        .transpose()?;
    
    // Create configuration
    let config = HaplotypeConfig::new(mean_read_length);
    
    // Call the Rust graph building function
    let result = build_phasing_graph(
        &read_pair_map,
        &allele_depth_map,
        intrinsic_ad_map.as_ref(),
        &config,
    )?;
    
    // Convert the result back to Python format
    export_graph_result_to_python(py, result)
}

/// Convert Python read pair map data to Rust ReadPairMap
fn convert_read_pair_map(py: Python<'_>, data: &PyDict) -> PyResult<ReadPairMap> {
    // This is a placeholder - you'll need to implement based on your actual data format
    // The Python data structure needs to be converted to ReadPairMap
    todo!("Implement read pair map conversion from Python data")
}

/// Convert Python allele depth data to Rust AlleleDepthMap
fn convert_allele_depth_map(data: &PyDict) -> PyResult<AlleleDepthMap> {
    let mut map = AlleleDepthMap::new();
    
    // Iterate through chromosomes
    for (chrom_obj, positions_dict) in data.iter() {
        let chrom: String = chrom_obj.extract()?;
        let positions: &PyDict = positions_dict.downcast()?;
        
        // Iterate through positions
        for (pos_obj, depth_data) in positions.iter() {
            let pos: u32 = pos_obj.extract()?;
            
            // Convert depth data - expecting a dict with keys 0-5
            let depth_dict: &PyDict = depth_data.downcast()?;
            let mut depth_array = [0u16; 6];
            
            for i in 0..6 {
                if let Ok(Some(val)) = depth_dict.get_item(i as i32) {
                    depth_array[i] = val.extract()?;
                }
            }
            
            map.insert(&chrom, pos, depth_array);
        }
    }
    
    Ok(map)
}

/// Export the Rust graph result to Python format compatible with graph-tool
fn export_graph_result_to_python(
    py: Python<'_>,
    result: PhasingGraphResult,
) -> PyResult<PyObject> {
    let dict = PyDict::new(py);
    
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
    let edge_array = PyArray2::from_vec2(py, &edges)?;
    let weight_array = PyArray1::from_vec(py, weights);
    
    dict.set_item("edges", edge_array)?;
    dict.set_item("edge_weights", weight_array)?;
    
    // 2. Export node names in the correct order (index i contains qname for node i)
    // With direct correspondence, we need to get qnames from read_pair_map
    let mut node_names = vec![String::new(); result.graph.node_count()];
    // TODO: This would need the read_pair_map passed in to get qnames by index
    // For now, we'll create placeholder names based on index
    for i in 0..result.graph.node_count() {
        node_names[i] = format!("read_pair_{}", i);
    }
    let node_list = PyList::new(py, node_names);
    dict.set_item("node_names", node_list)?;
    
    // 3. Export weight matrix as numpy array - much simpler with ndarray!
    match result.weight_matrix() {
        Some(matrix) => {
            // Convert f16 matrix to f32 for Python compatibility
            let f32_matrix: ndarray::Array2<f32> = matrix.mapv(|x| x.to_f32());
            let weight_matrix_array = f32_matrix.to_pyarray(py);
            dict.set_item("weight_matrix", weight_matrix_array)?;
        }
        None => {
            // Create empty f32 matrix for Python
            let empty_matrix = ndarray::Array2::<f32>::zeros((0, 0));
            let weight_matrix_array = empty_matrix.to_pyarray(py);
            dict.set_item("weight_matrix", weight_matrix_array)?;
        }
    }
    
    // 4. Note: qname_to_node mapping removed due to direct correspondence optimization
    // Direct mapping: qname_idx = NodeIndex.index()
    
    // 5. Export read_hap_vectors
    let hap_vectors_dict = PyDict::new(py);
    for (qname, vec) in &result.read_hap_vectors {
        let py_array = PyArray1::from_vec(py, vec.clone());
        hap_vectors_dict.set_item(qname, py_array)?;
    }
    dict.set_item("read_hap_vectors", hap_vectors_dict)?;
    
    // 6. Export read_error_vectors
    let error_vectors_dict = PyDict::new(py);
    for (qname, vec) in &result.read_error_vectors {
        let py_array = PyArray1::from_vec(py, vec.clone());
        error_vectors_dict.set_item(qname, py_array)?;
    }
    dict.set_item("read_error_vectors", error_vectors_dict)?;
    
    // 7. Export read_ref_pos_dict
    let ref_pos_dict = PyDict::new(py);
    for (qname, (start, end)) in &result.read_ref_pos_dict {
        let tuple = PyTuple::new(py, &[start, end]);
        ref_pos_dict.set_item(qname, tuple)?;
    }
    dict.set_item("read_ref_pos_dict", ref_pos_dict)?;
    
    // 8. Note: low_qual_qnames tracking removed (handled in BAM reading stage)
    let low_qual_list = PyList::empty(py);
    dict.set_item("low_qual_qnames", low_qual_list)?;
    
    // 9. Add metadata
    dict.set_item("num_vertices", result.graph.node_count())?;
    dict.set_item("num_edges", result.graph.edge_count())?;
    
    Ok(dict.into())
}