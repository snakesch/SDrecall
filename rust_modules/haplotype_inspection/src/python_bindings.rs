/// Python bindings for haplotype inspection module
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PySet};
use std::collections::{HashMap, HashSet};
use crate::identify_misaligned_haps::inspect_haplotypes;

/// Main Python-facing function for haplotype inspection
/// Replaces Python's inspect_by_haplotypes function
#[pyfunction]
#[pyo3(signature = (
    bam_path,
    intrinsic_bam_path,
    hap_qname_info,
    qname_hap_info,
    qname_to_node,
    total_lowqual_qnames,
    compare_haplotype_meta_tab = "",
    mean_read_length = 148.0,
    recall_mq_cutoff = 20,
    basequal_median_cutoff = 10
))]
pub fn inspect_haplotypes_rust(
    _py: Python,
    bam_path: String,
    intrinsic_bam_path: String,
    hap_qname_info: &Bound<'_, PyDict>,
    qname_hap_info: &Bound<'_, PyDict>,
    qname_to_node: &Bound<'_, PyDict>,
    total_lowqual_qnames: &Bound<'_, PyAny>,
    compare_haplotype_meta_tab: &str,
    mean_read_length: f64,
    recall_mq_cutoff: u8,
    basequal_median_cutoff: u8,
) -> PyResult<(Vec<String>, Vec<String>)> {

    // Convert hap_qname_info: Dict[int, List[str]] → HashMap<i32, Vec<String>>
    let mut hap_qname_info_rs: HashMap<i32, Vec<String>> = HashMap::new();
    for (key, value) in hap_qname_info.iter() {
        let hap_id: i32 = key.extract()?;
        let qnames: Vec<String> = value.extract()?;
        hap_qname_info_rs.insert(hap_id, qnames);
    }

    // Convert qname_hap_info: Dict[int, int] (vertex_idx → hap_id) → HashMap<i32, i32>
    let mut qname_hap_info_rs: HashMap<i32, i32> = HashMap::new();
    for (key, value) in qname_hap_info.iter() {
        let vert_idx: i32 = key.extract()?;
        let hap_id: i32 = value.extract()?;
        qname_hap_info_rs.insert(vert_idx, hap_id);
    }

    // Convert qname_to_node: Dict[str, int] → HashMap<String, u32>
    let mut qname_to_node_rs: HashMap<String, u32> = HashMap::new();
    for (key, value) in qname_to_node.iter() {
        let qname: String = key.extract()?;
        let node_idx: u32 = value.extract()?;
        qname_to_node_rs.insert(qname, node_idx);
    }

    // Convert total_lowqual_qnames: Set[str] or List[str] → HashSet<String>
    let lowqual_qnames_rs: HashSet<String> = if let Ok(py_set) = total_lowqual_qnames.downcast::<PySet>() {
        py_set.iter().map(|item| item.extract::<String>()).collect::<PyResult<HashSet<String>>>()?
    } else if let Ok(py_list) = total_lowqual_qnames.downcast::<PyList>() {
        py_list.iter().map(|item| item.extract::<String>()).collect::<PyResult<HashSet<String>>>()?
    } else {
        // Try as any iterable
        let iter = total_lowqual_qnames.iter()?;
        iter.map(|item| item?.extract::<String>()).collect::<PyResult<HashSet<String>>>()?
    };

    // Call Rust inspection logic
    let (correct_qnames, mismap_qnames) = inspect_haplotypes(
        &bam_path,
        &intrinsic_bam_path,
        &hap_qname_info_rs,
        &qname_hap_info_rs,
        &qname_to_node_rs,
        &lowqual_qnames_rs,
        compare_haplotype_meta_tab,
        mean_read_length,
        recall_mq_cutoff,
        basequal_median_cutoff,
    ).map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
        format!("Rust haplotype inspection failed: {}", e)
    ))?;

    // Convert HashSet<String> → Vec<String> for Python
    let correct: Vec<String> = correct_qnames.into_iter().collect();
    let mismap: Vec<String> = mismap_qnames.into_iter().collect();

    Ok((correct, mismap))
}
