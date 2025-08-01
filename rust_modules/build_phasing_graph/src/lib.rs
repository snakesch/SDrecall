pub mod structs;
pub mod bam_reading;
pub mod graph_builder;
pub mod haplotype_determination;
pub mod python_bindings;

use pyo3::prelude::*;

/// PyO3 module initialization
#[pymodule]
fn build_phasing_graph_rs(_py: Python, m: &PyModule) -> PyResult<()> {
    // Initialize logging
    pyo3_log::init();
    
    // Add the main function
    m.add_function(wrap_pyfunction!(python_bindings::build_phasing_graph_rust, m)?)?;
    
    Ok(())
}
