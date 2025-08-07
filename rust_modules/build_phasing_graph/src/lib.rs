pub mod structs;
pub mod bam_reading;
pub mod graph_builder;
pub mod haplotype_determination;
pub mod python_bindings;

use pyo3::prelude::*;

/// PyO3 module initialization
#[pymodule]
fn build_phasing_graph_rs(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Initialize pyo3-log bridge between Rust log crate and Python logging
    // This enables Rust log messages to be forwarded to Python's logging system
    // and respect Python's logging level configuration
    pyo3_log::init();
    
    // Add the main function
    m.add_function(wrap_pyfunction!(python_bindings::build_phasing_graph_rust, m)?)?;
    
    Ok(())
}
