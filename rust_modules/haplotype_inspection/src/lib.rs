pub mod structs;
pub mod bam_lappers;
pub mod pairwise_read_inspection;
pub mod identify_misaligned_haps;
pub mod bilc_solver;
pub mod python_bindings;

use pyo3::prelude::*;

/// PyO3 module initialization
#[pymodule]
fn haplotype_inspection(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Initialize pyo3-log bridge between Rust log crate and Python logging
    // This enables Rust log messages to be forwarded to Python's logging system
    // and respect Python's logging level configuration
    pyo3_log::init();

    // Add the main function
    m.add_function(wrap_pyfunction!(python_bindings::inspect_haplotypes_rust, m)?)?;

    Ok(())
}
