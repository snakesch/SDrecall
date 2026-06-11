//! Haplotype phasing — pure-Rust port of `fp_control/phasing.py` + the Greedy-Clique-Expansion
//! core in `fp_control/gce_algorithm.py` (T3 of the SDrecall Rust migration).
//!
//! Given a phasing graph (weight matrix + edge set) it clusters read-pairs into haplotypes,
//! producing the same two maps the Python pipeline emits:
//!   * `vertex -> hap_id`        (`qname_hap_info`, keyed by vertex index)
//!   * `hap_id -> {qname}`       (`hap_qname_info`)
//!
//! The haplotype labels are arbitrary, so the validated invariant is the **partition** of
//! qnames (compared up to relabeling against the dumped Python output).

pub mod kernels;
pub mod gce;
pub mod phasing;

pub use phasing::{phase, qname_partition, PhasingInput, Round};
