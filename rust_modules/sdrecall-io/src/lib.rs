//! `sdrecall-io` — shared file I/O for the SDrecall Rust pipeline.
//!
//! This is the second T0 foundation crate (the first is `sdrecall-utils`). It
//! holds the shared BAM / BED / VCF / GraphML / TSV read+write operations. Most
//! work runs in process via `rust-htslib`, `bedrs`, `petgraph-graphml`, and
//! `quick-xml`; checked `samtools` and `bcftools` leaf operations remain where
//! they are the deliberate production implementation. It replaces the
//! scattered helpers in `src/utils.py`, `src/insert_size.py`, and
//! `fp_control/bam_ncls.py`.
//!
//! ## The one versatile BAM reader (DUP-1)
//!
//! [`bam::is_read_noisy`] is the single canonical noisy-read predicate (port of
//! `bam_ncls.py:115-193`) and [`bam::build_bam_index`] is the single reader that
//! collates by qname, drops noisy groups and builds the per-chrom point-query
//! index — consolidating the three pre-existing Rust copies
//! (`build_phasing_graph`, `haplotype_inspection`, `read_extraction`). No
//! parallel near-duplicate reader is introduced.
//!
//! ## Scope
//!
//! Implemented + tested: the geometry-typed BED ops ([`bed`]), sorted-VCF I/O
//! ([`vcf`]), consolidated BAM reader + filter ([`bam`]), GraphML I/O
//! ([`graphml`]), insert-size stats ([`insert_size`]), and TSV I/O ([`tsv`]).
//! Production BAM merging deliberately remains in the orchestrator's checked
//! `samtools` wrapper; the unused placeholder Rust merge API has been removed.

pub mod bam;
pub mod bed;
pub mod graphml;
pub mod insert_size;
pub mod tsv;
pub mod vcf;

// Re-export the most-used items at the crate root so callers write
// `sdrecall_io::{read_bed, build_bam_index, ...}`.
pub use bam::{build_bam_index, is_read_noisy, remap_masked_bam_to_genomic, BamIndex, NoisyFilter};
pub use bed::{
    complement, intersect, merge_bed_files, read_bed, slop, sort_merge_bed, subtract, write_bed,
};
pub use graphml::{read_graphml, write_graphml};
pub use insert_size::{get_insert_size_distribution, FragStats};
pub use tsv::{read_tsv, write_tsv};
pub use vcf::{concat_sort_vcfs, read_vcf, write_vcf};
