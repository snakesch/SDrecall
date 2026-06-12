//! `sdrecall-io` — in-process file I/O for the SDrecall Rust pipeline.
//!
//! This is the second T0 foundation crate (the first is `sdrecall-utils`). It
//! holds every BAM / BED / VCF / GraphML / TSV read+write the pipeline performs,
//! done **in process** via `rust-htslib`, `bedrs`, `petgraph-graphml` and
//! `quick-xml` — there is **no shelling out** to `samtools` / `bcftools` /
//! `bedtools`. It replaces the scattered subprocess helpers in `src/utils.py`,
//! `src/insert_size.py` and `fp_control/bam_ncls.py`.
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
//! ## Scope today (T0)
//!
//! Implemented + tested: the geometry-typed BED ops ([`bed`]), the sorted-VCF
//! merge ([`vcf`]), the consolidated BAM reader + filter ([`bam`]), GraphML I/O
//! ([`graphml`]), insert-size stats ([`insert_size`]) and TSV I/O ([`tsv`]).
//! Still stubbed with `TODO(T9)` and a typed [`sdrecall_utils::SdError`] (so the
//! API surface compiles, but with no current consumer): BAM merge
//! ([`bam::merge_bams`]) — its SQ-line reconcile is deferred to T9 (DESIGN §6).

pub mod bam;
pub mod bed;
pub mod graphml;
pub mod insert_size;
pub mod tsv;
pub mod vcf;

// Re-export the most-used items at the crate root so callers write
// `sdrecall_io::{read_bed, build_bam_index, ...}`.
pub use bam::{build_bam_index, is_read_noisy, merge_bams, BamIndex, NoisyFilter};
pub use bed::{
    complement, intersect, merge_bed_files, read_bed, slop, sort_merge_bed, write_bed,
};
pub use graphml::{read_graphml, write_graphml};
pub use insert_size::{get_insert_size_distribution, FragStats};
pub use tsv::{read_tsv, write_tsv};
pub use vcf::{concat_sort_vcfs, read_vcf, write_vcf};
