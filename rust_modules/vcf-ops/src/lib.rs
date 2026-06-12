//! `vcf-ops` (SDrecall T5) — sorted-VCF priority merge + inhouse-common binomial.
//!
//! Ports `src/merge_variants_with_priority.py::merge_with_priority` and
//! `identify_common_vars.py::annotate_inhouse_common`. Both Python files carry a
//! **byte-for-byte identical** sorted two-pointer co-iteration skeleton; per the
//! #1 coding rule that skeleton lives **once** in [`coiterate`] and is
//! parameterized by a [`coiterate::LocusOp`] callback. The two GT-correction
//! ladders + the binomial test are the only per-path logic; the four Python GT
//! threshold ladders collapse into one [`gt_rules::should_force_hom`] evaluator
//! over `static` rule tables.
//!
//! All VCF I/O is in-process via `rust-htslib`. The single deliberately-external
//! piece is indel left-normalization ([`norm`]), a `bcftools norm` leaf subprocess
//! — an external-ALGORITHM choice (design §6), not a missing-dep fallback.
//!
//! ## DESIGN deviations (recorded)
//!
//! - **rayon-over-contigs deferred.** The design (§5) calls for `contigs.par_iter()`.
//!   `rust_htslib::bcf::Record` holds an `Rc<HeaderView>` and carries
//!   `unsafe impl Send`; sharing translated records (which share the writer's `Rc`
//!   header) across rayon threads would race the non-atomic `Rc` refcount. The
//!   per-contig pass is therefore **sequential**. This is correct and still
//!   delivers the headline win — the Python `multiprocessing` pickling boundary
//!   (the documented bottleneck) is eliminated because records stay in-memory; the
//!   contig-axis parallelism was a secondary speedup. A future safe version can
//!   give each worker its own header-cloned reader.
//! - **Read-all-then-group instead of indexed `fetch`.** The design maps the
//!   Python `vcf.fetch(region)` to `IndexedReader::fetch`. Since the inputs are
//!   already globally sorted by `sort_vcf`, we stream each VCF once and bucket by
//!   contig in memory — no `.tbi` dependency, and records are translated into the
//!   output header up front so the engine output is write-ready.

mod coiterate;
mod gt_rules;
mod header_merge;
mod inhouse_common;
mod norm;
mod priority_merge;
mod record_io;
mod vcf_group;

// The one engine + its primitives.
pub use coiterate::{coiterate_sorted_vcfs, CoiterSets, LocusKey, LocusOp};
pub use gt_rules::{
    should_force_hom, GtRule, SampleStats, MATCHED_PAIR_LADDER, QUERY_ONLY_LADDER, REF_ONLY_LADDER,
};

// Public entry points (files in → files out).
pub use inhouse_common::{
    annotate_inhouse_common, determine_common, is_inhouse_contig, InhouseParams,
};
pub use priority_merge::{is_main_contig, merge_with_priority, MergeParams};

// Re-export the normalization leaf so the bin / harnesses can pre-normalize.
pub use norm::{norm_dedup_sort, sort_vcf};
