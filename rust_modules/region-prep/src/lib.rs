//! `region-prep` (SDrecall T7) — per-RG fc/nfc realignment BEDs via strand-aware
//! relative-coordinate segment projection + padding + merge.
//!
//! Rust port of `realign_recall/prepare_masked_align_region.py`. The **only
//! genuine compute** is the strand-aware projection in
//! [`project::extract_and_pad_segments`]; the FC-target build, table read and
//! tag filtering map 1:1 onto `sdrecall-io` / `sdrecall-utils` units (the
//! no-duplicate-helpers rule). See `docs/analysis/tasks/T7_region_prep.md`.
//!
//! ## Units (one versatile unit per job)
//!
//! - [`all_region_bed`] — the 7-col `all_homo_regions` reader + per-subgroup row
//!   splitter (the single tag-parse / row-select unit).
//! - [`project`] — THE compute unit: pure, borrow-in / owned-out, no I/O.
//! - [`per_rg`] — orchestration: FC target once + rayon over subgroups.
//!
//! ## Verbatim-port parity (`T7_region_prep.md` §6 + two findings from the
//! HG002 RG0 differential)
//!
//! R1 OR predicate, R2 unstranded final merge, R3 mtime+coverage reuse, R5
//! basename-derived label, R6 `NFC_PAD`(600)/`FC_SLOP`(500) distinct, R7 required
//! `.fai`, R8 empty-Vec on disjoint FC. Two parity facts the differential surfaced
//! (the original §6 assumptions were wrong on both):
//! - **R2-bis** — pybedtools `.merge()` FUSES book-ended intervals (canonical
//!   `bedtools merge -d 0`), so the durable merges use [`per_rg`]'s local
//!   bookended merge, NOT `sdrecall-io::sort_merge_bed` (which is strict-overlap).
//! - **R9** — the Python FC bed handed to the projection is 3-column (no strand),
//!   so the opposite-strand (reverse-complement) branch runs for virtually every
//!   NFC row; [`per_rg`] feeds the FC interval as `Strand::Unknown` to match.

pub mod all_region_bed;
pub mod per_rg;
pub mod project;

pub use all_region_bed::{
    fc_rows, read_all_region_bed, split_subgroup, AllRegionRow, RgTag, COL_SENTINEL,
};
pub use per_rg::{
    fai_path_for, prepare_masked_align_region_per_rg, read_fai, RgSubgroupRecord, FC_SLOP,
};
pub use project::{extract_and_pad_segments, NfcInterval, NFC_PAD};
