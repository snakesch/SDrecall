//! Shared geometry vocabulary: the genomic-interval value type, its strand-less
//! map key, and the two integer newtypes (`HapId`/`QnameIdx`) used across stages.
//!
//! These are the *value types* of the T0 interface contract (§4 of
//! `T0_foundation_crates.DESIGN.md`). They are deliberately cheap and
//! dependency-light: `GenomicInterval` owns its `chrom` `String` so an interval
//! can outlive any single BAM header and move freely across crate boundaries,
//! while `Strand`, `HapId` and `QnameIdx` are `Copy` (≤4 bytes, no alloc).
//!
//! The heavy interval set-ops engine (bedrs) lives in `sdrecall-io`; this module
//! only holds the value type plus the handful of pure predicates (`overlaps`,
//! `contains_point`, `len`, `region_key`) that every stage needs without pulling
//! in `rust-htslib`. The point-query index (`rust_lapper::Lapper`) likewise lives
//! in io — three roles, one value type, no overlap (the DUP-3 collapse).

use serde::{Deserialize, Serialize};

/// Haplotype label. `i32` matches the Python hap labels, including the `-1`
/// "no-haplotype" sentinel. `Copy` so it is passed by value everywhere.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug, Serialize, Deserialize)]
pub struct HapId(pub i32);

/// Dense `0..N` query-name index produced by the BAM reader. `u32` so it fits
/// the `Lapper<u32, QnameIdx>` value slot; `Copy` so it is passed by value.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug, Serialize, Deserialize)]
pub struct QnameIdx(pub u32);

/// Strand of a genomic feature. `Unknown` is the default (a BED record with no
/// strand column, or an interval where strand is irrelevant). A 1-byte `Copy`
/// enum — no allocation, free to pass by value.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default, Serialize, Deserialize)]
pub enum Strand {
    /// No strand information (BED with <6 columns, or strand-agnostic op).
    #[default]
    Unknown,
    /// `+` strand.
    Forward,
    /// `-` strand.
    Reverse,
}

/// A half-open genomic interval `[start, end)` on `chrom`, with optional strand.
///
/// `chrom` is an owned `String` (not a `&str` / tid): intervals outlive any single
/// BAM header and are moved across crates, so owning the name keeps them
/// self-contained. `start`/`end` are `i64` to match htslib's `pos` and the
/// half-open BED convention.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
pub struct GenomicInterval {
    /// Contig / chromosome name.
    pub chrom: String,
    /// 0-based inclusive start (half-open `[start, end)`).
    pub start: i64,
    /// 0-based exclusive end.
    pub end: i64,
    /// Strand (defaults to [`Strand::Unknown`]).
    pub strand: Strand,
}

impl GenomicInterval {
    /// Construct an interval with [`Strand::Unknown`] — the common BED3 case.
    pub fn new(chrom: impl Into<String>, start: i64, end: i64) -> Self {
        Self {
            chrom: chrom.into(),
            start,
            end,
            strand: Strand::Unknown,
        }
    }

    /// Construct an interval with an explicit strand (BED6).
    pub fn with_strand(chrom: impl Into<String>, start: i64, end: i64, strand: Strand) -> Self {
        Self {
            chrom: chrom.into(),
            start,
            end,
            strand,
        }
    }

    /// True iff the two intervals are on the same contig and their half-open
    /// spans intersect. Half-open semantics: `[10,20)` and `[20,30)` do **not**
    /// overlap (they only touch). Strand is ignored — overlap is a positional
    /// relation. Read-only borrow, no allocation.
    pub fn overlaps(&self, other: &Self) -> bool {
        self.chrom == other.chrom && self.start < other.end && other.start < self.end
    }

    /// True iff `pos` lies within `[start, end)` on `chrom`. The end is exclusive,
    /// matching BED half-open semantics.
    pub fn contains_point(&self, chrom: &str, pos: i64) -> bool {
        self.chrom == chrom && pos >= self.start && pos < self.end
    }

    /// Span length `end - start`. May be `0` for an empty interval and is never
    /// negative for a well-formed interval.
    pub fn len(&self) -> i64 {
        self.end - self.start
    }

    /// True iff the interval is empty (`start >= end`). Provided so clippy does
    /// not flag the inherent `len` method as missing an `is_empty` companion.
    pub fn is_empty(&self) -> bool {
        self.start >= self.end
    }

    /// The strand-insensitive `(chrom, start, end)` key — the Rust analog of the
    /// Python tuple used to key per-region dicts. Clones the contig name into an
    /// owned key (intervals are borrowed, keys are stored).
    pub fn region_key(&self) -> RegionKey {
        RegionKey {
            chrom: self.chrom.clone(),
            start: self.start,
            end: self.end,
        }
    }
}

/// Strand-insensitive hashable map key — the Rust analog of the `(chrom, start,
/// end)` Python tuple used to key per-region dictionaries.
///
/// A distinct type from [`GenomicInterval`] on purpose: it drops `strand` (so two
/// intervals differing only in strand collapse to one key), is smaller, and makes
/// "this is a lookup key, not a feature" intent-revealing at call sites.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
pub struct RegionKey {
    /// Contig / chromosome name.
    pub chrom: String,
    /// 0-based inclusive start.
    pub start: i64,
    /// 0-based exclusive end.
    pub end: i64,
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── overlaps ─────────────────────────────────────────────────────────────

    #[test]
    fn overlapping_intervals_on_same_chrom() {
        let a = GenomicInterval::new("chr1", 10, 30);
        let b = GenomicInterval::new("chr1", 20, 40);
        assert!(a.overlaps(&b));
        assert!(b.overlaps(&a)); // symmetric
    }

    #[test]
    fn touching_intervals_do_not_overlap_half_open() {
        // [10,20) and [20,30) share only the point 20, which is excluded from
        // the first interval — half-open, so NO overlap.
        let a = GenomicInterval::new("chr1", 10, 20);
        let b = GenomicInterval::new("chr1", 20, 30);
        assert!(!a.overlaps(&b));
        assert!(!b.overlaps(&a));
    }

    #[test]
    fn nested_interval_overlaps() {
        let outer = GenomicInterval::new("chr1", 10, 100);
        let inner = GenomicInterval::new("chr1", 40, 50);
        assert!(outer.overlaps(&inner));
        assert!(inner.overlaps(&outer));
    }

    #[test]
    fn identical_intervals_overlap() {
        let a = GenomicInterval::new("chr1", 10, 20);
        let b = GenomicInterval::new("chr1", 10, 20);
        assert!(a.overlaps(&b));
    }

    #[test]
    fn different_chrom_never_overlaps() {
        let a = GenomicInterval::new("chr1", 10, 30);
        let b = GenomicInterval::new("chr2", 10, 30);
        assert!(!a.overlaps(&b));
    }

    #[test]
    fn strand_is_ignored_for_overlap() {
        let a = GenomicInterval::with_strand("chr1", 10, 30, Strand::Forward);
        let b = GenomicInterval::with_strand("chr1", 20, 40, Strand::Reverse);
        assert!(a.overlaps(&b));
    }

    #[test]
    fn one_base_overlap_is_detected() {
        // [10,21) and [20,30) share base 20 → overlap.
        let a = GenomicInterval::new("chr1", 10, 21);
        let b = GenomicInterval::new("chr1", 20, 30);
        assert!(a.overlaps(&b));
    }

    // ── contains_point ───────────────────────────────────────────────────────

    #[test]
    fn contains_point_inclusive_start_exclusive_end() {
        let iv = GenomicInterval::new("chr1", 10, 20);
        assert!(iv.contains_point("chr1", 10)); // start is inclusive
        assert!(iv.contains_point("chr1", 19));
        assert!(!iv.contains_point("chr1", 20)); // end is exclusive
        assert!(!iv.contains_point("chr1", 9));
    }

    #[test]
    fn contains_point_wrong_chrom() {
        let iv = GenomicInterval::new("chr1", 10, 20);
        assert!(!iv.contains_point("chr2", 15));
    }

    // ── len / is_empty ───────────────────────────────────────────────────────

    #[test]
    fn len_is_end_minus_start() {
        assert_eq!(GenomicInterval::new("chr1", 10, 20).len(), 10);
        assert_eq!(GenomicInterval::new("chr1", 5, 5).len(), 0);
    }

    #[test]
    fn is_empty_when_start_ge_end() {
        assert!(GenomicInterval::new("chr1", 5, 5).is_empty());
        assert!(!GenomicInterval::new("chr1", 5, 6).is_empty());
    }

    // ── region_key ───────────────────────────────────────────────────────────

    #[test]
    fn region_key_drops_strand() {
        let fwd = GenomicInterval::with_strand("chr1", 10, 20, Strand::Forward);
        let rev = GenomicInterval::with_strand("chr1", 10, 20, Strand::Reverse);
        // Same position, different strand → same key (strand-insensitive).
        assert_eq!(fwd.region_key(), rev.region_key());
        assert_eq!(
            fwd.region_key(),
            RegionKey {
                chrom: "chr1".to_string(),
                start: 10,
                end: 20
            }
        );
    }

    #[test]
    fn region_key_distinguishes_position() {
        let a = GenomicInterval::new("chr1", 10, 20);
        let b = GenomicInterval::new("chr1", 10, 21);
        assert_ne!(a.region_key(), b.region_key());
    }

    // ── newtypes ─────────────────────────────────────────────────────────────

    #[test]
    fn hapid_sentinel_and_ordering() {
        let none = HapId(-1);
        let zero = HapId(0);
        assert!(none < zero);
        assert_eq!(none, HapId(-1));
    }

    #[test]
    fn strand_default_is_unknown() {
        assert_eq!(Strand::default(), Strand::Unknown);
    }

    #[test]
    fn qname_idx_is_copy() {
        let q = QnameIdx(7);
        let q2 = q; // Copy, not move
        assert_eq!(q, q2);
    }
}
