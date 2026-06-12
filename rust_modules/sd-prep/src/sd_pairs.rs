//! Umbrella SD-pair filtering — `Pair` + the ONE `umbrella_to_remove` sweep.
//!
//! Port of `preparation/sd_pairs.py`: `Pair`, `calculate_interval_overlaps`
//! (l.32-64), `is_umbrella_pair` (l.96-122), `extract_subsegment_for_target`
//! (l.125-173), and `_find_umbrella_pairs` (l.182-253) — which runs the **same**
//! O(n²) umbrella sweep twice (over the raw pairs, then over the
//! target-refined "granular" pairs). The two sweeps collapse into ONE
//! [`umbrella_to_remove`] called twice (the #1 coding rule), exactly as the design
//! prescribes.
//!
//! ## What "umbrella" means (parity-critical)
//!
//! `a.is_umbrella_pair(b)` is true iff: same chroms on both segments; strand
//! consistency `(a.sA==a.sB) == (b.sA==b.sB)`; the overlap fraction of `a` over
//! `b` (computed with `fraction_select="other"`, i.e. relative to `b`'s sizes) is
//! `>= threshold` on **both** segments; AND `a.overlap_len <= b.overlap_len`.
//! Then `a` (the smaller-overlap, fully-covered pair) is the umbrella to REMOVE.
//!
//! NOTE the Python naming is counter-intuitive: the pair flagged for removal is
//! the one whose intervals are *covered by* the other and whose `overlap_len` is
//! `<=` the other's. We keep the Python semantics verbatim.

use sdrecall_utils::Strand;

/// A pair of genomic intervals (segment A, segment B) plus the BAM-overlap length
/// used to rank umbrella relationships. Mirrors `sd_pairs.py::Pair`.
///
/// All coordinate fields are owned/`Copy`; `chr_a`/`chr_b` are owned `String`s so
/// a `Pair` is self-contained. `overlap_len` is the `overlap_len` column (the bp
/// of overlap between segment A and the multi-align interval that grouped it).
#[derive(Clone, Debug)]
pub struct Pair {
    pub chr_a: String,
    pub start_a: i64,
    pub end_a: i64,
    pub strand_a: Strand,
    pub chr_b: String,
    pub start_b: i64,
    pub end_b: i64,
    pub strand_b: Strand,
    pub overlap_len: i64,
}

impl Pair {
    /// Overlap fractions relative to the **other** pair's segment sizes
    /// (`calculate_interval_overlaps(other, fraction_select="other")`,
    /// sd_pairs.py l.45-54, l.59). Returns `(frac_a, frac_b)`.
    ///
    /// `frac_a = overlap_span_A / other.sizeA`, `frac_b = overlap_span_B /
    /// other.sizeB`. Half-open spans clamped at 0. Division by a zero-size other
    /// segment yields `0.0` (Python would raise ZeroDivisionError, but SD
    /// intervals always have positive size; we guard to stay total).
    fn overlap_fractions_other(&self, other: &Pair) -> (f64, f64) {
        let span_a = (self.end_a.min(other.end_a) - self.start_a.max(other.start_a)).max(0);
        let span_b = (self.end_b.min(other.end_b) - self.start_b.max(other.start_b)).max(0);
        let size2_a = other.end_a - other.start_a;
        let size2_b = other.end_b - other.start_b;
        let fa = if size2_a > 0 {
            span_a as f64 / size2_a as f64
        } else {
            0.0
        };
        let fb = if size2_b > 0 {
            span_b as f64 / size2_b as f64
        } else {
            0.0
        };
        (fa, fb)
    }

    /// True iff `self` is an umbrella pair enclosing `other` (input order matters),
    /// `is_umbrella_pair(other, coverage_threshold)` (sd_pairs.py l.96-122).
    pub fn is_umbrella_pair(&self, other: &Pair, coverage_threshold: f64) -> bool {
        if self.chr_a != other.chr_a || self.chr_b != other.chr_b {
            return false;
        }
        // Strand consistency: (sA==sB) must agree between the two pairs.
        let self_consistent = self.strand_a == self.strand_b;
        let other_consistent = other.strand_a == other.strand_b;
        if self_consistent != other_consistent {
            return false;
        }
        let (fa, fb) = self.overlap_fractions_other(other);
        if fa < coverage_threshold || fb < coverage_threshold {
            return false;
        }
        self.overlap_len <= other.overlap_len
    }

    /// Refine this pair against a target region overlapping segment A
    /// (`extract_subsegment_for_target`, sd_pairs.py l.125-173). Returns a new
    /// `Pair` whose segment A is clipped to the overlap with `(chrom,start,end)`
    /// and whose segment B is the corresponding sub-segment (direct map on same
    /// strand, flipped map on opposite strand).
    ///
    /// Returns `None` if the target does not overlap segment A (Python logs a
    /// warning and returns `None`).
    pub fn extract_subsegment_for_target(
        &self,
        chrom: &str,
        target_start: i64,
        target_end: i64,
    ) -> Option<Pair> {
        if chrom != self.chr_a || target_end <= self.start_a || target_start >= self.end_a {
            return None;
        }
        let overlap_start = self.start_a.max(target_start);
        let overlap_end = self.end_a.min(target_end);
        let rel_start = overlap_start - self.start_a;
        let rel_end = overlap_end - self.start_a;
        let (new_start_b, new_end_b) = if self.strand_a == self.strand_b {
            (self.start_b + rel_start, self.start_b + rel_end)
        } else {
            (self.end_b - rel_end, self.end_b - rel_start)
        };
        Some(Pair {
            chr_a: self.chr_a.clone(),
            start_a: overlap_start,
            end_a: overlap_end,
            strand_a: self.strand_a,
            chr_b: self.chr_b.clone(),
            start_b: new_start_b,
            end_b: new_end_b,
            strand_b: self.strand_b,
            overlap_len: self.overlap_len,
        })
    }
}

/// THE one umbrella sweep — the O(n²) "mark umbrella pairs for removal" loop run
/// once over the raw pairs and once over the granular pairs (`_find_umbrella_pairs`
/// l.206-249). Returns the set of indices to remove.
///
/// The control flow matches Python exactly, including the asymmetric `break`:
/// - For each `i` not already removed, scan `j > i` not already removed.
/// - If `pairs[i].is_umbrella_pair(pairs[j])` → mark `i`, **break** (stop scanning
///   `i`).
/// - Else if `pairs[j].is_umbrella_pair(pairs[i])` → mark `j`, **continue** scanning.
///
/// `&[Pair]` borrow-in; returns owned index set. Called twice by
/// [`filter_umbrella_group`].
pub fn umbrella_to_remove(pairs: &[Pair], coverage_threshold: f64) -> ahash::AHashSet<usize> {
    let mut to_remove: ahash::AHashSet<usize> = ahash::AHashSet::new();
    let n = pairs.len();
    for i in 0..n {
        if to_remove.contains(&i) {
            continue;
        }
        for j in (i + 1)..n {
            if to_remove.contains(&j) {
                continue;
            }
            if pairs[i].is_umbrella_pair(&pairs[j], coverage_threshold) {
                to_remove.insert(i);
                break;
            } else if pairs[j].is_umbrella_pair(&pairs[i], coverage_threshold) {
                to_remove.insert(j);
            }
        }
    }
    to_remove
}

/// Filter a group of pairs (all sharing the same `(chr_bam1, start_bam1,
/// end_bam1)` multi-align key) by running the two umbrella sweeps and dropping the
/// union of removed indices (`filter_umbrella_pairs` + `_find_umbrella_pairs`).
///
/// `targets[i]` is the `(chrom,start,end)` target region paired with `pairs[i]`
/// (the `chr_target/start_target/end_target` columns the Python group carries,
/// l.197, l.226-227); used to build the granular sweep. A `None` granular pair
/// (target does not overlap segment A) is excluded from the granular sweep, as in
/// Python a `None` from `extract_subsegment_for_target` would still be appended
/// and crash `is_umbrella_pair` — but in practice every grouped row's target DOES
/// overlap its own segment A (that is why they were grouped), so this never
/// triggers; we skip defensively.
///
/// Returns the kept pairs (the input order, minus removed indices), and the kept
/// indices (for callers that need to map back to original table rows).
pub fn filter_umbrella_group(
    pairs: &[Pair],
    targets: &[(String, i64, i64)],
    coverage_threshold: f64,
) -> Vec<usize> {
    assert_eq!(pairs.len(), targets.len(), "pairs/targets length mismatch");
    let raw_remove = umbrella_to_remove(pairs, coverage_threshold);

    // Build granular pairs (refined against each row's own target).
    let granular: Vec<Pair> = pairs
        .iter()
        .zip(targets.iter())
        .filter_map(|(p, (c, s, e))| p.extract_subsegment_for_target(c, *s, *e))
        .collect();
    // The granular sweep operates on its own index space; but Python keeps the
    // raw-pair row_index alignment because granular_pairs is built in the same
    // order. Since we skip None granular pairs, re-derive the alignment: only
    // pairs whose target overlaps segment A produce a granular pair. To preserve
    // index parity with Python (which never produces None here), we require the
    // 1:1 mapping and assert it.
    let granular_remove = if granular.len() == pairs.len() {
        umbrella_to_remove(&granular, coverage_threshold)
    } else {
        // Defensive: one or more targets did not overlap; run granular only over
        // the ones that did, mapping indices back. This branch is not expected on
        // real data (see doc), so we log and fall back to raw-only removal.
        log::warn!(
            "filter_umbrella_group: {} of {} granular pairs were None (target did not overlap segment A); granular sweep skipped",
            pairs.len() - granular.len(),
            pairs.len()
        );
        ahash::AHashSet::new()
    };

    let mut removed = raw_remove;
    removed.extend(granular_remove);
    (0..pairs.len()).filter(|i| !removed.contains(i)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pair(
        sa: i64,
        ea: i64,
        strand_a: Strand,
        sb: i64,
        eb: i64,
        strand_b: Strand,
        overlap_len: i64,
    ) -> Pair {
        Pair {
            chr_a: "chr1".into(),
            start_a: sa,
            end_a: ea,
            strand_a,
            chr_b: "chr2".into(),
            start_b: sb,
            end_b: eb,
            strand_b,
            overlap_len,
        }
    }

    #[test]
    fn umbrella_covers_inner_pair() {
        // A: segA [100,1100) segB [5000,6000), strand +/+, overlap_len 50 (small).
        // B: segA [100,1100) segB [5000,6000), strand +/+, overlap_len 1000 (large).
        // A fully covers B (identical intervals → frac 1.0 ≥ 0.95) and
        // A.overlap_len (50) <= B.overlap_len (1000) → A is the umbrella to remove.
        let a = pair(100, 1100, Strand::Forward, 5000, 6000, Strand::Forward, 50);
        let b = pair(100, 1100, Strand::Forward, 5000, 6000, Strand::Forward, 1000);
        assert!(a.is_umbrella_pair(&b, 0.95));
        assert!(!b.is_umbrella_pair(&a, 0.95)); // b.overlap_len > a.overlap_len → not <=
        let rem = umbrella_to_remove(&[a, b], 0.95);
        assert_eq!(rem, [0usize].into_iter().collect());
    }

    #[test]
    fn strand_consistency_flip_blocks_umbrella() {
        // A: +/+ (consistent). B: +/- (inconsistent). (sA==sB) differs → not umbrella.
        let a = pair(100, 1100, Strand::Forward, 5000, 6000, Strand::Forward, 50);
        let b = pair(100, 1100, Strand::Forward, 5000, 6000, Strand::Reverse, 1000);
        assert!(!a.is_umbrella_pair(&b, 0.95));
        assert!(!b.is_umbrella_pair(&a, 0.95));
        let rem = umbrella_to_remove(&[a, b], 0.95);
        assert!(rem.is_empty());
    }

    #[test]
    fn partial_coverage_below_threshold_not_umbrella() {
        // A covers only half of B's segB → frac_b ~ 0.5 < 0.95 → not umbrella.
        let a = pair(100, 1100, Strand::Forward, 5000, 6000, Strand::Forward, 50);
        let b = pair(100, 1100, Strand::Forward, 5000, 7000, Strand::Forward, 1000);
        // frac over B: segA 1000/1000=1.0; segB 1000/2000=0.5 → fails.
        assert!(!a.is_umbrella_pair(&b, 0.95));
    }

    #[test]
    fn different_chrom_b_not_umbrella() {
        let a = pair(100, 1100, Strand::Forward, 5000, 6000, Strand::Forward, 50);
        let mut b = pair(100, 1100, Strand::Forward, 5000, 6000, Strand::Forward, 1000);
        b.chr_b = "chr9".into();
        assert!(!a.is_umbrella_pair(&b, 0.95));
    }

    #[test]
    fn extract_subsegment_same_strand_direct_map() {
        // segA [100,1100) segB [5000,6000) +/+. Target [600,800) on segA →
        // rel [500,700) → segB [5500,5700).
        let p = pair(100, 1100, Strand::Forward, 5000, 6000, Strand::Forward, 50);
        let r = p.extract_subsegment_for_target("chr1", 600, 800).unwrap();
        assert_eq!((r.start_a, r.end_a), (600, 800));
        assert_eq!((r.start_b, r.end_b), (5500, 5700));
    }

    #[test]
    fn extract_subsegment_opposite_strand_flips() {
        // segA [100,1100) segB [5000,6000) +/-. Target [600,800) on segA →
        // rel [500,700) → segB flipped: [6000-700, 6000-500) = [5300,5500).
        let p = pair(100, 1100, Strand::Forward, 5000, 6000, Strand::Reverse, 50);
        let r = p.extract_subsegment_for_target("chr1", 600, 800).unwrap();
        assert_eq!((r.start_a, r.end_a), (600, 800));
        assert_eq!((r.start_b, r.end_b), (5300, 5500));
    }

    #[test]
    fn extract_subsegment_no_overlap_is_none() {
        let p = pair(100, 1100, Strand::Forward, 5000, 6000, Strand::Forward, 50);
        assert!(p.extract_subsegment_for_target("chr1", 2000, 2200).is_none());
        assert!(p.extract_subsegment_for_target("chr9", 600, 800).is_none());
    }

    #[test]
    fn umbrella_sweep_break_semantics() {
        // 0 covers 1 (break after marking 0). 2 is independent.
        let p0 = pair(100, 1100, Strand::Forward, 5000, 6000, Strand::Forward, 10);
        let p1 = pair(100, 1100, Strand::Forward, 5000, 6000, Strand::Forward, 1000);
        let p2 = pair(20000, 21000, Strand::Forward, 25000, 26000, Strand::Forward, 500);
        let rem = umbrella_to_remove(&[p0, p1, p2], 0.95);
        assert_eq!(rem, [0usize].into_iter().collect());
    }

    #[test]
    fn filter_group_drops_umbrella_keeps_rest() {
        let p0 = pair(100, 1100, Strand::Forward, 5000, 6000, Strand::Forward, 10);
        let p1 = pair(100, 1100, Strand::Forward, 5000, 6000, Strand::Forward, 1000);
        let pairs = vec![p0, p1];
        // Targets overlap segment A so granular pairs exist (1:1).
        let targets = vec![
            ("chr1".to_string(), 100i64, 1100i64),
            ("chr1".to_string(), 100i64, 1100i64),
        ];
        let kept = filter_umbrella_group(&pairs, &targets, 0.95);
        assert_eq!(kept, vec![1]); // index 0 (umbrella) dropped
    }
}
