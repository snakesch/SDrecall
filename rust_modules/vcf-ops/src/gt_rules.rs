//! The ONE versatile GT-correction unit.
//!
//! The Python carries **four** separate first-hit-wins threshold ladders that all
//! read the same `(ref_dp, alt_dp, gq, num_hps)` quadruple and decide whether to
//! force `GT=(1,1)`:
//!
//! 1. `merge_variants_with_priority.py::modify_gt_based_on_ad_gq` (L87-110) — the
//!    matched-pair ladder (reads from the REFERENCE record).
//! 2. `merge_with_priority` query-only finalizer (L598-606).
//! 3. `merge_with_priority` ref-only finalizer (L650).
//! 4. (inhouse-common has no GT correction.)
//!
//! Per the #1 coding rule these collapse into one evaluator: a `&[GtRule]` rule
//! table (one [`GtRule`] per Python `if`-clause) evaluated against a
//! [`SampleStats`] quadruple, returning `true` on the first matching rule. The
//! call-sites pass three different `static` rule tables — no duplicated ladders.
//!
//! ## Exact ratio semantics (load-bearing parity detail)
//!
//! The Python ladders are NOT uniform in their denominator:
//!
//! - `modify_gt_based_on_ad_gq` uses `rdp = max(1, rref + ralt)` and tests
//!   `ralt / rdp` (L76, L87…). So a `0/0` AD gives ratio `0/1 = 0`.
//! - The query-only finalizer uses `alt_dp / (alt_dp + ref_dp)` with **no**
//!   zero-guard (L598) — Python would raise `ZeroDivisionError` on `0/0`, but in
//!   practice the `(alt_dp + ref_dp) >= 5` clause is ANDed first and short-circuits
//!   it. We mirror that: when total is 0 the ratio is treated as 0 and the
//!   `min_dp >= 5` clause already fails, so the rule never fires.
//! - The ref-only finalizer uses `alt_dp / total_dp` with `total_dp = max(1, …)`
//!   (L647) and tests `>= 0.7`.
//!
//! All three reduce to "ratio = alt_dp / max(1, ref_dp + alt_dp)", which equals
//! every Python expression on the inputs that actually reach the comparison
//! (the `>= 5` gate makes the zero-guard difference unobservable). [`SampleStats`]
//! exposes that single ratio.

/// Per-sample stats extracted once from a record's FORMAT fields, with the exact
/// Python missing-value defaults already applied (AD→`[0,0]`, GQ→`0`, HPSUP→`.`).
/// Plain `Copy` scalars so the htslib buffer borrow ends before any mutation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SampleStats {
    pub ref_dp: i32,
    pub alt_dp: i32,
    pub gq: i32,
    /// `len(HPSUP[0].split(";"))`; `0` when HPSUP is absent on the query-only path,
    /// but the matched/ref ladders use a default of `1` (Python `".".split(";")`
    /// has length 1). The extractor records the as-seen count; the default for a
    /// missing field is the caller's responsibility (see [`crate::priority_merge`]).
    pub num_hps: usize,
}

impl SampleStats {
    /// `alt_dp / max(1, ref_dp + alt_dp)` — the single ratio every Python ladder
    /// reduces to on inputs that reach the comparison (see module docs).
    #[inline]
    pub fn alt_ratio(&self) -> f64 {
        let total = self.ref_dp + self.alt_dp;
        let denom = if total == 0 { 1 } else { total };
        self.alt_dp as f64 / denom as f64
    }

    /// `ref_dp + alt_dp` — the un-guarded depth used by the `>= min_dp` clauses.
    #[inline]
    pub fn total_dp(&self) -> i32 {
        self.ref_dp + self.alt_dp
    }
}

/// One first-hit-wins clause of a GT-correction ladder. `Copy`, cheap; the rule
/// tables are tiny `static` arrays. Each field is an optional constraint ANDed
/// with the others; `None`/`false`/`0.0` means "no constraint from this field".
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GtRule {
    /// `num_hps >= min_num_hps` (0 = no constraint).
    pub min_num_hps: usize,
    /// `alt_ratio() >= min_ratio` (0.0 = no constraint).
    pub min_ratio: f64,
    /// `total_dp() >= min_dp` (0 = no constraint).
    pub min_dp: i32,
    /// `gq < max_gq` when `Some` (the ref-only / GQ-rescue clauses).
    pub max_gq: Option<i32>,
    /// `alt_dp >= ref_dp` when `true` (the `num_hps>=4 & alt>=ref` clause).
    pub alt_ge_ref: bool,
}

impl GtRule {
    /// Does this single clause fire on the given stats? All present constraints
    /// must hold (logical AND), mirroring one Python `if`.
    #[inline]
    fn matches(&self, s: &SampleStats) -> bool {
        if self.min_num_hps > 0 && s.num_hps < self.min_num_hps {
            return false;
        }
        if self.min_ratio > 0.0 && s.alt_ratio() < self.min_ratio {
            return false;
        }
        if self.min_dp > 0 && s.total_dp() < self.min_dp {
            return false;
        }
        if let Some(mg) = self.max_gq {
            if s.gq >= mg {
                return false;
            }
        }
        if self.alt_ge_ref && s.alt_dp < s.ref_dp {
            return false;
        }
        true
    }
}

/// THE evaluator: returns `true` if any rule in the table fires (first-hit-wins
/// short-circuits, like the Python `continue`-after-set ladder). `&` everywhere —
/// read-only, no allocation.
#[inline]
pub fn should_force_hom(stats: &SampleStats, rules: &[GtRule]) -> bool {
    rules.iter().any(|r| r.matches(stats))
}

// ── The three Python ladders as `static` rule tables ────────────────────────

/// `modify_gt_based_on_ad_gq` (L87-110) — the matched-pair ladder, read from the
/// REFERENCE record. Note: `ralt/rdp` uses `rdp = max(1, rref+ralt)`, captured by
/// [`SampleStats::alt_ratio`]. The final clause additionally needs `rdp > 5`
/// (strict), encoded as `min_dp: 6` (`>5` ⇔ `>=6` on integers).
pub static MATCHED_PAIR_LADDER: &[GtRule] = &[
    // ralt/rdp >= 0.9
    GtRule {
        min_num_hps: 0,
        min_ratio: 0.9,
        min_dp: 0,
        max_gq: None,
        alt_ge_ref: false,
    },
    // num_hps >= 2 and ralt/rdp >= 0.33
    GtRule {
        min_num_hps: 2,
        min_ratio: 0.33,
        min_dp: 0,
        max_gq: None,
        alt_ge_ref: false,
    },
    // num_hps >= 3 and ralt/rdp >= 0.30
    GtRule {
        min_num_hps: 3,
        min_ratio: 0.30,
        min_dp: 0,
        max_gq: None,
        alt_ge_ref: false,
    },
    // num_hps >= 4 and ralt/rdp >= 0.25
    GtRule {
        min_num_hps: 4,
        min_ratio: 0.25,
        min_dp: 0,
        max_gq: None,
        alt_ge_ref: false,
    },
    // rgq < 5 and ralt/rdp >= 0.5 and rdp > 5  (rdp>5 ⇔ total_dp>=6)
    GtRule {
        min_num_hps: 0,
        min_ratio: 0.5,
        min_dp: 6,
        max_gq: Some(5),
        alt_ge_ref: false,
    },
];

/// `merge_with_priority` query-only finalizer (L598-606). All three clauses AND
/// `(alt_dp + ref_dp) >= 5`.
pub static QUERY_ONLY_LADDER: &[GtRule] = &[
    // num_hps >= 2 and ratio >= 0.55 and dp >= 5
    GtRule {
        min_num_hps: 2,
        min_ratio: 0.55,
        min_dp: 5,
        max_gq: None,
        alt_ge_ref: false,
    },
    // num_hps >= 4 and alt_dp >= ref_dp and dp >= 5
    GtRule {
        min_num_hps: 4,
        min_ratio: 0.0,
        min_dp: 5,
        max_gq: None,
        alt_ge_ref: true,
    },
    // ratio >= 0.9 and dp >= 5
    GtRule {
        min_num_hps: 0,
        min_ratio: 0.9,
        min_dp: 5,
        max_gq: None,
        alt_ge_ref: false,
    },
];

/// `merge_with_priority` ref-only finalizer (L650): `gq < 5 and alt/total >= 0.7`.
/// `total_dp = max(1, …)` is captured by [`SampleStats::alt_ratio`].
pub static REF_ONLY_LADDER: &[GtRule] = &[GtRule {
    min_num_hps: 0,
    min_ratio: 0.7,
    min_dp: 0,
    max_gq: Some(5),
    alt_ge_ref: false,
}];

#[cfg(test)]
mod tests {
    use super::*;

    fn stats(ref_dp: i32, alt_dp: i32, gq: i32, num_hps: usize) -> SampleStats {
        SampleStats {
            ref_dp,
            alt_dp,
            gq,
            num_hps,
        }
    }

    #[test]
    fn alt_ratio_zero_guard() {
        // 0/0 → 0/1 = 0.0 (Python rdp = max(1, 0)).
        assert_eq!(stats(0, 0, 0, 1).alt_ratio(), 0.0);
        // 1/4 = 0.25
        assert_eq!(stats(3, 1, 0, 1).alt_ratio(), 0.25);
        // 9/10 = 0.9
        assert_eq!(stats(1, 9, 0, 1).alt_ratio(), 0.9);
    }

    // ── matched-pair ladder boundary points (design §7) ─────────────────────

    #[test]
    fn matched_ratio_090_boundary() {
        // ralt/rdp = 0.89 (89/100) → no rule fires; 0.90 (90/100) → first rule.
        assert!(!should_force_hom(
            &stats(11, 89, 30, 1),
            MATCHED_PAIR_LADDER
        ));
        assert!(should_force_hom(&stats(10, 90, 30, 1), MATCHED_PAIR_LADDER));
    }

    #[test]
    fn matched_numhps2_033_boundary() {
        // num_hps=2: ratio 0.32 → no; 0.33 → yes.
        assert!(!should_force_hom(
            &stats(68, 32, 30, 2),
            MATCHED_PAIR_LADDER
        ));
        assert!(should_force_hom(&stats(67, 33, 30, 2), MATCHED_PAIR_LADDER));
        // num_hps=1 at ratio 0.33 → no (needs >=2, and 0.33 < 0.9).
        assert!(!should_force_hom(
            &stats(67, 33, 30, 1),
            MATCHED_PAIR_LADDER
        ));
    }

    #[test]
    fn matched_numhps3_030_boundary() {
        assert!(!should_force_hom(
            &stats(71, 29, 30, 3),
            MATCHED_PAIR_LADDER
        ));
        assert!(should_force_hom(&stats(70, 30, 30, 3), MATCHED_PAIR_LADDER));
    }

    #[test]
    fn matched_numhps4_025_boundary() {
        assert!(!should_force_hom(
            &stats(76, 24, 30, 4),
            MATCHED_PAIR_LADDER
        ));
        assert!(should_force_hom(&stats(75, 25, 30, 4), MATCHED_PAIR_LADDER));
    }

    #[test]
    fn matched_gq_rescue_boundary() {
        // rgq<5 & ratio>=0.5 & rdp>5(>=6). num_hps=1 so earlier rules can't fire
        // (ratio 0.5 < 0.9). gq=4 fires, gq=5 does not.
        assert!(should_force_hom(&stats(3, 3, 4, 1), MATCHED_PAIR_LADDER)); // 3/6=0.5, dp6, gq4
        assert!(!should_force_hom(&stats(3, 3, 5, 1), MATCHED_PAIR_LADDER)); // gq=5 not <5
                                                                             // dp must be >5: 2/4=0.5 but dp=4 → no.
        assert!(!should_force_hom(&stats(2, 2, 4, 1), MATCHED_PAIR_LADDER));
        // ratio 0.5 at dp exactly 6 (alt 3 ref 3) already covered; dp 6 with ratio
        // 0.49 (alt slightly less) cannot occur on ints — use 5/12≈0.416 → no.
        assert!(!should_force_hom(&stats(7, 5, 4, 1), MATCHED_PAIR_LADDER));
    }

    // ── query-only ladder ───────────────────────────────────────────────────

    #[test]
    fn query_numhps2_055_dp5_boundary() {
        // num_hps=2, ratio 0.54 (no) vs 0.55 (yes), dp must be >=5.
        // ratio 0.54: pick alt/ref so total>=5 and ratio just below 0.55.
        // 0.5 (5/10) < 0.55 → no.
        assert!(!should_force_hom(&stats(5, 5, 30, 2), QUERY_ONLY_LADDER));
        // 0.6 (6/10) >= 0.55 → yes.
        assert!(should_force_hom(&stats(4, 6, 30, 2), QUERY_ONLY_LADDER));
        // dp gate: ratio 1.0 but dp 4 (0,4) and num_hps=2 → ratio>=0.9 rule needs dp>=5 too → no.
        assert!(!should_force_hom(&stats(0, 4, 30, 2), QUERY_ONLY_LADDER));
    }

    #[test]
    fn query_numhps4_alt_ge_ref_boundary() {
        // num_hps=4, alt<ref → no; alt>=ref & dp>=5 → yes.
        assert!(!should_force_hom(&stats(4, 3, 30, 4), QUERY_ONLY_LADDER)); // alt<ref, ratio 3/7<0.55
        assert!(should_force_hom(&stats(3, 4, 30, 4), QUERY_ONLY_LADDER)); // alt>=ref, dp7
                                                                           // alt>=ref but dp 4 (<5) → no.
        assert!(!should_force_hom(&stats(2, 2, 30, 4), QUERY_ONLY_LADDER));
    }

    #[test]
    fn query_ratio_090_dp5_boundary() {
        // num_hps low so only the ratio>=0.9 rule applies.
        assert!(!should_force_hom(&stats(2, 7, 30, 1), QUERY_ONLY_LADDER)); // 7/9≈0.78 <0.9
        assert!(should_force_hom(&stats(1, 9, 30, 1), QUERY_ONLY_LADDER)); // 9/10=0.9, dp10
                                                                           // ratio 1.0 but dp 4 → no.
        assert!(!should_force_hom(&stats(0, 4, 30, 1), QUERY_ONLY_LADDER));
    }

    // ── ref-only ladder ──────────────────────────────────────────────────────

    #[test]
    fn ref_only_gq_alt_total_boundary() {
        // gq<5 & alt/total>=0.7. alt/total 0.69 → no; 0.70 → yes.
        // 7/10 = 0.7 → yes; 69/100 → no.
        assert!(should_force_hom(&stats(3, 7, 4, 1), REF_ONLY_LADDER));
        assert!(!should_force_hom(&stats(31, 69, 4, 1), REF_ONLY_LADDER));
        // gq=5 (not <5) → no even at ratio 1.0.
        assert!(!should_force_hom(&stats(0, 10, 5, 1), REF_ONLY_LADDER));
    }
}
