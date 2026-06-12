//! THE one sorted-VCF two-pointer co-iteration engine.
//!
//! Both `merge_variants_with_priority.py` and `identify_common_vars.py` carry a
//! **byte-for-byte identical** `compare_variant_positions` /
//! `deal_with_same_loc_variants` / `process_region` skeleton (the only diff is the
//! per-locus callback and what happens to the three result sets). Per the #1
//! coding rule that skeleton is migrated **once** here and parameterized by a
//! [`LocusOp`] trait.
//!
//! ## Faithful port of the Python control flow
//!
//! `process_region` (merge L246-394 / inhouse L222-370) is a single-pass two-
//! pointer merge over two coordinate-sorted, single-contig record streams:
//!
//! - Order is by **0-based `start` then `stop`** (`compare_variant_positions`,
//!   merge L148-165), NOT by alleles. "Same location" (cmp==0) is a coarser bucket
//!   than "same variant" (allele-level `__eq__`, L34-39).
//! - When two records share a *location*, `deal_with_same_loc_variants`
//!   (L168-238) drains **both** streams of all co-located records, then does an
//!   O(n·m) allele-equality match: matched pairs → callback → `matched`; unmatched
//!   query → `query_only`; unmatched ref → `ref_only`.
//! - The Python `set` semantics dedup on the wrapper's `__hash__`/`__eq__`
//!   (chrom,pos,alleles). We mirror that with a [`LocusKey`] de-dup so a record
//!   added to the same set twice (the Python `set.add`) appears once.
//!
//! Because we operate on already-materialized, header-translated `Vec<bcf::Record>`
//! (one contig), the engine is index-free and the records it returns are ready to
//! write — replacing the Python `.pickable()` cross-process tuple round-trip.

use rust_htslib::bcf;
use std::cmp::Ordering;
use std::collections::HashSet;

/// The three result sets of one contig's co-iteration. Owned `Vec<bcf::Record>`
/// because the records are mutated per-set and later written (they already live
/// in the output header). Replaces the Python `(frozenset, frozenset, frozenset)`.
pub struct CoiterSets {
    pub matched: Vec<bcf::Record>,
    pub query_only: Vec<bcf::Record>,
    pub ref_only: Vec<bcf::Record>,
}

/// Allele-level identity key — the Python `VariantRecordWrapper.__hash__`/`__eq__`
/// on `(chrom, pos, alleles)` (merge L31-39). Used both for the same-location
/// allele match and for the set-dedup the Python `set` provides.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct LocusKey {
    pub rid: i64,
    pub pos: i64,
    pub alleles: Vec<Box<[u8]>>,
}

impl LocusKey {
    /// Build from a record. `rid` is the record's contig index; within one
    /// engine call all records share the output header so `rid` is comparable.
    pub fn of(rec: &bcf::Record) -> Self {
        LocusKey {
            rid: rec.rid().map(|r| r as i64).unwrap_or(-1),
            pos: rec.pos(),
            alleles: rec.alleles().iter().map(|a| (*a).into()).collect(),
        }
    }
}

/// `compare_variant_positions` (merge L148-165): total order on `(start, stop)`.
/// Returns `Ordering::Less` if `a` is strictly before `b` (Python `1`), `Greater`
/// if after (Python `2`), `Equal` if same location (Python `0`). Different contig
/// (Python `-1`) is impossible here — the engine is called per-contig — and we
/// debug-assert it. `start = pos()` (0-based), `stop = end()`.
fn position_cmp(a: &bcf::Record, b: &bcf::Record) -> Ordering {
    debug_assert_eq!(a.rid(), b.rid(), "coiterate is called per-contig");
    a.pos().cmp(&b.pos()).then_with(|| a.end().cmp(&b.end()))
}

/// The per-locus callback. Called once per (query, ref) pair sharing BOTH location
/// AND alleles (the Python `merging_func` / `process_func`). May mutate either
/// record in place; returns the record to keep in the `matched` set, or `None` to
/// drop the matched pair (no current path drops, but the option keeps the engine
/// generic).
pub trait LocusOp {
    fn on_match(&self, query: &mut bcf::Record, refr: &mut bcf::Record) -> Option<bcf::Record>;
}

/// A small set whose membership is by [`LocusKey`] but which preserves first-seen
/// insertion order (the Python iterates `frozenset`s, but for deterministic output
/// we keep order; the final write is sorted anyway). Dedups like a Python `set`.
struct KeyedRecords {
    seen: HashSet<LocusKey>,
    recs: Vec<bcf::Record>,
}

impl KeyedRecords {
    fn new() -> Self {
        KeyedRecords { seen: HashSet::new(), recs: Vec::new() }
    }
    /// `set.add(rec)` — insert only if its key is new.
    fn add(&mut self, rec: bcf::Record) {
        if self.seen.insert(LocusKey::of(&rec)) {
            self.recs.push(rec);
        }
    }
}

/// Resolve all records co-located with `head` from a stream, given the engine has
/// already taken `head` off the stream. Returns `(co_located, next_downstream)`
/// where `next_downstream` is the first record strictly after `head` (pushed back
/// for the outer loop), if any. Mirrors the two drain-loops in
/// `deal_with_same_loc_variants` (merge L186-217).
fn drain_co_located(
    head: bcf::Record,
    stream: &mut std::vec::IntoIter<bcf::Record>,
) -> (Vec<bcf::Record>, Option<bcf::Record>) {
    let mut co_located = vec![head];
    let mut downstream = None;
    for next in stream.by_ref() {
        match position_cmp(&co_located[0], &next) {
            Ordering::Equal => co_located.push(next),
            Ordering::Less => {
                downstream = Some(next);
                break;
            }
            Ordering::Greater => {
                // Python raises ValueError("records are not sorted by position").
                // Inputs are bcftools-sorted, so this is a programmer/data error.
                panic!("co-iteration input is not sorted by position");
            }
        }
    }
    (co_located, downstream)
}

/// Handle a same-location cluster: drain both streams, O(n·m) allele match, route
/// matched pairs through `op`, unmatched into their `*_only` sets. The drained
/// downstream records (one per stream, if any) are returned so the outer loop can
/// resume with them (the Python `buffer_var{1,2}`).
#[allow(clippy::too_many_arguments)]
fn deal_with_same_loc(
    q_head: bcf::Record,
    r_head: bcf::Record,
    q_stream: &mut std::vec::IntoIter<bcf::Record>,
    r_stream: &mut std::vec::IntoIter<bcf::Record>,
    matched: &mut KeyedRecords,
    query_only: &mut KeyedRecords,
    ref_only: &mut KeyedRecords,
    op: &dyn LocusOp,
) -> (Option<bcf::Record>, Option<bcf::Record>) {
    let (q_loc, q_down) = drain_co_located(q_head, q_stream);
    let (mut r_loc, r_down) = drain_co_located(r_head, r_stream);

    // Track which ref records got matched (by index) so the rest go to ref_only.
    let mut r_matched = vec![false; r_loc.len()];

    for mut q in q_loc {
        let q_key = LocusKey::of(&q);
        let mut found = false;
        for (i, r) in r_loc.iter_mut().enumerate() {
            if r_matched[i] {
                continue;
            }
            if LocusKey::of(r) == q_key {
                // matched (location AND alleles) → callback.
                if let Some(kept) = op.on_match(&mut q, r) {
                    matched.add(kept);
                }
                r_matched[i] = true;
                found = true;
                break;
            }
        }
        if !found {
            query_only.add(q);
        }
    }
    // Unmatched ref records → ref_only (Python L236).
    for (i, r) in r_loc.into_iter().enumerate() {
        if !r_matched[i] {
            ref_only.add(r);
        }
    }

    (q_down, r_down)
}

/// THE engine. Co-iterate two coordinate-sorted, single-contig record streams,
/// routing each locus through `op`. Records are consumed by value (already
/// translated into the output header by the caller).
///
/// `&dyn LocusOp` (not generic) keeps this module non-generic and compile-fast —
/// there are exactly two ops and the I/O dominates dynamic-dispatch cost.
pub fn coiterate_sorted_vcfs(
    query_recs: Vec<bcf::Record>,
    ref_recs: Vec<bcf::Record>,
    op: &dyn LocusOp,
) -> CoiterSets {
    let mut matched = KeyedRecords::new();
    let mut query_only = KeyedRecords::new();
    let mut ref_only = KeyedRecords::new();

    let mut q_stream = query_recs.into_iter();
    let mut r_stream = ref_recs.into_iter();

    // Edge cases: an empty stream sends the whole other stream to its *_only set
    // (Python process_region L268-272 / inhouse L245-251).
    let mut q_next = q_stream.next();
    let mut r_next = r_stream.next();
    if q_next.is_none() {
        if let Some(r) = r_next {
            ref_only.add(r);
        }
        for r in r_stream {
            ref_only.add(r);
        }
        return CoiterSets {
            matched: matched.recs,
            query_only: query_only.recs,
            ref_only: ref_only.recs,
        };
    }
    if r_next.is_none() {
        if let Some(q) = q_next {
            query_only.add(q);
        }
        for q in q_stream {
            query_only.add(q);
        }
        return CoiterSets {
            matched: matched.recs,
            query_only: query_only.recs,
            ref_only: ref_only.recs,
        };
    }

    // Main two-pointer loop. Invariant: q_next / r_next hold the current heads.
    loop {
        let (q, r) = match (q_next.take(), r_next.take()) {
            (Some(q), Some(r)) => (q, r),
            // One stream exhausted → flush the other entirely.
            (Some(q), None) => {
                query_only.add(q);
                for q in q_stream.by_ref() {
                    query_only.add(q);
                }
                break;
            }
            (None, Some(r)) => {
                ref_only.add(r);
                for r in r_stream.by_ref() {
                    ref_only.add(r);
                }
                break;
            }
            (None, None) => break,
        };

        match position_cmp(&q, &r) {
            Ordering::Less => {
                // query before ref → query is non-overlapping; advance query,
                // keep ref as the current head.
                query_only.add(q);
                q_next = q_stream.next();
                r_next = Some(r);
            }
            Ordering::Greater => {
                // ref before query → ref is non-overlapping; advance ref,
                // keep query as the current head.
                ref_only.add(r);
                q_next = Some(q);
                r_next = r_stream.next();
            }
            Ordering::Equal => {
                let (q_down, r_down) =
                    deal_with_same_loc(q, r, &mut q_stream, &mut r_stream, &mut matched, &mut query_only, &mut ref_only, op);
                q_next = q_down.or_else(|| q_stream.next());
                r_next = r_down.or_else(|| r_stream.next());
            }
        }
    }

    CoiterSets {
        matched: matched.recs,
        query_only: query_only.recs,
        ref_only: ref_only.recs,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bcf::record::GenotypeAllele;
    use rust_htslib::bcf::{Format, Writer};

    /// A `LocusOp` that keeps the QUERY record verbatim (no mutation) — lets us
    /// test the engine's set-routing in isolation.
    struct KeepQuery;
    impl LocusOp for KeepQuery {
        fn on_match(&self, q: &mut bcf::Record, _r: &mut bcf::Record) -> Option<bcf::Record> {
            Some(q.clone())
        }
    }

    /// Build a writer over a temp VCF whose header has one contig + GT, so we can
    /// mint records that all share one header (required for `LocusKey` comparison).
    fn writer() -> (Writer, tempfile::TempPath) {
        let tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        let mut header = bcf::Header::new();
        header.push_record(b"##contig=<ID=chr1,length=1000000>");
        header.push_record(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        header.push_sample(b"s1");
        let w = Writer::from_path(tmp.path(), &header, true, Format::Vcf).unwrap();
        (w, tmp.into_temp_path())
    }

    /// Mint a record at 0-based `pos` with `(ref, alt)` alleles and a 0/1 GT.
    fn rec(w: &Writer, pos: i64, ref_a: &str, alt: &str) -> bcf::Record {
        let mut r = w.empty_record();
        let rid = w.header().name2rid(b"chr1").unwrap();
        r.set_rid(Some(rid));
        r.set_pos(pos);
        r.set_alleles(&[ref_a.as_bytes(), alt.as_bytes()]).unwrap();
        r.push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
            .unwrap();
        r
    }

    fn keys(recs: &[bcf::Record]) -> Vec<(i64, String, String)> {
        let mut v: Vec<_> = recs
            .iter()
            .map(|r| {
                let a = r.alleles();
                (
                    r.pos(),
                    String::from_utf8_lossy(a[0]).to_string(),
                    String::from_utf8_lossy(a[1]).to_string(),
                )
            })
            .collect();
        v.sort();
        v
    }

    #[test]
    fn query_only_stream() {
        let (w, _t) = writer();
        let q = vec![rec(&w, 100, "A", "T"), rec(&w, 200, "C", "G")];
        let sets = coiterate_sorted_vcfs(q, vec![], &KeepQuery);
        assert_eq!(keys(&sets.query_only), vec![(100, "A".into(), "T".into()), (200, "C".into(), "G".into())]);
        assert!(sets.matched.is_empty() && sets.ref_only.is_empty());
    }

    #[test]
    fn ref_only_stream() {
        let (w, _t) = writer();
        let r = vec![rec(&w, 100, "A", "T"), rec(&w, 200, "C", "G")];
        let sets = coiterate_sorted_vcfs(vec![], r, &KeepQuery);
        assert_eq!(keys(&sets.ref_only), vec![(100, "A".into(), "T".into()), (200, "C".into(), "G".into())]);
        assert!(sets.matched.is_empty() && sets.query_only.is_empty());
    }

    #[test]
    fn same_locus_same_allele_matches() {
        let (w, _t) = writer();
        let q = vec![rec(&w, 100, "A", "T")];
        let r = vec![rec(&w, 100, "A", "T")];
        let sets = coiterate_sorted_vcfs(q, r, &KeepQuery);
        assert_eq!(keys(&sets.matched), vec![(100, "A".into(), "T".into())]);
        assert!(sets.query_only.is_empty() && sets.ref_only.is_empty());
    }

    #[test]
    fn same_locus_different_allele_splits() {
        // Same pos, same ref length (so start==stop), different ALT → NOT a match;
        // each goes to its own *_only set (the allele-level __eq__ is finer than
        // the location-level cmp==0).
        let (w, _t) = writer();
        let q = vec![rec(&w, 100, "A", "T")];
        let r = vec![rec(&w, 100, "A", "G")];
        let sets = coiterate_sorted_vcfs(q, r, &KeepQuery);
        assert!(sets.matched.is_empty());
        assert_eq!(keys(&sets.query_only), vec![(100, "A".into(), "T".into())]);
        assert_eq!(keys(&sets.ref_only), vec![(100, "A".into(), "G".into())]);
    }

    #[test]
    fn multiple_co_located_records_drain_and_match() {
        // Two query + two ref records at the same pos; one allele pair matches,
        // the rest split. Exercises the O(n·m) drain + match.
        let (w, _t) = writer();
        let q = vec![rec(&w, 100, "A", "T"), rec(&w, 100, "A", "C")];
        let r = vec![rec(&w, 100, "A", "C"), rec(&w, 100, "A", "G")];
        let sets = coiterate_sorted_vcfs(q, r, &KeepQuery);
        assert_eq!(keys(&sets.matched), vec![(100, "A".into(), "C".into())]);
        assert_eq!(keys(&sets.query_only), vec![(100, "A".into(), "T".into())]);
        assert_eq!(keys(&sets.ref_only), vec![(100, "A".into(), "G".into())]);
    }

    #[test]
    fn interleaved_non_overlapping() {
        // q: 100, 300 ; r: 200, 400 → all four non-overlapping, correct routing.
        let (w, _t) = writer();
        let q = vec![rec(&w, 100, "A", "T"), rec(&w, 300, "A", "T")];
        let r = vec![rec(&w, 200, "C", "G"), rec(&w, 400, "C", "G")];
        let sets = coiterate_sorted_vcfs(q, r, &KeepQuery);
        assert_eq!(keys(&sets.query_only), vec![(100, "A".into(), "T".into()), (300, "A".into(), "T".into())]);
        assert_eq!(keys(&sets.ref_only), vec![(200, "C".into(), "G".into()), (400, "C".into(), "G".into())]);
        assert!(sets.matched.is_empty());
    }

    #[test]
    fn overlapping_indels_ordered_by_stop() {
        // Two records sharing start=100 but different stop (REF length): a SNV
        // (stop 101) and a deletion (REF "ACGT", stop 104). They are at different
        // *locations* by the (start,stop) order, so a query SNV at 100 and a ref
        // deletion at 100 do NOT co-locate and both go to *_only.
        let (w, _t) = writer();
        let q = vec![rec(&w, 100, "A", "T")]; // stop 101
        let r = vec![rec(&w, 100, "ACGT", "A")]; // stop 104
        let sets = coiterate_sorted_vcfs(q, r, &KeepQuery);
        assert!(sets.matched.is_empty());
        assert_eq!(keys(&sets.query_only), vec![(100, "A".into(), "T".into())]);
        assert_eq!(keys(&sets.ref_only), vec![(100, "ACGT".into(), "A".into())]);
    }

    #[test]
    fn set_dedup_on_locuskey() {
        // A query record duplicated (same chrom,pos,alleles) is added once (the
        // Python `set` semantics).
        let (w, _t) = writer();
        let q = vec![rec(&w, 100, "A", "T"), rec(&w, 100, "A", "T")];
        let sets = coiterate_sorted_vcfs(q, vec![], &KeepQuery);
        // Both are co-located + identical; with no ref they flush to query_only,
        // and the dedup keeps one.
        assert_eq!(sets.query_only.len(), 1);
    }
}
