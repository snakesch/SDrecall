//! THE compute unit: strand-aware relative-coordinate segment projection.
//!
//! A faithful port of the pure body of `extract_and_pad_segments`
//! (`realign_recall/prepare_masked_align_region.py:10-121`). The file-write +
//! freshness-reuse tail (lines 123-141) is **not** here — it stays in the
//! orchestrator ([`crate::per_rg`]), because it is I/O bookkeeping, not compute.
//!
//! ## What it does
//!
//! 1. Intersect the single FC interval with the target → overlapping segments
//!    (`prepare_masked_align_region.py:45`).
//! 2. Pad each by `padding`, convert to coordinates relative to the FC interval's
//!    start, and merge overlapping padded segments (lines 53-76).
//! 3. For every NFC interval, slice out the parts that correspond to each merged
//!    relative segment, **strand-aware** (same strand → straight shift; opposite
//!    strand → reverse-complement flip), clamp to the NFC interval bounds, and
//!    emit (lines 81-121).
//!
//! The output is sorted+merged + written by the caller (it is the durable nfc
//! artifact). This function only builds the un-merged emitted intervals.
//!
//! ## Parity hazards handled here (see `T7_region_prep.md` §6)
//!
//! - **R1** — the line-91 overlap predicate is ported as the exact **OR**
//!   (`rel_start < rei || rel_end > rsi`), not a real AND-overlap.
//! - **R8** — Python `relative_segments[0]` (line 66) IndexErrors when the FC
//!   interval doesn't intersect the target; we return an **empty `Vec`** instead
//!   (a deliberate, documented divergence — more correct, no crash).

use sdrecall_io::intersect;
use sdrecall_utils::{GenomicInterval, Strand};

/// Minimal carrier for an NFC interval's projection inputs.
///
/// Borrows `chrom`/`name` `&str` straight from the source row (zero copy until
/// the output [`GenomicInterval`] is built). `rel_start_interval`/
/// `rel_end_interval` are the row's col4/col5 (the NFC interval projected into
/// the FC node's frame, read at `prepare_masked_align_region.py:84-85`).
#[derive(Clone, Copy, Debug)]
pub struct NfcInterval<'a> {
    /// Contig of the NFC interval (BED col1).
    pub chrom: &'a str,
    /// 0-based start (BED col2).
    pub start: i64,
    /// 0-based exclusive end (BED col3).
    pub end: i64,
    /// Strand (BED col6).
    pub strand: Strand,
    /// `interval.name` kept verbatim (Python line 120). Currently the FC/NFC tag.
    pub name: &'a str,
    /// `rel_start_interval` (col4) — NFC start in the FC node's coordinate frame.
    pub rel_start_interval: i64,
    /// `rel_end_interval` (col5) — NFC end in the FC node's coordinate frame.
    pub rel_end_interval: i64,
}

/// Default NFC padding (`prepare_masked_align_region.py:14`, the `extract_and_pad_segments`
/// `padding=600` default). Distinct from the FC-side slop (`FC_SLOP`); see R6.
pub const NFC_PAD: i64 = 600;

/// Intersect the FC interval with the target, pad each overlap by `padding`,
/// convert to coordinates relative to `fc.start`, and merge overlapping padded
/// segments. Returns the merged `(relative_start, relative_end)` pairs.
///
/// Ports `prepare_masked_align_region.py:45-76`. The absolute start/end (Python
/// tuple slots 2,3) are used here only as the merge sort key / overlap test and
/// never escape, so the returned type carries only the two relative fields (the
/// design's "don't widen the public type" note).
///
/// **R8:** when the FC interval doesn't intersect the target, `overlapping` is
/// empty and we return an empty `Vec` (Python `relative_segments[0]` would
/// `IndexError`).
fn relative_padded_merged(
    fc: &GenomicInterval,
    target: &[GenomicInterval],
    padding: i64,
) -> Vec<(i64, i64)> {
    let fc_slice = std::slice::from_ref(fc);
    let overlapping = intersect(fc_slice, target);
    let main_start = fc.start;

    // (relative_start, relative_end, abs_start, abs_end) — abs are sort/merge-only.
    let mut relative_segments: Vec<(i64, i64, i64, i64)> = Vec::with_capacity(overlapping.len());
    for seg in &overlapping {
        // Pad; clamp the low side to 0 (line 55). No high clamp here — slop on
        // the FC side and the NFC-bound clamp handle the high side later.
        let start = (seg.start - padding).max(0);
        let end = seg.end + padding;
        let relative_start = start - main_start;
        let relative_end = end - main_start;
        relative_segments.push((relative_start, relative_end, start, end));
    }

    // R8: empty intersection → empty result (no IndexError).
    if relative_segments.is_empty() {
        return Vec::new();
    }

    // Sort by the original absolute start (line 65), then fold-merge (67-75).
    relative_segments.sort_by_key(|x| x.2);
    let mut merged: Vec<(i64, i64, i64, i64)> = vec![relative_segments[0]];
    for &current in &relative_segments[1..] {
        let last = *merged.last().unwrap();
        if current.2 <= last.3 {
            // Overlap (abs start ≤ last abs end). Parity quirk (line 72): the
            // merged relative-start is carried from the FIRST segment of the run
            // (last.0), the merged relative-end is recomputed as merged_end - main_start.
            let merged_end = last.3.max(current.3);
            *merged.last_mut().unwrap() = (last.0, merged_end - main_start, last.2, merged_end);
        } else {
            merged.push(current);
        }
    }

    merged.iter().map(|m| (m.0, m.1)).collect()
}

/// Compute the NFC-local `(rel_small_start, rel_small_end)` for one clamped
/// `(rel_start, rel_end)` segment — the collapse of Python's two duplicated
/// strand branches (lines 101-108) into ONE `match`.
///
/// - same strand (line 103-104): `(rel_start - rsi, rel_end - rsi)`.
/// - opposite strand (line 107-108): the reverse-complement flip
///   `(rei - rel_end, rei - rel_start)`.
fn nfc_local_coords(
    same_strand: bool,
    rel_start: i64,
    rel_end: i64,
    rsi: i64,
    rei: i64,
) -> (i64, i64) {
    match same_strand {
        true => (rel_start - rsi, rel_end - rsi),
        false => (rei - rel_end, rei - rel_start),
    }
}

/// Strand-aware relative-coordinate projection + pad + merge. PURE: no I/O, no
/// mtime, no temp files.
///
/// Ports `extract_and_pad_segments` (`prepare_masked_align_region.py:10-121`).
/// The output is the un-merged emitted NFC segments; the caller sorts+merges and
/// writes them.
///
/// `&GenomicInterval` / `&[..]`: inputs are read-only. Owned `Vec` out: the
/// emitted intervals have no borrow tie to the inputs.
pub fn extract_and_pad_segments(
    fc_interval: &GenomicInterval,
    nfc_intervals: &[NfcInterval],
    target: &[GenomicInterval],
    padding: i64,
) -> Vec<GenomicInterval> {
    let merged_segments = relative_padded_merged(fc_interval, target, padding);
    // R8 short-circuit: nothing projected.
    if merged_segments.is_empty() {
        return Vec::new();
    }
    let main_strand = fc_interval.strand;

    // Lower-bound pre-alloc (≥ one segment per NFC interval).
    let mut out: Vec<GenomicInterval> = Vec::with_capacity(nfc_intervals.len());

    for interval in nfc_intervals {
        let interval_strand = interval.strand;
        let rsi = interval.rel_start_interval;
        let rei = interval.rel_end_interval;
        let same_strand = interval_strand == main_strand;

        for &(rel_start, rel_end) in &merged_segments {
            // R1: the line-91 OR predicate, ported verbatim (almost always true).
            if !(rel_start < rei || rel_end > rsi) {
                continue;
            }
            // Clamp the merged relative segment to [rsi, rei) (lines 89-90).
            let rel_start = rel_start.max(rsi);
            let rel_end = rel_end.min(rei);
            // Guard (line 98).
            if rel_start >= rel_end {
                continue;
            }

            let (rel_small_start, rel_small_end) =
                nfc_local_coords(same_strand, rel_start, rel_end, rsi, rei);

            // Absolute coords on the NFC interval (lines 111-112), clamped to its
            // bounds (lines 115-116).
            let abs_start = (interval.start + rel_small_start).max(interval.start);
            let abs_end = (interval.start + rel_small_end).min(interval.end);

            if abs_start < abs_end {
                out.push(GenomicInterval::with_strand(
                    interval.chrom,
                    abs_start,
                    abs_end,
                    interval_strand,
                ));
                // `name` (interval.name, Python line 120) is not retained: the
                // output BED is BED3-after-merge (R2), so the name is dropped by
                // the unstranded merge anyway. Kept on `NfcInterval` for parity
                // documentation / future BED6 consumers.
                let _ = interval.name;
            }
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fc(start: i64, end: i64, s: Strand) -> GenomicInterval {
        GenomicInterval::with_strand("chr1", start, end, s)
    }
    fn tgt(start: i64, end: i64) -> GenomicInterval {
        GenomicInterval::new("chr1", start, end)
    }
    fn nfc<'a>(start: i64, end: i64, s: Strand, rsi: i64, rei: i64) -> NfcInterval<'a> {
        NfcInterval {
            chrom: "chr1",
            start,
            end,
            strand: s,
            name: "NFC:RG0_0",
            rel_start_interval: rsi,
            rel_end_interval: rei,
        }
    }

    // ── same-strand projection ──────────────────────────────────────────────
    #[test]
    fn proj_same_strand() {
        // FC chr1:1000-2000 +; target chr1:1200-1400; padding 100.
        // overlap = [1200,1400). padded abs = [max(1200-100,0), 1400+100) = [1100,1500).
        // relative (to main_start 1000) = [100, 500).
        // NFC chr1:5000-6000 +, col4/col5 = 0..1000 (rsi=0, rei=1000).
        // R1 OR true; clamp [100,500) to [0,1000) = [100,500); same strand:
        // rel_small = (100-0, 500-0) = (100,500); abs = (5000+100, 5000+500) = (5100,5500).
        let fc_iv = fc(1000, 2000, Strand::Forward);
        let target = [tgt(1200, 1400)];
        let nfcs = [nfc(5000, 6000, Strand::Forward, 0, 1000)];
        let out = extract_and_pad_segments(&fc_iv, &nfcs, &target, 100);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].chrom, "chr1");
        assert_eq!(out[0].start, 5100);
        assert_eq!(out[0].end, 5500);
        assert_eq!(out[0].strand, Strand::Forward);
    }

    // ── opposite-strand projection (reverse-complement flip) ────────────────
    #[test]
    fn proj_opposite_strand() {
        // Same FC/target → merged relative [100,500). NFC chr1:5000-6000 -, rsi=0,rei=1000.
        // opposite strand: rel_small = (rei-rel_end, rei-rel_start) = (1000-500, 1000-100) = (500,900).
        // abs = (5000+500, 5000+900) = (5500,5900).
        let fc_iv = fc(1000, 2000, Strand::Forward);
        let target = [tgt(1200, 1400)];
        let nfcs = [nfc(5000, 6000, Strand::Reverse, 0, 1000)];
        let out = extract_and_pad_segments(&fc_iv, &nfcs, &target, 100);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].start, 5500);
        assert_eq!(out[0].end, 5900);
        assert_eq!(out[0].strand, Strand::Reverse);
    }

    // ── padding clamp at contig 0 ───────────────────────────────────────────
    #[test]
    fn pad_clamp_at_zero() {
        // FC chr1:0-2000 +; target chr1:50-100; padding 600.
        // overlap [50,100); padded start = max(50-600,0) = 0 (NOT -550); end = 700.
        // relative = [0 - 0, 700 - 0) = [0,700). main_start is 0 so relative_start = 0.
        let fc_iv = fc(0, 2000, Strand::Forward);
        let target = [tgt(50, 100)];
        // NFC large frame so nothing is clamped away.
        let nfcs = [nfc(5000, 6000, Strand::Forward, 0, 2000)];
        let out = extract_and_pad_segments(&fc_iv, &nfcs, &target, 600);
        assert_eq!(out.len(), 1);
        // relative [0,700) projected onto NFC [5000..]: same strand → [5000, 5700).
        assert_eq!(out[0].start, 5000);
        assert_eq!(out[0].end, 5700);
    }

    // ── merge of adjacent/overlapping padded relative segments ──────────────
    #[test]
    fn merge_adjacent_relative() {
        // FC chr1:1000-3000 +; two target overlaps whose padded abs ranges touch.
        // target A [1200,1300] padded (pad 100) → abs [1100,1400); B [1350,1450] → abs [1250,1550).
        // B.abs_start 1250 <= A.abs_end 1400 → merge; merged_end = max(1400,1550) = 1550.
        // merged relative carried = (A.rel_start=100, 1550-1000=550) → [100,550).
        let fc_iv = fc(1000, 3000, Strand::Forward);
        let target = [tgt(1200, 1300), tgt(1350, 1450)];
        let nfcs = [nfc(5000, 7000, Strand::Forward, 0, 2000)];
        let out = extract_and_pad_segments(&fc_iv, &nfcs, &target, 100);
        // one merged relative segment [100,550) → same strand projects to [5100,5550).
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].start, 5100);
        assert_eq!(out[0].end, 5550);
    }

    // ── bookended projected output stays SEPARATE after sort_merge (d=0) ─────
    // NOTE: the crate's durable-output merge (`per_rg::merge_bookended`) FUSES
    // book-ended projected segments — pybedtools `.merge()` parity (R2-bis),
    // proven by the HG002 differential. That merge + its book-ended behaviour are
    // unit-tested in `per_rg::tests`; the projection here only produces the
    // un-merged segments, so there is no merge assertion at this layer.

    // ── R1: OR predicate true but clamp drops it to empty ───────────────────
    #[test]
    fn overlap_predicate_or_then_clamp_drops() {
        // Construct a merged relative segment fully ABOVE [rsi,rei): rel = [1500,2000),
        // NFC frame rsi=0,rei=1000. R1 OR: (1500<1000)=false || (2000>0)=true → passes.
        // clamp: rel_start=max(1500,0)=1500, rel_end=min(2000,1000)=1000 → 1500>=1000 → drop.
        // To get rel [1500,2000): FC chr1:1000-4000; target [2500,3000); pad 0 →
        // abs [2500,3000) relative [1500,2000).
        let fc_iv = fc(1000, 4000, Strand::Forward);
        let target = [tgt(2500, 3000)];
        let nfcs = [nfc(5000, 6000, Strand::Forward, 0, 1000)];
        let out = extract_and_pad_segments(&fc_iv, &nfcs, &target, 0);
        assert!(out.is_empty(), "expected empty, got {out:?}");
    }

    // ── interval-bound clamp (abs_end = min(interval.end, ..)) ───────────────
    #[test]
    fn interval_bound_clamp() {
        // Projection that would exceed interval.end. FC chr1:1000-2000 +; target
        // [1200,1900); pad 0 → abs [1200,1900) relative [200,900).
        // NFC chr1:5000-5500 (len 500), rsi=0,rei=1000. same strand:
        // rel_small=(200,900); abs=(5200,5900) → clamp abs_end=min(5500,5900)=5500.
        let fc_iv = fc(1000, 2000, Strand::Forward);
        let target = [tgt(1200, 1900)];
        let nfcs = [nfc(5000, 5500, Strand::Forward, 0, 1000)];
        let out = extract_and_pad_segments(&fc_iv, &nfcs, &target, 0);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].start, 5200);
        assert_eq!(out[0].end, 5500); // clamped
    }

    // ── R8: FC disjoint from target → empty Vec (no panic) ──────────────────
    #[test]
    fn empty_no_intersect() {
        let fc_iv = fc(1000, 2000, Strand::Forward);
        let target = [tgt(8000, 9000)]; // disjoint
        let nfcs = [nfc(5000, 6000, Strand::Forward, 0, 1000)];
        let out = extract_and_pad_segments(&fc_iv, &nfcs, &target, 100);
        assert!(out.is_empty());
    }

    #[test]
    fn empty_no_nfc_intervals() {
        let fc_iv = fc(1000, 2000, Strand::Forward);
        let target = [tgt(1200, 1400)];
        let out = extract_and_pad_segments(&fc_iv, &[], &target, 100);
        assert!(out.is_empty());
    }
}
