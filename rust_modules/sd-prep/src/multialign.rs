//! Multi-align depth — the ONE per-base coverage sweep + the 4-way AND filter.
//!
//! Ports `preparation/pick_multialign_regions.py::pick_multialigned_regions`
//! (l.7-63) + `inferred_depths.py::calculate_inferred_coverage` (l.67-136).
//!
//! ## What is ported vs stubbed
//!
//! - [`depth_sweep`] — the per-base coverage kernel: given a set of read intervals
//!   (already filtered), produce per-position `(chrom, pos, depth)` rows EXACTLY
//!   matching `BedTool.genome_coverage(bg=True)` then the row-expansion
//!   `for i in range(start+1, end+1)` (inferred_depths.py l.124-127). This is the
//!   ONE coverage unit, reused for all 4 passes (raw / high-MQ / XA / XS).
//!   UNIT-TESTED.
//! - [`multialign_filter_mask`] — the 4-way AND filter
//!   `raw>=min_depth & (XA_frac>=frac | XS_frac>=frac) & high_MQ<=hq_depth`
//!   (pick_multialign_regions.py l.43-48) over the merged per-position table.
//!   UNIT-TESTED.
//! - [`DepthPass`] — the 4 pass variants (Raw / HighMq / Xa / Xs), each a
//!   `(min_mapq, tag predicate)` tuple. ONE enum drives the read predicate.
//! - **`inferred_coverage` (BAM read) — STUBBED** (`TODO(T8)`): the
//!   `bam.fetch(region)` + `filter_and_process_read` per-read filter
//!   (inferred_depths.py l.10-29, l.91-110). The read predicate is fully specified
//!   here but the rust-htslib `IndexedReader::fetch` loop is left as a TODO so the
//!   driver wiring is explicit and testable later. The per-read filter semantics
//!   (the `AS-XS<=5` XS rule, the flag filters) are documented on
//!   [`read_passes`].

use ahash::AHashMap;

/// A depth pass: the MAPQ floor and which tag (if any) a read must carry.
///
/// The 4 production passes (pick_multialign_regions.py l.21-24):
/// - `Raw`     — `min_mapq=0`, no tag filter.
/// - `HighMq`  — `min_mapq=MQ_threshold` (e.g. 41), no tag filter.
/// - `Xa`      — `min_mapq=0`, read must have an `XA` tag.
/// - `Xs`      — `min_mapq=0`, read must satisfy the XS rule (`AS-XS<=5`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum DepthPass {
    Raw,
    HighMq,
    Xa,
    Xs,
}

/// Per-base coverage from a set of read intervals — the EXACT analog of
/// `BedTool.genome_coverage(bg=True)` followed by the per-position row expansion
/// `for i in range(interval.start + 1, interval.end + 1)` (inferred_depths.py
/// l.124-127).
///
/// `reads` are half-open `[start, end)` 0-based read intervals per contig (already
/// passing the per-read filter). The bedtools genomecov `bg` output reports, per
/// maximal constant-depth run `[s, e)`, the depth `d`; Python then emits one row
/// per position `pos = s+1 .. e` (1-based) with depth `d`. We compute the same via
/// a sweep-line over interval endpoints and emit one `(chrom, pos_1based, depth)`
/// row per covered base with `depth > 0` (bedtools `bg` omits zero-depth runs).
///
/// Output is sorted by `(chrom, pos)` and contains only `depth > 0` positions.
/// `&AHashMap<chrom, Vec<(start,end)>>` borrow-in; owned rows out.
pub fn depth_sweep(reads: &AHashMap<String, Vec<(i64, i64)>>) -> Vec<(String, i64, i64)> {
    let mut out: Vec<(String, i64, i64)> = Vec::new();
    // Deterministic chrom order.
    let mut chroms: Vec<&String> = reads.keys().collect();
    chroms.sort();
    for chrom in chroms {
        let intervals = &reads[chrom];
        if intervals.is_empty() {
            continue;
        }
        // Sweep-line: +1 at start, -1 at end.
        let mut events: Vec<(i64, i64)> = Vec::with_capacity(intervals.len() * 2);
        for &(s, e) in intervals {
            if e <= s {
                continue;
            }
            events.push((s, 1));
            events.push((e, -1));
        }
        events.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        // Walk runs of constant depth. Between consecutive distinct event
        // positions p0..p1 the depth is the running sum; emit per-base rows
        // (1-based pos = base0 + 1, matching range(start+1, end+1)).
        let mut depth: i64 = 0;
        let mut i = 0;
        let n = events.len();
        while i < n {
            let pos = events[i].0;
            // Apply all events at this position.
            while i < n && events[i].0 == pos {
                depth += events[i].1;
                i += 1;
            }
            // Depth `depth` now holds from `pos` until the next event position.
            if depth > 0 && i < n {
                let next = events[i].0;
                // Emit per-base rows for [pos, next): 1-based positions pos+1..=next.
                for base0 in pos..next {
                    out.push((chrom.clone(), base0 + 1, depth));
                }
            }
        }
    }
    out
}

/// The merged per-position depth table for one base: the 4 pass depths.
/// `raw`/`high_mq`/`xa`/`xs` mirror the renamed columns the Python merge produces
/// (pick_multialign_regions.py l.29-34), with missing passes filled to 0.
#[derive(Clone, Copy, Debug, Default)]
pub struct MergedDepth {
    pub raw: i64,
    pub high_mq: i64,
    pub xa: i64,
    pub xs: i64,
}

/// The 4-way AND filter (pick_multialign_regions.py l.43-48): a base is kept iff
///
/// `raw >= min_depth && (xa/raw >= frac || xs/raw >= frac) && high_mq <= hq_depth`.
///
/// Returns the per-base boolean mask aligned to `table`. Division uses the raw
/// depth; `raw == 0` rows fail `raw >= min_depth` (min_depth ≥ 1 in production), so
/// the fraction is only evaluated when `raw > 0` — guarding the divide.
pub fn multialign_filter_mask(
    table: &[MergedDepth],
    min_depth: i64,
    hq_depth: i64,
    frac: f64,
) -> Vec<bool> {
    table
        .iter()
        .map(|d| {
            let min_ok = d.raw >= min_depth;
            if !min_ok {
                return false;
            }
            let xa_frac = d.xa as f64 / d.raw as f64;
            let xs_frac = d.xs as f64 / d.raw as f64;
            let multi_ok = xa_frac >= frac || xs_frac >= frac;
            let not_enough_evidence = d.high_mq <= hq_depth;
            multi_ok && not_enough_evidence
        })
        .collect()
}

/// Per-read accept predicate for a pass — the analog of
/// `inferred_depths.py::filter_and_process_read` (l.10-29).
///
/// A read passes iff: `mapq >= min_mapq`, not unmapped/dup/secondary/supplementary/
/// qcfail, AND the pass's tag condition:
/// - `Raw`/`HighMq`: no tag condition.
/// - `Xa`: read has an `XA` tag.
/// - `Xs`: read has both `AS` and `XS` tags AND `AS - XS <= 5` (the
///   `read.get_tag("AS") - read.get_tag("XS") <= 5` rule, l.18). If `AS`/`XS` are
///   absent the Python `read.has_tag("XS")` branch is `False` → read rejected.
///
/// `flags` is the SAM flag bitfield; `mapq` the mapping quality; the four
/// `Option<i64>` are the read's `XA`(presence-only, modeled as `has_xa`), `AS`,
/// `XS` tag values. Pure (no htslib) so it is unit-tested directly; the BAM loop
/// (STUBBED) feeds it.
pub fn read_passes(
    pass: DepthPass,
    flags: u16,
    mapq: u8,
    min_mapq: u8,
    has_xa: bool,
    as_tag: Option<i64>,
    xs_tag: Option<i64>,
) -> bool {
    const UNMAPPED: u16 = 0x4;
    const SECONDARY: u16 = 0x100;
    const QCFAIL: u16 = 0x200;
    const DUP: u16 = 0x400;
    const SUPPLEMENTARY: u16 = 0x800;

    let flag_ok = (flags & UNMAPPED) == 0
        && (flags & DUP) == 0
        && (flags & SECONDARY) == 0
        && (flags & SUPPLEMENTARY) == 0
        && (flags & QCFAIL) == 0;
    if !flag_ok || mapq < min_mapq {
        return false;
    }
    match pass {
        DepthPass::Raw | DepthPass::HighMq => true,
        DepthPass::Xa => has_xa,
        DepthPass::Xs => match (as_tag, xs_tag) {
            (Some(a), Some(x)) => a - x <= 5,
            _ => false,
        },
    }
}

/// The pass's MAPQ floor. `HighMq` uses `mq_threshold`; the rest use 0.
pub fn pass_min_mapq(pass: DepthPass, mq_threshold: u8) -> u8 {
    match pass {
        DepthPass::HighMq => mq_threshold,
        _ => 0,
    }
}

use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{self, Read as _};
use sdrecall_utils::{GenomicInterval, Result, SdError};
use std::path::Path;

/// Read a signed-integer aux tag (`AS`/`XS`), or `None` if absent / not integral.
fn aux_i64(rec: &bam::Record, tag: &[u8]) -> Option<i64> {
    match rec.aux(tag) {
        Ok(bam::record::Aux::I8(v)) => Some(v as i64),
        Ok(bam::record::Aux::U8(v)) => Some(v as i64),
        Ok(bam::record::Aux::I16(v)) => Some(v as i64),
        Ok(bam::record::Aux::U16(v)) => Some(v as i64),
        Ok(bam::record::Aux::I32(v)) => Some(v as i64),
        Ok(bam::record::Aux::U32(v)) => Some(v as i64),
        _ => None,
    }
}

/// One depth pass over the BAM — the analog of `calculate_inferred_coverage`
/// (inferred_depths.py l.67-136). Fetches each merged `target` region, applies the
/// per-read filter ([`read_passes`]) for `pass`, collects `(reference_start,
/// reference_end)` per read, then runs [`depth_sweep`] → per-position
/// `(chrom, pos, depth)` rows.
///
/// `bam_path` is opened with an index; `mq_threshold` is the `HighMq` MAPQ floor.
/// `&[GenomicInterval]` targets are the already sorted+merged target regions
/// (Python fetches per `chrom:start-end`). The ONE coverage kernel, parameterized
/// by `DepthPass` — no per-pass helper (the #1 coding rule).
pub fn inferred_coverage(
    bam_path: &Path,
    pass: DepthPass,
    targets: &[GenomicInterval],
    mq_threshold: u8,
    threads: u8,
) -> Result<Vec<(String, i64, i64)>> {
    let mut reader = bam::IndexedReader::from_path(bam_path)
        .map_err(|e| SdError::Htslib(format!("open {}: {e}", bam_path.display())))?;
    if threads > 1 {
        let _ = reader.set_threads(threads as usize);
    }
    let header = reader.header().to_owned();
    let min_mapq = pass_min_mapq(pass, mq_threshold);

    let mut reads: AHashMap<String, Vec<(i64, i64)>> = AHashMap::new();
    let mut rec = bam::Record::new();
    for t in targets {
        // fetch by (tid, start, end); 0-based half-open like pysam fetch.
        let tid = match header.tid(t.chrom.as_bytes()) {
            Some(tid) => tid as i32,
            None => continue, // contig absent from BAM header → no coverage
        };
        reader.fetch((tid, t.start, t.end)).map_err(|e| {
            SdError::Htslib(format!("fetch {}:{}-{}: {e}", t.chrom, t.start, t.end))
        })?;
        while let Some(r) = reader.read(&mut rec) {
            r.map_err(|e| SdError::Htslib(format!("read record: {e}")))?;
            let has_xa = rec.aux(b"XA").is_ok();
            let keep = read_passes(
                pass,
                rec.flags(),
                rec.mapq(),
                min_mapq,
                has_xa,
                aux_i64(&rec, b"AS"),
                aux_i64(&rec, b"XS"),
            );
            if !keep {
                continue;
            }
            let start = rec.reference_start();
            let end = rec.reference_end();
            if end <= start {
                continue;
            }
            reads.entry(t.chrom.clone()).or_default().push((start, end));
        }
    }
    Ok(depth_sweep(&reads))
}

/// Merge the 4 per-position depth tables on `(chrom, pos)` into [`MergedDepth`]
/// rows, exactly as the Python `raw.merge(high_mq).merge(xa).merge(xs).fillna(0)`
/// (pick_multialign_regions.py l.34-35): the `raw` pass defines the row universe
/// (a left-merge), missing passes fill 0. Returns the merged table aligned to a
/// sorted `(chrom, pos)` key vector.
fn merge_depth_tables(
    raw: &[(String, i64, i64)],
    high_mq: &[(String, i64, i64)],
    xa: &[(String, i64, i64)],
    xs: &[(String, i64, i64)],
) -> Vec<((String, i64), MergedDepth)> {
    use std::collections::BTreeMap;
    let index = |rows: &[(String, i64, i64)]| -> AHashMap<(String, i64), i64> {
        rows.iter().map(|(c, p, d)| ((c.clone(), *p), *d)).collect()
    };
    let hi = index(high_mq);
    let xa_i = index(xa);
    let xs_i = index(xs);

    // raw defines the universe (Python left-merges onto raw). Dedup raw keys and
    // keep deterministic (chrom,pos) order via BTreeMap.
    let mut out: BTreeMap<(String, i64), MergedDepth> = BTreeMap::new();
    for (c, p, d) in raw {
        let key = (c.clone(), *p);
        let entry = out.entry(key.clone()).or_default();
        entry.raw = *d;
        entry.high_mq = hi.get(&key).copied().unwrap_or(0);
        entry.xa = xa_i.get(&key).copied().unwrap_or(0);
        entry.xs = xs_i.get(&key).copied().unwrap_or(0);
    }
    out.into_iter().collect()
}

/// The full multi-align region pick — port of `pick_multialigned_regions`
/// (pick_multialign_regions.py l.7-63). Runs the 4 depth passes (raw / high-MQ /
/// XA / XS) in parallel, merges them, applies the 4-way AND filter
/// ([`multialign_filter_mask`]), converts kept positions to `start=pos-1,end=pos`
/// BED intervals, sort+merges them, and intersects with `target` (sort+merge again).
///
/// Returns the final multi-align BED intervals. `bam_path` + `target` borrow in;
/// the 4 passes share the BAM read-only via rayon. Mirrors the Python AND logic
/// `raw>=min_depth & (xa_frac>=frac | xs_frac>=frac) & high_mq<=hq_depth`.
pub fn pick_multialigned_regions(
    bam_path: &Path,
    target: &[GenomicInterval],
    mq_threshold: u8,
    high_quality_depth: i64,
    minimum_depth: i64,
    multialign_frac: f64,
    threads: u8,
) -> Result<Vec<GenomicInterval>> {
    use rayon::prelude::*;
    // Python merges the target BED before fetching (inferred_depths l.93).
    let merged_target = sdrecall_io::sort_merge_bed(target, false);

    // 4 passes in parallel — one closure per DepthPass (no per-pass helper).
    // Python uses Pool(min(4, threads)); a scoped rayon pool bounds the fan-out.
    let n_workers = (threads as usize).clamp(1, 4);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_workers)
        .build()
        .map_err(|e| SdError::Compute(format!("rayon pool build failed: {e}")))?;
    let passes = [
        DepthPass::Raw,
        DepthPass::HighMq,
        DepthPass::Xa,
        DepthPass::Xs,
    ];
    let results: Vec<Result<Vec<(String, i64, i64)>>> = pool.install(|| {
        passes
            .par_iter()
            .map(|&p| inferred_coverage(bam_path, p, &merged_target, mq_threshold, 1))
            .collect()
    });
    let mut tables = Vec::with_capacity(4);
    for r in results {
        tables.push(r?);
    }
    let merged = merge_depth_tables(&tables[0], &tables[1], &tables[2], &tables[3]);

    // Apply the AND filter and emit BED intervals (start=pos-1, end=pos).
    let table: Vec<MergedDepth> = merged.iter().map(|(_, d)| *d).collect();
    let mask = multialign_filter_mask(&table, minimum_depth, high_quality_depth, multialign_frac);
    let mut bed: Vec<GenomicInterval> = Vec::new();
    for ((c, pos), keep) in merged.iter().map(|(k, _)| k).zip(mask.iter()) {
        if *keep {
            bed.push(GenomicInterval::new(c.clone(), pos - 1, *pos));
        }
    }
    let multi = sdrecall_io::sort_merge_bed(&bed, false);

    // Intersect with target, sort+merge (Python l.55-57).
    let isect = sdrecall_io::intersect(&multi, &merged_target);
    Ok(sdrecall_io::sort_merge_bed(&isect, false))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reads(pairs: &[(&str, i64, i64)]) -> AHashMap<String, Vec<(i64, i64)>> {
        let mut m: AHashMap<String, Vec<(i64, i64)>> = AHashMap::new();
        for &(c, s, e) in pairs {
            m.entry(c.to_string()).or_default().push((s, e));
        }
        m
    }

    #[test]
    fn depth_single_read() {
        // One read [10,13) → bases 11,12,13 (1-based) each depth 1.
        let r = reads(&[("chr1", 10, 13)]);
        let out = depth_sweep(&r);
        assert_eq!(
            out,
            vec![
                ("chr1".to_string(), 11, 1),
                ("chr1".to_string(), 12, 1),
                ("chr1".to_string(), 13, 1),
            ]
        );
    }

    #[test]
    fn depth_overlapping_reads_sum() {
        // [10,15) and [12,18). bg runs: [10,12)=1, [12,15)=2, [15,18)=1.
        // 1-based rows: 11=1,12=1, 13=2,14=2,15=2, 16=1,17=1,18=1.
        let r = reads(&[("chr1", 10, 15), ("chr1", 12, 18)]);
        let out = depth_sweep(&r);
        let expected = vec![
            ("chr1".to_string(), 11, 1),
            ("chr1".to_string(), 12, 1),
            ("chr1".to_string(), 13, 2),
            ("chr1".to_string(), 14, 2),
            ("chr1".to_string(), 15, 2),
            ("chr1".to_string(), 16, 1),
            ("chr1".to_string(), 17, 1),
            ("chr1".to_string(), 18, 1),
        ];
        assert_eq!(out, expected);
    }

    #[test]
    fn depth_zero_gaps_omitted() {
        // [10,12) and [20,22): gap [12,20) has depth 0 and is omitted (bg behavior).
        let r = reads(&[("chr1", 10, 12), ("chr1", 20, 22)]);
        let out = depth_sweep(&r);
        assert_eq!(
            out,
            vec![
                ("chr1".to_string(), 11, 1),
                ("chr1".to_string(), 12, 1),
                ("chr1".to_string(), 21, 1),
                ("chr1".to_string(), 22, 1),
            ]
        );
    }

    #[test]
    fn filter_mask_and_logic() {
        // min_depth=3, hq_depth=10, frac=0.5.
        let table = vec![
            // pass: raw 5, xa 3 (0.6>=0.5), xs 0, high_mq 2 (<=10) → keep
            MergedDepth {
                raw: 5,
                xa: 3,
                xs: 0,
                high_mq: 2,
            },
            // fail min_depth: raw 2
            MergedDepth {
                raw: 2,
                xa: 2,
                xs: 2,
                high_mq: 0,
            },
            // fail multi: xa 1/5=0.2, xs 1/5=0.2 both <0.5
            MergedDepth {
                raw: 5,
                xa: 1,
                xs: 1,
                high_mq: 0,
            },
            // fail high_mq: high_mq 11 > 10
            MergedDepth {
                raw: 5,
                xa: 5,
                xs: 0,
                high_mq: 11,
            },
            // keep via XS branch: xs 4/5=0.8>=0.5
            MergedDepth {
                raw: 5,
                xa: 0,
                xs: 4,
                high_mq: 3,
            },
        ];
        let mask = multialign_filter_mask(&table, 3, 10, 0.5);
        assert_eq!(mask, vec![true, false, false, false, true]);
    }

    #[test]
    fn read_passes_flag_and_mapq() {
        // proper mapped read, mapq 60, Raw pass, min_mapq 0 → pass.
        assert!(read_passes(DepthPass::Raw, 0x2, 60, 0, false, None, None));
        // unmapped → fail.
        assert!(!read_passes(DepthPass::Raw, 0x4, 60, 0, false, None, None));
        // mapq below floor for HighMq.
        assert!(!read_passes(
            DepthPass::HighMq,
            0x2,
            40,
            41,
            false,
            None,
            None
        ));
        assert!(read_passes(
            DepthPass::HighMq,
            0x2,
            41,
            41,
            false,
            None,
            None
        ));
    }

    #[test]
    fn read_passes_tag_logic() {
        // Xa pass needs XA tag.
        assert!(read_passes(DepthPass::Xa, 0x2, 0, 0, true, None, None));
        assert!(!read_passes(DepthPass::Xa, 0x2, 0, 0, false, None, None));
        // Xs pass: AS-XS <= 5 keeps; > 5 rejects; missing tags reject.
        assert!(read_passes(
            DepthPass::Xs,
            0x2,
            0,
            0,
            false,
            Some(50),
            Some(48)
        )); // 2<=5
        assert!(!read_passes(
            DepthPass::Xs,
            0x2,
            0,
            0,
            false,
            Some(50),
            Some(40)
        )); // 10>5
        assert!(!read_passes(
            DepthPass::Xs,
            0x2,
            0,
            0,
            false,
            None,
            Some(48)
        ));
        assert!(!read_passes(
            DepthPass::Xs,
            0x2,
            0,
            0,
            false,
            Some(50),
            None
        ));
    }

    #[test]
    fn pass_min_mapq_only_highmq_uses_threshold() {
        assert_eq!(pass_min_mapq(DepthPass::HighMq, 41), 41);
        assert_eq!(pass_min_mapq(DepthPass::Raw, 41), 0);
        assert_eq!(pass_min_mapq(DepthPass::Xa, 41), 0);
        assert_eq!(pass_min_mapq(DepthPass::Xs, 41), 0);
    }
}
