//! Insert-size (fragment) distribution — the ONE frag-stat unit.
//!
//! Ports `src/insert_size.py::get_insert_size_distribution`. Each of `num_runs`
//! passes streams the BAM and randomly accepts proper-pair reads (flag `0x2` set,
//! none of `0x90C` set — i.e. not unmapped / secondary / supplementary), whose
//! mate is on the same contig and whose `|TLEN| ≤ max_template_length`, with
//! per-read probability `sample_prob`, stopping at `num_samples`. The accepted
//! `|TLEN|` values are trimmed at the 99th percentile (NumPy *linear*
//! interpolation, replicated exactly — not "element at 0.99·n"), then mean /
//! median / population-std (`np.std` ddof=0) are computed and averaged across the
//! runs that produced data.
//!
//! ## Parity notes (DESIGN §6)
//!
//! - **Monte-Carlo sampling.** Python uses `random.random()`; this uses the
//!   `rand` crate. Exact RNG parity with CPython is impossible, so the
//!   differential test asserts the mean within a tolerance band and the median
//!   exact *after* the percentile trim (the trim + statistics are deterministic
//!   given the accepted set). The deterministic core is exercised in the unit
//!   test by setting `sample_prob = 1.0` on a hand-built BAM.
//! - **`np.percentile(., 99)` linear method.** For sorted `x` of length `n`, the
//!   virtual rank is `h = (n-1)·0.99`; the percentile is
//!   `x[⌊h⌋] + (h-⌊h⌋)·(x[⌊h⌋+1] - x[⌊h⌋])`. Replicated in [`percentile_linear`].

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rust_htslib::bam::{self, Read};
use sdrecall_utils::{Result, SdError};
use std::path::Path;

/// Fragment-size summary statistics. `None` from
/// [`get_insert_size_distribution`] mirrors Python's `(None, None, None)` when a
/// BAM has too few usable pairs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FragStats {
    pub mean: f64,
    pub median: f64,
    pub std: f64,
}

// Flag bits (samtools / SAM spec).
const PROPER_PAIR: u16 = 0x2;
/// `0x90C` = supplementary (0x800) | secondary (0x100) | unmapped (0x4) — the
/// exact mask Python tests with `not (flag & 0x90C)`.
const EXCLUDE_MASK: u16 = 0x90C;

/// NumPy `np.percentile(sorted_data, q, method="linear")` for `q ∈ [0,100]`.
/// `data` must be sorted ascending and non-empty.
fn percentile_linear(sorted: &[f64], q: f64) -> f64 {
    let n = sorted.len();
    if n == 1 {
        return sorted[0];
    }
    let h = (n as f64 - 1.0) * (q / 100.0);
    let lo = h.floor() as usize;
    let frac = h - lo as f64;
    if lo + 1 < n {
        sorted[lo] + frac * (sorted[lo + 1] - sorted[lo])
    } else {
        sorted[lo]
    }
}

/// Mean of a non-empty slice.
fn mean(xs: &[f64]) -> f64 {
    xs.iter().sum::<f64>() / xs.len() as f64
}

/// `np.median`: mean of the two middle elements for even length. `xs` is sorted.
fn median_sorted(sorted: &[f64]) -> f64 {
    let n = sorted.len();
    let mid = n / 2;
    if n % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[mid]
    }
}

/// `np.std` with `ddof=0` (population standard deviation).
fn std_pop(xs: &[f64], m: f64) -> f64 {
    let var = xs.iter().map(|x| (x - m).powi(2)).sum::<f64>() / xs.len() as f64;
    var.sqrt()
}

/// One sampling pass: stream the BAM, accept reads per the Python predicate, stop
/// at `num_samples`. Returns the accepted `|TLEN|` list (unsorted, as Python
/// collects it). `accept` decides per read (the RNG hook — `sample_prob` test in
/// production, always-accept in the deterministic unit test).
fn collect_one_pass(
    bam_path: &Path,
    num_samples: usize,
    max_template_length: i64,
    mut accept: impl FnMut() -> bool,
) -> Result<Vec<f64>> {
    let mut reader = bam::Reader::from_path(bam_path)
        .map_err(|e| SdError::Htslib(format!("open {}: {e}", bam_path.display())))?;
    let mut insert_sizes: Vec<f64> = Vec::new();
    let mut rec = bam::Record::new();
    while let Some(res) = reader.read(&mut rec) {
        res.map_err(|e| SdError::Htslib(format!("read record: {e}")))?;
        let flag = rec.flags();
        if (flag & PROPER_PAIR) != 0
            && (flag & EXCLUDE_MASK) == 0
            && rec.tid() == rec.mtid()
            && rec.insert_size().abs() <= max_template_length
            && accept()
        {
            insert_sizes.push(rec.insert_size().unsigned_abs() as f64);
            if insert_sizes.len() >= num_samples {
                break;
            }
        }
    }
    Ok(insert_sizes)
}

/// Reduce one pass's accepted `|TLEN|` values to `[mean, median, std]` after the
/// 99th-percentile trim, or `None` if the pass yielded nothing usable. This is
/// the deterministic core (given the accepted set) — the only place the trim +
/// statistics live.
fn stats_for_pass(insert_sizes: &[f64], num_samples: usize) -> Option<FragStats> {
    if insert_sizes.is_empty() {
        return None;
    }
    let mut sorted = insert_sizes.to_vec();
    // `total_cmp` is a total order (never panics on a non-finite TLEN; any NaN
    // sorts to the end). TLEN-derived sizes are finite in practice — this is
    // defensive.
    sorted.sort_by(|a, b| a.total_cmp(b));
    let threshold = percentile_linear(&sorted, 99.0);

    // Python filters the ORIGINAL (insertion-order) list then truncates to
    // num_samples. The percentile/threshold is order-independent, and mean/std
    // are too; median needs sorting, which we do below. Filtering the sorted
    // copy yields the identical multiset, so the statistics match.
    let mut filtered: Vec<f64> = sorted.into_iter().filter(|&s| s <= threshold).collect();
    filtered.truncate(num_samples);
    if filtered.is_empty() {
        return None;
    }
    // `filtered` is still sorted (we filtered a sorted vec, truncation keeps the
    // prefix sorted) → median directly.
    let m = mean(&filtered);
    Some(FragStats {
        mean: m,
        median: median_sorted(&filtered),
        std: std_pop(&filtered, m),
    })
}

/// Average the per-run `FragStats` (Python's `np.mean(np.stack(all_stats),
/// axis=0)`), or `None` if no run produced stats.
fn average_runs(runs: &[FragStats]) -> Option<FragStats> {
    if runs.is_empty() {
        return None;
    }
    let n = runs.len() as f64;
    Some(FragStats {
        mean: runs.iter().map(|s| s.mean).sum::<f64>() / n,
        median: runs.iter().map(|s| s.median).sum::<f64>() / n,
        std: runs.iter().map(|s| s.std).sum::<f64>() / n,
    })
}

/// Estimate the fragment-size distribution of a BAM — port of
/// `get_insert_size_distribution` with the Python defaults (`num_runs=10`,
/// `num_samples=2000`, `max_template_length=8000`, `sample_prob=0.001`).
///
/// Returns `Ok(None)` (Python's `(None, None, None)`) when no run accumulates any
/// usable proper-pair reads. Thin caller of [`insert_size_distribution_params`].
pub fn get_insert_size_distribution(bam: &Path) -> Result<Option<FragStats>> {
    insert_size_distribution_params(bam, 10, 2000, 8000, 0.001)
}

/// Parameterized core (the ONE frag-stat unit): runs `num_runs` Monte-Carlo
/// sampling passes and averages their stats. The public
/// [`get_insert_size_distribution`] is the Python-default caller; tests call this
/// with `sample_prob = 1.0` for a deterministic accepted set.
pub fn insert_size_distribution_params(
    bam: &Path,
    num_runs: usize,
    num_samples: usize,
    max_template_length: i64,
    sample_prob: f64,
) -> Result<Option<FragStats>> {
    let mut rng = StdRng::seed_from_u64(42);
    let mut all_stats: Vec<FragStats> = Vec::new();
    for _ in 0..num_runs {
        let insert_sizes = collect_one_pass(bam, num_samples, max_template_length, || {
            sample_prob >= 1.0 || rng.gen::<f64>() < sample_prob
        })?;
        if let Some(s) = stats_for_pass(&insert_sizes, num_samples) {
            all_stats.push(s);
        }
    }
    Ok(average_runs(&all_stats))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::header::HeaderRecord;
    use rust_htslib::bam::record::{Cigar, CigarString};
    use rust_htslib::bam::{Format, Header, Writer};

    const PAIRED: u16 = 0x1;
    const PROPER: u16 = 0x2;
    const UNMAPPED: u16 = 0x4;
    const READ1: u16 = 0x40;
    const READ2: u16 = 0x80;

    fn header() -> Header {
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 100_000);
        let mut sq2 = HeaderRecord::new(b"SQ");
        sq2.push_tag(b"SN", "chr2");
        sq2.push_tag(b"LN", 100_000);
        let mut h = Header::new();
        h.push_record(&sq);
        h.push_record(&sq2);
        h
    }

    /// One read with the given flags, contig (`tid`), mate contig (`mtid`) and
    /// signed TLEN.
    fn read(flags: u16, tid: i32, mtid: i32, tlen: i64) -> bam::Record {
        let span = 100u32;
        let seq = vec![b'A'; span as usize];
        let quals = vec![40u8; span as usize];
        let cs = CigarString(vec![Cigar::Match(span)]);
        let mut r = bam::Record::new();
        r.set(b"q1", Some(&cs), &seq, &quals);
        r.set_flags(flags);
        r.set_mapq(60);
        r.set_tid(tid);
        r.set_pos(1000);
        r.set_mtid(mtid);
        r.set_mpos(1000 + tlen);
        r.set_insert_size(tlen);
        r
    }

    fn write_bam(records: &[bam::Record]) -> tempfile::NamedTempFile {
        let tmp = tempfile::Builder::new().suffix(".bam").tempfile().unwrap();
        {
            let mut w = Writer::from_path(tmp.path(), &header(), Format::Bam).unwrap();
            for r in records {
                w.write(r).unwrap();
            }
        }
        tmp
    }

    // ── percentile parity with numpy ─────────────────────────────────────────

    #[test]
    fn percentile_linear_matches_numpy() {
        // np.percentile([1,2,3,4], 99) = 3.97 (h = 3*0.99 = 2.97 → 3 + 0.97*(4-3))
        let d = vec![1.0, 2.0, 3.0, 4.0];
        assert!((percentile_linear(&d, 99.0) - 3.97).abs() < 1e-9);
        // np.percentile([10,20,30], 50) = 20 (h = 2*0.5 = 1)
        let d = vec![10.0, 20.0, 30.0];
        assert!((percentile_linear(&d, 50.0) - 20.0).abs() < 1e-9);
        // single element
        assert_eq!(percentile_linear(&[5.0], 99.0), 5.0);
        // np.percentile([1..=100], 99) = 99.01
        let d: Vec<f64> = (1..=100).map(|x| x as f64).collect();
        assert!((percentile_linear(&d, 99.0) - 99.01).abs() < 1e-9);
    }

    #[test]
    fn median_and_std_match_numpy() {
        // np.median([1,2,3,4]) = 2.5 ; np.std([1,2,3,4]) (ddof=0) = 1.1180339...
        let d = vec![1.0, 2.0, 3.0, 4.0];
        assert_eq!(median_sorted(&d), 2.5);
        let m = mean(&d);
        assert!((std_pop(&d, m) - 1.118_033_988_749_895).abs() < 1e-12);
    }

    // ── deterministic end-to-end (sample_prob = 1.0) ─────────────────────────

    #[test]
    fn deterministic_stats_with_full_sampling() {
        // 5 proper pairs with |TLEN| = 100,200,300,400,5000. The 5000 is an
        // outlier above the 99th pct → trimmed. Build > a couple records.
        let recs = vec![
            read(PAIRED | PROPER | READ1, 0, 0, 100),
            read(PAIRED | PROPER | READ2, 0, 0, 200),
            read(PAIRED | PROPER | READ1, 0, 0, 300),
            read(PAIRED | PROPER | READ2, 0, 0, 400),
            read(PAIRED | PROPER | READ1, 0, 0, 5000),
        ];
        let bam = write_bam(&recs);
        // sample_prob = 1.0 → accept every qualifying read; num_runs = 1 → exact.
        let stats = insert_size_distribution_params(bam.path(), 1, 2000, 8000, 1.0)
            .unwrap()
            .expect("some stats");

        // threshold = percentile([100,200,300,400,5000], 99)
        //           = h = 4*0.99 = 3.96 → 400 + 0.96*(5000-400) = 4816.0
        // 5000 > 4816 → trimmed. Remaining: [100,200,300,400].
        // mean = 250, median = 250, std(ddof=0) = sqrt(12500) = 111.8033...
        assert!((stats.mean - 250.0).abs() < 1e-9, "mean {}", stats.mean);
        assert!(
            (stats.median - 250.0).abs() < 1e-9,
            "median {}",
            stats.median
        );
        assert!(
            (stats.std - 111.803_398_874_989_48).abs() < 1e-9,
            "std {}",
            stats.std
        );
    }

    #[test]
    fn excludes_improper_and_cross_contig_reads() {
        // Only one read qualifies (proper, same contig). The others are excluded:
        // not-proper, unmapped, mate-on-other-contig, |TLEN| over the cap.
        let recs = vec![
            read(PAIRED | PROPER | READ1, 0, 0, 300), // qualifies → |TLEN| 300
            read(PAIRED | READ1, 0, 0, 200),          // not proper (no 0x2)
            read(PAIRED | PROPER | UNMAPPED, 0, 0, 200), // unmapped → masked
            read(PAIRED | PROPER | READ1, 0, 1, 200), // mate other contig
            read(PAIRED | PROPER | READ1, 0, 0, 9000), // |TLEN| > 8000
        ];
        let bam = write_bam(&recs);
        let stats = insert_size_distribution_params(bam.path(), 1, 2000, 8000, 1.0)
            .unwrap()
            .expect("one usable read");
        // single value 300 → mean = median = 300, std = 0
        assert_eq!(stats.mean, 300.0);
        assert_eq!(stats.median, 300.0);
        assert_eq!(stats.std, 0.0);
    }

    #[test]
    fn no_usable_reads_returns_none() {
        // all reads not-proper → nothing accepted → None
        let recs = vec![read(PAIRED | READ1, 0, 0, 300)];
        let bam = write_bam(&recs);
        let got = insert_size_distribution_params(bam.path(), 3, 2000, 8000, 1.0).unwrap();
        assert!(got.is_none());
    }
}
