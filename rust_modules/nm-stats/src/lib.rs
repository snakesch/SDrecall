//! `nm-stats` — Rust port of `realign_recall/cal_edge_NM_values.py`
//! (`calculate_NM_distribution_poisson`).
//!
//! Streams a BAM's records, keeps proper-pair / primary / non-dup / non-qcfail /
//! no-`XA` / MAPQ-60 reads, computes per-read `NM − max_indel_gap` ("scattered
//! edit distance"), takes the mean, then linear-searches the smallest integer
//! `cutoff` with `Poisson(mean).cdf(cutoff) ≥ 1 − conf_level`, floored at 4.
//!
//! This crate is the **template** every later lib+bin stage copies: a streaming
//! record scan with **zero per-record allocation**, one versatile filter+map
//! unit ([`passing_scatter_dist`]), a thin numeric tail ([`poisson_cutoff`]), and
//! a typed error via [`sdrecall_utils::SdError`].

use rust_htslib::bam::{self, record::Aux, record::Cigar, Read};
use sdrecall_utils::{Result, SdError};
use statrs::distribution::{DiscreteCDF, Poisson};
use std::path::Path;

/// Result of the NM-distribution estimate: `(cutoff, mean)`, mirroring the
/// Python return tuple. A plain struct keeps the bin / Python-dump format
/// (`cutoff\tmean`) obvious.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NmCutoff {
    pub cutoff: u64,
    pub mean: f64,
}

/// Read the `NM` (edit-distance) aux tag as an `i64`, accepting any integer
/// width (BWA/minimap may emit signed `i`). A missing or non-integer `NM` on a
/// read we accepted is a hard error — mirroring Python's fatal `KeyError`, and
/// per the dependency-availability rule (fail clearly, no skip-fallback).
fn nm_tag(rec: &bam::Record) -> Result<i64> {
    match rec.aux(b"NM") {
        Ok(Aux::U8(v)) => Ok(v as i64),
        Ok(Aux::U16(v)) => Ok(v as i64),
        Ok(Aux::U32(v)) => Ok(v as i64),
        Ok(Aux::I8(v)) => Ok(v as i64),
        Ok(Aux::I16(v)) => Ok(v as i64),
        Ok(Aux::I32(v)) => Ok(v as i64),
        Ok(other) => Err(SdError::Htslib(format!(
            "NM tag has a non-integer type: {other:?}"
        ))),
        Err(e) => Err(SdError::Htslib(format!("NM tag missing/unreadable: {e}"))),
    }
}

/// Largest single indel gap = `max(len)` over `I`/`D` CIGAR ops with `len > 1`
/// (Python op codes 1=Ins, 2=Del; the strict `> 1` is load-bearing — a `1I`/`1D`
/// does not count). `0` if there is no such op. Inlined into the one place CIGAR
/// is walked rather than exposed as a near-duplicate helper.
fn max_indel_gap(rec: &bam::Record) -> u32 {
    let mut max_gap = 0u32;
    for op in rec.cigar().iter() {
        if let Cigar::Ins(l) | Cigar::Del(l) = op {
            if *l > 1 && *l > max_gap {
                max_gap = *l;
            }
        }
    }
    max_gap
}

/// THE unit: the L15–26 filter + map for a single record, fused so the CIGAR/tag
/// access happens once.
///
/// Returns `Ok(None)` if the read is filtered out (normal), `Ok(Some(d))` with
/// `d = NM − max_indel_gap` if it passes, and `Err` if a passing read is missing
/// its `NM` tag.
///
/// `&Record`: read-only borrow, nothing escapes — so the driver can reuse one
/// `Record` buffer across the whole scan and there is never a clone in the hot
/// path. The cheap bit-flag / MAPQ checks run before the pricier `aux`/CIGAR
/// work so `||` short-circuits on the majority of rejected reads.
pub fn passing_scatter_dist(rec: &bam::Record) -> Result<Option<f64>> {
    if !rec.is_proper_pair()
        || rec.is_secondary()
        || rec.is_supplementary()
        || rec.is_duplicate()
        || rec.is_quality_check_failed()
        || rec.mapq() != 60
    {
        return Ok(None);
    }
    // "no XA tag present": htslib returns Err when the tag is absent.
    if rec.aux(b"XA").is_ok() {
        return Ok(None);
    }
    let nm = nm_tag(rec)?;
    let gap = max_indel_gap(rec) as i64;
    Ok(Some((nm - gap) as f64))
}

/// The pure numeric tail: smallest integer `cutoff` with
/// `Poisson(mean).cdf(cutoff) ≥ 1 − conf_level`, floored at 4 (Python L44).
///
/// Split out only because it is independently unit-testable against hand-computed
/// CDF thresholds and has zero I/O. The `max(_, 4)` floor lives here (one place).
///
/// `mean ≤ 0` is the degenerate case: scipy's `poisson.cdf(0, 0) = 1.0` makes the
/// raw cutoff `0`, which the floor lifts to `4`; statrs rejects `lambda ≤ 0`, so
/// we short-circuit to the same answer.
pub fn poisson_cutoff(mean: f64, conf_level: f64) -> Result<u64> {
    // Validate inputs before the search loop: a `conf_level` outside [0, 1] makes
    // `target = 1 - conf_level` exceed 1.0, which the Poisson CDF never reaches →
    // the `while` below would loop forever. A non-finite mean is likewise rejected.
    if !(0.0..=1.0).contains(&conf_level) {
        return Err(SdError::Compute(format!(
            "poisson_cutoff: conf_level {conf_level} is out of range [0, 1]"
        )));
    }
    if !mean.is_finite() {
        return Err(SdError::Compute(format!(
            "poisson_cutoff: mean {mean} is not finite"
        )));
    }
    if mean <= 0.0 {
        return Ok(4);
    }
    let pois = Poisson::new(mean)
        .map_err(|e| SdError::Compute(format!("Poisson::new(mean={mean}): {e}")))?;
    let target = 1.0 - conf_level;
    let mut cutoff: u64 = 0;
    while pois.cdf(cutoff) < target {
        cutoff += 1;
    }
    Ok(cutoff.max(4))
}

/// Driver / orchestration: stream the BAM, collect up to `sample_size` scattered
/// edit distances, take the mean, and resolve the Poisson cutoff.
///
/// `&Path` (open-only, no ownership); scalars by value. `threads` is forwarded to
/// htslib's BGZF decompression pool (the scan is I/O-bound; record-level rayon
/// would not help and would break the early `sample_size` break).
pub fn nm_distribution_poisson(
    bam_path: &Path,
    conf_level: f64,
    sample_size: usize,
    threads: u8,
) -> Result<NmCutoff> {
    let mut reader = bam::Reader::from_path(bam_path)
        .map_err(|e| SdError::Htslib(format!("open {}: {e}", bam_path.display())))?;
    // Best-effort decompression threads; single-threaded htslib default is fine.
    let _ = reader.set_threads(threads as usize);

    let mut rec = bam::Record::new(); // ONE buffer, reused every iteration
    let mut acc = Vec::<f64>::with_capacity(sample_size.min(1 << 20));
    let mut n = 0usize;

    while let Some(read_result) = reader.read(&mut rec) {
        read_result.map_err(|e| SdError::Htslib(format!("read record: {e}")))?;
        if let Some(d) = passing_scatter_dist(&rec)? {
            acc.push(d);
            n += 1;
            if n >= sample_size {
                break;
            }
        }
    }

    if n == 0 {
        // Python computes np.mean([]) = NaN; we fail loudly instead.
        return Err(SdError::InsufficientPairs(0));
    }
    if n < sample_size {
        log::warn!(
            "BAM {} only had {n} usable reads (< sample_size {sample_size})",
            bam_path.display()
        );
    }

    let mean = acc.iter().sum::<f64>() / n as f64;
    let cutoff = poisson_cutoff(mean, conf_level)?;
    log::info!(
        "NM distribution for {}: mean={mean}, {}-percentile cutoff={cutoff}",
        bam_path.display(),
        1.0 - conf_level
    );
    Ok(NmCutoff { cutoff, mean })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::header::HeaderRecord;
    use rust_htslib::bam::record::CigarString;
    use rust_htslib::bam::{Format, Header, Writer};

    // SAM flag bits.
    const PAIRED: u16 = 0x1;
    const PROPER_PAIR: u16 = 0x2;
    const SECONDARY: u16 = 0x100;
    const QCFAIL: u16 = 0x200;
    const DUP: u16 = 0x400;
    const SUPPLEMENTARY: u16 = 0x800;
    const PASS_FLAGS: u16 = PAIRED | PROPER_PAIR;

    /// Build a synthetic record with the given flags / MAPQ / CIGAR / NM / XA.
    fn make_read(
        flags: u16,
        mapq: u8,
        cigar: Vec<Cigar>,
        nm: Option<i32>,
        xa: bool,
    ) -> bam::Record {
        // query length = sum of query-consuming ops (M/=/X/I/S); D/N/H/P don't consume.
        let qlen: usize = cigar
            .iter()
            .map(|c| match c {
                Cigar::Match(l)
                | Cigar::Ins(l)
                | Cigar::SoftClip(l)
                | Cigar::Equal(l)
                | Cigar::Diff(l) => *l as usize,
                _ => 0,
            })
            .sum();
        let seq = vec![b'A'; qlen];
        let qual = vec![30u8; qlen];
        let cs = CigarString(cigar);

        let mut rec = bam::Record::new();
        rec.set(b"read1", Some(&cs), &seq, &qual);
        rec.set_flags(flags);
        rec.set_mapq(mapq);
        rec.set_tid(0);
        rec.set_pos(1000);
        if let Some(v) = nm {
            rec.push_aux(b"NM", Aux::I32(v)).unwrap();
        }
        if xa {
            rec.push_aux(b"XA", Aux::String("chr1,+100,50M,0;"))
                .unwrap();
        }
        rec
    }

    // ── passing_scatter_dist ────────────────────────────────────────────────

    #[test]
    fn pass_computes_nm_minus_max_gap() {
        // 10=2I3=4D5= → max indel gap = max(2,4) = 4; NM 5 → 5 - 4 = 1.0
        let rec = make_read(
            PASS_FLAGS,
            60,
            vec![
                Cigar::Equal(10),
                Cigar::Ins(2),
                Cigar::Equal(3),
                Cigar::Del(4),
                Cigar::Equal(5),
            ],
            Some(5),
            false,
        );
        assert_eq!(passing_scatter_dist(&rec).unwrap(), Some(1.0));
    }

    #[test]
    fn indel_length_one_does_not_count_as_gap() {
        // 1I/1D have len == 1 (strict > 1 excludes them) → gap 0; NM 3 → 3.0
        let rec = make_read(
            PASS_FLAGS,
            60,
            vec![
                Cigar::Equal(10),
                Cigar::Ins(1),
                Cigar::Equal(5),
                Cigar::Del(1),
                Cigar::Equal(5),
            ],
            Some(3),
            false,
        );
        assert_eq!(passing_scatter_dist(&rec).unwrap(), Some(3.0));
    }

    #[test]
    fn each_filter_rejects() {
        let cig = || vec![Cigar::Equal(10)];
        // not proper pair (only PAIRED set)
        assert_eq!(
            passing_scatter_dist(&make_read(PAIRED, 60, cig(), Some(2), false)).unwrap(),
            None
        );
        // secondary / supplementary / duplicate / qcfail
        for extra in [SECONDARY, SUPPLEMENTARY, DUP, QCFAIL] {
            assert_eq!(
                passing_scatter_dist(&make_read(PASS_FLAGS | extra, 60, cig(), Some(2), false))
                    .unwrap(),
                None
            );
        }
        // MAPQ != 60
        assert_eq!(
            passing_scatter_dist(&make_read(PASS_FLAGS, 59, cig(), Some(2), false)).unwrap(),
            None
        );
        // has XA tag
        assert_eq!(
            passing_scatter_dist(&make_read(PASS_FLAGS, 60, cig(), Some(2), true)).unwrap(),
            None
        );
    }

    #[test]
    fn missing_nm_on_passing_read_is_error() {
        let rec = make_read(PASS_FLAGS, 60, vec![Cigar::Equal(10)], None, false);
        let err = passing_scatter_dist(&rec).unwrap_err();
        assert!(matches!(err, SdError::Htslib(_)), "got {err:?}");
    }

    // ── poisson_cutoff ──────────────────────────────────────────────────────

    #[test]
    fn poisson_cutoff_mean_one() {
        // mean 1, conf 0.01: cdf(3)=0.9810 < 0.99, cdf(4)=0.9963 ≥ 0.99 → 4
        assert_eq!(poisson_cutoff(1.0, 0.01).unwrap(), 4);
    }

    #[test]
    fn poisson_cutoff_floor_at_four() {
        // small mean → raw cutoff 1, floored up to 4
        assert_eq!(poisson_cutoff(0.05, 0.01).unwrap(), 4);
        // degenerate mean 0 → 4 (matches scipy poisson.cdf(0,0)=1 then floor)
        assert_eq!(poisson_cutoff(0.0, 0.01).unwrap(), 4);
    }

    #[test]
    fn poisson_cutoff_large_mean_not_lowered_by_floor() {
        // mean 8, conf 0.001 → cutoff well above 4 (floor must not reduce it)
        let c = poisson_cutoff(8.0, 0.001).unwrap();
        assert!(c > 4, "expected cutoff > 4, got {c}");
    }

    #[test]
    fn poisson_cutoff_rejects_bad_conf_level() {
        // Out-of-range conf_level would make `target > 1` and loop forever.
        assert!(poisson_cutoff(5.0, -0.1).is_err());
        assert!(poisson_cutoff(5.0, 1.5).is_err());
        // The boundaries are valid.
        assert!(poisson_cutoff(5.0, 0.0).is_ok());
        assert!(poisson_cutoff(5.0, 1.0).is_ok());
    }

    #[test]
    fn poisson_cutoff_rejects_nonfinite_mean() {
        assert!(poisson_cutoff(f64::NAN, 0.01).is_err());
        assert!(poisson_cutoff(f64::INFINITY, 0.01).is_err());
    }

    // ── nm_distribution_poisson (end-to-end on a synthetic BAM) ──────────────

    fn write_bam(path: &Path, records: &[bam::Record]) {
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 100_000);
        let mut header = Header::new();
        header.push_record(&sq);
        let mut w = Writer::from_path(path, &header, Format::Bam).unwrap();
        for r in records {
            w.write(r).unwrap();
        }
    }

    #[test]
    fn driver_end_to_end() {
        let scatter1 = vec![
            Cigar::Equal(10),
            Cigar::Ins(2),
            Cigar::Equal(3),
            Cigar::Del(4),
            Cigar::Equal(5),
        ];
        let records = vec![
            make_read(PASS_FLAGS, 60, scatter1.clone(), Some(5), false), // pass → 1.0
            make_read(PASS_FLAGS, 60, scatter1.clone(), Some(5), false), // pass → 1.0
            make_read(PASS_FLAGS, 60, scatter1.clone(), Some(5), false), // pass → 1.0
            make_read(PASS_FLAGS, 30, vec![Cigar::Equal(10)], Some(2), false), // filtered (mapq)
        ];
        let tmp = tempfile::Builder::new().suffix(".bam").tempfile().unwrap();
        write_bam(tmp.path(), &records);

        let out = nm_distribution_poisson(tmp.path(), 0.01, 1_000_000, 1).unwrap();
        // 3 passing reads each with scatter 1.0 → mean 1.0 → cutoff 4
        assert_eq!(
            out,
            NmCutoff {
                cutoff: 4,
                mean: 1.0
            }
        );
    }

    #[test]
    fn driver_errors_on_no_usable_reads() {
        let records = vec![make_read(
            PAIRED,
            60,
            vec![Cigar::Equal(10)],
            Some(2),
            false,
        )];
        let tmp = tempfile::Builder::new().suffix(".bam").tempfile().unwrap();
        write_bam(tmp.path(), &records);
        let err = nm_distribution_poisson(tmp.path(), 0.01, 1_000_000, 1).unwrap_err();
        assert!(matches!(err, SdError::InsufficientPairs(0)), "got {err:?}");
    }
}
