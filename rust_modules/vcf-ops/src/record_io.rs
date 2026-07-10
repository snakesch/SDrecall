//! Shared FORMAT/INFO readers + FILTER-tag ordering — the plumbing both callbacks
//! and the finalizers reuse (no duplicated per-field extraction).
//!
//! All readers apply the **exact Python missing-value defaults**: `AD → [0,0]`,
//! `GQ → 0`, `HPSUP → "."` (so `num_hps` defaults to 1 via `".".split(";")`), and
//! for the query-only finalizer a genuinely-absent HPSUP yields `num_hps = 0`
//! (Python L573 `len(...) if len(hps)>0 else 0`).

use crate::gt_rules::SampleStats;
use rust_htslib::bcf::{self, record::GenotypeAllele};
use sdrecall_utils::{Result, SdError};

/// Read `AD` for `sample` as `(ref_dp, alt_dp)`, Python default `[0, 0]` when the
/// field is missing. A record with >2 AD values keeps the first two (Python
/// L580-584 "not biallelic" branch); after `bcftools norm -m -both` records are
/// split so AD is `[ref, alt]`, but we stay defensive.
pub fn ad_ref_alt(rec: &bcf::Record, sample: usize) -> (i32, i32) {
    match rec.format(b"AD").integer() {
        Ok(buf) => {
            let per_sample = buf.get(sample);
            match per_sample {
                Some(ad) if ad.len() >= 2 => (clamp_missing(ad[0]), clamp_missing(ad[1])),
                Some(ad) if ad.len() == 1 => (clamp_missing(ad[0]), 0),
                _ => (0, 0),
            }
        }
        Err(_) => (0, 0),
    }
}

/// Read `GQ` for `sample`, Python default `0` when missing (also None→0).
pub fn gq(rec: &bcf::Record, sample: usize) -> i32 {
    match rec.format(b"GQ").integer() {
        Ok(buf) => match buf.get(sample) {
            Some(g) if !g.is_empty() => clamp_missing(g[0]),
            _ => 0,
        },
        Err(_) => 0,
    }
}

/// `len(HPSUP[0].split(";"))` for `sample`. Returns:
/// - `Some(n)` with `n = count of ';'-separated tokens` when HPSUP is present,
/// - `None` when the field is absent (the caller decides the default: the
///   matched/ref ladders use Python `".".split(";")` = 1; the query-only finalizer
///   uses 0).
pub fn num_hps_opt(rec: &bcf::Record, sample: usize) -> Option<usize> {
    let buf = rec.format(b"HPSUP").string().ok()?;
    let per_sample = buf.get(sample)?;
    // htslib returns one byte-slice per sample; the value itself may carry the
    // ';'-joined haplotype list. Python's `hps[0]` is that string; split on ';'.
    if per_sample.is_empty() {
        return None;
    }
    Some(per_sample.split(|&c| c == b';').count())
}

/// Build [`SampleStats`] for `sample` with the matched/ref-ladder HPSUP default of
/// `1` (Python `".".split(";")` length). Used by [`crate::priority_merge`]'s
/// matched-pair and ref-only paths.
pub fn sample_stats_default1(rec: &bcf::Record, sample: usize) -> SampleStats {
    let (ref_dp, alt_dp) = ad_ref_alt(rec, sample);
    SampleStats {
        ref_dp,
        alt_dp,
        gq: gq(rec, sample),
        num_hps: num_hps_opt(rec, sample).unwrap_or(1),
    }
}

/// Build [`SampleStats`] with the query-only finalizer's HPSUP default of `0`
/// (Python L573).
pub fn sample_stats_default0(rec: &bcf::Record, sample: usize) -> SampleStats {
    let (ref_dp, alt_dp) = ad_ref_alt(rec, sample);
    SampleStats {
        ref_dp,
        alt_dp,
        gq: gq(rec, sample),
        num_hps: num_hps_opt(rec, sample).unwrap_or(0),
    }
}

/// htslib uses INT32_MIN-ish sentinels for missing scalar integers. Map any value
/// `< -1` (the htslib `bcf_int32_missing`/`vector_end` region) to 0 — matching the
/// Python `0 if x is None else x` guard. Real depths/GQ are non-negative.
#[inline]
fn clamp_missing(v: i32) -> i32 {
    if v < 0 {
        0
    } else {
        v
    }
}

/// Read the `AC[0]` / `AN` INFO integers (Python L46 takes `AC[0]`). Returns
/// `None` if either is absent or NA — mirroring `determine_common_per_pysam_record`
/// returning `False` (not common) on missing AC/AN (L40-53).
pub fn ac_an(rec: &bcf::Record) -> Option<(u64, u64)> {
    let ac = {
        let buf = rec.info(b"AC").integer().ok()??;
        let v = *buf.first()?;
        if v < 0 {
            return None;
        }
        v as u64
    };
    let an = {
        let buf = rec.info(b"AN").integer().ok()??;
        let v = *buf.first()?;
        if v < 0 {
            return None;
        }
        v as u64
    };
    Some((ac, an))
}

/// Whether the record already carries the `SDrecall` FILTER (Python
/// `"SDrecall" in qrecord.filter`, inhouse L108).
pub fn has_filter(rec: &bcf::Record, tag: &str) -> bool {
    rec.has_filter(tag.as_bytes())
}

/// Force `GT=(1,1)` (homozygous-alt) for the single SDrecall sample. SDrecall is
/// always single-sample diploid (design §6); `push_genotypes` overwrites sample 0.
pub fn set_gt_hom_alt(rec: &mut bcf::Record) -> Result<()> {
    rec.push_genotypes(&[GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)])
        .map_err(|e| SdError::Vcf(format!("set GT=(1,1): {e}")))
}

/// Read the current FILTER tag NAMES of a record, in stored order (the Python
/// `list(record.filter)` order). `PASS` (empty filter) yields an empty Vec.
pub fn filter_names(rec: &bcf::Record) -> Vec<Vec<u8>> {
    let hdr = rec.header();
    rec.filters().map(|id| hdr.id_to_name(id)).collect()
}

/// Replicate the Python FILTER reordering `[f for f in filters if f != tag] + [tag]`
/// (merge L538-541/L559-562): move `tag` to the END, append it if absent. Returns
/// the new ordered name list. `existing` are the current names (from
/// [`filter_names`]); `tag` is the source/priority tag to push last.
pub fn reorder_filter_tag(existing: &[Vec<u8>], tag: &str) -> Vec<Vec<u8>> {
    let tag_bytes = tag.as_bytes();
    let mut out: Vec<Vec<u8>> = existing
        .iter()
        .filter(|f| f.as_slice() != tag_bytes)
        .cloned()
        .collect();
    out.push(tag_bytes.to_vec());
    out
}

/// Apply an ordered FILTER name list to a record via one `set_filters` call,
/// preserving exact order (the parity-critical write — design §6 "set_filters
/// ordering parity"). An empty list sets `PASS`. Names absent from the header are
/// a hard error (the tags are added to the header before writing).
pub fn apply_filters(rec: &mut bcf::Record, names: &[Vec<u8>]) -> Result<()> {
    let id_refs: Vec<&[u8]> = names.iter().map(|n| n.as_slice()).collect();
    rec.set_filters(&id_refs)
        .map_err(|e| SdError::Vcf(format!("set_filters {names:?}: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bcf::record::GenotypeAllele;
    use rust_htslib::bcf::{Format, Writer};

    fn writer() -> (Writer, tempfile::TempPath) {
        let tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        let mut header = bcf::Header::new();
        header.push_record(b"##contig=<ID=chr1,length=1000000>");
        header.push_record(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        header
            .push_record(b"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">");
        header.push_record(b"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"GQ\">");
        header.push_record(b"##FORMAT=<ID=HPSUP,Number=A,Type=String,Description=\"HP support\">");
        header.push_sample(b"s1");
        let w = Writer::from_path(tmp.path(), &header, true, Format::Vcf).unwrap();
        (w, tmp.into_temp_path())
    }

    fn base_rec(w: &Writer) -> bcf::Record {
        let mut r = w.empty_record();
        let rid = w.header().name2rid(b"chr1").unwrap();
        r.set_rid(Some(rid));
        r.set_pos(100);
        r.set_alleles(&[b"A", b"T"]).unwrap();
        r.push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
            .unwrap();
        r
    }

    #[test]
    fn ad_gq_read_with_values() {
        let (w, _t) = writer();
        let mut r = base_rec(&w);
        r.push_format_integer(b"AD", &[3, 7]).unwrap();
        r.push_format_integer(b"GQ", &[42]).unwrap();
        assert_eq!(ad_ref_alt(&r, 0), (3, 7));
        assert_eq!(gq(&r, 0), 42);
    }

    #[test]
    fn ad_gq_default_when_absent() {
        // No AD/GQ → Python defaults [0,0] / 0.
        let (w, _t) = writer();
        let r = base_rec(&w);
        assert_eq!(ad_ref_alt(&r, 0), (0, 0));
        assert_eq!(gq(&r, 0), 0);
    }

    #[test]
    fn num_hps_counts_semicolon_tokens() {
        let (w, _t) = writer();
        let mut r = base_rec(&w);
        r.push_format_string(b"HPSUP", &[b"chunk1;chunk2;chunk3".as_slice()])
            .unwrap();
        assert_eq!(num_hps_opt(&r, 0), Some(3));
        // matched/ref default is 1 when absent; query-only default is 0.
        let r2 = base_rec(&w);
        assert_eq!(num_hps_opt(&r2, 0), None);
        assert_eq!(sample_stats_default1(&r2, 0).num_hps, 1);
        assert_eq!(sample_stats_default0(&r2, 0).num_hps, 0);
    }

    #[test]
    fn sample_stats_full() {
        let (w, _t) = writer();
        let mut r = base_rec(&w);
        r.push_format_integer(b"AD", &[4, 6]).unwrap();
        r.push_format_integer(b"GQ", &[20]).unwrap();
        r.push_format_string(b"HPSUP", &[b"a;b".as_slice()])
            .unwrap();
        let s = sample_stats_default1(&r, 0);
        assert_eq!((s.ref_dp, s.alt_dp, s.gq, s.num_hps), (4, 6, 20, 2));
        assert_eq!(s.alt_ratio(), 0.6);
    }

    #[test]
    fn reorder_filter_moves_tag_to_end() {
        // pure-Vec test of the ordering helper (no header needed).
        let existing = vec![b"RG0".to_vec(), b"RAW".to_vec()];
        let out = reorder_filter_tag(&existing, "RAW");
        assert_eq!(out, vec![b"RG0".to_vec(), b"RAW".to_vec()]);
        let out2 = reorder_filter_tag(&existing, "CLEAN");
        assert_eq!(
            out2,
            vec![b"RG0".to_vec(), b"RAW".to_vec(), b"CLEAN".to_vec()]
        );
    }
}
