//! Inhouse-common path — ports `identify_common_vars.py::annotate_inhouse_common`.
//!
//! The matched-pair callback ([`InhouseOp`]) runs the binomial test on the COHORT
//! record and, if the variant is common AND the QUERY record carries the
//! `SDrecall` filter, adds the `INHOUSE_COMMON` filter to the kept QUERY record
//! (inhouse `process_target_records` L100-111). Crucially the orchestrator writes
//! only `matched` + `query_only` — cohort-only records are dropped (inhouse
//! L426-469 has no cohort-only writer).

use crate::coiterate::{coiterate_sorted_vcfs, CoiterSets, LocusOp};
use crate::header_merge::build_query_header;
use crate::norm::sort_vcf;
use crate::record_io::{ac_an, apply_filters, filter_names, has_filter};
use rust_htslib::bcf;
use sdrecall_utils::{Result, SdError};
use statrs::distribution::{Binomial, DiscreteCDF};
use std::path::Path;

/// The inhouse contig filter: regex `^chr[0-9MTXY]+$` (inhouse L418). `chr`-prefixed
/// only — keep SEPARATE from the merge `main_contigs` set (design §6).
pub fn is_inhouse_contig(name: &str) -> bool {
    let rest = match name.strip_prefix("chr") {
        Some(r) => r,
        None => return false,
    };
    !rest.is_empty()
        && rest
            .bytes()
            .all(|b| b.is_ascii_digit() || matches!(b, b'M' | b'T' | b'X' | b'Y'))
}

/// Parameters for [`annotate_inhouse_common`].
pub struct InhouseParams<'a> {
    pub query_vcf: &'a Path,
    pub cohort_vcf: &'a Path,
    pub output_vcf: &'a Path,
    pub ref_genome: &'a Path,
    /// `added_filter` (Python default "INHOUSE_COMMON").
    pub added_filter: &'a str,
    pub inhouse_common_cutoff: f64,
    pub conf_level: f64,
    pub threads: u8,
    pub tmp_dir: &'a Path,
}

/// The binomial common-variant test (inhouse `determine_common_per_pysam_record`
/// L18-59): `binom.cdf(AC, AN, cutoff) > conf_level`. scipy `binom.cdf(k,n,p)` maps
/// to `Binomial::new(p, n).cdf(k)` — note the (p, n) arg order. Missing/NA AC or AN
/// → `false` (not common). `AN == 0` makes `Binomial::new` fail (n must be ≥ 0 with
/// the trivial dist at 0; scipy `binom.cdf(k,0,p)=1.0` for k≥0) — handled explicitly.
pub fn determine_common(rec: &bcf::Record, cutoff: f64, conf_level: f64) -> Result<bool> {
    let (ac, an) = match ac_an(rec) {
        Some(v) => v,
        None => return Ok(false), // missing AC/AN → not common
    };
    // scipy: binom.cdf(k, 0, p) = 1.0 for any k >= 0 → stat_power=1.0 > conf_level.
    if an == 0 {
        return Ok(1.0 > conf_level);
    }
    let dist = Binomial::new(cutoff, an)
        .map_err(|e| SdError::Vcf(format!("Binomial::new(p={cutoff}, n={an}): {e}")))?;
    let stat_power = dist.cdf(ac);
    Ok(stat_power > conf_level)
}

/// The matched-pair callback. `query` is mutated (the `INHOUSE_COMMON` filter is
/// added when common + has `SDrecall`); `cohort` is read-only data. Keeps the
/// QUERY record.
struct InhouseOp<'a> {
    added_filter: &'a str,
    cutoff: f64,
    conf_level: f64,
}

impl LocusOp for InhouseOp<'_> {
    fn on_match(&self, query: &mut bcf::Record, cohort: &mut bcf::Record) -> Option<bcf::Record> {
        let is_common = determine_common(cohort, self.cutoff, self.conf_level).unwrap_or(false);
        if is_common && has_filter(query, "SDrecall") {
            // filter.add(INHOUSE_COMMON): append if not already present.
            let mut names = filter_names(query);
            if !names.iter().any(|n| n == self.added_filter.as_bytes()) {
                names.push(self.added_filter.as_bytes().to_vec());
                let _ = apply_filters(query, &names);
            }
        }
        Some(query.clone())
    }
}

/// Orchestrate the inhouse-common annotation (Python `annotate_inhouse_common`).
/// Writes only matched + query-only records (cohort-only dropped), then a final
/// `sort_vcf` normalization of the output.
pub fn annotate_inhouse_common(p: InhouseParams<'_>) -> Result<()> {
    use crate::vcf_group::{group_by_contig, read_translate_group};

    let sorted_query = tmp_path(p.tmp_dir, "vcfops.ic.query.sorted.vcf.gz");
    let sorted_cohort = tmp_path(p.tmp_dir, "vcfops.ic.cohort.sorted.vcf.gz");
    sort_vcf(p.query_vcf, p.ref_genome, &sorted_query, p.threads)?;
    sort_vcf(p.cohort_vcf, p.ref_genome, &sorted_cohort, p.threads)?;

    // Output header = query header + INHOUSE_COMMON FILTER (inhouse L420 + L397).
    let query_reader = bcf::Reader::from_path(&sorted_query)
        .map_err(|e| SdError::Vcf(format!("open {}: {e}", sorted_query.display())))?;
    let out_header = build_query_header(&query_reader, &[p.added_filter]);
    drop(query_reader);

    let tmp_out = tmp_path(p.tmp_dir, "vcfops.ic.tmp.vcf.gz");
    let mut writer = bcf::Writer::from_path(&tmp_out, &out_header, false, bcf::Format::Vcf)
        .map_err(|e| SdError::Vcf(format!("create {}: {e}", tmp_out.display())))?;

    let q_recs = read_translate_group(&sorted_query, &mut writer)?;
    let c_recs = read_translate_group(&sorted_cohort, &mut writer)?;
    let q_by_contig = group_by_contig(q_recs);
    let c_by_contig = group_by_contig(c_recs);

    let mut rids: Vec<u32> = q_by_contig
        .keys()
        .chain(c_by_contig.keys())
        .copied()
        .collect();
    rids.sort_unstable();
    rids.dedup();

    let op = InhouseOp {
        added_filter: p.added_filter,
        cutoff: p.inhouse_common_cutoff,
        conf_level: p.conf_level,
    };

    for rid in rids {
        let name = match writer.header().rid2name(rid).ok() {
            Some(n) => String::from_utf8_lossy(n).to_string(),
            None => continue,
        };
        if !is_inhouse_contig(&name) {
            continue;
        }
        let q = q_by_contig.get(&rid).cloned().unwrap_or_default();
        let c = c_by_contig.get(&rid).cloned().unwrap_or_default();
        let CoiterSets {
            matched,
            query_only,
            ref_only: _cohort_only,
        } = coiterate_sorted_vcfs(q, c, &op)?;
        // Write matched + query-only ONLY (cohort-only dropped — inhouse L426-469).
        for rec in matched.iter().chain(query_only.iter()) {
            writer
                .write(rec)
                .map_err(|e| SdError::Vcf(format!("write record: {e}")))?;
        }
    }
    drop(writer);

    // Final sort_vcf of the output (inhouse L482-484).
    sort_vcf(&tmp_out, p.ref_genome, p.output_vcf, p.threads)?;
    let _ = std::fs::remove_file(&tmp_out);
    let _ = std::fs::remove_file(&sorted_query);
    let _ = std::fs::remove_file(&sorted_cohort);
    Ok(())
}

/// A process-unique temp path under `tmp_dir`: `{pid}.{seq}.{name}` (see the
/// twin in `priority_merge`). The atomic `seq` makes concurrent calls within one
/// process collision-free. The caller cleans up.
fn tmp_path(tmp_dir: &Path, name: &str) -> std::path::PathBuf {
    use std::sync::atomic::{AtomicU64, Ordering};
    static SEQ: AtomicU64 = AtomicU64::new(0);
    let seq = SEQ.fetch_add(1, Ordering::Relaxed);
    tmp_dir.join(format!("{}.{}.{}", std::process::id(), seq, name))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inhouse_contig_filter() {
        for c in ["chr1", "chr22", "chrX", "chrY", "chrM", "chrMT", "chr10"] {
            assert!(is_inhouse_contig(c), "{c} should match the inhouse regex");
        }
        // No `chr` prefix → excluded (unlike main_contigs).
        for c in [
            "1",
            "X",
            "MT",
            "chrUn_KI270302v1",
            "chr1_random",
            "chr",
            "GL000220.1",
        ] {
            assert!(
                !is_inhouse_contig(c),
                "{c} should NOT match the inhouse regex"
            );
        }
    }

    // ── binomial cutoff boundary (design §7) ────────────────────────────────
    //
    // scipy reference values (binom.cdf(k, n, p) for n=100, p=0.01), verified by
    // running `scipy.stats.binom.cdf` directly:
    //   cdf(4,100,0.01) = 0.99656768  → not > 0.999 → not common
    //   cdf(5,100,0.01) = 0.99946547  → > 0.999     → COMMON
    // The cdf crosses the 0.999 conf_level between k=4 and k=5. statrs matches
    // scipy to ~1e-8 (both use the regularized incomplete beta).

    #[test]
    fn binomial_cdf_crosses_0999_between_ac4_and_ac5_for_an100() {
        let dist = Binomial::new(0.01, 100).unwrap();
        let c4 = dist.cdf(4);
        let c5 = dist.cdf(5);
        assert!(c4 <= 0.999, "cdf(4)={c4} should be <= 0.999");
        assert!(c5 > 0.999, "cdf(5)={c5} should be > 0.999");
        // statrs == scipy to ~1e-7 at both boundary points
        assert!((c4 - 0.99656768).abs() < 1e-7, "cdf(4)={c4}");
        assert!((c5 - 0.99946547).abs() < 1e-7, "cdf(5)={c5}");
    }

    #[test]
    fn binomial_arg_order_is_p_then_n() {
        // If the (p, n) order were swapped, Binomial::new(100, 0.01) would error
        // (p must be in [0,1]); confirm the correct order constructs fine and an
        // obviously-common count is flagged.
        assert!(Binomial::new(0.01, 100).is_ok());
        assert!(Binomial::new(100.0, 1).is_err()); // p=100 invalid → proves order
    }

    // ── record-level determine_common (AC/AN from INFO) ──────────────────────

    use rust_htslib::bcf::record::GenotypeAllele;
    use rust_htslib::bcf::{Format, Writer};

    fn writer_with_ac_an() -> (Writer, tempfile::TempPath) {
        let tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        let mut header = bcf::Header::new();
        header.push_record(b"##contig=<ID=chr1,length=1000000>");
        header.push_record(b"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">");
        header.push_record(b"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">");
        header.push_record(b"##FILTER=<ID=SDrecall,Description=\"sdrecall\">");
        header.push_record(b"##FILTER=<ID=INHOUSE_COMMON,Description=\"common\">");
        header.push_record(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        header.push_sample(b"s1");
        let w = Writer::from_path(tmp.path(), &header, true, Format::Vcf).unwrap();
        (w, tmp.into_temp_path())
    }

    fn rec_ac_an(w: &Writer, ac: i32, an: i32) -> bcf::Record {
        let mut r = w.empty_record();
        let rid = w.header().name2rid(b"chr1").unwrap();
        r.set_rid(Some(rid));
        r.set_pos(100);
        r.set_alleles(&[b"A", b"T"]).unwrap();
        r.push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
            .unwrap();
        r.push_info_integer(b"AC", &[ac]).unwrap();
        r.push_info_integer(b"AN", &[an]).unwrap();
        r
    }

    #[test]
    fn determine_common_record_boundary() {
        let (w, _t) = writer_with_ac_an();
        // AN=100, cutoff=0.01, conf=0.999: AC=4 not common, AC=5 common.
        let r4 = rec_ac_an(&w, 4, 100);
        let r5 = rec_ac_an(&w, 5, 100);
        assert!(!determine_common(&r4, 0.01, 0.999).unwrap());
        assert!(determine_common(&r5, 0.01, 0.999).unwrap());
    }

    #[test]
    fn determine_common_missing_ac_an_is_false() {
        // No AC/AN INFO → not common (Python returns False on missing).
        let (w, _t) = writer_with_ac_an();
        let mut r = w.empty_record();
        let rid = w.header().name2rid(b"chr1").unwrap();
        r.set_rid(Some(rid));
        r.set_pos(100);
        r.set_alleles(&[b"A", b"T"]).unwrap();
        r.push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
            .unwrap();
        assert!(!determine_common(&r, 0.01, 0.999).unwrap());
    }

    #[test]
    fn inhouse_op_requires_sdrecall_filter() {
        // common cohort variant (AC=10,AN=100) but query lacks SDrecall → no tag;
        // with SDrecall → INHOUSE_COMMON added.
        let (w, _t) = writer_with_ac_an();
        let op = InhouseOp {
            added_filter: "INHOUSE_COMMON",
            cutoff: 0.01,
            conf_level: 0.999,
        };

        // query without SDrecall filter
        let mut q_plain = rec_ac_an(&w, 0, 2);
        let mut cohort = rec_ac_an(&w, 10, 100);
        op.on_match(&mut q_plain, &mut cohort);
        assert!(
            !has_filter(&q_plain, "INHOUSE_COMMON"),
            "no tag without SDrecall"
        );

        // query WITH SDrecall filter
        let mut q_sd = rec_ac_an(&w, 0, 2);
        q_sd.set_filters(&[b"SDrecall".as_slice()]).unwrap();
        let mut cohort2 = rec_ac_an(&w, 10, 100);
        op.on_match(&mut q_sd, &mut cohort2);
        assert!(
            has_filter(&q_sd, "INHOUSE_COMMON"),
            "tag added with SDrecall + common"
        );
    }
}
