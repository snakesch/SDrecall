//! Priority-merge path — ports `merge_variants_with_priority.py::merge_with_priority`.
//!
//! The matched-pair callback ([`MergeOp`]) tags the kept (reference-derived) record
//! with the source tags and GT-corrects it; the three per-set finalizers
//! (`finalize_matched` / `finalize_query_only` / `finalize_ref_only`) replicate the
//! Python write-time transforms (merge L536-657), each keyed by which set the
//! record landed in. The GT ladders are the shared [`crate::gt_rules`] tables.

use crate::coiterate::{coiterate_sorted_vcfs, CoiterSets, LocusOp};
use crate::gt_rules::{should_force_hom, MATCHED_PAIR_LADDER, QUERY_ONLY_LADDER, REF_ONLY_LADDER};
use crate::header_merge::build_merged_header;
use crate::norm::{norm_dedup_sort, sort_vcf};
use crate::record_io::{
    apply_filters, filter_names, reorder_filter_tag, sample_stats_default0, sample_stats_default1,
    set_gt_hom_alt,
};
use crate::vcf_group::{group_by_contig, read_translate_group};
use rust_htslib::bcf;
use sdrecall_utils::{Result, SdError};
use std::path::Path;

/// The `main_contigs` set (merge L513): `chr1..22,X,Y,M` plus the no-`chr` aliases.
/// Keep SEPARATE from the inhouse regex (design §6 — do NOT unify).
pub fn is_main_contig(name: &str) -> bool {
    let stripped = name.strip_prefix("chr").unwrap_or(name);
    match stripped {
        "X" | "Y" | "M" | "MT" => true,
        n => n
            .parse::<u8>()
            .map(|v| (1..=22).contains(&v))
            .unwrap_or(false),
    }
}

/// Parameters for [`merge_with_priority`] (borrowed paths, scalars by value).
pub struct MergeParams<'a> {
    pub query_vcf: &'a Path,
    pub reference_vcf: &'a Path,
    pub output_vcf: &'a Path,
    pub ref_genome: &'a Path,
    /// `added_filter` (Python "MISALIGNED") tagged onto query-only records.
    pub added_filter: Option<&'a str>,
    /// `qv_tag` (Python "RAW") — source tag for query.
    pub qv_tag: Option<&'a str>,
    /// `rv_tag` (Python "CLEAN") — source tag for reference.
    pub rv_tag: Option<&'a str>,
    /// Python `modify_gt`; the production call passes `False`.
    pub modify_gt: bool,
    pub threads: u8,
    pub tmp_dir: &'a Path,
}

/// The matched-pair `LocusOp` (Python `merge_same_variant_rec`, L115-135). Adds
/// `qv_tag`/`rv_tag` to the REFERENCE record's FILTER (Python adds them via
/// `filter.add`, order fixed up later by the finalizer), GT-corrects it from the
/// reference record's own AD/GQ/HPSUP when `modify_gt`, and keeps the reference
/// record.
struct MergeOp<'a> {
    qv_tag: Option<&'a str>,
    rv_tag: Option<&'a str>,
    modify_gt: bool,
}

impl LocusOp for MergeOp<'_> {
    fn on_match(&self, _query: &mut bcf::Record, refr: &mut bcf::Record) -> Option<bcf::Record> {
        // Add source tags to the kept (reference) record. `filter.add` is a set-add;
        // the finalizer re-orders so the tags end up last regardless.
        let mut names = filter_names(refr);
        for tag in [self.qv_tag, self.rv_tag].into_iter().flatten() {
            if !names.iter().any(|n| n == tag.as_bytes()) {
                names.push(tag.as_bytes().to_vec());
            }
        }
        if apply_filters(refr, &names).is_err() {
            return Some(refr.clone());
        }

        // GT correction (Python `modify_gt_based_on_ad_gq`) reads the reference
        // record's stats; only when GT is not already (1,1) (the Python `continue`).
        if self.modify_gt && !is_hom_alt(refr) {
            let stats = sample_stats_default1(refr, 0);
            if should_force_hom(&stats, MATCHED_PAIR_LADDER) {
                let _ = set_gt_hom_alt(refr);
            }
        }
        Some(refr.clone())
    }
}

/// Whether sample 0's GT is `(1,1)` (Python `samples[s]['GT'] == (1,1)`).
fn is_hom_alt(rec: &bcf::Record) -> bool {
    use rust_htslib::bcf::record::GenotypeAllele;
    match rec.genotypes() {
        Ok(gts) => {
            if rec.sample_count() == 0 {
                return false;
            }
            let gt = gts.get(0);
            let alleles: &[GenotypeAllele] = &gt;
            alleles.len() == 2
                && alleles
                    .iter()
                    .all(|a| matches!(a, GenotypeAllele::Unphased(1) | GenotypeAllele::Phased(1)))
        }
        Err(_) => false,
    }
}

/// Orchestrate the full priority merge (Python `merge_with_priority`, files in →
/// file out). Sorts both inputs (leaf bcftools), co-iterates per contig, applies
/// the three finalizers, writes, then `bcftools norm -d exact | sort`.
pub fn merge_with_priority(p: MergeParams<'_>) -> Result<()> {
    // 1. Sort/normalize both inputs (the one external leaf).
    let sorted_query = tmp_path(p.tmp_dir, "vcfops.query.sorted.vcf.gz");
    let sorted_ref = tmp_path(p.tmp_dir, "vcfops.ref.sorted.vcf.gz");
    sort_vcf(p.query_vcf, p.ref_genome, &sorted_query, p.threads)?;
    sort_vcf(p.reference_vcf, p.ref_genome, &sorted_ref, p.threads)?;

    // 2. Build the merged output header (ref template + query FILTER + new tags).
    let query_reader = bcf::Reader::from_path(&sorted_query)
        .map_err(|e| SdError::Vcf(format!("open {}: {e}", sorted_query.display())))?;
    let ref_reader = bcf::Reader::from_path(&sorted_ref)
        .map_err(|e| SdError::Vcf(format!("open {}: {e}", sorted_ref.display())))?;
    let extra: Vec<&str> = [p.added_filter, p.qv_tag, p.rv_tag]
        .into_iter()
        .flatten()
        .collect();
    let out_header = build_merged_header(&ref_reader, &query_reader, &extra);
    drop(query_reader);
    drop(ref_reader);

    // 3. Writer over a temp file (its HeaderView is the translation target).
    let tmp_out = tmp_path(p.tmp_dir, "vcfops.merged.tmp.vcf.gz");
    let mut writer = bcf::Writer::from_path(&tmp_out, &out_header, false, bcf::Format::Vcf)
        .map_err(|e| SdError::Vcf(format!("create {}: {e}", tmp_out.display())))?;

    // 4. Read + translate all records into the output header, group by contig.
    let q_recs = read_translate_group(&sorted_query, &mut writer)?;
    let r_recs = read_translate_group(&sorted_ref, &mut writer)?;
    let q_by_contig = group_by_contig(q_recs);
    let r_by_contig = group_by_contig(r_recs);

    // Contig order: union, restricted to main contigs, in sorted rid order.
    let mut rids: Vec<u32> = q_by_contig
        .keys()
        .chain(r_by_contig.keys())
        .copied()
        .collect();
    rids.sort_unstable();
    rids.dedup();

    let op = MergeOp {
        qv_tag: p.qv_tag,
        rv_tag: p.rv_tag,
        modify_gt: p.modify_gt,
    };

    // 5. Per-contig co-iteration + finalizers (sequential; see DESIGN note in
    //    lib.rs on why rayon-over-contigs is deferred for Rc<HeaderView> soundness).
    for rid in rids {
        let name = writer.header().rid2name(rid).ok();
        let name = match name {
            Some(n) => String::from_utf8_lossy(n).to_string(),
            None => continue,
        };
        if !is_main_contig(&name) {
            continue;
        }
        let q = q_by_contig.get(&rid).cloned().unwrap_or_default();
        let r = r_by_contig.get(&rid).cloned().unwrap_or_default();
        let sets = coiterate_sorted_vcfs(q, r, &op)?;
        write_merge_sets(&mut writer, sets, &p)?;
    }
    drop(writer);

    // 6. Final dedup + sort (Python L671-676).
    norm_dedup_sort(&tmp_out, p.output_vcf, p.threads)?;
    let _ = std::fs::remove_file(&tmp_out);
    let _ = std::fs::remove_file(&sorted_query);
    let _ = std::fs::remove_file(&sorted_ref);
    Ok(())
}

/// Apply the three per-set finalizers and write each record (merge L536-657).
fn write_merge_sets(writer: &mut bcf::Writer, sets: CoiterSets, p: &MergeParams<'_>) -> Result<()> {
    let CoiterSets {
        matched,
        mut query_only,
        mut ref_only,
    } = sets;

    // matched (L536-556): reorder rv_tag then qv_tag to the end; write as-is.
    for mut rec in matched {
        finalize_filter_order(&mut rec, &[p.rv_tag, p.qv_tag])?;
        writer.write(&rec).map_err(werr)?;
    }
    // query-only (L557-612): qv_tag then added_filter to end; query-only GT ladder.
    for rec in query_only.iter_mut() {
        finalize_filter_order(rec, &[p.qv_tag, p.added_filter])?;
        if p.modify_gt && !is_hom_alt(rec) {
            let stats = sample_stats_default0(rec, 0);
            if should_force_hom(&stats, QUERY_ONLY_LADDER) {
                set_gt_hom_alt(rec)?;
            }
        }
    }
    for rec in query_only {
        writer.write(&rec).map_err(werr)?;
    }
    // ref-only (L613-657): rv_tag to end; ref-only GT ladder.
    for rec in ref_only.iter_mut() {
        finalize_filter_order(rec, &[p.rv_tag])?;
        if p.modify_gt && !is_hom_alt(rec) {
            let stats = sample_stats_default1(rec, 0);
            if should_force_hom(&stats, REF_ONLY_LADDER) {
                set_gt_hom_alt(rec)?;
            }
        }
    }
    for rec in ref_only {
        writer.write(&rec).map_err(werr)?;
    }
    Ok(())
}

/// Reorder the given source/priority tags to the END of the record's FILTER list,
/// in the order given (Python `[f for f in filters if f != tag] + [tag]` applied
/// successively, merge L538-541/L559-562). `None` tags are skipped.
fn finalize_filter_order(rec: &mut bcf::Record, tags: &[Option<&str>]) -> Result<()> {
    let mut names = filter_names(rec);
    for tag in tags.iter().flatten() {
        names = reorder_filter_tag(&names, tag);
    }
    apply_filters(rec, &names)
}

fn werr(e: rust_htslib::errors::Error) -> SdError {
    SdError::Vcf(format!("write record: {e}"))
}

/// A process-unique temp path under `tmp_dir`: `{pid}.{seq}.{name}`. The pid
/// disambiguates concurrent SDrecall processes sharing a tmp dir; the atomic
/// `seq` disambiguates concurrent calls WITHIN one process, so a future
/// rayon-over-contigs cannot collide on a fixed `name`. The caller cleans up.
fn tmp_path(tmp_dir: &Path, name: &str) -> std::path::PathBuf {
    use std::sync::atomic::{AtomicU64, Ordering};
    static SEQ: AtomicU64 = AtomicU64::new(0);
    let seq = SEQ.fetch_add(1, Ordering::Relaxed);
    tmp_dir.join(format!("{}.{}.{}", std::process::id(), seq, name))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bcf::record::GenotypeAllele;
    use rust_htslib::bcf::{Format, Writer};

    #[test]
    fn main_contig_set() {
        for c in [
            "chr1", "chr22", "chrX", "chrY", "chrM", "1", "22", "X", "Y", "MT",
        ] {
            assert!(is_main_contig(c), "{c} should be a main contig");
        }
        for c in [
            "chr23",
            "chrUn_KI270302v1",
            "chr1_KI270706v1_random",
            "GL000220.1",
            "0",
            "chr0",
        ] {
            assert!(!is_main_contig(c), "{c} should NOT be a main contig");
        }
    }

    /// Header with two FILTER tags (RG0 + the source tags) so we can test the
    /// "move source tag to END" reordering parity.
    fn writer_with_filters() -> (Writer, tempfile::TempPath) {
        let tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        let mut header = bcf::Header::new();
        header.push_record(b"##contig=<ID=chr1,length=1000000>");
        header.push_record(b"##FILTER=<ID=RG0,Description=\"rg0\">");
        header.push_record(b"##FILTER=<ID=RAW,Description=\"raw\">");
        header.push_record(b"##FILTER=<ID=CLEAN,Description=\"clean\">");
        header.push_record(b"##FILTER=<ID=MISALIGNED,Description=\"mis\">");
        header.push_record(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        header.push_sample(b"s1");
        let w = Writer::from_path(tmp.path(), &header, true, Format::Vcf).unwrap();
        (w, tmp.into_temp_path())
    }

    fn rec_with_filters(w: &Writer, filters: &[&str]) -> bcf::Record {
        let mut r = w.empty_record();
        let rid = w.header().name2rid(b"chr1").unwrap();
        r.set_rid(Some(rid));
        r.set_pos(100);
        r.set_alleles(&[b"A", b"T"]).unwrap();
        r.push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
            .unwrap();
        let ids: Vec<&[u8]> = filters.iter().map(|f| f.as_bytes()).collect();
        r.set_filters(&ids).unwrap();
        r
    }

    #[test]
    fn filter_tag_reordered_to_end() {
        // Start with [RG0, RAW]; finalize with [CLEAN, RAW] → RAW already present
        // gets moved to the end after CLEAN, and CLEAN appended:
        // [RG0] + CLEAN + RAW = [RG0, CLEAN, RAW].
        let (w, _t) = writer_with_filters();
        let mut r = rec_with_filters(&w, &["RG0", "RAW"]);
        finalize_filter_order(&mut r, &[Some("CLEAN"), Some("RAW")]).unwrap();
        let names: Vec<String> = filter_names(&r)
            .iter()
            .map(|n| String::from_utf8_lossy(n).to_string())
            .collect();
        assert_eq!(names, vec!["RG0", "CLEAN", "RAW"]);
    }

    #[test]
    fn filter_tag_appended_when_absent() {
        // [RG0] finalized with [RAW] → [RG0, RAW] (RAW appended at end).
        let (w, _t) = writer_with_filters();
        let mut r = rec_with_filters(&w, &["RG0"]);
        finalize_filter_order(&mut r, &[Some("RAW")]).unwrap();
        let names: Vec<String> = filter_names(&r)
            .iter()
            .map(|n| String::from_utf8_lossy(n).to_string())
            .collect();
        assert_eq!(names, vec!["RG0", "RAW"]);
    }
}
