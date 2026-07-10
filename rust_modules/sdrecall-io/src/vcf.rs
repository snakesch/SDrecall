//! VCF/BCF read + sorted-merge — the ONE vcf unit.
//!
//! Ports `src/utils.py::combine_vcfs` (`bcftools concat -a -d exact | sort`).
//! `concat_sort_vcfs` does a k-way merge over already coordinate-sorted inputs
//! (each island VCF is sorted) and exact-dedups on `(CHROM, POS, REF, ALT)` — a
//! streaming O(N) merge rather than load-all-then-sort.
//!
//! Note on signatures: `rust_htslib::bcf::Writer::from_path` takes an owned
//! `bcf::Header` (not a `&HeaderView`), so [`write_vcf`] takes `&bcf::Header`
//! — a small deviation from the DESIGN's `&HeaderView`, recorded here. Records
//! read from an input are bound to that input's header, so they are written via
//! `Writer::translate` into the output header.

use rust_htslib::bcf::{self, Read};
use sdrecall_utils::{Result, SdError};
use std::cmp::Ordering;
use std::path::{Path, PathBuf};

/// Open a VCF/BCF reader (a sorted cursor; the caller streams records). The
/// inputs are assumed coordinate-sorted (as every SDrecall island VCF is).
pub fn read_vcf(path: &Path) -> Result<bcf::Reader> {
    bcf::Reader::from_path(path).map_err(|e| SdError::Vcf(format!("open {}: {e}", path.display())))
}

/// Write records to a VCF/BCF. `header` is the output header (build it from a
/// template via `bcf::Header::from_template(reader.header())`). Output format is
/// inferred from the path extension (`.bcf` → BCF, else bgzipped VCF). The caller
/// supplies records already translated to `header`.
pub fn write_vcf(
    path: &Path,
    header: &bcf::Header,
    recs: impl IntoIterator<Item = bcf::Record>,
) -> Result<()> {
    let is_bcf = path
        .extension()
        .map(|e| e.eq_ignore_ascii_case("bcf"))
        .unwrap_or(false);
    let format = if is_bcf {
        bcf::Format::Bcf
    } else {
        bcf::Format::Vcf
    };
    // uncompressed=false → bgzip for VCF (.vcf.gz), compressed BCF.
    let mut writer = bcf::Writer::from_path(path, header, false, format)
        .map_err(|e| SdError::Vcf(format!("create {}: {e}", path.display())))?;
    for rec in recs {
        writer
            .write(&rec)
            .map_err(|e| SdError::Vcf(format!("write record: {e}")))?;
    }
    Ok(())
}

/// The exact-dedup / sort key: `(rid, pos, ref_allele, alt_alleles)`. `rid` is
/// translated to the contig NAME so the key is comparable across inputs whose
/// headers may order contigs differently. `bcftools sort` orders by the header's
/// contig order; here we sort by contig order as it appears in the *output*
/// header (built from the first input), then by pos, then alleles — matching
/// `concat | sort` for inputs that share a header (the SDrecall island case).
#[derive(Clone, PartialEq, Eq)]
struct VariantKey {
    rid: i64, // output-header contig index (ordering axis)
    pos: i64,
    alleles: Vec<Vec<u8>>,
}

impl PartialOrd for VariantKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for VariantKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.rid
            .cmp(&other.rid)
            .then(self.pos.cmp(&other.pos))
            .then(self.alleles.cmp(&other.alleles))
    }
}

/// Concatenate + coordinate-sort several VCFs into one output, optionally
/// exact-deduplicating on `(CHROM, POS, REF, ALT)` (`bcftools concat -a -d
/// exact`). The output header is built from the first input. Records are
/// gathered, keyed, sorted by `(contig-order, pos, alleles)`, deduped, and
/// written.
///
/// The inputs are individually sorted, so this is effectively a sort over the
/// concatenation; for the typical small per-island VCF set the simple
/// gather-sort is O(N log N) and clear. `threads` forwards to the BGZF pools.
pub fn concat_sort_vcfs(
    inputs: &[&Path],
    out: &Path,
    dedup_exact: bool,
    threads: u8,
) -> Result<()> {
    if inputs.is_empty() {
        return Err(SdError::Vcf("concat_sort_vcfs: no inputs".to_string()));
    }
    remove_vcf_indexes(out)?;

    // Output header from the first input.
    let first = read_vcf(inputs[0])?;
    let out_header = bcf::Header::from_template(first.header());
    // A writer we will translate every record into (gives us its HeaderView for rid lookup).
    let mut writer = {
        let is_bcf = out
            .extension()
            .map(|e| e.eq_ignore_ascii_case("bcf"))
            .unwrap_or(false);
        let format = if is_bcf {
            bcf::Format::Bcf
        } else {
            bcf::Format::Vcf
        };
        bcf::Writer::from_path(out, &out_header, false, format)
            .map_err(|e| SdError::Vcf(format!("create {}: {e}", out.display())))?
    };
    let _ = writer.set_threads(threads as usize);

    // Gather all records, translated into the output header, with a sort key.
    let mut keyed: Vec<(VariantKey, bcf::Record)> = Vec::new();
    for input in inputs {
        let mut reader = read_vcf(input)?;
        let _ = reader.set_threads(threads as usize);
        for res in reader.records() {
            let mut rec = res.map_err(|e| SdError::Vcf(format!("read record: {e}")))?;
            // Translate the record from its source header into the writer header.
            writer.translate(&mut rec);
            let pos = rec.pos();
            let alleles: Vec<Vec<u8>> = rec.alleles().iter().map(|a| a.to_vec()).collect();
            // After translate, rec.rid() indexes the OUTPUT header's contigs.
            let rid = rec.rid().map(|r| r as i64).unwrap_or(-1);
            keyed.push((VariantKey { rid, pos, alleles }, rec));
        }
    }

    // Stable sort by key (coordinate order via output-header contig index).
    keyed.sort_by(|a, b| a.0.cmp(&b.0));

    let mut written = 0usize;
    let mut prev_key: Option<VariantKey> = None;
    for (key, rec) in &keyed {
        if dedup_exact {
            if let Some(pk) = &prev_key {
                if pk == key {
                    continue; // exact duplicate on (CHROM,POS,REF,ALT)
                }
            }
        }
        writer
            .write(rec)
            .map_err(|e| SdError::Vcf(format!("write record: {e}")))?;
        written += 1;
        prev_key = Some(key.clone());
    }

    log::info!(
        "concat_sort_vcfs: {} inputs → {} records written ({} total, dedup_exact={})",
        inputs.len(),
        written,
        keyed.len(),
        dedup_exact
    );
    // Ensure all BGZF blocks and the terminator are flushed before an external
    // bcftools process opens the file for indexing.
    drop(writer);
    if is_indexable_vcf(out) {
        index_vcf(out)?;
    }
    Ok(())
}

fn is_indexable_vcf(path: &Path) -> bool {
    let name = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    name.ends_with(".vcf.gz") || name.ends_with(".bcf")
}

fn index_vcf(path: &Path) -> Result<()> {
    let status = std::process::Command::new("bcftools")
        .arg("index")
        .arg("-f")
        .arg(path)
        .status()
        .map_err(|e| SdError::Vcf(format!("spawn bcftools index {}: {e}", path.display())))?;
    if !status.success() {
        return Err(SdError::Vcf(format!(
            "bcftools index failed (exit {:?}) for {}",
            status.code(),
            path.display()
        )));
    }
    Ok(())
}

fn remove_vcf_indexes(vcf: &Path) -> Result<()> {
    for path in vcf_index_paths(vcf) {
        match std::fs::remove_file(&path) {
            Ok(()) => {}
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => {
                return Err(SdError::Vcf(format!(
                    "remove stale VCF index {}: {e}",
                    path.display()
                )));
            }
        }
    }
    Ok(())
}

fn vcf_index_paths(vcf: &Path) -> [PathBuf; 2] {
    [
        append_path_suffix(vcf, ".csi"),
        append_path_suffix(vcf, ".tbi"),
    ]
}

fn append_path_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut s = path.as_os_str().to_os_string();
    s.push(suffix);
    PathBuf::from(s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bcf::record::GenotypeAllele;

    /// Build a tiny VCF with a single contig and the given variants
    /// `(pos0based, ref, alt)`. Returns the temp file (kept alive by the caller).
    fn write_test_vcf(variants: &[(i64, &str, &str)]) -> tempfile::NamedTempFile {
        let tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        let mut header = bcf::Header::new();
        header.push_record(b"##contig=<ID=chr1,length=100000>");
        header.push_record(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        header.push_sample(b"sample1");
        {
            let mut w =
                bcf::Writer::from_path(tmp.path(), &header, true, bcf::Format::Vcf).unwrap();
            for (pos, r, a) in variants {
                let mut rec = w.empty_record();
                let rid = w.header().name2rid(b"chr1").unwrap();
                rec.set_rid(Some(rid));
                rec.set_pos(*pos);
                let alleles = [r.as_bytes(), a.as_bytes()];
                rec.set_alleles(&alleles).unwrap();
                rec.push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
                    .unwrap();
                w.write(&rec).unwrap();
            }
        }
        tmp
    }

    fn read_variants(path: &Path) -> Vec<(i64, String, String)> {
        let mut reader = read_vcf(path).unwrap();
        let mut out = Vec::new();
        for res in reader.records() {
            let rec = res.unwrap();
            let alleles = rec.alleles();
            out.push((
                rec.pos(),
                String::from_utf8_lossy(alleles[0]).to_string(),
                String::from_utf8_lossy(alleles[1]).to_string(),
            ));
        }
        out
    }

    #[test]
    fn only_bgzipped_vcf_and_bcf_outputs_are_indexed() {
        assert!(is_indexable_vcf(Path::new("a.vcf.gz")));
        assert!(is_indexable_vcf(Path::new("a.bcf")));
        assert!(!is_indexable_vcf(Path::new("a.vcf")));
    }

    #[test]
    fn round_trip_single_vcf() {
        let v = write_test_vcf(&[(100, "A", "T"), (200, "C", "G")]);
        let got = read_variants(v.path());
        assert_eq!(
            got,
            vec![
                (100, "A".to_string(), "T".to_string()),
                (200, "C".to_string(), "G".to_string())
            ]
        );
    }

    #[test]
    fn concat_sorts_across_inputs() {
        // Input A has pos 300 and 100; input B has 200 — concat+sort → 100,200,300.
        let a = write_test_vcf(&[(300, "A", "T"), (100, "A", "T")]);
        let b = write_test_vcf(&[(200, "C", "G")]);
        let out = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        concat_sort_vcfs(&[a.path(), b.path()], out.path(), false, 1).unwrap();
        let got = read_variants(out.path());
        let positions: Vec<i64> = got.iter().map(|(p, _, _)| *p).collect();
        assert_eq!(positions, vec![100, 200, 300]);
    }

    #[test]
    fn concat_dedup_exact_removes_identical_records() {
        // Same variant (100,A,T) appears in both inputs → one copy with dedup.
        let a = write_test_vcf(&[(100, "A", "T"), (150, "G", "C")]);
        let b = write_test_vcf(&[(100, "A", "T")]);
        let out = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        concat_sort_vcfs(&[a.path(), b.path()], out.path(), true, 1).unwrap();
        let got = read_variants(out.path());
        assert_eq!(
            got,
            vec![
                (100, "A".to_string(), "T".to_string()),
                (150, "G".to_string(), "C".to_string())
            ]
        );
    }

    #[test]
    fn concat_without_dedup_keeps_duplicates() {
        let a = write_test_vcf(&[(100, "A", "T")]);
        let b = write_test_vcf(&[(100, "A", "T")]);
        let out = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        concat_sort_vcfs(&[a.path(), b.path()], out.path(), false, 1).unwrap();
        let got = read_variants(out.path());
        assert_eq!(got.len(), 2);
    }

    #[test]
    fn concat_bgzipped_vcf_is_closed_before_indexing() {
        if std::process::Command::new("bcftools")
            .arg("--version")
            .output()
            .is_err()
        {
            return;
        }

        let a = write_test_vcf(&[(100, "A", "T")]);
        let b = write_test_vcf(&[(200, "C", "G")]);
        let out = tempfile::Builder::new()
            .suffix(".vcf.gz")
            .tempfile()
            .unwrap();

        concat_sort_vcfs(&[a.path(), b.path()], out.path(), false, 1).unwrap();

        let got = read_variants(out.path());
        let positions: Vec<i64> = got.iter().map(|(p, _, _)| *p).collect();
        assert_eq!(positions, vec![100, 200]);
        assert!(vcf_index_paths(out.path()).iter().any(|p| p.is_file()));
    }

    #[test]
    fn concat_no_inputs_errors() {
        let out = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        let err = concat_sort_vcfs(&[], out.path(), false, 1).unwrap_err();
        assert!(matches!(err, SdError::Vcf(_)), "got {err:?}");
    }
}
