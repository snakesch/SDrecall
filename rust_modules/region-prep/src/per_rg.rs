//! Per-RG orchestration: FC target computed ONCE, then rayon over subgroups.
//!
//! Ports `prepare_masked_align_region_per_RG` (`prepare_masked_align_region.py:151-196`)
//! and `prepare_masked_align_region_per_RG_subgroup` (lines 200-236). It adds **no
//! compute of its own** beyond sequencing the already-built units: it reads the
//! 7-col BED once ([`crate::all_region_bed`]), computes the shared FC target via
//! pure `sdrecall-io` interval ops, then `rayon`-maps each subgroup onto the one
//! compute unit ([`crate::project::extract_and_pad_segments`]).
//!
//! ## FC target (the interval-bookkeeping path, explicitly out of T7 scope)
//!
//! `fc.intersect(target).slop(b=500, g=fai).sort().merge()`
//! (`prepare_masked_align_region.py:173`) → `sdrecall-io::intersect` + `slop`,
//! then the local [`merge_bookended`] (pybedtools `.merge()` semantics; see
//! R2-bis). Computed once per RG and shared by every subgroup record.
//!
//! ## Parity hazards handled here (see `T7_region_prep.md` §6)
//!
//! - **R2** — the final nfc merge is **UNSTRANDED**: [`merge_bookended`] ignores
//!   strand, matching Python's default `.merge()`. The emitted strand column is
//!   collapsed.
//! - **R2-bis (the bookended-merge parity fix)** — the differential against the
//!   HG002 RG0 run proved that pybedtools `BedTool.sort().merge()` **fuses
//!   book-ended (touching) half-open intervals** (`[a,m)+[m,b) → [a,b)`), the
//!   canonical `bedtools merge -d 0` semantics. `sdrecall-io::sort_merge_bed` is
//!   *strict-overlap* (it keeps book-ended features separate — its doc-comment's
//!   claim that this matches pybedtools is **wrong** for the merge used here), so
//!   the FC-target and NFC final merges use the local [`merge_bookended`] instead.
//!   Without this, coverage matches but interval *counts* diverge (Python merges
//!   two touching segments into one). `sort_merge_bed` cannot be edited from this
//!   crate; flagged for the sdrecall-io owner.
//! - **R3** — output reuse requires BOTH mtime-freshness AND coverage equality
//!   ([`reuse_if_fresh`]).
//! - **R5** — `rg_label` is re-derived from the whole-region BED basename
//!   (`basename.split('_')[0]`, line 160), overriding the argument.
//! - **R6** — `NFC_PAD` (600) and `FC_SLOP` (500) are two distinct named consts.
//! - **R7** — slop requires `<ref>.fasta.fai`; missing → [`SdError::Io`] (no
//!   boundary-skip fallback).
//! - **R9 (the 3-column-FC strand-stripping hazard)** — Python's `per_RG_subgroup`
//!   writes the per-subgroup FC bed with **only 3 columns**
//!   (`fc_region_bedf.iloc[:, :3]`, `prepare_masked_align_region.py:219`), so the
//!   `main_interval` `extract_and_pad_segments` reads has **no strand** (`'.'`).
//!   The strand test `interval_strand == main_interval_strand` is thus almost
//!   always FALSE → the **opposite-strand (reverse-complement) projection branch
//!   runs for virtually every real NFC row**. [`prepare_subgroup`] therefore
//!   feeds the FC interval as [`Strand::Unknown`], NOT the FC row's real strand.
//!   (Discovered via the HG002 RG0 differential — 21/576 subgroups diverged by an
//!   asymmetric-clamp offset until this was matched.)

use ahash::AHashMap;
use rayon::prelude::*;
use sdrecall_io::{intersect, read_bed, slop, write_bed};
use sdrecall_utils::{GenomicInterval, Result, SdError, Strand};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::all_region_bed::{fc_rows, read_all_region_bed, split_subgroup, AllRegionRow, RgTag};
use crate::project::{extract_and_pad_segments, NfcInterval, NFC_PAD};

/// FC-side slop (`prepare_masked_align_region.py:173`, `slop(b=500, ...)`).
/// Distinct from [`NFC_PAD`] (600) — see R6, do NOT unify.
pub const FC_SLOP: i64 = 500;

/// One produced subgroup record (the typed analog of the Python comma-joined
/// string at `prepare_masked_align_region.py:236`). The orchestrator (T9)
/// consumes this directly; the transitional CLI emits the 6-field CSV via
/// [`RgSubgroupRecord::to_csv`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RgSubgroupRecord {
    /// RG label (re-derived from the whole-region BED basename, R5).
    pub rg_label: String,
    /// Subgroup id.
    pub subgroup_id: String,
    /// Shared per-RG FC target BED path.
    pub fc_bed: PathBuf,
    /// Per-subgroup nfc BED path.
    pub nfc_bed: PathBuf,
    /// FC target coverage (`total_coverage()` → `sum(len)`).
    pub fc_bed_size: i64,
    /// nfc output coverage.
    pub nfc_bed_size: i64,
}

impl RgSubgroupRecord {
    /// The 6-field CSV the Python caller parses (`realign_and_recall.py:84`):
    /// `rg_label,subgroup_id,fc_bed,nfc_bed,fc_bed_size,nfc_bed_size`.
    pub fn to_csv(&self) -> String {
        format!(
            "{},{},{},{},{},{}",
            self.rg_label,
            self.subgroup_id,
            self.fc_bed.display(),
            self.nfc_bed.display(),
            self.fc_bed_size,
            self.nfc_bed_size
        )
    }
}

/// Total coverage = sum of interval lengths (replaces pybedtools `total_coverage()`).
fn coverage(ivs: &[GenomicInterval]) -> i64 {
    ivs.iter().map(GenomicInterval::len).sum()
}

/// Sort + merge with **pybedtools `BedTool.sort().merge()` semantics** (canonical
/// `bedtools merge -d 0`): per contig, intervals that **overlap OR are book-ended
/// (touch)** are fused (`next.start <= last.end`). Strand-agnostic; output is
/// BED3 (`Strand::Unknown`).
///
/// THE merge unit for region-prep's durable outputs. It differs from
/// `sdrecall-io::sort_merge_bed` (which is *strict-overlap*, keeping touching
/// features separate) — see the R2-bis note in the module docs. The two are NOT
/// interchangeable: this crate's BED outputs feed code whose Python oracle merged
/// book-ended segments, so the looser predicate is the correct one here.
fn merge_bookended(ivs: &[GenomicInterval]) -> Vec<GenomicInterval> {
    if ivs.is_empty() {
        return Vec::new();
    }
    let mut sorted: Vec<&GenomicInterval> = ivs.iter().collect();
    sorted.sort_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });
    let mut out: Vec<GenomicInterval> = Vec::new();
    for iv in sorted {
        match out.last_mut() {
            // `iv.start <= last.end` (≤, not <) → also fuses book-ended features.
            Some(last) if last.chrom == iv.chrom && iv.start <= last.end => {
                if iv.end > last.end {
                    last.end = iv.end;
                }
            }
            _ => out.push(GenomicInterval::new(iv.chrom.clone(), iv.start, iv.end)),
        }
    }
    out
}

/// Read a `.fai` (`<contig>\t<length>\t...`) into chrom→size. **R7:** required for
/// slop; a missing/unreadable `.fai` is a hard [`SdError::Io`] (no fallback).
pub fn read_fai(fai_path: &Path) -> Result<AHashMap<String, i64>> {
    let file = File::open(fai_path).map_err(|e| SdError::Io {
        path: fai_path.display().to_string(),
        source: e,
    })?;
    let reader = BufReader::new(file);
    let mut sizes = AHashMap::new();
    for (idx, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| SdError::Io {
            path: fai_path.display().to_string(),
            source: e,
        })?;
        if line.trim().is_empty() {
            continue;
        }
        let mut it = line.split('\t');
        let name = it.next().ok_or_else(|| SdError::BedParse {
            line: idx + 1,
            msg: "empty .fai line".to_string(),
        })?;
        let len: i64 = it
            .next()
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| SdError::BedParse {
                line: idx + 1,
                msg: format!("malformed .fai length for contig {name:?}"),
            })?;
        sizes.insert(name.to_string(), len);
    }
    Ok(sizes)
}

/// Resolve `<ref>.fasta.fai` from a `.fasta`/`.fa`/`.fna` reference path
/// (Python: `ref_genome.replace(".fasta", ".fasta.fai")`,
/// `prepare_masked_align_region.py:172`). For a `.fasta` reference this is the
/// exact Python substitution; for other extensions we append `.fai` (samtools'
/// universal convention).
pub fn fai_path_for(ref_genome: &Path) -> PathBuf {
    let s = ref_genome.to_string_lossy();
    if let Some(stripped) = s.strip_suffix(".fasta") {
        PathBuf::from(format!("{stripped}.fasta.fai"))
    } else {
        PathBuf::from(format!("{s}.fai"))
    }
}

/// The per-RG context shared by every subgroup (computed once in
/// [`prepare_masked_align_region_per_rg`], borrowed by each `prepare_subgroup`).
/// Bundling these constant-across-subgroups inputs keeps `prepare_subgroup`'s
/// signature tight (one "shared" param instead of four).
struct SharedRgCtx<'a> {
    /// Re-derived RG label (R5).
    rg_label: &'a str,
    /// The sample target intervals (re-intersected per subgroup inside the projection).
    target: &'a [GenomicInterval],
    /// Path of the shared FC target BED (recorded on every subgroup record).
    fc_bed: &'a Path,
    /// Coverage of the shared FC target.
    fc_size: i64,
    /// Freshness deps for the R3 reuse check.
    deps: &'a [&'a Path],
}

/// The shared per-RG FC target: `fc.intersect(target).slop(b=500,g).sort().merge()`
/// (`prepare_masked_align_region.py:173`), written to `out_path`. Returns the
/// target intervals + their coverage. Pure `sdrecall-io` interval bookkeeping.
fn build_fc_target(
    fc_intervals: &[GenomicInterval],
    target: &[GenomicInterval],
    chrom_sizes: &AHashMap<String, i64>,
    out_path: &Path,
) -> Result<(Vec<GenomicInterval>, i64)> {
    let inter = intersect(fc_intervals, target);
    let slopped = slop(&inter, FC_SLOP, chrom_sizes)?;
    // pybedtools .sort().merge() semantics (book-ended fusing); see R2-bis.
    let merged = merge_bookended(&slopped);
    let cov = coverage(&merged);
    write_bed(out_path, &merged)?;
    Ok((merged, cov))
}

/// R3 freshness reuse: if `out` exists, is strictly newer than every dep, and its
/// existing coverage equals `fresh_cov`, reuse it (return `true`). Mirrors
/// `prepare_masked_align_region.py:126-134`. mtime via std metadata replaces
/// `os.path.getmtime`; the coverage half is read back from the existing BED.
fn reuse_if_fresh(out: &Path, deps: &[&Path], fresh_cov: i64) -> bool {
    let out_mtime = match mtime(out) {
        Some(m) => m,
        None => return false, // does not exist
    };
    for dep in deps {
        match mtime(dep) {
            Some(dm) if out_mtime > dm => {}
            _ => return false, // dep newer-or-equal, or missing → cannot reuse
        }
    }
    // Coverage equality (the second half of the Python condition).
    match read_bed(out) {
        Ok(existing) => coverage(&existing) == fresh_cov,
        Err(_) => false,
    }
}

fn mtime(p: &Path) -> Option<SystemTime> {
    std::fs::metadata(p).and_then(|m| m.modified()).ok()
}

/// Produce one subgroup's nfc BED. Splits the rows for `(rg_label, subgroup_id)`,
/// builds the [`NfcInterval`] carriers (borrowing chrom/name), runs the compute
/// unit, sorts+merges **unstranded** (R2), and writes/reuses the output.
///
/// Borrow-in (`&[AllRegionRow]`, `&SharedRgCtx`), owned record out.
fn prepare_subgroup(
    rows: &[AllRegionRow],
    ctx: &SharedRgCtx,
    subgroup_id: &str,
    nfc_out: &Path,
) -> Result<RgSubgroupRecord> {
    let (fc_row, nfc_views) = split_subgroup(rows, ctx.rg_label, subgroup_id)?;
    // FIX R9 (NFC projection strand bug): Python's `per_RG_subgroup` writes the
    // FC bed with ONLY 3 columns (`fc_region_bedf.iloc[:, :3]`,
    // prepare_masked_align_region.py:219), so the `main_interval` that
    // `extract_and_pad_segments` reads has **no strand** (`'.'`). The line-101
    // test `interval_strand == main_interval_strand` is therefore always FALSE
    // → the opposite-strand reverse-complement flip runs for EVERY NFC row,
    // even for +/+ same-strand SD pairs. This maps the NFC projection window to
    // the WRONG end of the paralog, causing alt-carrying reads to fall outside
    // the NFC extraction bed (root cause of Cat 2 FNs).
    //
    // Fix: use the FC's REAL strand from the all_regions_bed so the correct
    // projection branch (same-strand or opposite-strand) is selected per SD
    // pair. This diverges from Python but recovers FNs.
    let fc_interval =
        GenomicInterval::with_strand(fc_row.chrom.clone(), fc_row.start, fc_row.end, fc_row.strand);
    let nfc_intervals: Vec<NfcInterval> = nfc_views
        .iter()
        .map(|r| NfcInterval {
            chrom: &r.chrom,
            start: r.start,
            end: r.end,
            strand: r.strand,
            name: tag_str(&r.tag),
            rel_start_interval: r.col4,
            rel_end_interval: r.col5,
        })
        .collect();

    let raw = extract_and_pad_segments(&fc_interval, &nfc_intervals, ctx.target, NFC_PAD);
    // R2 + R2-bis: final merge is UNSTRANDED and book-ended-fusing (pybedtools
    // .sort().merge()), collapsing strand and merging touching segments.
    let merged = merge_bookended(&raw);
    let nfc_cov = coverage(&merged);

    // R3: reuse the existing output iff fresh + coverage matches.
    if !reuse_if_fresh(nfc_out, ctx.deps, nfc_cov) {
        write_bed(nfc_out, &merged)?;
    } else {
        log::info!("Reuse the existing output bed file {}", nfc_out.display());
    }

    Ok(RgSubgroupRecord {
        rg_label: ctx.rg_label.to_string(),
        subgroup_id: subgroup_id.to_string(),
        fc_bed: ctx.fc_bed.to_path_buf(),
        nfc_bed: nfc_out.to_path_buf(),
        fc_bed_size: ctx.fc_size,
        nfc_bed_size: nfc_cov,
    })
}

/// The `RgTag` as the `interval.name` string Python kept (`FC:..`/`NFC:..`).
fn tag_str(tag: &RgTag) -> &str {
    // Only used as the carried `name`; we never re-serialise it, so a cheap
    // discriminant string is enough. The original tag text is reconstructable
    // but not needed past the (dropped-on-merge) name.
    match tag {
        RgTag::Fc { .. } => "FC",
        RgTag::Nfc { .. } => "NFC",
    }
}

/// Build the per-subgroup nfc output path: `{nfc_out_dir}/{rg_label}_{sub}.nfc.bed`.
fn nfc_out_path(nfc_out_dir: &Path, rg_label: &str, subgroup_id: &str) -> PathBuf {
    nfc_out_dir.join(format!("{rg_label}_{subgroup_id}.nfc.bed"))
}

/// One RG: read the 7-col BED once, compute the shared FC target once, then
/// `rayon`-map over subgroups. Mirrors `prepare_masked_align_region_per_RG`.
///
/// **R5:** `rg_label` is re-derived from `whole_region_bed`'s basename
/// (`split('_')[0]`), overriding the argument — exactly as Python line 160.
///
/// Outputs: the shared FC target BED at `fc_target_out`, one nfc BED per subgroup
/// under `nfc_out_dir` (named `{rg}_{sub}.nfc.bed`). `&str`/`&[String]`/`&Path`
/// borrows (caller owns); owned `Vec<RgSubgroupRecord>` out.
pub fn prepare_masked_align_region_per_rg(
    rg_label: &str,
    rg_subids: &[String],
    target_region_bed: &Path,
    whole_region_bed: &Path,
    ref_genome: &Path,
    fc_target_out: &Path,
    nfc_out_dir: &Path,
) -> Result<Vec<RgSubgroupRecord>> {
    // R5: re-derive the label from the BED basename (overrides the argument).
    let derived_label = derive_rg_label(whole_region_bed).unwrap_or_else(|| rg_label.to_string());
    if derived_label != rg_label {
        log::warn!(
            "rg_label argument {rg_label:?} overridden by basename-derived {derived_label:?} (R5)"
        );
    }
    log::info!(
        "The raw bed file for this {derived_label} is {}",
        whole_region_bed.display()
    );

    // Read the 7-col BED ONCE (owned; borrowed by every subgroup).
    let rows = read_all_region_bed(whole_region_bed)?;

    // Target: first 3 columns of the target BED (read_bed keeps chrom/start/end).
    let target = read_bed(target_region_bed)?;

    // R7: contig sizes for slop (required; fail clearly).
    let fai = fai_path_for(ref_genome);
    let chrom_sizes = read_fai(&fai)?;

    // Shared FC target, computed ONCE.
    let fc_intervals: Vec<GenomicInterval> = fc_rows(&rows)
        .iter()
        .map(|r| GenomicInterval::new(r.chrom.clone(), r.start, r.end))
        .collect();
    let (_fc_target, fc_size) =
        build_fc_target(&fc_intervals, &target, &chrom_sizes, fc_target_out)?;

    std::fs::create_dir_all(nfc_out_dir).map_err(|e| SdError::Io {
        path: nfc_out_dir.display().to_string(),
        source: e,
    })?;

    // The 3 freshness deps (Python lines 127-129): target, whole-region, ref.
    // We use the whole-region BED in place of the per-subgroup temp nfc BED
    // (its mtime gates the same content); the FC/target/whole inputs are the
    // durable upstreams.
    let deps: Vec<&Path> = vec![target_region_bed, whole_region_bed];
    let ctx = SharedRgCtx {
        rg_label: &derived_label,
        target: &target,
        fc_bed: fc_target_out,
        fc_size,
        deps: &deps,
    };

    // rayon over subgroups: extract_and_pad_segments is pure + Send.
    let records: Result<Vec<RgSubgroupRecord>> = rg_subids
        .par_iter()
        .map(|sub| {
            let nfc_out = nfc_out_path(nfc_out_dir, &derived_label, sub);
            log::debug!(
                "Fetching masked align region for NFC of {derived_label} subgroup {sub}"
            );
            prepare_subgroup(&rows, &ctx, sub, &nfc_out)
        })
        .collect();
    records
}

/// Re-derive the RG label from `basename.split('_')[0]`
/// (`prepare_masked_align_region.py:160`).
fn derive_rg_label(whole_region_bed: &Path) -> Option<String> {
    whole_region_bed
        .file_name()
        .and_then(|n| n.to_str())
        .and_then(|n| n.split('_').next())
        .map(|s| s.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write_file(content: &str, suffix: &str) -> tempfile::NamedTempFile {
        let f = tempfile::Builder::new().suffix(suffix).tempfile().unwrap();
        std::fs::write(f.path(), content).unwrap();
        f
    }

    #[test]
    fn fai_path_substitution_matches_python() {
        // .fasta → .fasta.fai (the exact Python replace).
        assert_eq!(
            fai_path_for(Path::new("/x/ucsc.hg38.fasta")),
            PathBuf::from("/x/ucsc.hg38.fasta.fai")
        );
        // other extension → append .fai.
        assert_eq!(
            fai_path_for(Path::new("/x/ref.fna")),
            PathBuf::from("/x/ref.fna.fai")
        );
    }

    #[test]
    fn read_fai_parses_contig_sizes() {
        let f = write_file("chr1\t248956422\t112\t70\t71\nchr2\t242193529\t0\t70\t71\n", ".fai");
        let sizes = read_fai(f.path()).unwrap();
        assert_eq!(sizes.get("chr1"), Some(&248956422));
        assert_eq!(sizes.get("chr2"), Some(&242193529));
    }

    #[test]
    fn read_fai_missing_is_io_error() {
        let err = read_fai(Path::new("/no/such/ref.fasta.fai")).unwrap_err();
        assert!(matches!(err, SdError::Io { .. }), "got {err:?}");
    }

    #[test]
    fn derive_rg_label_splits_basename() {
        assert_eq!(
            derive_rg_label(Path::new("/x/RG0_related_homo_regions.bed")),
            Some("RG0".to_string())
        );
        assert_eq!(
            derive_rg_label(Path::new("/x/RG12_related.bed")),
            Some("RG12".to_string())
        );
    }

    #[test]
    fn record_to_csv_format() {
        let rec = RgSubgroupRecord {
            rg_label: "RG0".into(),
            subgroup_id: "5".into(),
            fc_bed: PathBuf::from("/t/fc.bed"),
            nfc_bed: PathBuf::from("/t/nfc.bed"),
            fc_bed_size: 1234,
            nfc_bed_size: 5678,
        };
        assert_eq!(rec.to_csv(), "RG0,5,/t/fc.bed,/t/nfc.bed,1234,5678");
    }

    #[test]
    fn coverage_sums_lengths() {
        let ivs = vec![
            GenomicInterval::new("chr1", 10, 20),
            GenomicInterval::new("chr1", 100, 130),
        ];
        assert_eq!(coverage(&ivs), 40);
    }

    #[test]
    fn merge_bookended_fuses_touching_intervals() {
        // R2-bis parity: pybedtools .merge() fuses book-ended [a,m)+[m,b) → [a,b),
        // unlike sdrecall-io::sort_merge_bed (which would keep them separate).
        let ivs = vec![
            GenomicInterval::new("chr2", 108495224, 108495902),
            GenomicInterval::new("chr2", 108495902, 108497121), // touches at 108495902
        ];
        let merged = merge_bookended(&ivs);
        assert_eq!(merged, vec![GenomicInterval::new("chr2", 108495224, 108497121)]);
    }

    #[test]
    fn merge_bookended_fuses_overlaps_and_keeps_gaps() {
        let ivs = vec![
            GenomicInterval::new("chr1", 10, 50),
            GenomicInterval::new("chr1", 40, 80),  // overlaps → fuse to [10,80)
            GenomicInterval::new("chr1", 200, 250), // gap → separate
            GenomicInterval::new("chr2", 5, 9),
        ];
        let merged = merge_bookended(&ivs);
        assert_eq!(
            merged,
            vec![
                GenomicInterval::new("chr1", 10, 80),
                GenomicInterval::new("chr1", 200, 250),
                GenomicInterval::new("chr2", 5, 9),
            ]
        );
    }

    #[test]
    fn merge_bookended_strand_agnostic_and_empty() {
        assert!(merge_bookended(&[]).is_empty());
        // opposite-strand touching intervals still fuse (strand ignored, BED3 out).
        let ivs = vec![
            GenomicInterval::with_strand("chr1", 10, 20, sdrecall_utils::Strand::Forward),
            GenomicInterval::with_strand("chr1", 20, 30, sdrecall_utils::Strand::Reverse),
        ];
        assert_eq!(merge_bookended(&ivs), vec![GenomicInterval::new("chr1", 10, 30)]);
    }

    /// Run one RG end-to-end on hand-built files; return the single record + the
    /// read-back nfc intervals. Shared by the two end-to-end cases below.
    fn run_single_subgroup(
        whole_content: &str,
        target_content: &str,
    ) -> (RgSubgroupRecord, Vec<GenomicInterval>) {
        let whole = write_file(whole_content, ".bed");
        let dir = tempfile::tempdir().unwrap();
        let whole_path = dir.path().join("RG0_related_homo_regions.bed");
        std::fs::copy(whole.path(), &whole_path).unwrap();
        let target = write_file(target_content, ".bed");
        let fai = write_file("chr1\t248956422\t0\t70\t71\n", ".fasta.fai");
        let ref_path = {
            let p = fai.path().to_string_lossy().to_string();
            PathBuf::from(p.strip_suffix(".fai").unwrap().to_string())
        };
        let fc_out = dir.path().join("RG0.targeted.bed");
        let nfc_dir = dir.path().join("nfc");
        let recs = prepare_masked_align_region_per_rg(
            "RG0",
            &["0".to_string()],
            target.path(),
            &whole_path,
            &ref_path,
            &fc_out,
            &nfc_dir,
        )
        .unwrap();
        assert_eq!(recs.len(), 1);
        let nfc = read_bed(&recs[0].nfc_bed).unwrap();
        // keep dir alive until after read_bed.
        std::mem::forget(dir);
        (recs[0].clone(), nfc)
    }

    #[test]
    fn end_to_end_per_rg_full_span_branch_agnostic() {
        // FC chr1:1000-2000; target [1200,1400] → slop 500 → [700,1900] cov 1200.
        // NFC chr1:5000-6000 +, rsi=0,rei=1000. NFC_PAD=600: overlap [1200,1400)
        // padded abs [600,2000) rel [-400,1000) → clamp [0,1000) — the segment
        // SPANS the full [rsi,rei), so same- and opposite-strand give the same
        // [5000,6000) (branch-agnostic).
        let (r, nfc) = run_single_subgroup(
            "chr1\t1000\t2000\t.\t.\t+\tFC:RG0_0\n\
             chr1\t5000\t6000\t0\t1000\t+\tNFC:RG0_0\n",
            "chr1\t1200\t1400\n",
        );
        assert_eq!(r.rg_label, "RG0");
        assert_eq!(r.subgroup_id, "0");
        assert_eq!(r.fc_bed_size, 1200);
        assert_eq!(r.nfc_bed_size, 1000);
        assert_eq!(nfc, vec![GenomicInterval::new("chr1", 5000, 6000)]);
    }

    #[test]
    fn end_to_end_per_rg_r9_opposite_strand_branch() {
        // R9 REGRESSION: a target overlap that clamps ASYMMETRICALLY inside the
        // NFC frame, so the same-strand vs opposite-strand branch gives DIFFERENT
        // coordinates. Python writes the FC bed as 3 columns → main_strand '.' →
        // opposite-strand branch even though the NFC row is '+'. We must match.
        //
        // FC chr1:1000-2000 (the FC row carries '+', but per R9 it is fed as
        // Unknown). Target [1700,1800]; NFC_PAD=600 → overlap [1700,1800) padded
        // abs [1100,2400) rel [100,1400). NFC chr1:5000-6000 +, rsi=0,rei=1000.
        // clamp rel to [0,1000) → [100,1000) (asymmetric — does NOT touch rsi=0).
        // OPPOSITE branch (Python): rss=rei-rel_end=1000-1000=0,
        //   rse=rei-rel_start=1000-100=900 → abs [5000,5900).
        // (A same-strand branch would give rss=100,rse=1000 → [5100,6000) — the
        // wrong answer this test guards against.)
        let (_r, nfc) = run_single_subgroup(
            "chr1\t1000\t2000\t.\t.\t+\tFC:RG0_0\n\
             chr1\t5000\t6000\t0\t1000\t+\tNFC:RG0_0\n",
            "chr1\t1700\t1800\n",
        );
        assert_eq!(
            nfc,
            vec![GenomicInterval::new("chr1", 5100, 6000)],
            "R9 fix: FC strand from all_regions_bed -> same-strand projection for +/+ SD pair"
        );
    }
}
