//! The ONE versatile BAM reader + noisy-read predicate — the DUP-1 collapse.
//!
//! Consolidates the three pre-existing Rust copies (`build_phasing_graph/
//! src/bam_reading.rs`, `haplotype_inspection/src/bam_lappers.rs`,
//! `read_extraction/src/lib.rs`) into a single canonical home:
//!
//! - [`is_read_noisy`] — the single noisy-read predicate, a 1:1 port of
//!   `fp_control/bam_ncls.py:115-193` (the canonical Python version, which has
//!   *more* checks than the older Rust copy in `bam_lappers.rs`: qcfail,
//!   `reference_end`/`query_sequence` presence, and mate-same-contig).
//! - [`build_bam_index`] — collate-by-qname → drop noisy groups → build the
//!   per-chrom `Lapper` point-query index + the read / qname maps (port of
//!   `migrate_bam_to_ncls`, `_collate_bam_file`, `_process_qname_group`).
//!
//! HYG-6 fix: the `samtools collate` child's stdout pipe is read via `/dev/fd/N`
//! and the fd handle is dropped (not `mem::forget`-leaked) once htslib has dup'd
//! it; the child is reaped with `wait()`. No per-island fd leak.

use ahash::{AHashMap, AHashSet};
use rust_htslib::bam::{self, HeaderView, Read, Record};
use rust_lapper::{Interval, Lapper};
use rustc_hash::FxHashMap;
use sdrecall_utils::{QnameIdx, Result, SdError};
#[cfg(unix)]
use std::os::unix::io::AsRawFd;
use std::path::Path;
use std::process::{Child, Command, Stdio};
use tempfile::NamedTempFile;

/// Tunable thresholds for [`is_read_noisy`]. `Default` mirrors the Python
/// defaults (`migrate_bam_to_ncls`: mapq 10, basequal median 15, paired, filter).
#[derive(Clone, Copy, Debug)]
pub struct NoisyFilter {
    /// Minimum MAPQ; reads below are noisy.
    pub mapq: u8,
    /// Median base-quality threshold (`<=` is noisy when `filter_noisy`).
    pub basequal_median: u8,
    /// Paired-end mode (enables span / mate-contig checks).
    pub paired: bool,
    /// Enable the quality / soft-clip checks (Python `filter_noisy`).
    pub filter_noisy: bool,
}

impl Default for NoisyFilter {
    fn default() -> Self {
        Self {
            mapq: 10,
            basequal_median: 15,
            paired: true,
            filter_noisy: true,
        }
    }
}

/// Median of base qualities, matching `np.median`: for an even count the mean of
/// the two middle values (NOT integer-truncated), so the `median <= cutoff`
/// boundary matches Python. Returns `f32`.
fn fast_median(q: &[u8]) -> f32 {
    if q.is_empty() {
        return 0.0;
    }
    let mut sorted = q.to_vec();
    sorted.sort_unstable();
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] as f32 + sorted[mid] as f32) / 2.0
    } else {
        sorted[mid] as f32
    }
}

fn ref_name<'a>(rec: &Record, hv: &'a HeaderView) -> Option<&'a [u8]> {
    let tid = rec.tid();
    if tid >= 0 {
        Some(hv.tid2name(tid as u32))
    } else {
        None
    }
}

fn mate_ref_name<'a>(rec: &Record, hv: &'a HeaderView) -> Option<&'a [u8]> {
    let mtid = rec.mtid();
    if mtid >= 0 {
        Some(hv.tid2name(mtid as u32))
    } else {
        None
    }
}

/// THE noisy-read predicate — 1:1 port of `bam_ncls.py:115-193`.
///
/// Returns `true` if the read should be dropped as noisy. Mirrors the Python
/// branch order exactly (so debug parity is line-traceable):
/// - secondary/supplementary → **not** noisy (skipped, `false`) [py 126-128]
/// - unmapped / qcfail / MAPQ<filter / no ref_end / no query_seq → noisy [py 131-149]
/// - paired: span<75 → noisy [py 153-156]; mate on a different contig → noisy
///   [py 159-161]; not-proper-pair is only a warning (NOT a filter) [py 163-164]
/// - single-end: duplicate → noisy [py 167-169]
/// - if `filter_noisy` and qualities present: median baseQ ≤ cut → noisy
///   [py 172-177]; #(Q<cut) ≥ 75 → noisy [py 179-182]; total soft-clip ≥ 75 →
///   noisy [py 184-190]
///
/// Read-only borrow; nothing escapes, so the driver reuses one `Record` buffer.
pub fn is_read_noisy(rec: &Record, hv: &HeaderView, f: &NoisyFilter) -> bool {
    // py 126-128: only primary alignments are evaluated; sec/supp are not noisy.
    if rec.is_secondary() || rec.is_supplementary() {
        return false;
    }

    // py 131-149: common structural checks (both paired and single).
    if rec.is_unmapped() {
        return true;
    }
    if rec.is_quality_check_failed() {
        return true;
    }
    if rec.mapq() < f.mapq {
        return true;
    }
    // py 143-145: reference_end is None — in rust-htslib, an aligned primary read
    // always has a defined cigar().end_pos(); the analogue is "no mapped query",
    // i.e. an empty CIGAR. py 147-149: query_sequence is None == empty seq.
    if rec.cigar_len() == 0 {
        return true;
    }
    if rec.seq_len() == 0 {
        return true;
    }

    if f.paired {
        // py 153-156: alignment span < 75.
        let aln_len = rec.cigar().end_pos() - rec.pos();
        if aln_len < 75 {
            return true;
        }
        // py 159-161: mate must be on the same contig.
        match (ref_name(rec, hv), mate_ref_name(rec, hv)) {
            (Some(r), Some(m)) if r == m => {}
            // Python compares reference_name vs next_reference_name; if either is
            // None or they differ, the strings are unequal → noisy.
            _ => return true,
        }
        // py 163-164: not-a-proper-pair is logged but NOT filtered. No-op here.
    } else {
        // py 167-169: single-end — drop marked duplicates.
        if rec.is_duplicate() {
            return true;
        }
    }

    // py 171-190: quality / soft-clip checks (gated on filter_noisy).
    if f.filter_noisy {
        let q = rec.qual();
        if !q.is_empty() {
            // py 175-177
            if fast_median(q) <= f.basequal_median as f32 {
                return true;
            }
            // py 179-182
            let num_low = q.iter().filter(|&&b| b < f.basequal_median).count();
            if num_low >= 75 {
                return true;
            }
            // py 184-190: total soft-clip length (CIGAR op 4) ≥ 75.
            let softclip: u32 = rec
                .cigar()
                .iter()
                .filter_map(|c| match c {
                    bam::record::Cigar::SoftClip(l) => Some(*l),
                    _ => None,
                })
                .sum();
            if softclip >= 75 {
                return true;
            }
        }
    }

    false
}

/// The per-island in-memory BAM index — the owned maps the haplotype-inspection
/// path stores for an island's lifetime. Integer-keyed maps use `FxHashMap`,
/// string-keyed maps use `ahash` (project convention).
#[derive(Debug, Default)]
pub struct BamIndex {
    /// Per-chrom point-query interval tree: `chrom → Lapper<u32, QnameIdx>`.
    pub lapper: FxHashMap<String, Lapper<u32, QnameIdx>>,
    /// Retained reads, keyed by dense qname index.
    pub reads: FxHashMap<QnameIdx, Vec<Record>>,
    /// `qname_idx → qname`.
    pub qname: FxHashMap<QnameIdx, String>,
    /// `qname → qname_idx`.
    pub qname_idx: AHashMap<String, QnameIdx>,
    /// Qnames dropped as noisy.
    pub noisy: AHashSet<String>,
}

/// Whether `samtools` is on PATH (only needed for paired-end collation).
fn samtools_available() -> bool {
    Command::new("samtools")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

/// Owns the `samtools collate` child so it is reaped on drop (HYG-6: no leaked
/// pipe fd, no zombie). The stdout pipe fd was already dup'd by htslib's
/// `hts_open`, so we close our copy immediately by dropping the `ChildStdout`
/// inside `spawn_collate_pipe` — only the `Child` lifetime is tracked here.
struct CollateChild {
    child: Child,
}

impl Drop for CollateChild {
    fn drop(&mut self) {
        // Best-effort reap; if the reader consumed all records the child will
        // have already exited. We do not panic in Drop.
        match self.child.wait() {
            Ok(status) if !status.success() => {
                log::warn!("samtools collate exited with {status}");
            }
            Err(e) => log::warn!("waiting on samtools collate child failed: {e}"),
            _ => {}
        }
    }
}

/// Spawn `samtools collate -f` with stdout piped, open a BAM reader on the pipe
/// via `/dev/fd/N`. Returns `None` if samtools is unavailable. The returned
/// `CollateChild` must be kept alive until all records are read, then dropped
/// (which reaps the child). HYG-6: the `ChildStdout` is dropped here once htslib
/// has dup'd the fd — NOT `mem::forget`-ed.
#[cfg(unix)]
fn spawn_collate_pipe(
    bam_path: &Path,
    threads: u8,
) -> Result<Option<(bam::Reader, CollateChild)>> {
    if !samtools_available() {
        return Ok(None);
    }
    let mut child = Command::new("samtools")
        .args([
            "collate",
            "-f",
            "-@",
            &threads.to_string(),
            &bam_path.display().to_string(),
            "-o",
            "-",
        ])
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .map_err(|e| SdError::Io {
            path: "samtools collate".to_string(),
            source: e,
        })?;

    let stdout = child
        .stdout
        .take()
        .ok_or_else(|| SdError::Htslib("failed to capture samtools stdout".to_string()))?;
    let fd = stdout.as_raw_fd();
    let fd_path = format!("/dev/fd/{fd}");

    match bam::Reader::from_path(&fd_path) {
        Ok(reader) => {
            // htslib's hts_open dup'd the fd; drop our copy now (HYG-6: explicit
            // drop, no mem::forget). The pipe stays open via the child's own end.
            drop(stdout);
            Ok(Some((reader, CollateChild { child })))
        }
        Err(e) => {
            let _ = child.kill();
            let _ = child.wait();
            Err(SdError::Htslib(format!("open collate pipe {fd_path}: {e}")))
        }
    }
}

/// Fallback: collate to a temp BAM (kept alive via the returned handle).
fn collate_to_tempfile(bam_path: &Path, threads: u8) -> Result<Option<NamedTempFile>> {
    if !samtools_available() {
        return Ok(None);
    }
    let tf = NamedTempFile::with_suffix(".bam").map_err(|e| SdError::Io {
        path: "collate temp".to_string(),
        source: e,
    })?;
    let out = Command::new("samtools")
        .args([
            "collate",
            "-f",
            "-@",
            &threads.to_string(),
            &bam_path.display().to_string(),
            "-o",
            &tf.path().display().to_string(),
        ])
        .output()
        .map_err(|e| SdError::Io {
            path: "samtools collate".to_string(),
            source: e,
        })?;
    if !out.status.success() {
        log::warn!(
            "samtools collate to temp failed: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        return Ok(None);
    }
    Ok(Some(tf))
}

/// Build the per-island BAM index — port of `migrate_bam_to_ncls`.
///
/// Paired-end data is collated by qname (pipe → temp-file → in-memory fallback)
/// so same-qname reads are adjacent and can be flushed group-by-group. A qname
/// group is dropped if ANY read is noisy ([`is_read_noisy`]); paired groups also
/// require both R1 and R2. Returns owned maps.
///
/// `&Path` open-only; `threads` forwards to htslib's BGZF pool. Records are
/// cloned into `reads` exactly once per retained read (the index outlives the
/// reader, so owned `Record`s are genuinely needed here).
pub fn build_bam_index(bam_path: &Path, f: &NoisyFilter, threads: u8) -> Result<BamIndex> {
    // ── obtain a (possibly collated) reader ──────────────────────────────────
    enum Source {
        Collated(bam::Reader, Option<CollateChild>, Option<NamedTempFile>),
        Plain(bam::Reader),
    }

    let source = if f.paired {
        #[cfg(unix)]
        let piped = spawn_collate_pipe(bam_path, threads)?;
        #[cfg(not(unix))]
        let piped: Option<(bam::Reader, CollateChild)> = None;

        match piped {
            Some((reader, child)) => Source::Collated(reader, Some(child), None),
            None => match collate_to_tempfile(bam_path, threads)? {
                Some(tf) => {
                    let reader = bam::Reader::from_path(tf.path())
                        .map_err(|e| SdError::Htslib(format!("open collated temp: {e}")))?;
                    Source::Collated(reader, None, Some(tf))
                }
                None => Source::Plain(
                    bam::Reader::from_path(bam_path)
                        .map_err(|e| SdError::Htslib(format!("open {}: {e}", bam_path.display())))?,
                ),
            },
        }
    } else {
        Source::Plain(
            bam::Reader::from_path(bam_path)
                .map_err(|e| SdError::Htslib(format!("open {}: {e}", bam_path.display())))?,
        )
    };

    let collated = matches!(source, Source::Collated(..));
    let (mut reader, _child, _tf) = match source {
        Source::Collated(r, c, t) => (r, c, t),
        Source::Plain(r) => (r, None, None),
    };
    let _ = reader.set_threads(threads as usize);

    let header = reader.header().clone();
    let chroms: Vec<String> = header
        .target_names()
        .iter()
        .map(|n| String::from_utf8_lossy(n).to_string())
        .collect();

    // ── accumulators ─────────────────────────────────────────────────────────
    let mut reads: FxHashMap<QnameIdx, Vec<Record>> = FxHashMap::default();
    let mut intervals: AHashMap<String, Vec<(i64, i64, QnameIdx)>> =
        chroms.iter().map(|c| (c.clone(), Vec::new())).collect();
    let mut qname_idx: AHashMap<String, QnameIdx> = AHashMap::new();
    let mut qname: FxHashMap<QnameIdx, String> = FxHashMap::default();
    let mut noisy: AHashSet<String> = AHashSet::new();
    let mut total_qnames: AHashSet<String> = AHashSet::new();
    let mut counter: u32 = 0;

    // One closure that processes a finished qname group — the single per-group
    // unit (port of _process_qname_group), shared by the collated and the
    // in-memory paths so there is no duplicated grouping logic.
    let flush = |qn: &str,
                 group: &[Record],
                     counter: &mut u32,
                     reads: &mut FxHashMap<QnameIdx, Vec<Record>>,
                     intervals: &mut AHashMap<String, Vec<(i64, i64, QnameIdx)>>,
                     qname_idx: &mut AHashMap<String, QnameIdx>,
                     qname: &mut FxHashMap<QnameIdx, String>,
                     noisy: &mut AHashSet<String>,
                     total_qnames: &mut AHashSet<String>| {
        total_qnames.insert(qn.to_string());
        if group.iter().any(|r| is_read_noisy(r, &header, f)) {
            noisy.insert(qn.to_string());
            return;
        }
        if f.paired {
            let has_r1 = group.iter().any(|r| r.is_first_in_template());
            let has_r2 = group.iter().any(|r| r.is_last_in_template());
            if !has_r1 || !has_r2 {
                return;
            }
        }
        let idx = *qname_idx.entry(qn.to_string()).or_insert_with(|| {
            let i = QnameIdx(*counter);
            qname.insert(i, qn.to_string());
            *counter += 1;
            i
        });
        reads.entry(idx).or_default().extend(group.iter().cloned());
        for r in group {
            if let Some(name) = ref_name(r, &header) {
                let chrom = String::from_utf8_lossy(name).to_string();
                if let Some(v) = intervals.get_mut(&chrom) {
                    v.push((r.pos(), r.cigar().end_pos(), idx));
                }
            }
        }
    };

    // ── streaming pass ───────────────────────────────────────────────────────
    if collated || !f.paired {
        // Collated (or single-end): same-qname reads are adjacent — buffer one
        // group at a time. For single-end each record is effectively its own
        // group, which the adjacency buffer handles transparently.
        let mut cur_qname: Option<String> = None;
        let mut group: Vec<Record> = Vec::new();
        let mut rec = Record::new();
        while let Some(res) = reader.read(&mut rec) {
            res.map_err(|e| SdError::Htslib(format!("read record: {e}")))?;
            if rec.is_secondary() || rec.is_supplementary() || rec.is_duplicate() {
                continue;
            }
            let qn = std::str::from_utf8(rec.qname())
                .map_err(|_| SdError::NonUtf8ReadName)?
                .to_string();
            if cur_qname.as_deref() == Some(&qn) {
                group.push(rec.clone());
            } else {
                if let Some(prev) = cur_qname.take() {
                    if !group.is_empty() {
                        flush(&prev, &group, &mut counter, &mut reads, &mut intervals,
                              &mut qname_idx, &mut qname, &mut noisy, &mut total_qnames);
                    }
                }
                cur_qname = Some(qn);
                group.clear();
                group.push(rec.clone());
            }
        }
        if let Some(prev) = cur_qname {
            if !group.is_empty() {
                flush(&prev, &group, &mut counter, &mut reads, &mut intervals,
                      &mut qname_idx, &mut qname, &mut noisy, &mut total_qnames);
            }
        }
    } else {
        // Paired but no samtools: load all primary reads grouped by qname.
        log::warn!("build_bam_index: paired data without samtools — in-memory grouping fallback");
        let mut by_qname: AHashMap<String, Vec<Record>> = AHashMap::new();
        let mut rec = Record::new();
        while let Some(res) = reader.read(&mut rec) {
            res.map_err(|e| SdError::Htslib(format!("read record: {e}")))?;
            if rec.is_secondary() || rec.is_supplementary() || rec.is_duplicate() {
                continue;
            }
            let qn = std::str::from_utf8(rec.qname())
                .map_err(|_| SdError::NonUtf8ReadName)?
                .to_string();
            by_qname.entry(qn).or_default().push(rec.clone());
        }
        for (qn, group) in &by_qname {
            flush(qn, group, &mut counter, &mut reads, &mut intervals,
                  &mut qname_idx, &mut qname, &mut noisy, &mut total_qnames);
        }
    }

    // Reaping the collate child happens when `_child` drops at end of scope.
    drop(_child);
    drop(_tf);

    // ── build the per-chrom Lapper index ─────────────────────────────────────
    let mut lapper: FxHashMap<String, Lapper<u32, QnameIdx>> = FxHashMap::default();
    for chrom in &chroms {
        let ivs = &intervals[chrom];
        if ivs.is_empty() {
            continue;
        }
        let mut lap_ivs: Vec<Interval<u32, QnameIdx>> = ivs
            .iter()
            .map(|(s, e, idx)| Interval {
                start: *s as u32,
                stop: *e as u32,
                val: *idx,
            })
            .collect();
        lap_ivs.sort_by_key(|iv| iv.start);
        lapper.insert(chrom.clone(), Lapper::new(lap_ivs));
    }

    // Python returns None (→ here a typed error) for paired data with ≤2 qnames.
    if f.paired && total_qnames.len() <= 2 {
        return Err(SdError::InsufficientPairs(total_qnames.len()));
    }

    log::info!(
        "build_bam_index({}): {} qnames retained, {} noisy, {} total",
        bam_path.display(),
        qname_idx.len(),
        noisy.len(),
        total_qnames.len()
    );

    Ok(BamIndex {
        lapper,
        reads,
        qname,
        qname_idx,
        noisy,
    })
}

/// Merge several BAMs into one sorted, indexed output.
///
/// STUB — `TODO(T9)`. No current consumer needs this until the orchestrator (T9),
/// and the SQ-line reconciliation is genuinely non-trivial in rust-htslib, so per
/// the "do not rat-hole" guidance it is left stubbed with the exact semantics
/// documented below rather than half-built.
///
/// ## SQ-reconcile semantics to port at T9 (from `shell_utils.sh::modify_bam_sq_lines`)
///
/// The Python path rebuilds the merged header before merging: it takes input
/// `bam_list[0]`'s header, **strips every `@SQ` and `@PG` line**, then appends
/// fresh `@SQ` lines generated from the reference FASTA index — one
/// `@SQ\tSN:<name>\tLN:<len>` per `.fai` line, **in `.fai` order** (running
/// `samtools faidx` first if the index is stale). `samtools merge -h <header>`
/// then emits records under this reference-ordered SQ dictionary.
///
/// The hard part in rust-htslib: the reference `.fai` SQ order generally differs
/// from (and is a superset of) each input BAM's SQ order, so **every record's
/// `tid` and `mtid` must be remapped from the input's SQ index to the new SQ
/// index, matched by contig NAME** (records whose contig is absent from the
/// reference are a hard error — they would mis-map). rust-htslib has no
/// header-SQ-replace + tid-translate primitive, so this means: build the target
/// `Header` from the `.fai`, construct a `name → new_tid` map, then for each input
/// stream copy each `Record`, rewrite `tid`/`mtid` via the map (and re-validate
/// `pos`), write to one `Writer`, finally coordinate-sort (k-way merge over the
/// already-sorted inputs) and `bam::index::build`. `_ref_fasta` supplies the
/// `.fai`; `_threads` forwards to the BGZF pools.
pub fn merge_bams(_inputs: &[&Path], _out: &Path, _ref_fasta: &Path, _threads: u8) -> Result<()> {
    // TODO(T9): port modify_bam_sq_lines — replace header @SQ/@PG with the
    // reference .fai SQ lines (in .fai order) and remap every record tid/mtid by
    // contig name into the new SQ index, then k-way coordinate-sort + index.
    Err(SdError::Htslib(
        "merge_bams not yet implemented: SQ-line reconcile (modify_bam_sq_lines) deferred to T9; \
         no current consumer (DESIGN §6)"
            .to_string(),
    ))
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
    const SECONDARY: u16 = 0x100;
    const QCFAIL: u16 = 0x200;
    const DUP: u16 = 0x400;
    const SUPPLEMENTARY: u16 = 0x800;

    fn test_header() -> Header {
        let mut h = Header::new();
        for (sn, ln) in [("chr1", 100_000), ("chr2", 100_000)] {
            let mut sq = HeaderRecord::new(b"SQ");
            sq.push_tag(b"SN", sn);
            sq.push_tag(b"LN", ln);
            h.push_record(&sq);
        }
        h
    }

    fn header_view() -> HeaderView {
        HeaderView::from_header(&test_header())
    }

    /// Build a record. `span` controls the Match length (and so the aligned
    /// span); `qual` is the per-base quality (uniform). `tid`/`mtid` set the
    /// contig / mate contig.
    fn make(
        flags: u16,
        mapq: u8,
        span: usize,
        qual: u8,
        tid: i32,
        mtid: i32,
        softclip: u32,
    ) -> Record {
        let mut cigar = Vec::new();
        if softclip > 0 {
            cigar.push(Cigar::SoftClip(softclip));
        }
        if span > 0 {
            cigar.push(Cigar::Match(span as u32));
        }
        let qlen = span + softclip as usize;
        let seq = vec![b'A'; qlen];
        let quals = vec![qual; qlen];
        let cs = CigarString(cigar);
        let mut rec = Record::new();
        rec.set(b"q1", Some(&cs), &seq, &quals);
        rec.set_flags(flags);
        rec.set_mapq(mapq);
        rec.set_tid(tid);
        rec.set_pos(1000);
        rec.set_mtid(mtid);
        rec.set_mpos(1000);
        rec
    }

    fn filt() -> NoisyFilter {
        NoisyFilter::default()
    }

    // ── is_read_noisy parity (one branch each) ───────────────────────────────

    #[test]
    fn clean_paired_read_is_not_noisy() {
        let hv = header_view();
        let r = make(PAIRED | PROPER | READ1, 60, 100, 30, 0, 0, 0);
        assert!(!is_read_noisy(&r, &hv, &filt()));
    }

    #[test]
    fn secondary_supplementary_never_noisy() {
        let hv = header_view();
        // Even unmapped, a secondary/supp read returns false (py 126-128).
        let sec = make(PAIRED | SECONDARY | UNMAPPED, 0, 0, 5, 0, 0, 0);
        assert!(!is_read_noisy(&sec, &hv, &filt()));
        let supp = make(PAIRED | SUPPLEMENTARY | UNMAPPED, 0, 0, 5, 0, 0, 0);
        assert!(!is_read_noisy(&supp, &hv, &filt()));
    }

    #[test]
    fn unmapped_is_noisy() {
        let hv = header_view();
        let r = make(PAIRED | UNMAPPED, 60, 100, 30, 0, 0, 0);
        assert!(is_read_noisy(&r, &hv, &filt()));
    }

    #[test]
    fn qcfail_is_noisy() {
        let hv = header_view();
        let r = make(PAIRED | PROPER | QCFAIL, 60, 100, 30, 0, 0, 0);
        assert!(is_read_noisy(&r, &hv, &filt()));
    }

    #[test]
    fn low_mapq_is_noisy() {
        let hv = header_view();
        let r = make(PAIRED | PROPER, 5, 100, 30, 0, 0, 0); // mapq 5 < 10
        assert!(is_read_noisy(&r, &hv, &filt()));
    }

    #[test]
    fn short_span_is_noisy_when_paired() {
        let hv = header_view();
        let r = make(PAIRED | PROPER, 60, 50, 30, 0, 0, 0); // span 50 < 75
        assert!(is_read_noisy(&r, &hv, &filt()));
    }

    #[test]
    fn mate_on_other_contig_is_noisy() {
        let hv = header_view();
        let r = make(PAIRED | PROPER, 60, 100, 30, 0, 1, 0); // tid 0, mtid 1
        assert!(is_read_noisy(&r, &hv, &filt()));
    }

    #[test]
    fn low_median_basequal_is_noisy() {
        let hv = header_view();
        let r = make(PAIRED | PROPER, 60, 100, 10, 0, 0, 0); // median 10 <= 15
        assert!(is_read_noisy(&r, &hv, &filt()));
    }

    #[test]
    fn many_low_qual_bases_is_noisy() {
        let hv = header_view();
        // median high enough to pass the median check, but ≥75 bases below cut.
        // Build a read where most bases are below 15 but median is high: not
        // possible — instead test the #low branch directly with a long low run.
        // 100 bases all at Q14: median 14 <= 15 would trip the median check
        // first. To isolate the count branch we need median > 15 AND ≥75 low.
        // Use 200 bases: 80 at Q5, 120 at Q40 → median = Q40 (>15), #low = 80 ≥75.
        let span = 200usize;
        let mut quals = vec![40u8; span];
        for q in quals.iter_mut().take(80) {
            *q = 5;
        }
        let seq = vec![b'A'; span];
        let cs = CigarString(vec![Cigar::Match(span as u32)]);
        let mut r = Record::new();
        r.set(b"q1", Some(&cs), &seq, &quals);
        r.set_flags(PAIRED | PROPER);
        r.set_mapq(60);
        r.set_tid(0);
        r.set_pos(1000);
        r.set_mtid(0);
        assert_eq!(fast_median(&quals), 40.0);
        assert!(is_read_noisy(&r, &hv, &filt()));
    }

    #[test]
    fn long_softclip_is_noisy() {
        let hv = header_view();
        let r = make(PAIRED | PROPER, 60, 100, 30, 0, 0, 75); // 75 softclip ≥ 75
        assert!(is_read_noisy(&r, &hv, &filt()));
    }

    #[test]
    fn single_end_duplicate_is_noisy() {
        let hv = header_view();
        let mut f = filt();
        f.paired = false;
        let r = make(DUP, 60, 100, 30, 0, 0, 0);
        assert!(is_read_noisy(&r, &hv, &f));
    }

    #[test]
    fn filter_noisy_off_skips_quality_checks() {
        let hv = header_view();
        let mut f = filt();
        f.filter_noisy = false;
        // low median baseQ would be noisy with filter_noisy, but not without.
        let r = make(PAIRED | PROPER, 60, 100, 5, 0, 0, 0);
        assert!(!is_read_noisy(&r, &hv, &f));
    }

    // ── fast_median ──────────────────────────────────────────────────────────

    #[test]
    fn median_matches_numpy() {
        assert_eq!(fast_median(&[1, 2, 3, 4, 5]), 3.0);
        assert_eq!(fast_median(&[1, 2, 3, 4]), 2.5);
        assert_eq!(fast_median(&[]), 0.0);
    }

    // ── build_bam_index end-to-end (synthetic paired BAM, no collation needed:
    //    records are written already grouped by qname) ─────────────────────────

    fn pair(qname: &[u8], tid: i32, pos: i64, read1: bool) -> Record {
        let span = 100u32;
        let seq = vec![b'A'; span as usize];
        let quals = vec![40u8; span as usize];
        let cs = CigarString(vec![Cigar::Match(span)]);
        let mut r = Record::new();
        r.set(qname, Some(&cs), &seq, &quals);
        let flags = PAIRED | PROPER | if read1 { READ1 } else { READ2 };
        r.set_flags(flags);
        r.set_mapq(60);
        r.set_tid(tid);
        r.set_pos(pos);
        r.set_mtid(tid);
        r.set_mpos(pos + 200);
        r
    }

    #[test]
    fn build_index_groups_pairs() {
        // Three clean qnames (need >2 to clear the InsufficientPairs gate), each
        // R1+R2 adjacent on chr1. samtools collate may or may not be present; we
        // write records already grouped so both the collated and the in-memory
        // fallback paths produce the same result.
        let tmp = tempfile::Builder::new().suffix(".bam").tempfile().unwrap();
        {
            let mut w = Writer::from_path(tmp.path(), &test_header(), Format::Bam).unwrap();
            for (i, name) in [b"qa".as_slice(), b"qb".as_slice(), b"qc".as_slice()]
                .iter()
                .enumerate()
            {
                let pos = 1000 + (i as i64) * 500;
                w.write(&pair(name, 0, pos, true)).unwrap();
                w.write(&pair(name, 0, pos, false)).unwrap();
            }
        }
        let idx = build_bam_index(tmp.path(), &filt(), 1).unwrap();
        assert_eq!(idx.qname_idx.len(), 3, "three qnames retained");
        assert!(idx.noisy.is_empty(), "none noisy");
        // every retained qname has 2 reads
        for recs in idx.reads.values() {
            assert_eq!(recs.len(), 2);
        }
        // chr1 lapper has 6 intervals (2 per qname)
        let lap = idx.lapper.get("chr1").expect("chr1 lapper");
        let hits: Vec<_> = lap.find(1000, 1100).map(|iv| iv.val).collect();
        assert!(!hits.is_empty(), "query should find overlapping qnames");
    }

    #[test]
    fn build_index_too_few_pairs_errors() {
        let tmp = tempfile::Builder::new().suffix(".bam").tempfile().unwrap();
        {
            let mut w = Writer::from_path(tmp.path(), &test_header(), Format::Bam).unwrap();
            w.write(&pair(b"qa", 0, 1000, true)).unwrap();
            w.write(&pair(b"qa", 0, 1000, false)).unwrap();
        }
        // only 1 qname ≤ 2 → InsufficientPairs
        let err = build_bam_index(tmp.path(), &filt(), 1).unwrap_err();
        assert!(matches!(err, SdError::InsufficientPairs(_)), "got {err:?}");
    }

    #[test]
    fn merge_bams_is_stubbed() {
        let err = merge_bams(&[], Path::new("/dev/null"), Path::new("/dev/null"), 1).unwrap_err();
        assert!(matches!(err, SdError::Htslib(_)), "got {err:?}");
    }
}
