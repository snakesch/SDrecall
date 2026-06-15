//! Masked-genome construction — N-end masking + N-bridge contig merge.
//!
//! Ports `preparation/genome.py::Genome.mask` (l.42-152). The masked FASTA md5 is
//! a pass criterion, so the coordinate-bearing logic must be exact.
//!
//! ## Port status
//!
//! - [`apply_end_mask`] — the per-contig N-end mask: replace the first and last
//!   1000 bp with `N`, keeping `[1000:-1000]` of the slice (genome.py l.80-87).
//!   The exact slice arithmetic (`mask + seq[1000:-1000] + mask`) is load-bearing.
//!   Ported + UNIT-TESTED on byte sequences.
//! - [`merge_nearby_contigs`] — the same-chrom N-bridge merge: contigs within
//!   `max_gap` bp are concatenated with an **exact** N-bridge of length
//!   `next_start - (cur_start + len(cur_seq))` (genome.py l.89-152), preserving the
//!   `start + local_pos` coordinate mapping. Ported + UNIT-TESTED.
//! - **`mask_genome` (FASTA read/write) — STUBBED** (`TODO(T8)`): the
//!   `Fasta[chrom][start:stop]` slice + `SeqIO.write` + md5-gated update
//!   (genome.py l.51-61, l.154-167) via `bio::io::fasta::IndexedReader::fetch` +
//!   `bio::io::fasta::Writer` (the `bio` dep IS wired into the build and compiles;
//!   only this orchestration glue is deferred). The contig-merge + end-mask logic
//!   (the only coordinate-bearing parts) are ported + tested.

/// One masked contig: its chromosome, genomic start (the `{chrom}:{start}` FASTA
/// id), and sequence bytes. `start + local_position` maps a base back to the
/// genome (the invariant `modify_masked_genome_coords` relies on).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MaskedContig {
    pub chrom: String,
    pub start: i64,
    pub seq: Vec<u8>,
}

/// Apply the 1000-bp N-end mask to a slice (`Genome._apply_mask`, genome.py
/// l.80-87): `"N"*1000 + seq[1000:-1000] + "N"*1000`.
///
/// If the slice is `<= 2000 bp` the Python slice `seq[1000:-1000]` is empty (or
/// the start/stop cross), yielding `"N"*2000` of total length `min(len, 2000)`...
/// actually Python's `mask + ""[..] + mask` is always exactly `2000 + max(0,
/// len-2000)` long because `seq[1000:-1000]` keeps the middle. We reproduce the
/// Python slice semantics precisely: the kept middle is `seq[1000 .. len-1000]`
/// when `len > 2000`, else empty, and the result is `N*1000 + middle + N*1000`.
///
/// `&[u8]` borrow-in; owned `Vec<u8>` out. The result length equals `max(2000,
/// len)` only when `len >= 2000`; for `len < 2000` the result is `2*1000 +
/// max(0, len-2000)` = `2000` (the middle is empty) — matching Python where
/// `seq[1000:-1000]` on a short string is `""`.
pub fn apply_end_mask(seq: &[u8]) -> Vec<u8> {
    let n = seq.len();
    let mask = vec![b'N'; 1000];
    let mut out = Vec::with_capacity(n.max(2000));
    out.extend_from_slice(&mask);
    if n > 2000 {
        out.extend_from_slice(&seq[1000..n - 1000]);
    }
    // else: Python seq[1000:-1000] is empty for len <= 2000.
    out.extend_from_slice(&mask);
    out
}

/// Merge same-chromosome contigs within `max_gap` bp by inserting an exact-length
/// N-bridge (`Genome._merge_nearby_contigs`, genome.py l.89-152).
///
/// The bridge length between an accumulating contig and the next is
/// `next_start - (cur_start + cur_seq.len())`. Behavior by gap:
/// - `bridge < 0` (overlap, should not happen post sort+merge+slop): emit the
///   current contig and start a fresh one (Python logs a warning).
/// - `0 <= bridge <= max_gap`: concatenate `cur + N*bridge + next`.
/// - `bridge > max_gap`: emit current, start fresh.
///
/// Contigs are grouped by chrom (chroms processed in sorted order) and within a
/// chrom sorted by genomic start. `&[MaskedContig]` borrow-in; owned out.
pub fn merge_nearby_contigs(contigs: &[MaskedContig], max_gap: i64) -> Vec<MaskedContig> {
    use std::collections::BTreeMap;
    // Group by chrom (BTreeMap → sorted chrom order, matching Python's
    // `sorted(by_chrom.keys())`).
    let mut by_chrom: BTreeMap<String, Vec<&MaskedContig>> = BTreeMap::new();
    for c in contigs {
        by_chrom.entry(c.chrom.clone()).or_default().push(c);
    }

    let mut merged: Vec<MaskedContig> = Vec::new();
    for (chrom, mut group) in by_chrom {
        group.sort_by_key(|c| c.start);
        let mut cur_start = group[0].start;
        let mut cur_seq = group[0].seq.clone();

        for next in &group[1..] {
            let next_start = next.start;
            let cur_end = cur_start + cur_seq.len() as i64;
            let bridge = next_start - cur_end;
            if bridge < 0 {
                log::warn!(
                    "Contigs {chrom}:{cur_start} and {chrom}:{next_start} overlap by {}bp -- emitting separately",
                    -bridge
                );
                merged.push(MaskedContig {
                    chrom: chrom.clone(),
                    start: cur_start,
                    seq: std::mem::take(&mut cur_seq),
                });
                cur_start = next_start;
                cur_seq = next.seq.clone();
            } else if bridge <= max_gap {
                cur_seq.extend(std::iter::repeat_n(b'N', bridge as usize));
                cur_seq.extend_from_slice(&next.seq);
            } else {
                merged.push(MaskedContig {
                    chrom: chrom.clone(),
                    start: cur_start,
                    seq: std::mem::take(&mut cur_seq),
                });
                cur_start = next_start;
                cur_seq = next.seq.clone();
            }
        }
        merged.push(MaskedContig {
            chrom: chrom.clone(),
            start: cur_start,
            seq: cur_seq,
        });
    }
    log::info!(
        "Contig merging (max_gap={max_gap}bp): {} contigs -> {} merged contigs",
        contigs.len(),
        merged.len()
    );
    merged
}

/// Build the masked genome FASTA for one RG, the in-process port of
/// `Genome.mask` (genome.py l.42-152) + `_mask_intervals`/`_write_masked_genome`.
///
/// Pipeline (parity-critical — the output FASTA md5 is a pass criterion):
/// 1. sort+merge the `query_bed`, slop by `avg_frag + 2*std_frag + 1000` clamped to
///    contig sizes, sort+merge again (genome.py l.68-73). One slop+merge unit via
///    [`sdrecall_io::sort_merge_bed`] + [`sdrecall_io::slop`].
/// 2. For each padded interval `[start, stop)` fetch `ref[chrom][start:stop]` via
///    `bio::io::fasta::IndexedReader` (faidx random access) → a contig with id
///    `{chrom}:{start}` (genome.py l.77-78).
/// 3. [`apply_end_mask`] each contig (`N*1000 + seq[1000:-1000] + N*1000`).
/// 4. if `merge_gap > 0`, [`merge_nearby_contigs`] (same-chrom N-bridge merge).
/// 5. write a 60-char-wrapped FASTA (Biopython `SeqIO.write(..., "fasta")` default)
///    with `desc=None` so headers are `>{chrom}:{start}` exactly. md5-gate the
///    update (only replace `out` when its md5 differs — `update_plain_file_on_md5`).
///
/// `contig_sizes` is the reference `.fai` (chrom → length) used for slop clamping
/// (Python passes the `.contigsize.genome` derived from the `.fai`). `&Path`s
/// borrow; the output FASTA is written at `out`.
pub fn mask_genome(
    query_bed: &[sdrecall_utils::GenomicInterval],
    ref_fa: &std::path::Path,
    out: &std::path::Path,
    avg_frag: f64,
    std_frag: f64,
    merge_gap: i64,
) -> sdrecall_utils::Result<()> {
    use sdrecall_utils::SdError;

    // ---- 1. sort+merge → slop → sort+merge (genome.py l.68-73) ----
    // Python pads by int(avg + 2*std + 1000) on both sides; slop clamps to contig
    // bounds (start>=0, end<=contig_len). The .fai gives contig sizes.
    let fai_path = format!("{}.fai", ref_fa.display());
    let contig_sizes = read_fai_sizes(std::path::Path::new(&fai_path))?;
    let pad = (avg_frag + 2.0 * std_frag + 1000.0) as i64;

    let merged = sdrecall_io::sort_merge_bed(query_bed, false);
    // slop grows by `pad` on both sides, clamped to [0, contig_len].
    let slopped = sdrecall_io::slop(&merged, pad, &contig_sizes)?;
    let padded = sdrecall_io::sort_merge_bed(&slopped, false);

    // ---- 2-3. fetch each interval, end-mask, build contigs ----
    let mut reader = bio::io::fasta::IndexedReader::from_file(&ref_fa).map_err(|e| {
        SdError::Io {
            path: ref_fa.display().to_string(),
            source: std::io::Error::other(e.to_string()),
        }
    })?;
    let mut contigs: Vec<MaskedContig> = Vec::with_capacity(padded.len());
    let mut buf: Vec<u8> = Vec::new();
    for iv in &padded {
        // Python uses pyfaidx slice ref[chrom][start:stop] (0-based half-open).
        reader
            .fetch(&iv.chrom, iv.start.max(0) as u64, iv.end.max(0) as u64)
            .map_err(|e| SdError::Io {
                path: format!("{}:{}-{}", iv.chrom, iv.start, iv.end),
                source: e,
            })?;
        buf.clear();
        reader.read(&mut buf).map_err(|e| SdError::Io {
            path: format!("{}:{}-{}", iv.chrom, iv.start, iv.end),
            source: e,
        })?;
        let masked = apply_end_mask(&buf);
        contigs.push(MaskedContig {
            chrom: iv.chrom.clone(),
            start: iv.start,
            seq: masked,
        });
    }

    // ---- 4. merge nearby same-chrom contigs ----
    if merge_gap > 0 {
        contigs = merge_nearby_contigs(&contigs, merge_gap);
    }

    // ---- 5. write 60-wrapped FASTA to a temp, then md5-gate ----
    let tmp = format!("{}.{}.tmp.fasta", out.display(), std::process::id());
    write_wrapped_fasta(std::path::Path::new(&tmp), &contigs)?;
    md5_gated_replace(std::path::Path::new(&tmp), out)?;
    // Build the .fai for the masked genome so downstream intrinsic alignment can
    // open it (Python re-initialises Genome(masked) which runs `samtools faidx`).
    build_fai(out)?;
    Ok(())
}

/// Parse a reference `.fai` into a chrom→length map (the slop clamp bounds).
fn read_fai_sizes(fai: &std::path::Path) -> sdrecall_utils::Result<ahash::AHashMap<String, i64>> {
    use sdrecall_utils::SdError;
    let text = std::fs::read_to_string(fai).map_err(|e| SdError::Io {
        path: fai.display().to_string(),
        source: e,
    })?;
    let mut m = ahash::AHashMap::new();
    for line in text.lines() {
        let mut f = line.split('\t');
        if let (Some(name), Some(len)) = (f.next(), f.next()) {
            if let Ok(len) = len.parse::<i64>() {
                m.insert(name.to_string(), len);
            }
        }
    }
    Ok(m)
}

/// Write contigs as a 60-char-wrapped FASTA with header `>{chrom}:{start}` and no
/// description — byte-identical to Biopython `SeqIO.write(records, "fasta")` where
/// each record has `id="{chrom}:{start}"`, `description=""`. The md5 criterion
/// depends on the 60-char wrap + trailing newline per line.
fn write_wrapped_fasta(
    path: &std::path::Path,
    contigs: &[MaskedContig],
) -> sdrecall_utils::Result<()> {
    use sdrecall_utils::SdError;
    use std::io::Write;
    let f = std::fs::File::create(path).map_err(|e| SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    let mut w = std::io::BufWriter::new(f);
    let io_err = |e: std::io::Error| SdError::Io {
        path: path.display().to_string(),
        source: e,
    };
    for c in contigs {
        writeln!(w, ">{}:{}", c.chrom, c.start).map_err(io_err)?;
        for chunk in c.seq.chunks(60) {
            w.write_all(chunk).map_err(io_err)?;
            w.write_all(b"\n").map_err(io_err)?;
        }
    }
    w.flush().map_err(io_err)?;
    Ok(())
}

/// md5-gated replace (`update_plain_file_on_md5`, src/utils.py): only move `tmp`
/// onto `out` when their md5 differs (or `out` is absent). Avoids touching the
/// file mtime when content is unchanged (the freshness check downstream relies on
/// this). The tmp is removed either way.
fn md5_gated_replace(
    tmp: &std::path::Path,
    out: &std::path::Path,
) -> sdrecall_utils::Result<()> {
    use sdrecall_utils::SdError;
    let new_md5 = file_md5(tmp)?;
    let same = out
        .exists()
        .then(|| file_md5(out).ok())
        .flatten()
        .map(|old| old == new_md5)
        .unwrap_or(false);
    if same {
        let _ = std::fs::remove_file(tmp);
    } else {
        std::fs::rename(tmp, out).map_err(|e| SdError::Io {
            path: out.display().to_string(),
            source: e,
        })?;
    }
    Ok(())
}

/// Streaming md5 hex digest of a file (no full-file buffering).
fn file_md5(path: &std::path::Path) -> sdrecall_utils::Result<String> {
    use sdrecall_utils::SdError;
    use std::io::Read;
    let mut f = std::fs::File::open(path).map_err(|e| SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    let mut ctx = md5::Context::new();
    let mut buf = [0u8; 65536];
    loop {
        let n = f.read(&mut buf).map_err(|e| SdError::Io {
            path: path.display().to_string(),
            source: e,
        })?;
        if n == 0 {
            break;
        }
        ctx.consume(&buf[..n]);
    }
    Ok(format!("{:x}", ctx.compute()))
}

/// Build the `.fai` for a written FASTA via samtools faidx (the same external tool
/// Python's `Genome._create_index` shells; sdrecall-io itself shells samtools, so
/// this is consistent tool orchestration, not a re-implemented algorithm).
fn build_fai(fasta: &std::path::Path) -> sdrecall_utils::Result<()> {
    use sdrecall_utils::SdError;
    let out = std::process::Command::new("samtools")
        .args(["faidx", &fasta.display().to_string()])
        .output()
        .map_err(|e| SdError::Io {
            path: "samtools faidx".to_string(),
            source: e,
        })?;
    if !out.status.success() {
        return Err(SdError::Compute(format!(
            "samtools faidx {} failed: {}",
            fasta.display(),
            String::from_utf8_lossy(&out.stderr)
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn contig(chrom: &str, start: i64, seq: &str) -> MaskedContig {
        MaskedContig {
            chrom: chrom.into(),
            start,
            seq: seq.as_bytes().to_vec(),
        }
    }

    #[test]
    fn end_mask_long_sequence() {
        // length 2500: keep [1000:1500] of the middle, N*1000 either side → 2500.
        let seq: Vec<u8> = (0..2500).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        let masked = apply_end_mask(&seq);
        assert_eq!(masked.len(), 2500);
        assert!(masked[..1000].iter().all(|&b| b == b'N'));
        assert!(masked[1500..].iter().all(|&b| b == b'N'));
        // middle preserved == original[1000..1500].
        assert_eq!(&masked[1000..1500], &seq[1000..1500]);
    }

    #[test]
    fn end_mask_short_sequence_all_n() {
        // length 1500 (<2000): middle empty → result is 2000 Ns.
        let seq = vec![b'A'; 1500];
        let masked = apply_end_mask(&seq);
        assert_eq!(masked.len(), 2000);
        assert!(masked.iter().all(|&b| b == b'N'));
    }

    #[test]
    fn merge_bridges_within_gap() {
        // contig1 @100 len 10, contig2 @130 len 5. cur_end=110, bridge=130-110=20.
        // max_gap 1000 → merge: seq1 + N*20 + seq2, length 10+20+5=35.
        let c1 = contig("chr1", 100, "AAAAAAAAAA"); // 10
        let c2 = contig("chr1", 130, "CCCCC"); // 5
        let merged = merge_nearby_contigs(&[c1, c2], 1000);
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].start, 100);
        assert_eq!(merged[0].seq.len(), 35);
        assert_eq!(&merged[0].seq[10..30], &vec![b'N'; 20][..]);
        assert_eq!(&merged[0].seq[30..], b"CCCCC");
    }

    #[test]
    fn merge_gap_too_large_keeps_separate() {
        // bridge 2000 > max_gap 1000 → two contigs.
        let c1 = contig("chr1", 100, "AAAAAAAAAA");
        let c2 = contig("chr1", 2110, "CCCCC");
        let merged = merge_nearby_contigs(&[c1, c2], 1000);
        assert_eq!(merged.len(), 2);
        assert_eq!(merged[0].start, 100);
        assert_eq!(merged[1].start, 2110);
    }

    #[test]
    fn merge_different_chroms_independent() {
        let c1 = contig("chr1", 100, "AAAAA");
        let c2 = contig("chr2", 110, "CCCCC");
        let merged = merge_nearby_contigs(&[c1, c2], 1000);
        assert_eq!(merged.len(), 2);
        // sorted chrom order.
        assert_eq!(merged[0].chrom, "chr1");
        assert_eq!(merged[1].chrom, "chr2");
    }

    #[test]
    fn merge_overlap_emits_separately() {
        // contig1 @100 len 50 (end 150), contig2 @120 (bridge = 120-150 = -30 < 0).
        let c1 = contig("chr1", 100, &"A".repeat(50));
        let c2 = contig("chr1", 120, "CCCCC");
        let merged = merge_nearby_contigs(&[c1, c2], 1000);
        assert_eq!(merged.len(), 2);
    }

    #[test]
    fn mask_genome_end_mask_and_wrap_bytes() {
        use sdrecall_utils::GenomicInterval;
        use std::io::Write;
        // Build a tiny 1-contig reference: chrZ, length 4000 of repeated ACGT.
        let dir = tempfile::tempdir().unwrap();
        let ref_fa = dir.path().join("ref.fasta");
        let seq: Vec<u8> = (0..4000).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        {
            let mut f = std::fs::File::create(&ref_fa).unwrap();
            writeln!(f, ">chrZ").unwrap();
            for chunk in seq.chunks(60) {
                f.write_all(chunk).unwrap();
                f.write_all(b"\n").unwrap();
            }
        }
        // Need a .fai; build via samtools (test env has it). Skip test if absent.
        let ok = std::process::Command::new("samtools")
            .args(["faidx", ref_fa.to_str().unwrap()])
            .status()
            .map(|s| s.success())
            .unwrap_or(false);
        if !ok {
            eprintln!("samtools unavailable; skipping mask_genome byte test");
            return;
        }

        // Query a single interval [1000, 3000) (size 2000). With pad=0 (avg/std 0,
        // +1000 base pad) the slop would grow it; to assert exact bytes we use a
        // merge_gap of 0 and pad-free check by querying the produced FASTA structure
        // instead: just assert the contig is fully N-masked where expected.
        // Use avg=0,std=0 → pad = 1000; interval [1000,3000) → slop [0,4000) clamped.
        let q = vec![GenomicInterval::with_strand("chrZ", 1000, 3000, sdrecall_utils::Strand::Forward)];
        let out = dir.path().join("masked.fasta");
        mask_genome(&q, &ref_fa, &out, 0.0, 0.0, 1000).unwrap();

        let content = std::fs::read_to_string(&out).unwrap();
        // One contig, header is >chrZ:{slop_start}. slop start = max(1000-1000,0)=0.
        assert!(content.starts_with(">chrZ:0\n"), "header was: {}", &content[..20]);
        // Lines after header are 60-wide. The first 1000 bases are N (end mask).
        let body: String = content.lines().skip(1).collect();
        assert!(body.as_bytes()[..1000].iter().all(|&b| b == b'N'), "first 1000 should be N");
        // The masked contig length = slop span = [0,4000) clamped → 4000.
        assert_eq!(body.len(), 4000, "masked contig length");
        assert!(body.as_bytes()[3000..].iter().all(|&b| b == b'N'), "last 1000 should be N");
        // Re-running with identical content must not change the file (md5-gate).
        let m1 = file_md5(&out).unwrap();
        mask_genome(&q, &ref_fa, &out, 0.0, 0.0, 1000).unwrap();
        let m2 = file_md5(&out).unwrap();
        assert_eq!(m1, m2, "md5-gate keeps content stable on re-run");
    }

    #[test]
    fn merge_three_contigs_chain() {
        // @100 len10 (end110); @120 len10 bridge10 -> merge (end now 100..130 len30);
        // @140 len10 bridge=140-(100+30)=10 -> merge. Total len 10+10+10+10+10=50.
        let c1 = contig("chr1", 100, &"A".repeat(10));
        let c2 = contig("chr1", 120, &"C".repeat(10));
        let c3 = contig("chr1", 140, &"G".repeat(10));
        let merged = merge_nearby_contigs(&[c1, c2, c3], 1000);
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].seq.len(), 50);
    }
}
