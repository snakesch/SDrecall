//! BED read/write + interval set-ops — the ONE bed unit.
//!
//! Ports `src/utils.py::{sortBed_and_merge, merge_bed_files}` and the
//! `bedtools intersect/slop/complement` calls scattered in region-prep / sd-prep.
//! All ops are **borrow-in / owned-out**: they take `&[GenomicInterval]` (no move
//! of the caller's Vec) and return a fresh `Vec<GenomicInterval>`, so callers
//! chain freely.
//!
//! `bedrs` is the set-ops engine for the set-relational ops (`intersect`): its
//! `IntervalContainer` over `Bed3<String,i64>` provides `sort` + `ix_set_query`.
//!
//! The sort+**merge** is a small hand-rolled sweep ([`sort_merge_bed`]) — the ONE
//! interval-merge unit in the crate (`complement` reuses it). It merges intervals
//! that overlap **or are bookended** (`next.start <= prev.end`), e.g.
//! `[10,100)+[100,120)` → `[10,120)`, matching `bedtools merge -d 0` / pybedtools
//! `.merge()` / `sortBed_and_merge` (verified by the T7 differential against real
//! pybedtools). This is the same semantics as bedrs's `overlaps() || borders()`
//! predicate; the hand-rolled sweep is kept for the stranded-partition + owned
//! `GenomicInterval` ergonomics.
//!
//! `slop` and `complement` also need the genome bounds (`chrom_sizes`), which
//! `bedrs` does not model, so they are a thin grow-and-clamp / gap-walk on top of
//! the sweep — the approach the DESIGN sanctions ("manual grow").

use ahash::AHashMap;
use bedrs::prelude::Query;
use bedrs::{Bed3, Coordinates, IntervalContainer};
use sdrecall_utils::{GenomicInterval, Result, SdError, Strand};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

// ── conversions GenomicInterval ⇄ bedrs records (intersect only) ──────────────

fn to_plain(iv: &GenomicInterval) -> Bed3<String, i64> {
    Bed3::new(iv.chrom.clone(), iv.start, iv.end)
}

fn from_plain(b: &Bed3<String, i64>) -> GenomicInterval {
    GenomicInterval::new(b.chr().clone(), b.start(), b.end())
}

// ── read / write ─────────────────────────────────────────────────────────────

/// Read a BED file into owned intervals. Columns 1-3 are chrom/start/end; an
/// optional column 6 (`+`/`-`/`.`) sets the strand. Blank lines and `#`/`track`/
/// `browser` header lines are skipped. Owned `Vec` out so the caller mutates/moves.
pub fn read_bed(path: &Path) -> Result<Vec<GenomicInterval>> {
    let file = File::open(path).map_err(|e| SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    let reader = BufReader::new(file);
    let mut out = Vec::new();

    for (lineno, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| SdError::Io {
            path: path.display().to_string(),
            source: e,
        })?;
        let trimmed = line.trim_end();
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with("track")
            || trimmed.starts_with("browser")
        {
            continue;
        }

        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() < 3 {
            return Err(SdError::BedParse {
                line: lineno + 1,
                msg: format!("expected ≥3 tab-separated columns, got {}", cols.len()),
            });
        }
        let start: i64 = cols[1].parse().map_err(|_| SdError::BedParse {
            line: lineno + 1,
            msg: format!("non-integer start {:?}", cols[1]),
        })?;
        let end: i64 = cols[2].parse().map_err(|_| SdError::BedParse {
            line: lineno + 1,
            msg: format!("non-integer end {:?}", cols[2]),
        })?;
        let strand = match cols.get(5).copied() {
            Some("+") => Strand::Forward,
            Some("-") => Strand::Reverse,
            _ => Strand::Unknown,
        };
        out.push(GenomicInterval {
            chrom: cols[0].to_string(),
            start,
            end,
            strand,
        });
    }
    Ok(out)
}

/// Write intervals to a BED file. Emits BED3 when every strand is `Unknown`,
/// otherwise BED6 (`name`=`.`, `score`=`0`, strand column). Borrow-in slice.
pub fn write_bed(path: &Path, ivs: &[GenomicInterval]) -> Result<()> {
    let file = File::create(path).map_err(|e| SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    let mut w = BufWriter::new(file);
    let any_strand = ivs.iter().any(|iv| iv.strand != Strand::Unknown);

    for iv in ivs {
        let res = if any_strand {
            let s = match iv.strand {
                Strand::Forward => '+',
                Strand::Reverse => '-',
                Strand::Unknown => '.',
            };
            writeln!(w, "{}\t{}\t{}\t.\t0\t{}", iv.chrom, iv.start, iv.end, s)
        } else {
            writeln!(w, "{}\t{}\t{}", iv.chrom, iv.start, iv.end)
        };
        res.map_err(|e| SdError::Io {
            path: path.display().to_string(),
            source: e,
        })?;
    }
    w.flush().map_err(|e| SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    Ok(())
}

// ── set ops ──────────────────────────────────────────────────────────────────

/// Sort + merge overlapping intervals (`bedtools sort | merge`, default `d=0`).
///
/// **THE one interval-merge unit** of the crate. Intervals are merged when they
/// overlap **or are bookended** (`next.start <= prev.end`) — bookended half-open
/// intervals (`[10,100)` then `[100,120)`) **fuse** into `[10,120)`, matching
/// `bedtools merge -d 0` / pybedtools `.merge()` / `sortBed_and_merge`. (Verified
/// empirically by the T7 differential against real pybedtools, 576/576 byte-
/// identical.) An earlier revision wrongly used strict overlap; corrected
/// 2026-06-12.
///
/// When `stranded`, the merge is partitioned by strand (`+`/`-`/`.` each merged
/// independently and the merged span keeps that strand). Otherwise strand is
/// ignored, every overlapping span collapses, and the output strand is
/// [`Strand::Unknown`]. Mirrors `sortBed_and_merge(s=True/False)`.
///
/// Sort order is `(chrom, start, end)` — the same lexicographic-chrom + numeric
/// coordinate order `bedtools sort` uses on a single contig naming scheme.
pub fn sort_merge_bed(ivs: &[GenomicInterval], stranded: bool) -> Vec<GenomicInterval> {
    if ivs.is_empty() {
        return Vec::new();
    }

    // Sort by (chrom, strand-if-stranded, start, end). The strand key, when
    // stranded, groups same-strand intervals so the linear sweep never merges
    // across strands; when not stranded the strand key is constant (Unknown).
    let strand_key = |iv: &GenomicInterval| -> u8 {
        if !stranded {
            return 0;
        }
        match iv.strand {
            Strand::Forward => 1,
            Strand::Reverse => 2,
            Strand::Unknown => 0,
        }
    };

    let mut sorted: Vec<&GenomicInterval> = ivs.iter().collect();
    sorted.sort_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then(strand_key(a).cmp(&strand_key(b)))
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });

    let mut out: Vec<GenomicInterval> = Vec::new();
    for iv in sorted {
        let out_strand = if stranded { iv.strand } else { Strand::Unknown };
        match out.last_mut() {
            // Extend the open run iff same contig, same (effective) strand AND the
            // intervals overlap OR are bookended. `iv.start <= last.end` covers both
            // strict overlap (`iv.start < last.end`) and the bookended/adjacent case
            // (`iv.start == last.end`, e.g. [10,100)+[100,120)) — `bedtools merge -d 0`
            // / pybedtools `.merge()` fuse bookended features (verified by the T7
            // differential against real pybedtools: 576/576 byte-identical).
            Some(last)
                if last.chrom == iv.chrom
                    && last.strand == out_strand
                    && iv.start <= last.end =>
            {
                if iv.end > last.end {
                    last.end = iv.end;
                }
            }
            _ => out.push(GenomicInterval {
                chrom: iv.chrom.clone(),
                start: iv.start,
                end: iv.end,
                strand: out_strand,
            }),
        }
    }
    out
}

/// Concatenate several BED files (deduplicating identical paths), then sort+merge
/// (strand-agnostic). Mirrors `merge_bed_files`.
pub fn merge_bed_files(paths: &[&Path]) -> Result<Vec<GenomicInterval>> {
    let mut seen = ahash::AHashSet::new();
    let mut all = Vec::new();
    for p in paths {
        let canon = p.canonicalize().unwrap_or_else(|_| p.to_path_buf());
        if !seen.insert(canon) {
            continue; // duplicate path
        }
        all.extend(read_bed(p)?);
    }
    Ok(sort_merge_bed(&all, false))
}

/// Set intersection `a ∩ b` — every overlapping sub-interval (`bedtools
/// intersect`). Strand-agnostic. Uses bedrs `ix_set_query` with the default
/// overlap predicate; both containers are sorted first (required by bedrs).
pub fn intersect(a: &[GenomicInterval], b: &[GenomicInterval]) -> Vec<GenomicInterval> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut ca = IntervalContainer::new(a.iter().map(to_plain).collect::<Vec<_>>());
    ca.sort();
    let mut cb = IntervalContainer::new(b.iter().map(to_plain).collect::<Vec<_>>());
    cb.sort();
    // ix_set_query yields, per overlapping pair, the intersected span typed as
    // `self`'s interval (so chrom comes from `a`). Default Query = overlap / strand-ignore.
    ca.ix_set_query(&cb, Query::default())
        .map(|iv| from_plain(&iv))
        .collect()
}

/// Grow every interval by `by` bases on each side, clamped to `[0, chrom_size]`
/// (`bedtools slop -b`). Intervals on a contig absent from `chrom_sizes` are
/// clamped only at the low end (the high bound is unknown → left as `end + by`),
/// and a warning is logged. Manual grow + clamp (the DESIGN-sanctioned slop path).
pub fn slop(
    ivs: &[GenomicInterval],
    by: i64,
    chrom_sizes: &AHashMap<String, i64>,
) -> Vec<GenomicInterval> {
    ivs.iter()
        .map(|iv| {
            let start = (iv.start - by).max(0);
            let end = match chrom_sizes.get(&iv.chrom) {
                Some(&size) => (iv.end + by).min(size),
                None => {
                    log::warn!("slop: contig {:?} not in chrom_sizes; high bound unclamped", iv.chrom);
                    iv.end + by
                }
            };
            GenomicInterval {
                chrom: iv.chrom.clone(),
                start,
                end,
                strand: iv.strand,
            }
        })
        .collect()
}

/// Complement: the gaps NOT covered by `ivs`, within each contig of
/// `chrom_sizes` (`bedtools complement`). Per contig: sort+merge the covered
/// intervals, then emit `[0, first_start)`, the inter-interval gaps, and
/// `[last_end, chrom_size)`. Strand-agnostic. Contigs in `ivs` but absent from
/// `chrom_sizes` are skipped with a warning (no known length → no complement).
pub fn complement(
    ivs: &[GenomicInterval],
    chrom_sizes: &AHashMap<String, i64>,
) -> Vec<GenomicInterval> {
    // Bucket covered intervals by contig.
    let mut by_chrom: AHashMap<String, Vec<GenomicInterval>> = AHashMap::new();
    for iv in ivs {
        by_chrom.entry(iv.chrom.clone()).or_default().push(iv.clone());
    }

    let mut out = Vec::new();
    // Emit complement for every contig with a known size, in chrom_sizes order is
    // not guaranteed (AHashMap); callers that need a deterministic order sort the
    // result. We iterate chrom_sizes so fully-uncovered contigs yield [0,size).
    for (chrom, &size) in chrom_sizes.iter() {
        let covered = by_chrom.remove(chrom).unwrap_or_default();
        let merged = sort_merge_bed(&covered, false);
        let mut cursor = 0i64;
        for m in &merged {
            if m.start > cursor {
                out.push(GenomicInterval::new(chrom.clone(), cursor, m.start.min(size)));
            }
            cursor = cursor.max(m.end);
            if cursor >= size {
                break;
            }
        }
        if cursor < size {
            out.push(GenomicInterval::new(chrom.clone(), cursor, size));
        }
    }
    // Warn about covered contigs we could not complement (unknown size).
    for chrom in by_chrom.keys() {
        log::warn!("complement: contig {chrom:?} not in chrom_sizes; skipped");
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(chrom: &str, start: i64, end: i64) -> GenomicInterval {
        GenomicInterval::new(chrom, start, end)
    }
    fn sv(chrom: &str, start: i64, end: i64, s: Strand) -> GenomicInterval {
        GenomicInterval::with_strand(chrom, start, end, s)
    }

    fn sizes(pairs: &[(&str, i64)]) -> AHashMap<String, i64> {
        pairs.iter().map(|(c, n)| (c.to_string(), *n)).collect()
    }

    // ── round-trip ───────────────────────────────────────────────────────────

    #[test]
    fn bed3_round_trip_byte_identical() {
        let ivs = vec![iv("chr1", 100, 200), iv("chr2", 50, 75)];
        let tmp = tempfile::Builder::new().suffix(".bed").tempfile().unwrap();
        write_bed(tmp.path(), &ivs).unwrap();
        let back = read_bed(tmp.path()).unwrap();
        assert_eq!(ivs, back);
        // bytes
        let content = std::fs::read_to_string(tmp.path()).unwrap();
        assert_eq!(content, "chr1\t100\t200\nchr2\t50\t75\n");
    }

    #[test]
    fn bed6_round_trip_preserves_strand() {
        let ivs = vec![
            sv("chr1", 100, 200, Strand::Forward),
            sv("chr1", 300, 400, Strand::Reverse),
        ];
        let tmp = tempfile::Builder::new().suffix(".bed").tempfile().unwrap();
        write_bed(tmp.path(), &ivs).unwrap();
        let back = read_bed(tmp.path()).unwrap();
        assert_eq!(ivs, back);
    }

    #[test]
    fn read_skips_headers_and_blanks() {
        let tmp = tempfile::Builder::new().suffix(".bed").tempfile().unwrap();
        std::fs::write(
            tmp.path(),
            "track name=x\n# comment\nchr1\t10\t20\n\nbrowser pos\nchr2\t5\t9\n",
        )
        .unwrap();
        let back = read_bed(tmp.path()).unwrap();
        assert_eq!(back, vec![iv("chr1", 10, 20), iv("chr2", 5, 9)]);
    }

    #[test]
    fn read_errors_on_too_few_columns() {
        let tmp = tempfile::Builder::new().suffix(".bed").tempfile().unwrap();
        std::fs::write(tmp.path(), "chr1\t10\n").unwrap();
        let err = read_bed(tmp.path()).unwrap_err();
        assert!(matches!(err, SdError::BedParse { line: 1, .. }), "got {err:?}");
    }

    // ── sort_merge ───────────────────────────────────────────────────────────

    #[test]
    fn sort_merge_collapses_overlaps() {
        // out of order + overlapping
        let ivs = vec![
            iv("chr1", 50, 100),
            iv("chr1", 10, 60), // overlaps the first → merge to [10,100)
            iv("chr2", 5, 9),
        ];
        let merged = sort_merge_bed(&ivs, false);
        assert_eq!(merged, vec![iv("chr1", 10, 100), iv("chr2", 5, 9)]);
    }

    #[test]
    fn sort_merge_fuses_bookended_intervals_bedtools_d0() {
        // PARITY PIN (corrected 2026-06-12 via the T7 differential against real
        // pybedtools): bookended half-open intervals [10,100) and [100,120) ARE
        // fused by `bedtools merge -d 0` / pybedtools `.merge()` into [10,120).
        // (An earlier revision wrongly kept them separate.)
        let ivs = vec![iv("chr1", 10, 100), iv("chr1", 100, 120)];
        let merged = sort_merge_bed(&ivs, false);
        assert_eq!(merged, vec![iv("chr1", 10, 120)]);
    }

    #[test]
    fn sort_merge_strictly_overlapping_intervals_fuse() {
        // The companion to the bookended pin: [10,100) and [50,120) genuinely
        // overlap (50 < 100) → one [10,120).
        let ivs = vec![iv("chr1", 10, 100), iv("chr1", 50, 120)];
        let merged = sort_merge_bed(&ivs, false);
        assert_eq!(merged, vec![iv("chr1", 10, 120)]);
    }

    #[test]
    fn sort_merge_nested_interval_is_absorbed() {
        // A fully-contained interval must not shrink the enclosing one.
        let ivs = vec![iv("chr1", 10, 100), iv("chr1", 30, 50)];
        let merged = sort_merge_bed(&ivs, false);
        assert_eq!(merged, vec![iv("chr1", 10, 100)]);
    }

    #[test]
    fn sort_merge_stranded_keeps_strands_apart() {
        let ivs = vec![
            sv("chr1", 10, 50, Strand::Forward),
            sv("chr1", 40, 80, Strand::Reverse), // overlaps positionally but opposite strand
            sv("chr1", 45, 90, Strand::Forward), // overlaps the +; but [10,50)+ and [45,90)+ overlap → merge
        ];
        let mut merged = sort_merge_bed(&ivs, true);
        merged.sort_by(|a, b| (a.start, &a.chrom).cmp(&(b.start, &b.chrom)));
        // + strand: [10,50) and [45,90) merge → [10,90); - strand: [40,80)
        assert!(merged.contains(&sv("chr1", 10, 90, Strand::Forward)));
        assert!(merged.contains(&sv("chr1", 40, 80, Strand::Reverse)));
        assert_eq!(merged.len(), 2);
    }

    #[test]
    fn sort_merge_empty_is_empty() {
        assert!(sort_merge_bed(&[], false).is_empty());
    }

    // ── intersect ────────────────────────────────────────────────────────────

    #[test]
    fn intersect_yields_overlap_spans() {
        let a = vec![iv("chr1", 10, 50)];
        let b = vec![iv("chr1", 30, 70)];
        let r = intersect(&a, &b);
        assert_eq!(r, vec![iv("chr1", 30, 50)]);
    }

    #[test]
    fn intersect_no_overlap_empty() {
        let a = vec![iv("chr1", 10, 20)];
        let b = vec![iv("chr1", 30, 40)];
        assert!(intersect(&a, &b).is_empty());
        // different chrom
        let c = vec![iv("chr2", 10, 50)];
        assert!(intersect(&a, &c).is_empty());
    }

    // ── slop ─────────────────────────────────────────────────────────────────

    #[test]
    fn slop_grows_and_clamps() {
        let ivs = vec![iv("chr1", 100, 200), iv("chr1", 5, 10)];
        let cs = sizes(&[("chr1", 250)]);
        let r = slop(&ivs, 20, &cs);
        assert_eq!(r[0], iv("chr1", 80, 220));
        assert_eq!(r[1], iv("chr1", 0, 30)); // start clamped to 0
    }

    #[test]
    fn slop_clamps_high_to_chrom_size() {
        let ivs = vec![iv("chr1", 200, 245)];
        let cs = sizes(&[("chr1", 250)]);
        let r = slop(&ivs, 20, &cs);
        assert_eq!(r[0], iv("chr1", 180, 250)); // end clamped to 250
    }

    // ── complement ───────────────────────────────────────────────────────────

    #[test]
    fn complement_emits_gaps_and_edges() {
        let ivs = vec![iv("chr1", 20, 40), iv("chr1", 60, 80)];
        let cs = sizes(&[("chr1", 100)]);
        let mut r = complement(&ivs, &cs);
        r.sort_by_key(|x| x.start);
        // [0,20), [40,60), [80,100)
        assert_eq!(
            r,
            vec![iv("chr1", 0, 20), iv("chr1", 40, 60), iv("chr1", 80, 100)]
        );
    }

    #[test]
    fn complement_fully_uncovered_contig() {
        let ivs: Vec<GenomicInterval> = vec![];
        let cs = sizes(&[("chrX", 50)]);
        let r = complement(&ivs, &cs);
        assert_eq!(r, vec![iv("chrX", 0, 50)]);
    }

    #[test]
    fn complement_fully_covered_contig_is_empty() {
        let ivs = vec![iv("chr1", 0, 100)];
        let cs = sizes(&[("chr1", 100)]);
        assert!(complement(&ivs, &cs).is_empty());
    }

    // ── merge_bed_files ──────────────────────────────────────────────────────

    #[test]
    fn merge_bed_files_dedups_paths_and_merges() {
        let f1 = tempfile::Builder::new().suffix(".bed").tempfile().unwrap();
        write_bed(f1.path(), &[iv("chr1", 10, 50)]).unwrap();
        let f2 = tempfile::Builder::new().suffix(".bed").tempfile().unwrap();
        write_bed(f2.path(), &[iv("chr1", 40, 80)]).unwrap();
        // pass f1 twice — the dedup should ignore the repeat
        let merged = merge_bed_files(&[f1.path(), f1.path(), f2.path()]).unwrap();
        assert_eq!(merged, vec![iv("chr1", 10, 80)]);
    }
}
