//! The 7-column `all_homo_regions` BED reader + per-subgroup row splitter.
//!
//! This is the I/O-bookkeeping unit of the crate: it parses the headerless
//! 7-column BED produced by `preparation/build_beds_and_masked_genomes.py:174-189`
//! **once** per RG into owned [`AllRegionRow`]s, then [`split_subgroup`] hands out
//! borrowed views (the one FC row + the ≥1 NFC rows) for each subgroup. It is the
//! single tag-parsing / row-selection unit, replacing the three separate
//! `df.loc[...]` slices at `prepare_masked_align_region.py:163/180/181`.
//!
//! ## Column meaning (build_beds:178/180 FC, 189 NFC)
//!
//! | col | 0 chrom | 1 start | 2 end | 3 col4 | 4 col5 | 5 strand | 6 tag |
//! |-----|---------|---------|-------|--------|--------|----------|-------|
//! | FC  | chrom   | start   | end   | `"."`  | `"."`  | +/-      | `FC:{label}_{idx}`  |
//! | NFC | chrom   | start   | end   | int    | int    | +/-      | `NFC:{label}_{idx}` |
//!
//! For NFC rows col4/col5 are the NFC interval projected into the FC node's own
//! coordinate frame — already-written integers we only *read*
//! (`prepare_masked_align_region.py:84-85`). FC rows carry `"."` in col4/col5,
//! stored as the [`COL_SENTINEL`] (they are never read on the FC path).

use sdrecall_utils::{Result, SdError, Strand};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Sentinel for the `"."` col4/col5 of FC rows. FC rows never have their
/// col4/col5 read (only NFC rows feed `rel_start_interval`/`rel_end_interval`),
/// so the value is purely a "not a number" marker.
pub const COL_SENTINEL: i64 = i64::MIN;

/// A parsed `FC:`/`NFC:` tag. Parsed once (`prepare_masked_align_region.py:163`
/// regex `^FC.*` plus the exact `FC:{rg}_{sub}` / `NFC:{rg}_{sub}` matches on
/// lines 180-181 collapse into this one enum).
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum RgTag {
    /// Functional / target row (`FC:{label}_{sub}`).
    Fc { label: String, sub: String },
    /// Non-functional / counterpart row (`NFC:{label}_{sub}`).
    Nfc { label: String, sub: String },
}

impl RgTag {
    /// Parse a `FC:{label}_{sub}` / `NFC:{label}_{sub}` tag string. `label` is the
    /// text before the first `_` in the `{label}_{sub}` body; `sub` is the rest
    /// (so a multi-`_` sub is kept whole — matches Python's exact-string compare,
    /// which never re-splits the matched tag).
    fn parse(tag: &str) -> Option<RgTag> {
        let (kind, body) = tag.split_once(':')?;
        let (label, sub) = body.split_once('_')?;
        let label = label.to_string();
        let sub = sub.to_string();
        match kind {
            "FC" => Some(RgTag::Fc { label, sub }),
            "NFC" => Some(RgTag::Nfc { label, sub }),
            _ => None,
        }
    }
}

/// One parsed row of the 7-col all_homo_regions BED. Owns its `chrom` `String`
/// once (read once per RG); col4/col5 are pre-parsed `i64`
/// ([`COL_SENTINEL`] for the `"."` FC rows).
#[derive(Clone, Debug)]
pub struct AllRegionRow {
    /// Contig / chromosome name.
    pub chrom: String,
    /// 0-based start.
    pub start: i64,
    /// 0-based exclusive end.
    pub end: i64,
    /// FC rows: [`COL_SENTINEL`]; NFC rows: `rel_start_interval` (col4).
    pub col4: i64,
    /// FC rows: [`COL_SENTINEL`]; NFC rows: `rel_end_interval` (col5).
    pub col5: i64,
    /// Strand (col6).
    pub strand: Strand,
    /// Parsed tag (col7).
    pub tag: RgTag,
}

/// Parse the `"."`-or-integer col4/col5 field. `"."` → [`COL_SENTINEL`], matching
/// the FC-row convention; everything else must parse as `i64` (mirrors Python's
/// `int(interval[3])`, which would `ValueError` on a malformed NFC col).
fn parse_col(s: &str, line: usize, which: &str) -> Result<i64> {
    if s == "." {
        return Ok(COL_SENTINEL);
    }
    s.parse::<i64>().map_err(|_| SdError::BedParse {
        line,
        msg: format!("non-integer {which} column {s:?}"),
    })
}

/// Read the whole-region 7-col BED ONCE per RG.
///
/// Owned `Vec` because every subgroup filters it (single alloc, many borrows).
/// Strand is taken verbatim from col6 (`+`/`-`/`.`); rows whose tag is neither
/// `FC:`/`NFC:` are a hard parse error (the producer only writes those two).
pub fn read_all_region_bed(path: &Path) -> Result<Vec<AllRegionRow>> {
    let file = File::open(path).map_err(|e| SdError::Io {
        path: path.display().to_string(),
        source: e,
    })?;
    let reader = BufReader::new(file);
    let mut out = Vec::new();

    for (idx, line) in reader.lines().enumerate() {
        let lineno = idx + 1;
        let line = line.map_err(|e| SdError::Io {
            path: path.display().to_string(),
            source: e,
        })?;
        let trimmed = line.trim_end();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() < 7 {
            return Err(SdError::BedParse {
                line: lineno,
                msg: format!("expected 7 tab-separated columns, got {}", cols.len()),
            });
        }
        let start: i64 = cols[1].parse().map_err(|_| SdError::BedParse {
            line: lineno,
            msg: format!("non-integer start {:?}", cols[1]),
        })?;
        let end: i64 = cols[2].parse().map_err(|_| SdError::BedParse {
            line: lineno,
            msg: format!("non-integer end {:?}", cols[2]),
        })?;
        let col4 = parse_col(cols[3], lineno, "col4")?;
        let col5 = parse_col(cols[4], lineno, "col5")?;
        let strand = match cols[5] {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => Strand::Unknown,
        };
        let tag = RgTag::parse(cols[6]).ok_or_else(|| SdError::BedParse {
            line: lineno,
            msg: format!("tag {:?} is not FC:/NFC:", cols[6]),
        })?;
        out.push(AllRegionRow {
            chrom: cols[0].to_string(),
            start,
            end,
            col4,
            col5,
            strand,
            tag,
        });
    }
    Ok(out)
}

/// Borrow the rows; return the one FC row + the NFC rows for a `(label, sub)`
/// subgroup (no row data copied).
///
/// Mirrors the two Python asserts:
/// `len(fc_region_bedf) == 1` (`prepare_masked_align_region.py:213`) and
/// `len(nfc_region_bedf) >= 1` (line 214) — a violation is an [`SdError::Compute`]
/// (the typed analog of the Python `AssertionError`).
pub fn split_subgroup<'a>(
    rows: &'a [AllRegionRow],
    label: &str,
    sub: &str,
) -> Result<(&'a AllRegionRow, Vec<&'a AllRegionRow>)> {
    let mut fc: Vec<&AllRegionRow> = Vec::new();
    let mut nfc: Vec<&AllRegionRow> = Vec::new();
    for row in rows {
        match &row.tag {
            RgTag::Fc { label: l, sub: s } if l == label && s == sub => fc.push(row),
            RgTag::Nfc { label: l, sub: s } if l == label && s == sub => nfc.push(row),
            _ => {}
        }
    }
    if fc.len() != 1 {
        return Err(SdError::Compute(format!(
            "FC region for subgroup {label}_{sub} is not unique (found {})",
            fc.len()
        )));
    }
    if nfc.is_empty() {
        return Err(SdError::Compute(format!(
            "NFC region for subgroup {label}_{sub} does not exist"
        )));
    }
    Ok((fc[0], nfc))
}

/// All `FC:` rows across every subgroup of this RG — the input to the shared FC
/// target (`whole_region_bedf.loc[tag.str.contains("^FC.*")]`,
/// `prepare_masked_align_region.py:163`). Borrowed views, no copy.
pub fn fc_rows(rows: &[AllRegionRow]) -> Vec<&AllRegionRow> {
    rows.iter()
        .filter(|r| matches!(r.tag, RgTag::Fc { .. }))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(content: &str) -> tempfile::NamedTempFile {
        let f = tempfile::Builder::new().suffix(".bed").tempfile().unwrap();
        std::fs::write(f.path(), content).unwrap();
        f
    }

    #[test]
    fn parses_fc_and_nfc_rows() {
        let f = write(
            "chr2\t100\t200\t.\t.\t+\tFC:RG0_0\n\
             chr2\t500\t600\t0\t100\t-\tNFC:RG0_0\n",
        );
        let rows = read_all_region_bed(f.path()).unwrap();
        assert_eq!(rows.len(), 2);
        // FC row: col4/col5 are the sentinel.
        assert_eq!(rows[0].col4, COL_SENTINEL);
        assert_eq!(rows[0].col5, COL_SENTINEL);
        assert_eq!(rows[0].strand, Strand::Forward);
        assert_eq!(rows[0].tag, RgTag::Fc { label: "RG0".into(), sub: "0".into() });
        // NFC row: col4/col5 parsed.
        assert_eq!(rows[1].col4, 0);
        assert_eq!(rows[1].col5, 100);
        assert_eq!(rows[1].strand, Strand::Reverse);
        assert_eq!(rows[1].tag, RgTag::Nfc { label: "RG0".into(), sub: "0".into() });
    }

    #[test]
    fn split_subgroup_returns_fc_and_nfc() {
        let f = write(
            "chr2\t100\t200\t.\t.\t+\tFC:RG0_0\n\
             chr2\t500\t600\t0\t100\t-\tNFC:RG0_0\n\
             chr2\t700\t800\t0\t100\t+\tNFC:RG0_0\n\
             chr2\t900\t950\t.\t.\t+\tFC:RG0_1\n",
        );
        let rows = read_all_region_bed(f.path()).unwrap();
        let (fc, nfc) = split_subgroup(&rows, "RG0", "0").unwrap();
        assert_eq!(fc.start, 100);
        assert_eq!(nfc.len(), 2);
        // Different subgroup RG0_1 has an FC row but NO NFC row → error (the
        // len(nfc) >= 1 assert, Python:214).
        let err = split_subgroup(&rows, "RG0", "1").unwrap_err();
        assert!(matches!(err, SdError::Compute(_)), "got {err:?}");
    }

    #[test]
    fn split_subgroup_missing_nfc_is_error() {
        let f = write("chr2\t100\t200\t.\t.\t+\tFC:RG0_3\n");
        let rows = read_all_region_bed(f.path()).unwrap();
        let err = split_subgroup(&rows, "RG0", "3").unwrap_err();
        assert!(matches!(err, SdError::Compute(_)), "got {err:?}");
    }

    #[test]
    fn split_subgroup_duplicate_fc_is_error() {
        let f = write(
            "chr2\t100\t200\t.\t.\t+\tFC:RG0_0\n\
             chr2\t300\t400\t.\t.\t+\tFC:RG0_0\n\
             chr2\t500\t600\t0\t100\t+\tNFC:RG0_0\n",
        );
        let rows = read_all_region_bed(f.path()).unwrap();
        let err = split_subgroup(&rows, "RG0", "0").unwrap_err();
        assert!(matches!(err, SdError::Compute(_)), "got {err:?}");
    }

    #[test]
    fn fc_rows_filters_all_subgroups() {
        let f = write(
            "chr2\t100\t200\t.\t.\t+\tFC:RG0_0\n\
             chr2\t500\t600\t0\t100\t-\tNFC:RG0_0\n\
             chr2\t900\t950\t.\t.\t+\tFC:RG0_1\n",
        );
        let rows = read_all_region_bed(f.path()).unwrap();
        let fcs = fc_rows(&rows);
        assert_eq!(fcs.len(), 2);
        assert!(fcs.iter().all(|r| matches!(r.tag, RgTag::Fc { .. })));
    }

    #[test]
    fn too_few_columns_is_error() {
        let f = write("chr2\t100\t200\t.\t.\t+\n");
        let err = read_all_region_bed(f.path()).unwrap_err();
        assert!(matches!(err, SdError::BedParse { line: 1, .. }), "got {err:?}");
    }

    #[test]
    fn malformed_nfc_col_is_error() {
        let f = write("chr2\t100\t200\tNaN\t100\t+\tNFC:RG0_0\n");
        let err = read_all_region_bed(f.path()).unwrap_err();
        assert!(matches!(err, SdError::BedParse { line: 1, .. }), "got {err:?}");
    }
}
