//! TSV read/write — the ONE tsv unit.
//!
//! Serializes / deserializes tab-separated rows via `serde` on top of the `csv`
//! crate with the delimiter set to `b'\t'` (the DESIGN-sanctioned implementation).
//! Used for the realign meta table and similar typed dumps. A header row is
//! written/expected and maps to the struct field names.
//!
//! Both functions are generic over a serde row type `T`, so callers define a
//! `#[derive(Serialize, Deserialize)]` row struct once and reuse it for both
//! directions — no per-table parsing helper.

use sdrecall_utils::{Result, SdError};
use std::path::Path;

fn io_err(path: &Path, e: std::io::Error) -> SdError {
    SdError::Io {
        path: path.display().to_string(),
        source: e,
    }
}

/// Map a `csv::Error` to `SdError`, preserving the inner `io::Error` when the
/// failure was I/O (so the caller still sees the path) and stringifying parse /
/// deserialize failures.
fn csv_err(path: &Path, e: csv::Error) -> SdError {
    match e.into_kind() {
        csv::ErrorKind::Io(io) => io_err(path, io),
        other => SdError::Io {
            path: path.display().to_string(),
            source: std::io::Error::new(std::io::ErrorKind::InvalidData, format!("{other:?}")),
        },
    }
}

/// Read a TSV file into typed rows. The first line is the header and maps to the
/// struct's field names (serde `rename`/`alias` apply as usual). Owned `Vec` out.
pub fn read_tsv<T: serde::de::DeserializeOwned>(path: &Path) -> Result<Vec<T>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .map_err(|e| csv_err(path, e))?;
    let mut out = Vec::new();
    for rec in reader.deserialize() {
        out.push(rec.map_err(|e| csv_err(path, e))?);
    }
    Ok(out)
}

/// Write typed rows to a TSV file. The struct's field names become the header
/// row; one tab-separated line per row. Borrow-in slice, read-only emit.
pub fn write_tsv<T: serde::Serialize>(path: &Path, rows: &[T]) -> Result<()> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .map_err(|e| csv_err(path, e))?;
    for row in rows {
        writer.serialize(row).map_err(|e| csv_err(path, e))?;
    }
    writer.flush().map_err(|e| io_err(path, e))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde::{Deserialize, Serialize};

    #[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
    struct MetaRow {
        rg: String,
        chrom: String,
        start: i64,
        end: i64,
        score: f64,
        flag: bool,
    }

    #[test]
    fn round_trip_typed_rows() {
        let rows = vec![
            MetaRow {
                rg: "RG0".into(),
                chrom: "chr1".into(),
                start: 100,
                end: 200,
                score: 0.5,
                flag: true,
            },
            MetaRow {
                rg: "RG1".into(),
                chrom: "chr2".into(),
                start: 50,
                end: 75,
                score: -1.0,
                flag: false,
            },
        ];
        let tmp = tempfile::Builder::new().suffix(".tsv").tempfile().unwrap();
        write_tsv(tmp.path(), &rows).unwrap();
        let back: Vec<MetaRow> = read_tsv(tmp.path()).unwrap();
        assert_eq!(rows, back);
    }

    #[test]
    fn writes_header_and_tab_delimited() {
        let rows = vec![MetaRow {
            rg: "RG0".into(),
            chrom: "chr1".into(),
            start: 1,
            end: 2,
            score: 3.0,
            flag: true,
        }];
        let tmp = tempfile::Builder::new().suffix(".tsv").tempfile().unwrap();
        write_tsv(tmp.path(), &rows).unwrap();
        let content = std::fs::read_to_string(tmp.path()).unwrap();
        // header line is the field names, tab separated
        let mut lines = content.lines();
        assert_eq!(lines.next().unwrap(), "rg\tchrom\tstart\tend\tscore\tflag");
        assert_eq!(lines.next().unwrap(), "RG0\tchr1\t1\t2\t3.0\ttrue");
    }

    #[test]
    fn read_empty_with_only_header_is_empty() {
        let tmp = tempfile::Builder::new().suffix(".tsv").tempfile().unwrap();
        std::fs::write(tmp.path(), "rg\tchrom\tstart\tend\tscore\tflag\n").unwrap();
        let back: Vec<MetaRow> = read_tsv(tmp.path()).unwrap();
        assert!(back.is_empty());
    }

    #[test]
    fn read_missing_file_is_io_error() {
        let err = read_tsv::<MetaRow>(Path::new("/no/such/file.tsv")).unwrap_err();
        assert!(matches!(err, SdError::Io { .. }), "got {err:?}");
    }
}
