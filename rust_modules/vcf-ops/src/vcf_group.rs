//! Read-all-then-group plumbing shared by both orchestrators.
//!
//! Each input is already coordinate-sorted by `sort_vcf`. We read every record,
//! translate it into the output writer's header (so the engine and finalizers
//! operate on records that are ready to write — the Rust replacement for the
//! Python `.pickable()` cross-process tuple), and bucket by contig `rid` for the
//! per-contig two-pointer pass.

use rust_htslib::bcf::{self, Read};
use rustc_hash::FxHashMap;
use sdrecall_utils::{Result, SdError};
use std::path::Path;

/// Read every record from `path`, translate it into `writer`'s header, and return
/// them in file order (which is coordinate order, since the input is sorted).
pub fn read_translate_group(path: &Path, writer: &mut bcf::Writer) -> Result<Vec<bcf::Record>> {
    let mut reader = bcf::Reader::from_path(path)
        .map_err(|e| SdError::Vcf(format!("open {}: {e}", path.display())))?;
    let mut out = Vec::new();
    for res in reader.records() {
        let mut rec = res.map_err(|e| SdError::Vcf(format!("read record: {e}")))?;
        writer.translate(&mut rec);
        out.push(rec);
    }
    Ok(out)
}

/// Bucket records by contig `rid`, preserving each bucket's input order (already
/// coordinate-sorted). `FxHashMap` for the cheap `u32` rid keys.
pub fn group_by_contig(recs: Vec<bcf::Record>) -> FxHashMap<u32, Vec<bcf::Record>> {
    let mut map: FxHashMap<u32, Vec<bcf::Record>> = FxHashMap::default();
    for rec in recs {
        if let Some(rid) = rec.rid() {
            map.entry(rid).or_default().push(rec);
        }
    }
    map
}
