//! Header construction — ports `merge_vcf_headers` (merge L433-463) + the FILTER
//! tag injection both orchestrators do.
//!
//! `merge_vcf_headers(ref_header, query_header)` unions FILTER/INFO/FORMAT/contig
//! records (ref first). In this pipeline the query and reference VCFs both come
//! out of the same `sort_vcf`, so they share every INFO/FORMAT/contig definition;
//! the only lines that genuinely differ are **FILTER** tags (the per-RG source
//! tags). We therefore build from the reference template (all standard lines) and
//! append the query's FILTER definitions plus the new source/priority tags.
//! htslib's `bcf_hdr_append` (via [`bcf::Header::push_record`]) silently ignores a
//! duplicate ID, so appending an already-present FILTER is a no-op.

use rust_htslib::bcf::header::HeaderRecord;
use rust_htslib::bcf::{self, Read};

/// Append every FILTER definition from `src_reader`'s header into `header`,
/// reconstructing `##FILTER=<ID=...,Description="...">`. Duplicate IDs are ignored.
fn copy_filter_defs(header: &mut bcf::Header, src_reader: &bcf::Reader) {
    for rec in src_reader.header().header_records() {
        if let HeaderRecord::Filter { values, .. } = rec {
            let id = match values.get("ID") {
                Some(id) => id,
                None => continue,
            };
            if id == "PASS" {
                continue;
            }
            let desc = values
                .get("Description")
                .cloned()
                .unwrap_or_else(|| "\"\"".to_string());
            // `Description` comes back from htslib already quoted (e.g. `"text"`),
            // so emit it verbatim.
            let line = format!("##FILTER=<ID={id},Description={desc}>");
            header.push_record(line.as_bytes());
        }
    }
}

/// Append `##FILTER=<ID=tag,Description="…">` for each tag. Duplicate IDs ignored.
fn push_filter_defs(header: &mut bcf::Header, tags: &[&str]) {
    for tag in tags {
        let line =
            format!("##FILTER=<ID={tag},Description=\"Variant is likely to be {tag}\">");
        header.push_record(line.as_bytes());
    }
}

/// Build the merged output header for the priority-merge path: reference template
/// + query FILTER definitions + the `extra_filters` (source/priority/added tags).
pub fn build_merged_header(
    ref_reader: &bcf::Reader,
    query_reader: &bcf::Reader,
    extra_filters: &[&str],
) -> bcf::Header {
    let mut header = bcf::Header::from_template(ref_reader.header());
    copy_filter_defs(&mut header, query_reader);
    push_filter_defs(&mut header, extra_filters);
    header
}

/// Build the output header for the inhouse-common path: the query header (Python
/// writes with `header=bcf_query.header`, inhouse L420) + `extra_filters`.
pub fn build_query_header(query_reader: &bcf::Reader, extra_filters: &[&str]) -> bcf::Header {
    let mut header = bcf::Header::from_template(query_reader.header());
    push_filter_defs(&mut header, extra_filters);
    header
}
