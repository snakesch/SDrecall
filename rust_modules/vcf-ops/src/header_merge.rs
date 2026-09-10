//! Header construction for the priority-merge and inhouse-common paths.
//!
//! Records from both callers are translated into the output writer's header.
//! Because BCF stores FILTER/INFO/FORMAT/contig references as numeric IDs, every
//! definition used by either caller must be present before translation. The
//! merged header keeps the reference header (including its sample layout) as the
//! template, then appends definitions that exist only in the query header.

use rust_htslib::bcf::header::{HeaderRecord, HeaderView};
use rust_htslib::bcf::{self, Read};
use std::collections::HashSet;

#[derive(Debug, Eq, Hash, PartialEq)]
enum DefinitionKey {
    Filter(String),
    Info(String),
    Format(String),
    Contig(String),
    Structured(String, String),
    Generic(String, String),
}

fn record_id<'a>(values: impl IntoIterator<Item = (&'a String, &'a String)>) -> String {
    let mut fields: Vec<_> = values
        .into_iter()
        .filter(|(key, _)| key.as_str() != "IDX")
        .collect();
    if let Some((_, id)) = fields.iter().find(|(key, _)| key.as_str() == "ID") {
        return (*id).clone();
    }
    fields.sort_unstable_by(|(left, _), (right, _)| left.cmp(right));
    fields
        .into_iter()
        .map(|(key, value)| format!("{key}={value}"))
        .collect::<Vec<_>>()
        .join(",")
}

fn definition_key(record: &HeaderRecord) -> DefinitionKey {
    match record {
        HeaderRecord::Filter { values, .. } => DefinitionKey::Filter(record_id(values.iter())),
        HeaderRecord::Info { values, .. } => DefinitionKey::Info(record_id(values.iter())),
        HeaderRecord::Format { values, .. } => DefinitionKey::Format(record_id(values.iter())),
        HeaderRecord::Contig { values, .. } => DefinitionKey::Contig(record_id(values.iter())),
        HeaderRecord::Structured { key, values } => {
            DefinitionKey::Structured(key.clone(), record_id(values.iter()))
        }
        HeaderRecord::Generic { key, value } => {
            // A VCF can contain several distinct `source`/command lines, but only
            // one fileformat declaration. Preserve unique generic metadata while
            // always retaining the reference fileformat.
            let value = if key == "fileformat" {
                String::new()
            } else {
                value.clone()
            };
            DefinitionKey::Generic(key.clone(), value)
        }
    }
}

fn structured_line<'a>(
    key: &str,
    values: impl IntoIterator<Item = (&'a String, &'a String)>,
) -> String {
    let values = values
        .into_iter()
        // `header_records()` exposes htslib's numeric dictionary slot. It is
        // header-local state, not portable VCF metadata; the merged header must
        // assign a fresh IDX for definitions copied from the query.
        .filter(|(field, _)| field.as_str() != "IDX")
        .map(|(field, value)| format!("{field}={value}"))
        .collect::<Vec<_>>()
        .join(",");
    format!("##{key}=<{values}>")
}

fn record_line(record: &HeaderRecord) -> String {
    match record {
        HeaderRecord::Filter { key, values }
        | HeaderRecord::Info { key, values }
        | HeaderRecord::Format { key, values }
        | HeaderRecord::Contig { key, values }
        | HeaderRecord::Structured { key, values } => structured_line(key, values.iter()),
        HeaderRecord::Generic { key, value } => format!("##{key}={value}"),
    }
}

fn definition_keys(header: &HeaderView) -> HashSet<DefinitionKey> {
    header.header_records().iter().map(definition_key).collect()
}

/// Append query definitions missing from the reference-derived output header.
fn copy_missing_defs(
    header: &mut bcf::Header,
    source: &HeaderView,
    definitions: &mut HashSet<DefinitionKey>,
) {
    for record in source.header_records() {
        if definitions.insert(definition_key(&record)) {
            header.push_record(record_line(&record).as_bytes());
        }
    }
}

/// Append `##FILTER=<ID=tag,Description="...">` for previously unseen tags.
fn push_filter_defs(
    header: &mut bcf::Header,
    tags: &[&str],
    definitions: &mut HashSet<DefinitionKey>,
) {
    for tag in tags {
        if definitions.insert(DefinitionKey::Filter((*tag).to_owned())) {
            let line = format!("##FILTER=<ID={tag},Description=\"Variant is likely to be {tag}\">");
            header.push_record(line.as_bytes());
        }
    }
}

/// Build the priority-merge output header: reference template, missing query
/// definitions, then the source/priority FILTER tags.
pub fn build_merged_header(
    ref_reader: &bcf::Reader,
    query_reader: &bcf::Reader,
    extra_filters: &[&str],
) -> bcf::Header {
    let mut header = bcf::Header::from_template(ref_reader.header());
    let mut definitions = definition_keys(ref_reader.header());
    copy_missing_defs(&mut header, query_reader.header(), &mut definitions);
    push_filter_defs(&mut header, extra_filters, &mut definitions);
    header
}

/// Build the inhouse-common output header: query template plus extra FILTERs.
pub fn build_query_header(query_reader: &bcf::Reader, extra_filters: &[&str]) -> bcf::Header {
    let mut header = bcf::Header::from_template(query_reader.header());
    let mut definitions = definition_keys(query_reader.header());
    push_filter_defs(&mut header, extra_filters, &mut definitions);
    header
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bcf::record::GenotypeAllele;
    use rust_htslib::bcf::{Format, Reader, Writer};

    fn reference_vcf() -> tempfile::TempPath {
        let temp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        let mut header = bcf::Header::new();
        header.push_record(b"##source=reference-caller");
        header.push_record(b"##contig=<ID=chrRef,length=101>");
        header.push_record(b"##contig=<ID=chrShared,length=202>");
        header.push_record(b"##FILTER=<ID=RFAIL,Description=\"reference only\">");
        header.push_record(b"##FILTER=<ID=SHARED_FILTER,Description=\"reference shared filter\">");
        header
            .push_record(b"##INFO=<ID=RINFO,Number=1,Type=Integer,Description=\"reference only\">");
        header.push_record(
            b"##INFO=<ID=SHARED_INFO,Number=1,Type=Integer,Description=\"reference shared info\">",
        );
        header.push_record(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        header.push_record(
            b"##FORMAT=<ID=RFMT,Number=1,Type=Integer,Description=\"reference only\">",
        );
        header.push_record(
            b"##FORMAT=<ID=SHARED_FMT,Number=1,Type=Integer,Description=\"reference shared format\">",
        );
        header.push_record(b"##ALT=<ID=SHARED_ALT,Description=\"reference shared allele\">");
        header.push_sample(b"sample");

        let mut writer = Writer::from_path(temp.path(), &header, true, Format::Vcf).unwrap();
        let mut record = writer.empty_record();
        record.set_rid(Some(writer.header().name2rid(b"chrRef").unwrap()));
        record.set_pos(10);
        record.set_alleles(&[b"A", b"T"]).unwrap();
        record.set_filters(&[b"RFAIL".as_slice()]).unwrap();
        record.push_info_integer(b"RINFO", &[31]).unwrap();
        record.push_info_integer(b"SHARED_INFO", &[29]).unwrap();
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
            .unwrap();
        record.push_format_integer(b"RFMT", &[37]).unwrap();
        record.push_format_integer(b"SHARED_FMT", &[41]).unwrap();
        writer.write(&record).unwrap();
        drop(writer);
        temp.into_temp_path()
    }

    fn query_vcf() -> tempfile::TempPath {
        let temp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        let mut header = bcf::Header::new();
        header.push_record(b"##source=query-caller");
        // Deliberately reverse shared/unique definitions relative to the
        // reference so their internal numeric IDs mean different tag names.
        header.push_record(b"##contig=<ID=chrShared,length=999>");
        header.push_record(b"##contig=<ID=chrQuery,length=303>");
        header.push_record(b"##FILTER=<ID=SHARED_FILTER,Description=\"query shared filter\">");
        header.push_record(b"##FILTER=<ID=QFAIL,Description=\"query only\">");
        header.push_record(
            b"##INFO=<ID=SHARED_INFO,Number=1,Type=Integer,Description=\"query shared info\">",
        );
        header.push_record(b"##INFO=<ID=QINFO,Number=1,Type=Integer,Description=\"query only\">");
        header.push_record(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        header.push_record(
            b"##FORMAT=<ID=SHARED_FMT,Number=1,Type=Integer,Description=\"query shared format\">",
        );
        header.push_record(b"##FORMAT=<ID=QFMT,Number=1,Type=Integer,Description=\"query only\">");
        header.push_record(b"##ALT=<ID=QALT,Description=\"query symbolic allele\">");
        header.push_record(b"##ALT=<ID=SHARED_ALT,Description=\"query shared allele\">");
        header.push_sample(b"sample");

        let mut writer = Writer::from_path(temp.path(), &header, true, Format::Vcf).unwrap();
        let mut record = writer.empty_record();
        record.set_rid(Some(writer.header().name2rid(b"chrQuery").unwrap()));
        record.set_pos(20);
        record.set_alleles(&[b"C", b"<QALT>"]).unwrap();
        record
            .set_filters(&[b"SHARED_FILTER".as_slice(), b"QFAIL".as_slice()])
            .unwrap();
        record.push_info_integer(b"SHARED_INFO", &[11]).unwrap();
        record.push_info_integer(b"QINFO", &[17]).unwrap();
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)])
            .unwrap();
        record.push_format_integer(b"SHARED_FMT", &[13]).unwrap();
        record.push_format_integer(b"QFMT", &[23]).unwrap();
        writer.write(&record).unwrap();
        drop(writer);
        temp.into_temp_path()
    }

    fn definition_count(records: &[HeaderRecord], kind: &str, id: &str) -> usize {
        records
            .iter()
            .filter(|record| match record {
                HeaderRecord::Filter { values, .. } => {
                    kind == "FILTER" && values.get("ID").is_some_and(|value| value == id)
                }
                HeaderRecord::Info { values, .. } => {
                    kind == "INFO" && values.get("ID").is_some_and(|value| value == id)
                }
                HeaderRecord::Format { values, .. } => {
                    kind == "FORMAT" && values.get("ID").is_some_and(|value| value == id)
                }
                HeaderRecord::Contig { values, .. } => {
                    kind == "contig" && values.get("ID").is_some_and(|value| value == id)
                }
                HeaderRecord::Structured { key, values } => {
                    key == kind && values.get("ID").is_some_and(|value| value == id)
                }
                HeaderRecord::Generic { .. } => false,
            })
            .count()
    }

    #[test]
    fn heterogeneous_headers_translate_unique_tags_without_id_corruption() {
        let reference = reference_vcf();
        let query = query_vcf();
        let output = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();

        let mut ref_reader = Reader::from_path(&reference).unwrap();
        let mut query_reader = Reader::from_path(&query).unwrap();
        let header = build_merged_header(&ref_reader, &query_reader, &["EXTRA_FILTER", "QFAIL"]);
        let mut writer = Writer::from_path(output.path(), &header, true, Format::Vcf).unwrap();

        let mut query_record = query_reader.records().next().unwrap().unwrap();
        writer.translate(&mut query_record);
        writer.write(&query_record).unwrap();
        let mut ref_record = ref_reader.records().next().unwrap().unwrap();
        writer.translate(&mut ref_record);
        writer.write(&ref_record).unwrap();
        drop(writer);

        let mut merged = Reader::from_path(output.path()).unwrap();
        assert_eq!(merged.header().samples(), [b"sample".as_slice()]);
        let header_records = merged.header().header_records();
        for (kind, id) in [
            ("FILTER", "SHARED_FILTER"),
            ("FILTER", "QFAIL"),
            ("FILTER", "EXTRA_FILTER"),
            ("INFO", "SHARED_INFO"),
            ("INFO", "QINFO"),
            ("FORMAT", "SHARED_FMT"),
            ("FORMAT", "QFMT"),
            ("contig", "chrShared"),
            ("contig", "chrQuery"),
            ("ALT", "SHARED_ALT"),
            ("ALT", "QALT"),
        ] {
            assert_eq!(
                definition_count(&header_records, kind, id),
                1,
                "{kind}/{id}"
            );
        }

        let shared_info = header_records
            .iter()
            .find_map(|record| match record {
                HeaderRecord::Info { values, .. }
                    if values.get("ID").is_some_and(|value| value == "SHARED_INFO") =>
                {
                    values.get("Description")
                }
                _ => None,
            })
            .unwrap();
        assert_eq!(shared_info, "\"reference shared info\"");

        let sources: HashSet<_> = header_records
            .iter()
            .filter_map(|record| match record {
                HeaderRecord::Generic { key, value } if key == "source" => Some(value.as_str()),
                _ => None,
            })
            .collect();
        assert_eq!(sources, HashSet::from(["reference-caller", "query-caller"]));

        let mut saw_query = false;
        let mut saw_reference = false;
        for result in merged.records() {
            let record = result.unwrap();
            let rid = record.rid().unwrap();
            let contig = record.header().rid2name(rid).unwrap();
            match contig {
                b"chrQuery" => {
                    saw_query = true;
                    assert_eq!(record.alleles(), [b"C".as_slice(), b"<QALT>".as_slice()]);
                    assert!(record.has_filter(b"SHARED_FILTER".as_slice()));
                    assert!(record.has_filter(b"QFAIL".as_slice()));
                    assert_eq!(
                        record.info(b"SHARED_INFO").integer().unwrap().unwrap()[0],
                        11
                    );
                    assert_eq!(record.info(b"QINFO").integer().unwrap().unwrap()[0], 17);
                    assert_eq!(record.format(b"SHARED_FMT").integer().unwrap()[0][0], 13);
                    assert_eq!(record.format(b"QFMT").integer().unwrap()[0][0], 23);
                }
                b"chrRef" => {
                    saw_reference = true;
                    assert!(record.has_filter(b"RFAIL".as_slice()));
                    assert_eq!(record.info(b"RINFO").integer().unwrap().unwrap()[0], 31);
                    assert_eq!(
                        record.info(b"SHARED_INFO").integer().unwrap().unwrap()[0],
                        29
                    );
                    assert_eq!(record.format(b"RFMT").integer().unwrap()[0][0], 37);
                    assert_eq!(record.format(b"SHARED_FMT").integer().unwrap()[0][0], 41);
                }
                other => panic!("unexpected contig {}", String::from_utf8_lossy(other)),
            }
        }
        assert!(saw_query && saw_reference);
    }
}
