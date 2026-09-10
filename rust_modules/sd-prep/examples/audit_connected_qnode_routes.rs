//! Audit why specific FC pairs do or do not become `ConnectedQnodes` edges.
//!
//! This consumes the exact production `filtered_SD_binary_map.tsv`, applies the
//! production pruning rules, and evaluates every candidate pair in both
//! directions through the production Dijkstra, route projection, and minimap2
//! similarity code.

use sd_prep::graph_build::{build_multiplex_graph, NodeKey, SdPairRow};
use sd_prep::minimap::{align_similarity, Preset};
use sd_prep::traversal::{
    audit_qnode_route_in_pruned_graph, prune_graph, FragParams, QnodeRouteAudit,
};
use sdrecall_utils::Strand;
use std::collections::HashMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

fn parse_strand(value: &str) -> Strand {
    match value {
        "+" => Strand::Forward,
        "-" => Strand::Reverse,
        _ => Strand::Unknown,
    }
}

fn strand_text(strand: Strand) -> &'static str {
    match strand {
        Strand::Forward => "+",
        Strand::Reverse => "-",
        Strand::Unknown => ".",
    }
}

fn node_text(node: &NodeKey) -> String {
    format!(
        "{}:{}-{}({})",
        node.chrom,
        node.start,
        node.end,
        strand_text(node.strand)
    )
}

fn parse_display_interval(value: &str) -> NodeKey {
    let (without_strand, strand_part) = value
        .rsplit_once('(')
        .unwrap_or_else(|| panic!("interval lacks strand suffix: {value}"));
    let strand = strand_part
        .strip_suffix(')')
        .unwrap_or_else(|| panic!("invalid strand suffix: {value}"));
    let (chrom, coords) = without_strand
        .split_once(':')
        .unwrap_or_else(|| panic!("interval lacks chromosome separator: {value}"));
    let (start, end) = coords
        .split_once('-')
        .unwrap_or_else(|| panic!("interval lacks coordinate separator: {value}"));
    NodeKey::new(
        chrom,
        start.parse().expect("interval start"),
        end.parse().expect("interval end"),
        parse_strand(strand),
    )
}

fn read_sd_map(path: &Path) -> Vec<SdPairRow> {
    let text = std::fs::read_to_string(path).expect("read filtered SD map");
    let mut rows = Vec::new();
    for (line_index, line) in text.lines().enumerate() {
        if line.trim().is_empty() || (line_index == 0 && line.starts_with("chr_1")) {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        assert!(
            fields.len() >= 9,
            "filtered SD map row {} has only {} columns",
            line_index + 1,
            fields.len()
        );
        rows.push(SdPairRow {
            a: NodeKey::new(
                fields[0],
                fields[1].parse().expect("a start"),
                fields[2].parse().expect("a end"),
                parse_strand(fields[3]),
            ),
            b: NodeKey::new(
                fields[4],
                fields[5].parse().expect("b start"),
                fields[6].parse().expect("b end"),
                parse_strand(fields[7]),
            ),
            mismatch_rate: fields[8].parse().expect("mismatch rate"),
        });
    }
    rows
}

fn header_indexes(header: &str) -> HashMap<&str, usize> {
    header
        .split('\t')
        .enumerate()
        .map(|(index, name)| (name, index))
        .collect()
}

fn field<'a>(fields: &'a [&str], indexes: &HashMap<&str, usize>, name: &str) -> &'a str {
    fields[*indexes
        .get(name)
        .unwrap_or_else(|| panic!("candidate table lacks column {name}"))]
}

fn route_nodes(audit: &QnodeRouteAudit) -> String {
    if audit.route_edges.is_empty() {
        return String::new();
    }
    let mut value = node_text(&audit.route_edges[0].from);
    for edge in &audit.route_edges {
        let _ = write!(value, ">{}", node_text(&edge.to));
    }
    value
}

fn route_edge_details(audit: &QnodeRouteAudit) -> String {
    audit
        .route_edges
        .iter()
        .map(|edge| {
            let kind = match (edge.is_sd, edge.is_overlap) {
                (true, true) => "SD+PO",
                (true, false) => "SD",
                (false, true) => "PO",
                (false, false) => "NONE",
            };
            format!("{kind}:{:.8}", edge.weight)
        })
        .collect::<Vec<_>>()
        .join(";")
}

fn option_bool(value: Option<bool>) -> &'static str {
    match value {
        Some(true) => "yes",
        Some(false) => "no",
        None => "not_evaluable",
    }
}

fn option_float(value: Option<f64>) -> String {
    value.map_or_else(String::new, |number| format!("{number:.8}"))
}

fn candidate_window(audit: &QnodeRouteAudit) -> String {
    audit
        .candidate_window
        .as_ref()
        .map_or_else(String::new, |window| {
            format!(
                "{}:{}-{}({})",
                window.0,
                window.1,
                window.2,
                strand_text(window.3)
            )
        })
}

fn fetch_node_sequence(
    reader: &mut bio::io::fasta::IndexedReader<std::fs::File>,
    node: &NodeKey,
    reverse_complement: bool,
) -> Vec<u8> {
    reader
        .fetch(&node.chrom, node.start as u64, node.end as u64)
        .expect("fetch direct FC sequence");
    let mut sequence = Vec::new();
    reader.read(&mut sequence).expect("read direct FC sequence");
    if reverse_complement {
        bio::alphabets::dna::revcomp(&sequence)
    } else {
        sequence
    }
}

fn direct_full_fc_similarity(source: &NodeKey, target: &NodeKey, reference: &Path) -> f64 {
    let mut reader =
        bio::io::fasta::IndexedReader::from_file(&reference).expect("open reference FASTA");
    let source_sequence = fetch_node_sequence(&mut reader, source, false);
    let target_sequence = fetch_node_sequence(&mut reader, target, target.strand != source.strand);
    align_similarity(&target_sequence, &source_sequence, Preset::Asm10)
        .expect("direct full-FC minimap2 similarity")
}

fn genomic_gap_or_overlap(source: &NodeKey, target: &NodeKey) -> String {
    if source.chrom != target.chrom {
        return String::new();
    }
    let overlap = source.end.min(target.end) - source.start.max(target.start);
    if overlap > 0 {
        format!("overlap:{overlap}")
    } else {
        let gap = if source.end <= target.start {
            target.start - source.end
        } else {
            source.start - target.end
        };
        format!("gap:{gap}")
    }
}

struct RunConfig<'a> {
    assembly: &'a str,
    reference: &'a Path,
}

struct PairMeta<'a> {
    allele: &'a str,
    rg: &'a str,
    source_label: &'a str,
    target_label: &'a str,
    direction: &'a str,
    direct_similarity: f64,
    remote_nfc_label: &'a str,
    remote_nfc_interval: &'a str,
    remote_nfc_to_destination_similarity: Option<f64>,
}

fn write_audit_row(output: &mut String, meta: &PairMeta<'_>, audit: &QnodeRouteAudit) {
    let columns = [
        meta.allele.to_string(),
        meta.rg.to_string(),
        meta.direction.to_string(),
        meta.source_label.to_string(),
        node_text(&audit.source),
        meta.target_label.to_string(),
        node_text(&audit.target),
        genomic_gap_or_overlap(&audit.source, &audit.target),
        format!("{:.8}", meta.direct_similarity),
        meta.remote_nfc_label.to_string(),
        meta.remote_nfc_interval.to_string(),
        option_float(meta.remote_nfc_to_destination_similarity),
        if audit.source_present { "yes" } else { "no" }.to_string(),
        if audit.target_present { "yes" } else { "no" }.to_string(),
        option_bool(audit.same_component).to_string(),
        option_float(audit.route_cost),
        audit.route_edges.len().to_string(),
        route_nodes(audit),
        route_edge_details(audit),
        audit
            .last_edge_is_overlap
            .map_or_else(String::new, |value| {
                if value { "yes" } else { "no" }.to_string()
            }),
        audit.adjacent_overlap_pairs.to_string(),
        option_float(audit.sd_product),
        candidate_window(audit),
        option_float(audit.similarity),
        if audit.target_is_query_node {
            "yes"
        } else {
            "no"
        }
        .to_string(),
        if audit.creates_connected_qnode_edge {
            "yes"
        } else {
            "no"
        }
        .to_string(),
        audit.decision.to_string(),
    ];
    let _ = writeln!(output, "{}", columns.join("\t"));
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 10 {
        eprintln!(
            "usage: {} <assembly> <filtered_SD_binary_map.tsv> <reference.fa> \
             <avg_frag> <std_frag> <mean_read_length> <relationship_candidates.tsv> \
             <output.tsv> <threads>",
            args[0]
        );
        std::process::exit(2);
    }
    let assembly = &args[1];
    let sd_map = PathBuf::from(&args[2]);
    let reference = PathBuf::from(&args[3]);
    let avg_frag: f64 = args[4].parse().expect("avg_frag");
    let std_frag: f64 = args[5].parse().expect("std_frag");
    let mean_read_length: f64 = args[6].parse().expect("mean_read_length");
    let candidates_path = PathBuf::from(&args[7]);
    let output_path = PathBuf::from(&args[8]);
    let threads: usize = args[9].parse().expect("threads");
    let config = RunConfig {
        assembly,
        reference: &reference,
    };

    let rows = read_sd_map(&sd_map);
    let graph = build_multiplex_graph(&rows, threads);
    let pruned = prune_graph(
        &graph,
        &FragParams {
            avg_frag,
            std_frag,
            mean_read_length,
        },
    );
    let frag = FragParams {
        avg_frag,
        std_frag,
        mean_read_length,
    };

    let candidate_text = std::fs::read_to_string(&candidates_path).expect("read candidates");
    let mut lines = candidate_text.lines();
    let header = lines.next().expect("candidate header");
    let indexes = header_indexes(header);
    let mut output = String::from(
        "FN allele\tRG\tDirection\tSource FC\tSource interval\tTarget FC\tTarget interval\
         \tGenomic gap or overlap bp\tDirect full-FC minimap2 similarity (diagnostic only)\
         \tExact remote NFC\tExact remote NFC interval\
         \tRemote NFC to destination FC minimap2 similarity\
         \tSource present after pruning\tTarget present after pruning\tSame multiplex component\
         \tDijkstra route cost\tRoute edge count\tRoute nodes\tRoute edge types and weights\
         \tLast edge is physical overlap\tAdjacent physical-overlap edge pairs\tSD route product\
         \tProjected target window\tMinimap2 similarity\tTarget is query node\
         \tCreates ConnectedQnodes edge\tExact production decision\n",
    );

    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if field(&fields, &indexes, "Assembly") != config.assembly {
            continue;
        }
        let allele = field(&fields, &indexes, "FN allele");
        let rg = field(&fields, &indexes, "RG");
        let remote_schema = indexes.contains_key("Exact remote NFC interval");
        let realigned_fc_pair_schema = indexes.contains_key("FC A interval");
        let (fn_label, observed_label, fn_node, observed_node, remote_label, remote_interval) =
            if remote_schema {
                (
                    field(&fields, &indexes, "Owner FC"),
                    field(&fields, &indexes, "Destination FC"),
                    parse_display_interval(field(&fields, &indexes, "Owner FC interval")),
                    parse_display_interval(field(&fields, &indexes, "Destination FC interval")),
                    field(&fields, &indexes, "Exact remote NFC"),
                    field(&fields, &indexes, "Exact remote NFC interval"),
                )
            } else if realigned_fc_pair_schema {
                (
                    field(&fields, &indexes, "FC A"),
                    field(&fields, &indexes, "FC B"),
                    parse_display_interval(field(&fields, &indexes, "FC A interval")),
                    parse_display_interval(field(&fields, &indexes, "FC B interval")),
                    "",
                    "",
                )
            } else {
                (
                    field(&fields, &indexes, "FN-site FC"),
                    field(&fields, &indexes, "Observed remote FC"),
                    parse_display_interval(field(&fields, &indexes, "FN-site FC interval")),
                    parse_display_interval(field(&fields, &indexes, "Observed remote FC interval")),
                    "",
                    "",
                )
            };
        let remote_similarity = if remote_schema {
            let remote_node = parse_display_interval(remote_interval);
            Some(direct_full_fc_similarity(
                &remote_node,
                &observed_node,
                config.reference,
            ))
        } else {
            None
        };
        let forward_direction = if realigned_fc_pair_schema {
            "FC_A_to_FC_B"
        } else {
            "FN-site_to_observed"
        };
        let reverse_direction = if realigned_fc_pair_schema {
            "FC_B_to_FC_A"
        } else {
            "observed_to_FN-site"
        };

        let forward = audit_qnode_route_in_pruned_graph(
            &fn_node,
            &observed_node,
            true,
            &pruned,
            config.reference,
            &frag,
        )
        .expect("forward route audit");
        let forward_direct_similarity =
            direct_full_fc_similarity(&fn_node, &observed_node, config.reference);
        write_audit_row(
            &mut output,
            &PairMeta {
                allele,
                rg,
                source_label: fn_label,
                target_label: observed_label,
                direction: forward_direction,
                direct_similarity: forward_direct_similarity,
                remote_nfc_label: remote_label,
                remote_nfc_interval: remote_interval,
                remote_nfc_to_destination_similarity: remote_similarity,
            },
            &forward,
        );

        let reverse = audit_qnode_route_in_pruned_graph(
            &observed_node,
            &fn_node,
            true,
            &pruned,
            config.reference,
            &frag,
        )
        .expect("reverse route audit");
        let reverse_direct_similarity =
            direct_full_fc_similarity(&observed_node, &fn_node, config.reference);
        write_audit_row(
            &mut output,
            &PairMeta {
                allele,
                rg,
                source_label: observed_label,
                target_label: fn_label,
                direction: reverse_direction,
                direct_similarity: reverse_direct_similarity,
                remote_nfc_label: remote_label,
                remote_nfc_interval: remote_interval,
                remote_nfc_to_destination_similarity: remote_similarity,
            },
            &reverse,
        );
    }

    std::fs::write(&output_path, output).expect("write route audit");
    eprintln!(
        "assembly={} sd_rows={} graph_nodes={} graph_edges={} pruned_nodes={} \
         pruned_edges={} output={}",
        assembly,
        rows.len(),
        graph.node_count(),
        graph.edge_count(),
        pruned.node_count(),
        pruned.edge_count(),
        output_path.display()
    );
}
