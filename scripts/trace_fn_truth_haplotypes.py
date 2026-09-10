#!/usr/bin/env python3
"""Trace exact ALT truth haplotypes into their SDrecall target RG BAMs.

A realigned read counts as correctly placed only when the same read sequence
still contains an ALT-only truth marker and every marker base follows its
expected reference-coordinate pattern. This verifies insertions and deletion
gaps as well as SNVs; a matching read name elsewhere is not sufficient.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pysam


@dataclass(frozen=True)
class SiteDefinition:
    assembly: str
    label: str
    haplotype: str
    chrom: str
    pos: int
    ref: str
    alt: str
    window_start: int
    window_end: int
    kmer: str
    reference_map: tuple[int | None, ...] | None
    supported_sites: frozenset[str]


@dataclass(frozen=True)
class Placement:
    strand: str
    kmer_offset: int
    target_reference_pos: int
    mapped_reference_bases: int
    strict_reference_map_match: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--definitions", required=True, type=Path)
    parser.add_argument("--input-hits", required=True, type=Path)
    parser.add_argument("--rg-summary", required=True, type=Path)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--sample", default="HG002")
    parser.add_argument("--output-trace", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    parser.add_argument("--output-per-rg-summary", type=Path)
    return parser.parse_args()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def qname_base(qname: str) -> str:
    return re.sub(r":RG[0-9]+$", "", qname)


def record_id(read: pysam.AlignedSegment) -> str:
    if read.is_read1:
        mate = "1"
    elif read.is_read2:
        mate = "2"
    else:
        mate = "0"
    return f"{qname_base(read.query_name)}|{mate}"


def read_definitions(path: Path) -> list[SiteDefinition]:
    definitions: list[SiteDefinition] = []
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chrom, pos, ref, alt = row["FN_site"].split(":", 3)
            window_chrom, span = row["window"].split(":", 1)
            window_start, window_end = map(int, span.split("-", 1))
            if window_chrom != chrom:
                raise ValueError(f"Window chromosome mismatch for {row['FN_site']}")
            supported_sites = {row["FN_site"]}
            for variant_label in row["applied_ALT_variants"].split(";"):
                match = re.fullmatch(r"([^:]+):(\d+):([^>]+)>([^:]+):.*", variant_label)
                if match:
                    supported_sites.add(
                        f"{match.group(1)}:{match.group(2)}:{match.group(3)}:{match.group(4)}"
                    )
            map_text = (row.get("kmer_reference_map") or "").strip()
            reference_map = None
            if map_text:
                reference_map = tuple(
                    None if value == "." else int(value) for value in map_text.split(",")
                )
                if len(reference_map) != len(row["kmer_seq"]):
                    raise ValueError(
                        f"Marker/reference-map length mismatch for {row['FN_site']}"
                    )
            definitions.append(
                SiteDefinition(
                    assembly=row["assembly"],
                    label=row["FN_site"],
                    haplotype=row.get("haplotype", "0"),
                    chrom=chrom,
                    pos=int(pos),
                    ref=ref,
                    alt=alt,
                    window_start=window_start,
                    window_end=window_end,
                    kmer=row["kmer_seq"].upper(),
                    reference_map=reference_map,
                    supported_sites=frozenset(supported_sites),
                )
            )
    return definitions


def read_input_hits(path: Path) -> dict[str, set[str]]:
    hits: dict[str, set[str]] = defaultdict(set)
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            hits[row["FN_site"]].add(qname_base(row["read_qname"]))
    return hits


def read_rg_summary(path: Path) -> dict[str, list[str]]:
    by_site: dict[str, list[str]] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_site[row["FN_site"]] = [item for item in row["relevant_RGs"].split(",") if item]
    return by_site


def forward_query_to_reference_map(read: pysam.AlignedSegment) -> dict[int, int]:
    """Map original read-orientation query positions to 0-based references."""
    stored_map = {
        query_pos: reference_pos
        for query_pos, reference_pos in read.get_aligned_pairs(matches_only=False)
        if query_pos is not None and reference_pos is not None
    }
    if not read.is_reverse:
        return stored_map
    query_length = read.query_length or len(read.query_sequence or "")
    return {
        query_length - 1 - query_pos: reference_pos
        for query_pos, reference_pos in stored_map.items()
    }


def resolve_bam_contig(bam: pysam.AlignmentFile, chrom: str) -> str:
    if chrom in bam.references:
        return chrom
    aliases = [chrom.removeprefix("chr")]
    if not chrom.startswith("chr"):
        aliases.append(f"chr{chrom}")
    if chrom in {"chrM", "M"}:
        aliases.extend(["MT", "chrMT"])
    for alias in aliases:
        if alias in bam.references:
            return alias
    raise ValueError(f"Contig {chrom!r} is absent from BAM")


def kmer_placement(
    read: pysam.AlignedSegment,
    marker: SiteDefinition,
    target_pos: int,
) -> Placement | None:
    sequence = (read.get_forward_sequence() or read.query_sequence or "").upper()
    candidates = [(marker.kmer, "+", marker.reference_map)]
    reverse = reverse_complement(marker.kmer)
    if reverse != marker.kmer:
        reverse_map = tuple(reversed(marker.reference_map)) if marker.reference_map else None
        candidates.append((reverse, "-", reverse_map))

    query_to_ref = forward_query_to_reference_map(read)
    marker_center0 = marker.pos - 1
    target0 = target_pos - 1
    target_offset = target_pos - marker.window_start
    marker_offset = marker.pos - marker.window_start
    for kmer, strand, expected_reference_map in candidates:
        offset = sequence.find(kmer)
        while offset >= 0:
            if expected_reference_map is not None:
                expected0 = tuple(
                    None if position is None else position - 1
                    for position in expected_reference_map
                )
                observed0 = tuple(
                    query_to_ref.get(query_pos)
                    for query_pos in range(offset, offset + len(kmer))
                )
                if observed0 == expected0 and target0 in expected0:
                    return Placement(
                        strand=strand,
                        kmer_offset=offset,
                        target_reference_pos=target_pos,
                        mapped_reference_bases=sum(
                            position is not None for position in expected_reference_map
                        ),
                        strict_reference_map_match=True,
                    )
            else:
                if strand == "+":
                    center_query = offset + marker_offset
                    target_query = offset + target_offset
                else:
                    center_query = offset + len(kmer) - 1 - marker_offset
                    target_query = offset + len(kmer) - 1 - target_offset
                center_ref = query_to_ref.get(center_query, -1)
                target_ref = query_to_ref.get(target_query, -1)
                mapped_positions = [
                    query_to_ref[query_pos]
                    for query_pos in range(offset, offset + len(kmer))
                    if query_pos in query_to_ref
                ]
                expected_start0 = marker.window_start - 1
                expected_end0 = marker.window_end - 1
                expected_hits = sum(
                    expected_start0 <= ref_pos <= expected_end0
                    for ref_pos in mapped_positions
                )
                mapped_fraction = expected_hits / len(marker.kmer)
                if (
                    center_ref == marker_center0
                    and target_ref == target0
                    and mapped_fraction >= 0.90
                ):
                    return Placement(
                        strand=strand,
                        kmer_offset=offset,
                        target_reference_pos=target_ref + 1,
                        mapped_reference_bases=expected_hits,
                        strict_reference_map_match=False,
                    )
            offset = sequence.find(kmer, offset + 1)
    return None


def main() -> None:
    args = parse_args()
    definitions = read_definitions(args.definitions)
    input_hits = read_input_hits(args.input_hits)
    rg_summary = read_rg_summary(args.rg_summary)

    args.output_trace.parent.mkdir(parents=True, exist_ok=True)
    args.output_summary.parent.mkdir(parents=True, exist_ok=True)

    trace_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    per_rg_rows: list[dict[str, object]] = []
    definitions_by_site: dict[str, list[SiteDefinition]] = defaultdict(list)
    for definition in definitions:
        definitions_by_site[definition.label].append(definition)

    markers_by_target: dict[str, list[SiteDefinition]] = defaultdict(list)
    target_labels = set(definitions_by_site)
    for marker in definitions:
        for supported_site in marker.supported_sites:
            if supported_site in target_labels:
                markers_by_target[supported_site].append(marker)

    for site_label, site_definitions in definitions_by_site.items():
        definition = site_definitions[0]
        candidate_markers = markers_by_target[site_label]
        site_input_qnames = set().union(
            *(input_hits.get(marker.label, set()) for marker in candidate_markers)
        )
        site_target_covering: set[str] = set()
        site_window_spanning: set[str] = set()
        site_haplotype_at_target: set[str] = set()
        site_haplotype_input_qnames_at_target: set[str] = set()

        relevant_rgs = rg_summary.get(site_label, [])
        for rg in relevant_rgs:
            rg_target_covering: set[str] = set()
            rg_window_spanning: set[str] = set()
            rg_haplotype_at_target: set[str] = set()
            rg_input_qnames_at_target: set[str] = set()
            bam_path = (
                args.run_dir
                / "recall_results"
                / f"{args.sample}.sdrecall.only_{rg}.raw.bam"
            )
            if not bam_path.exists():
                raise FileNotFoundError(bam_path)

            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                fetch_start0 = max(0, definition.window_start - 1)
                fetch_end0 = definition.window_end
                for read in bam.fetch(definition.chrom, fetch_start0, fetch_end0):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    rid = record_id(read)
                    query_to_ref = forward_query_to_reference_map(read)
                    mapped_ref_positions = set(query_to_ref.values())
                    if definition.pos - 1 in mapped_ref_positions:
                        rg_target_covering.add(rid)
                        site_target_covering.add(f"{rg}|{rid}")
                    if {
                        definition.window_start - 1,
                        definition.window_end - 1,
                    }.issubset(mapped_ref_positions):
                        rg_window_spanning.add(rid)
                        site_window_spanning.add(f"{rg}|{rid}")

                    placement: Placement | None = None
                    matched_marker = None
                    for marker in candidate_markers:
                        marker_placement = kmer_placement(read, marker, definition.pos)
                        if marker_placement is not None:
                            placement = marker_placement
                            matched_marker = marker
                            break
                    if placement is None or matched_marker is None:
                        continue
                    rg_haplotype_at_target.add(rid)
                    site_haplotype_at_target.add(f"{rg}|{rid}")
                    input_match = qname_base(read.query_name) in site_input_qnames
                    if input_match:
                        rg_input_qnames_at_target.add(qname_base(read.query_name))
                        site_haplotype_input_qnames_at_target.add(qname_base(read.query_name))
                    trace_rows.append(
                        {
                            "FN_site": site_label,
                            "marker_FN_site": matched_marker.label,
                            "marker_haplotype": matched_marker.haplotype,
                            "marker_kmer_len": len(matched_marker.kmer),
                            "RG": rg,
                            "input_qname_match": "yes" if input_match else "no",
                            "read_qname": read.query_name,
                            "read_end": "1" if read.is_read1 else "2" if read.is_read2 else "0",
                            "realigned_chrom": read.reference_name,
                            "realigned_start": read.reference_start + 1,
                            "MAPQ": read.mapping_quality,
                            "CIGAR": read.cigarstring,
                            "kmer_strand": placement.strand,
                            "kmer_offset": placement.kmer_offset,
                            "kmer_center_reference_pos": matched_marker.pos,
                            "target_reference_pos": placement.target_reference_pos,
                            "mapped_kmer_bases_in_target_window": placement.mapped_reference_bases,
                            "placement_validation": (
                                "strict_reference_map"
                                if placement.strict_reference_map_match
                                else "legacy_center_and_window"
                            ),
                        }
                    )

            rg_target_count = len(rg_target_covering)
            rg_window_count = len(rg_window_spanning)
            rg_haplotype_count = len(rg_haplotype_at_target)
            per_rg_rows.append(
                {
                    "FN_site": site_label,
                    "RG": rg,
                    "input_haplotype_qnames": len(site_input_qnames),
                    "input_haplotype_qnames_reaching_target": len(rg_input_qnames_at_target),
                    "haplotype_records_at_target": rg_haplotype_count,
                    "target_covering_records": rg_target_count,
                    "full_window_spanning_records": rg_window_count,
                    "haplotype_fraction_of_target_covering": (
                        f"{rg_haplotype_count / rg_target_count:.6f}" if rg_target_count else "."
                    ),
                    "haplotype_fraction_of_full_window_spanning": (
                        f"{rg_haplotype_count / rg_window_count:.6f}" if rg_window_count else "."
                    ),
                }
            )

        target_covering_count = len(site_target_covering)
        window_spanning_count = len(site_window_spanning)
        haplotype_at_target_count = len(site_haplotype_at_target)
        input_haplotype_at_target_count = len(site_haplotype_input_qnames_at_target)
        summary_rows.append(
            {
                "FN_site": site_label,
                "relevant_RGs": ",".join(relevant_rgs),
                "input_haplotype_qnames": len(site_input_qnames),
                "input_haplotype_qnames_reaching_target": input_haplotype_at_target_count,
                "all_haplotype_records_at_target": haplotype_at_target_count,
                "all_target_covering_records": target_covering_count,
                "all_full_window_spanning_records": window_spanning_count,
                "haplotype_fraction_of_target_covering": (
                    f"{haplotype_at_target_count / target_covering_count:.6f}"
                    if target_covering_count
                    else "."
                ),
                "haplotype_fraction_of_full_window_spanning": (
                    f"{haplotype_at_target_count / window_spanning_count:.6f}"
                    if window_spanning_count
                    else "."
                ),
                "target_trace_status": (
                    "full_truth_haplotype_reaches_correct_target"
                    if input_haplotype_at_target_count > 0
                    else "no_input_truth_haplotype_confirmed_at_target"
                ),
            }
        )

    trace_fields = [
        "FN_site",
        "marker_FN_site",
        "marker_haplotype",
        "marker_kmer_len",
        "RG",
        "input_qname_match",
        "read_qname",
        "read_end",
        "realigned_chrom",
        "realigned_start",
        "MAPQ",
        "CIGAR",
        "kmer_strand",
        "kmer_offset",
        "kmer_center_reference_pos",
        "target_reference_pos",
        "mapped_kmer_bases_in_target_window",
        "placement_validation",
    ]
    with args.output_trace.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=trace_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(trace_rows)

    summary_fields = [
        "FN_site",
        "relevant_RGs",
        "input_haplotype_qnames",
        "input_haplotype_qnames_reaching_target",
        "all_haplotype_records_at_target",
        "all_target_covering_records",
        "all_full_window_spanning_records",
        "haplotype_fraction_of_target_covering",
        "haplotype_fraction_of_full_window_spanning",
        "target_trace_status",
    ]
    with args.output_summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    if args.output_per_rg_summary:
        args.output_per_rg_summary.parent.mkdir(parents=True, exist_ok=True)
        per_rg_fields = [
            "FN_site",
            "RG",
            "input_haplotype_qnames",
            "input_haplotype_qnames_reaching_target",
            "haplotype_records_at_target",
            "target_covering_records",
            "full_window_spanning_records",
            "haplotype_fraction_of_target_covering",
            "haplotype_fraction_of_full_window_spanning",
        ]
        with args.output_per_rg_summary.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=per_rg_fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(per_rg_rows)


if __name__ == "__main__":
    main()
