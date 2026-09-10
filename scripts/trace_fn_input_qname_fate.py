#!/usr/bin/env python3
"""Trace marker-positive input qnames through relevant pre-FP RG raw BAMs."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import pysam

from trace_fn_truth_haplotypes import (
    SiteDefinition,
    kmer_placement,
    qname_base,
    read_definitions,
    read_input_hits,
    read_rg_summary,
    record_id,
    resolve_bam_contig,
    reverse_complement,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--definitions", required=True, type=Path)
    parser.add_argument("--input-hits", required=True, type=Path)
    parser.add_argument("--rg-summary", required=True, type=Path)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--sample", default="HG002")
    parser.add_argument("--output-trace", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    return parser.parse_args()


def marker_in_read(read: pysam.AlignedSegment, marker: SiteDefinition) -> bool:
    sequence = (read.get_forward_sequence() or read.query_sequence or "").upper()
    return marker.kmer in sequence or reverse_complement(marker.kmer) in sequence


def main() -> None:
    args = parse_args()
    definitions = read_definitions(args.definitions)
    input_hits = read_input_hits(args.input_hits)
    rgs_by_site = read_rg_summary(args.rg_summary)

    definitions_by_site: dict[str, list[SiteDefinition]] = defaultdict(list)
    for definition in definitions:
        definitions_by_site[definition.label].append(definition)

    markers_by_target: dict[str, list[SiteDefinition]] = defaultdict(list)
    target_labels = set(definitions_by_site)
    for marker in definitions:
        for supported_site in marker.supported_sites:
            if supported_site in target_labels:
                markers_by_target[supported_site].append(marker)

    input_qnames_by_site: dict[str, set[str]] = {}
    for site in definitions_by_site:
        input_qnames_by_site[site] = set().union(
            *(input_hits.get(marker.label, set()) for marker in markers_by_target[site])
        )

    sites_by_rg: dict[str, set[str]] = defaultdict(set)
    for site, rgs in rgs_by_site.items():
        for rg in rgs:
            sites_by_rg[rg].add(site)

    qnames_present: dict[str, set[str]] = defaultdict(set)
    marker_anywhere: dict[str, set[str]] = defaultdict(set)
    marker_overlapping_target: dict[str, set[str]] = defaultdict(set)
    marker_at_target: dict[str, set[str]] = defaultdict(set)
    marker_loci: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    trace_rows: list[dict[str, object]] = []

    for rg, sites in sorted(sites_by_rg.items()):
        qname_to_sites: dict[str, set[str]] = defaultdict(set)
        for site in sites:
            for qname in input_qnames_by_site[site]:
                qname_to_sites[qname].add(site)

        bam_path = (
            args.run_dir
            / "recall_results"
            / f"{args.sample}.sdrecall.only_{rg}.raw.bam"
        )
        if not bam_path.exists():
            raise FileNotFoundError(bam_path)

        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            bam_contigs = {
                site: resolve_bam_contig(bam, definitions_by_site[site][0].chrom)
                for site in sites
            }
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                base_qname = qname_base(read.query_name)
                matched_sites = qname_to_sites.get(base_qname)
                if not matched_sites:
                    continue
                for site in matched_sites:
                    qnames_present[site].add(base_qname)
                    matched_marker = next(
                        (
                            marker
                            for marker in markers_by_target[site]
                            if marker_in_read(read, marker)
                        ),
                        None,
                    )
                    exact_marker = matched_marker is not None
                    strict_target = False
                    overlaps_target = False
                    if exact_marker and read.reference_name == bam_contigs[site]:
                        target0 = definitions_by_site[site][0].pos - 1
                        overlaps_target = (
                            read.reference_start <= target0 < read.reference_end
                        )
                        strict_target = (
                            kmer_placement(
                                read,
                                matched_marker,
                                definitions_by_site[site][0].pos,
                            )
                            is not None
                        )
                    if exact_marker:
                        marker_anywhere[site].add(base_qname)
                        marker_loci[site].append(
                            (
                                read.reference_name,
                                read.reference_start + 1,
                                read.reference_end,
                            )
                        )
                    if overlaps_target:
                        marker_overlapping_target[site].add(base_qname)
                    if strict_target:
                        marker_at_target[site].add(base_qname)
                    trace_rows.append(
                        {
                            "FN_site": site,
                            "RG": rg,
                            "read_qname": read.query_name,
                            "read_record_id": record_id(read),
                            "mapping_chrom": read.reference_name,
                            "mapping_start": read.reference_start + 1,
                            "mapping_end": read.reference_end,
                            "MAPQ": read.mapping_quality,
                            "CIGAR": read.cigarstring,
                            "exact_marker_in_record": "yes" if exact_marker else "no",
                            "marker_alignment_overlaps_target": (
                                "yes" if overlaps_target else "no"
                            ),
                            "strict_marker_at_target": "yes" if strict_target else "no",
                            "matched_marker_FN_site": (
                                matched_marker.label if matched_marker else "."
                            ),
                            "matched_marker_haplotype": (
                                matched_marker.haplotype if matched_marker else "."
                            ),
                        }
                    )

    trace_fields = [
        "FN_site",
        "RG",
        "read_qname",
        "read_record_id",
        "mapping_chrom",
        "mapping_start",
        "mapping_end",
        "MAPQ",
        "CIGAR",
        "exact_marker_in_record",
        "marker_alignment_overlaps_target",
        "strict_marker_at_target",
        "matched_marker_FN_site",
        "matched_marker_haplotype",
    ]
    args.output_trace.parent.mkdir(parents=True, exist_ok=True)
    with args.output_trace.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=trace_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(trace_rows)

    summary_fields = [
        "FN_site",
        "input_marker_qnames",
        "raw_qnames_present",
        "raw_exact_marker_qnames_anywhere",
        "raw_exact_marker_qnames_overlapping_target",
        "raw_exact_marker_qnames_at_strict_target",
        "raw_exact_marker_mapping_span",
        "pre_fp_qname_fate",
    ]
    args.output_summary.parent.mkdir(parents=True, exist_ok=True)
    with args.output_summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        for site in definitions_by_site:
            loci_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
            for chrom, start, end in marker_loci[site]:
                loci_by_chrom[chrom].append((start, end))
            mapping_span = ";".join(
                f"{chrom}:{min(start for start, _ in spans)}-"
                f"{max(end for _, end in spans)}(n={len(spans)})"
                for chrom, spans in sorted(loci_by_chrom.items())
            ) or "."
            if marker_at_target[site]:
                fate = "recruited_and_correctly_placed"
            elif marker_anywhere[site]:
                fate = "recruited_but_marker_not_at_strict_target"
            elif qnames_present[site]:
                fate = "qname_only_no_marker_bearing_record"
            else:
                fate = "input_marker_qnames_absent_from_relevant_raw_bams"
            writer.writerow(
                {
                    "FN_site": site,
                    "input_marker_qnames": len(input_qnames_by_site[site]),
                    "raw_qnames_present": len(qnames_present[site]),
                    "raw_exact_marker_qnames_anywhere": len(marker_anywhere[site]),
                    "raw_exact_marker_qnames_overlapping_target": len(
                        marker_overlapping_target[site]
                    ),
                    "raw_exact_marker_qnames_at_strict_target": len(
                        marker_at_target[site]
                    ),
                    "raw_exact_marker_mapping_span": mapping_span,
                    "pre_fp_qname_fate": fate,
                }
            )


if __name__ == "__main__":
    main()
