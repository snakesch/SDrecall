#!/usr/bin/env python3
"""Trace strict target haplotypes through island FP-control and clean filtering."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import pysam

from trace_fn_truth_haplotypes import (
    SiteDefinition,
    kmer_placement,
    qname_base,
    read_definitions,
    read_input_hits,
    resolve_bam_contig,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--definitions", required=True, type=Path)
    parser.add_argument("--input-hits", required=True, type=Path)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--sample", default="HG002")
    parser.add_argument("--output-trace", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    return parser.parse_args()


def hp_class(read: pysam.AlignedSegment) -> tuple[str, str]:
    hp = read.get_tag("HP") if read.has_tag("HP") else "."
    if "LOWQUAL" in hp:
        return hp, "lowqual"
    if "HIGHVD" in hp:
        return hp, "mismap"
    if hp != ".":
        return hp, "correct"
    return hp, "untagged"


def island_id(path: Path) -> int:
    match = re.fullmatch(r"island_(\d+)\.bed", path.name)
    if not match:
        raise ValueError(f"Unexpected island BED name: {path}")
    return int(match.group(1))


def main() -> None:
    args = parse_args()
    definitions = read_definitions(args.definitions)
    input_hits = read_input_hits(args.input_hits)

    definitions_by_site: dict[str, list[SiteDefinition]] = defaultdict(list)
    for definition in definitions:
        definitions_by_site[definition.label].append(definition)

    markers_by_target: dict[str, list[SiteDefinition]] = defaultdict(list)
    targets = set(definitions_by_site)
    for marker in definitions:
        for supported_site in marker.supported_sites:
            if supported_site in targets:
                markers_by_target[supported_site].append(marker)

    input_qnames_by_site: dict[str, set[str]] = {}
    for site in definitions_by_site:
        input_qnames_by_site[site] = set().union(
            *(input_hits.get(marker.label, set()) for marker in markers_by_target[site])
        )

    islands_by_site: dict[str, list[int]] = defaultdict(list)
    island_dir = args.run_dir / "intermediates" / "islands"
    island_beds = [
        path
        for path in island_dir.glob("island_*.bed")
        if re.fullmatch(r"island_\d+\.bed", path.name)
    ]
    for bed_path in sorted(island_beds, key=island_id):
        current_island = island_id(bed_path)
        with bed_path.open() as handle:
            intervals = [
                (fields[0], int(fields[1]), int(fields[2]))
                for line in handle
                if line.strip() and not line.startswith("#")
                for fields in [line.rstrip("\n").split("\t")]
            ]
        for site, site_definitions in definitions_by_site.items():
            definition = site_definitions[0]
            target0 = definition.pos - 1
            if any(
                chrom == definition.chrom and start0 <= target0 < end0
                for chrom, start0, end0 in intervals
            ):
                islands_by_site[site].append(current_island)

    trace_rows: list[dict[str, object]] = []
    classes_by_site: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    raw_records_by_site: dict[str, set[str]] = defaultdict(set)
    clean_qnames_by_site: dict[str, set[str]] = defaultdict(set)
    clean_records_by_site: dict[str, set[str]] = defaultdict(set)

    sites_by_island: dict[int, set[str]] = defaultdict(set)
    for site, island_ids in islands_by_site.items():
        for current_island in island_ids:
            sites_by_island[current_island].add(site)

    recall_dir = args.run_dir / "recall_results"
    for current_island, sites in sorted(sites_by_island.items()):
        raw_bam_path = recall_dir / (
            f"{args.sample}.pooled.raw.deduped.{current_island}.bam"
        )
        clean_bam_path = raw_bam_path.with_suffix(".clean.bam")
        if not raw_bam_path.exists():
            raise FileNotFoundError(raw_bam_path)

        with pysam.AlignmentFile(str(raw_bam_path), "rb") as raw_bam:
            for site in sites:
                definition = definitions_by_site[site][0]
                bam_chrom = resolve_bam_contig(raw_bam, definition.chrom)
                for read in raw_bam.fetch(
                    bam_chrom, definition.window_start - 1, definition.window_end
                ):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    base_qname = qname_base(read.query_name)
                    if base_qname not in input_qnames_by_site[site]:
                        continue
                    matched_marker = next(
                        (
                            marker
                            for marker in markers_by_target[site]
                            if kmer_placement(read, marker, definition.pos) is not None
                        ),
                        None,
                    )
                    if matched_marker is None:
                        continue
                    hp, classification = hp_class(read)
                    classes_by_site[site][classification].add(base_qname)
                    raw_records_by_site[site].add(
                        f"{current_island}|{read.query_name}|{read.flag}"
                    )
                    trace_rows.append(
                        {
                            "FN_site": site,
                            "island": current_island,
                            "stage": "island_raw",
                            "read_qname": read.query_name,
                            "mapping_chrom": read.reference_name,
                            "mapping_start": read.reference_start + 1,
                            "MAPQ": read.mapping_quality,
                            "CIGAR": read.cigarstring,
                            "HP": hp,
                            "fp_control_class": classification,
                            "matched_marker_FN_site": matched_marker.label,
                            "matched_marker_haplotype": matched_marker.haplotype,
                        }
                    )

        if not clean_bam_path.exists():
            continue
        with pysam.AlignmentFile(str(clean_bam_path), "rb") as clean_bam:
            for site in sites:
                definition = definitions_by_site[site][0]
                bam_chrom = resolve_bam_contig(clean_bam, definition.chrom)
                for read in clean_bam.fetch(
                    bam_chrom, definition.window_start - 1, definition.window_end
                ):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    base_qname = qname_base(read.query_name)
                    if base_qname not in input_qnames_by_site[site]:
                        continue
                    matched_marker = next(
                        (
                            marker
                            for marker in markers_by_target[site]
                            if kmer_placement(read, marker, definition.pos) is not None
                        ),
                        None,
                    )
                    if matched_marker is None:
                        continue
                    clean_qnames_by_site[site].add(base_qname)
                    clean_records_by_site[site].add(
                        f"{current_island}|{read.query_name}|{read.flag}"
                    )
                    hp, classification = hp_class(read)
                    trace_rows.append(
                        {
                            "FN_site": site,
                            "island": current_island,
                            "stage": "island_clean",
                            "read_qname": read.query_name,
                            "mapping_chrom": read.reference_name,
                            "mapping_start": read.reference_start + 1,
                            "MAPQ": read.mapping_quality,
                            "CIGAR": read.cigarstring,
                            "HP": hp,
                            "fp_control_class": classification,
                            "matched_marker_FN_site": matched_marker.label,
                            "matched_marker_haplotype": matched_marker.haplotype,
                        }
                    )

    trace_fields = [
        "FN_site",
        "island",
        "stage",
        "read_qname",
        "mapping_chrom",
        "mapping_start",
        "MAPQ",
        "CIGAR",
        "HP",
        "fp_control_class",
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
        "islands",
        "pre_fp_input_haplotype_qnames_at_target",
        "pre_fp_input_haplotype_records_at_target",
        "fp_correct_qnames",
        "fp_mismap_qnames",
        "fp_lowqual_qnames",
        "fp_untagged_qnames",
        "clean_input_haplotype_qnames_at_target",
        "clean_input_haplotype_records_at_target",
        "fp_control_fate",
    ]
    args.output_summary.parent.mkdir(parents=True, exist_ok=True)
    with args.output_summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        for site in definitions_by_site:
            correct = classes_by_site[site]["correct"]
            mismap = classes_by_site[site]["mismap"]
            lowqual = classes_by_site[site]["lowqual"]
            untagged = classes_by_site[site]["untagged"]
            raw_qnames = correct | mismap | lowqual | untagged
            clean_qnames = clean_qnames_by_site[site]
            if not raw_qnames:
                fate = "no_strict_pre_fp_target_haplotype"
            elif clean_qnames:
                fate = "full_haplotype_survives_clean_bam"
            elif correct:
                fate = "fp_correct_but_removed_by_clean_policy"
            elif mismap and not lowqual:
                fate = "removed_as_mismap"
            elif lowqual and not mismap:
                fate = "removed_as_lowqual"
            else:
                fate = "removed_by_mixed_fp_control_classes"
            writer.writerow(
                {
                    "FN_site": site,
                    "islands": ",".join(map(str, islands_by_site[site])) or ".",
                    "pre_fp_input_haplotype_qnames_at_target": len(raw_qnames),
                    "pre_fp_input_haplotype_records_at_target": len(
                        raw_records_by_site[site]
                    ),
                    "fp_correct_qnames": len(correct),
                    "fp_mismap_qnames": len(mismap),
                    "fp_lowqual_qnames": len(lowqual),
                    "fp_untagged_qnames": len(untagged),
                    "clean_input_haplotype_qnames_at_target": len(clean_qnames),
                    "clean_input_haplotype_records_at_target": len(
                        clean_records_by_site[site]
                    ),
                    "fp_control_fate": fate,
                }
            )


if __name__ == "__main__":
    main()
