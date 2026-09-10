#!/usr/bin/env python3
"""Validate exact truth-haplotype markers in an indexed BAM at their targets."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import pysam

from trace_fn_truth_haplotypes import (
    kmer_placement,
    qname_base,
    read_definitions,
    read_input_hits,
    record_id,
    resolve_bam_contig,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--definitions", required=True, type=Path)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--source-label", required=True)
    parser.add_argument("--input-hits", type=Path)
    parser.add_argument("--output-hits", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    definitions = read_definitions(args.definitions)
    input_hits = read_input_hits(args.input_hits) if args.input_hits else {}
    input_qnames_by_target: dict[str, set[str]] = {}
    target_labels = {marker.label for marker in definitions}
    for target in target_labels:
        input_qnames_by_target[target] = set().union(
            *(
                input_hits.get(marker.label, set())
                for marker in definitions
                if target in marker.supported_sites
            )
        )
    args.output_hits.parent.mkdir(parents=True, exist_ok=True)
    args.output_summary.parent.mkdir(parents=True, exist_ok=True)

    hit_rows: list[dict[str, object]] = []
    qnames_by_marker: dict[tuple[str, str, int], set[str]] = defaultdict(set)
    records_by_marker: dict[tuple[str, str, int], set[str]] = defaultdict(set)
    input_qnames_by_marker: dict[tuple[str, str, int], set[str]] = defaultdict(set)

    with pysam.AlignmentFile(str(args.bam), "rb") as bam:
        for marker in definitions:
            marker_key = (marker.label, marker.haplotype, len(marker.kmer))
            bam_chrom = resolve_bam_contig(bam, marker.chrom)
            fetch_start0 = max(0, marker.window_start - 1)
            fetch_end0 = marker.window_end
            for read in bam.fetch(bam_chrom, fetch_start0, fetch_end0):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                placement = kmer_placement(read, marker, marker.pos)
                if placement is None:
                    continue
                qnames_by_marker[marker_key].add(qname_base(read.query_name))
                records_by_marker[marker_key].add(record_id(read))
                input_match = qname_base(read.query_name) in input_qnames_by_target.get(
                    marker.label, set()
                )
                if input_match:
                    input_qnames_by_marker[marker_key].add(qname_base(read.query_name))
                hit_rows.append(
                    {
                        "FN_site": marker.label,
                        "haplotype": marker.haplotype,
                        "marker_len": len(marker.kmer),
                        "source": args.source_label,
                        "input_qname_match": "yes" if input_match else "no",
                        "read_qname": read.query_name,
                        "read_end": (
                            "1" if read.is_read1 else "2" if read.is_read2 else "0"
                        ),
                        "mapping_chrom": read.reference_name,
                        "mapping_start": read.reference_start + 1,
                        "MAPQ": read.mapping_quality,
                        "CIGAR": read.cigarstring,
                        "marker_strand": placement.strand,
                        "marker_offset": placement.kmer_offset,
                        "placement_validation": "strict_reference_map",
                    }
                )

    hit_fields = [
        "FN_site",
        "haplotype",
        "marker_len",
        "source",
        "input_qname_match",
        "read_qname",
        "read_end",
        "mapping_chrom",
        "mapping_start",
        "MAPQ",
        "CIGAR",
        "marker_strand",
        "marker_offset",
        "placement_validation",
    ]
    with args.output_hits.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=hit_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(hit_rows)

    summary_fields = [
        "FN_site",
        "haplotype",
        "marker_len",
        "source",
        "target_haplotype_qnames",
        "target_haplotype_records",
        "input_target_haplotype_qnames",
        "validation_status",
    ]
    with args.output_summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        for marker in definitions:
            marker_key = (marker.label, marker.haplotype, len(marker.kmer))
            qname_count = len(qnames_by_marker[marker_key])
            writer.writerow(
                {
                    "FN_site": marker.label,
                    "haplotype": marker.haplotype,
                    "marker_len": len(marker.kmer),
                    "source": args.source_label,
                    "target_haplotype_qnames": qname_count,
                    "target_haplotype_records": len(records_by_marker[marker_key]),
                    "input_target_haplotype_qnames": len(
                        input_qnames_by_marker[marker_key]
                    ),
                    "validation_status": (
                        "exact_marker_verified_at_target"
                        if qname_count > 0
                        else "no_exact_marker_at_target"
                    ),
                }
            )


if __name__ == "__main__":
    main()
