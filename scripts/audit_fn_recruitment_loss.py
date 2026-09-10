#!/usr/bin/env python3
"""Sub-classify input haplotypes absent from relevant pre-FP RG BAMs."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

from trace_fn_truth_haplotypes import qname_base, read_rg_summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pre-fp-summary", required=True, type=Path)
    parser.add_argument("--input-hits", required=True, type=Path)
    parser.add_argument("--rg-summary", required=True, type=Path)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--sample", default="HG002")
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def load_intervals(paths: list[Path]) -> dict[str, list[tuple[int, int]]]:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for path in paths:
        with path.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3:
                    intervals[fields[0]].append((int(fields[1]), int(fields[2])))
    return intervals


def point_is_covered(
    intervals: dict[str, list[tuple[int, int]]], chrom: str, pos1: int
) -> bool:
    pos0 = pos1 - 1
    return any(start0 <= pos0 < end0 for start0, end0 in intervals.get(chrom, []))


def main() -> None:
    args = parse_args()
    lost_sites = {
        row["FN_site"]
        for row in csv.DictReader(args.pre_fp_summary.open(), delimiter="\t")
        if row["pre_fp_qname_fate"]
        == "input_marker_qnames_absent_from_relevant_raw_bams"
    }
    rgs_by_site = read_rg_summary(args.rg_summary)
    hits_by_site: dict[str, list[dict[str, str]]] = defaultdict(list)
    with args.input_hits.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["FN_site"] in lost_sites:
                hits_by_site[row["FN_site"]].append(row)

    rows: list[dict[str, object]] = []
    for site in sorted(lost_sites):
        relevant_rgs = rgs_by_site[site]
        bed_paths: list[Path] = []
        for rg in relevant_rgs:
            rg_dir = args.run_dir / "realign_groups" / rg
            bed_paths.append(rg_dir / f"{rg}.fc_target.bed")
            bed_paths.extend(sorted(rg_dir.glob(f"{rg}_*.nfc.bed")))
        intervals = load_intervals(bed_paths)

        input_qnames = {qname_base(row["read_qname"]) for row in hits_by_site[site]}
        covered_qnames = {
            qname_base(row["read_qname"])
            for row in hits_by_site[site]
            if point_is_covered(
                intervals, row["mapping_chrom"], int(row["mapping_pos"])
            )
        }
        fastq_qnames: set[str] = set()
        for rg in relevant_rgs:
            for mate in ("r1", "r2"):
                fastq = (
                    args.run_dir
                    / "recall_results"
                    / f"{args.sample}.sdrecall.only_{rg}.{mate}.fastq"
                )
                with fastq.open() as handle:
                    for line_index, line in enumerate(handle):
                        if line_index % 4:
                            continue
                        qname = qname_base(
                            line[1:].split()[0].removesuffix("/1").removesuffix("/2")
                        )
                        if qname in input_qnames:
                            fastq_qnames.add(qname)

        if fastq_qnames:
            subcause = "fastq_present_but_missing_after_realign"
        elif covered_qnames:
            subcause = "fc_nfc_covered_but_missing_from_fastq"
        else:
            subcause = "source_locus_outside_relevant_fc_nfc"
        source_loci = sorted(
            {
                f"{row['mapping_chrom']}:{row['mapping_pos']}"
                for row in hits_by_site[site]
            }
        )
        rows.append(
            {
                "FN_site": site,
                "relevant_RGs": ",".join(relevant_rgs),
                "input_marker_qnames": len(input_qnames),
                "source_locus_covered_qnames": len(covered_qnames),
                "fastq_qnames": len(fastq_qnames),
                "source_mapping_span": (
                    f"{source_loci[0]}..{source_loci[-1]}" if source_loci else "."
                ),
                "recruitment_subcause": subcause,
            }
        )

    fields = [
        "FN_site",
        "relevant_RGs",
        "input_marker_qnames",
        "source_locus_covered_qnames",
        "fastq_qnames",
        "source_mapping_span",
        "recruitment_subcause",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
