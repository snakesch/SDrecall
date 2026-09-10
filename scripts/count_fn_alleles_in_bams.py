#!/usr/bin/env python3
"""Count exact REF/ALT alleles for FN sites in one or more indexed BAMs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pysam

from build_fn_truth_haplotypes import Variant, allele_call, read_sites
from trace_fn_truth_haplotypes import resolve_bam_contig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", required=True, type=Path)
    parser.add_argument(
        "--bam",
        action="append",
        required=True,
        help="Named BAM as LABEL=PATH; repeat for multiple stages",
    )
    parser.add_argument("--min-base-quality", type=int, default=0)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def parse_named_bam(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise ValueError(f"Expected LABEL=PATH for --bam, observed {value!r}")
    label, path = value.split("=", 1)
    return label, Path(path)


def main() -> None:
    args = parse_args()
    sites = read_sites(args.sites)
    rows: list[dict[str, object]] = []

    for stage, bam_path in map(parse_named_bam, args.bam):
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for site in sites:
                bam_chrom = resolve_bam_contig(bam, site.chrom)
                variant = Variant(
                    chrom=site.chrom,
                    pos=site.pos,
                    ref=site.ref,
                    alt=site.alt,
                    gt=(),
                    phased=False,
                )
                target0 = site.pos - 1
                target_covering = 0
                ref_records = 0
                alt_records = 0
                other_or_uncalled = 0
                ref_qnames: set[str] = set()
                alt_qnames: set[str] = set()
                fetch_end0 = target0 + max(1, len(site.ref))
                for read in bam.fetch(bam_chrom, target0, fetch_end0):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    mapped_positions = {
                        reference_pos
                        for query_pos, reference_pos in read.get_aligned_pairs(
                            matches_only=False
                        )
                        if query_pos is not None and reference_pos is not None
                    }
                    if target0 not in mapped_positions:
                        continue
                    target_covering += 1
                    call = allele_call(read, variant, args.min_base_quality)
                    if call == "ref":
                        ref_records += 1
                        ref_qnames.add(read.query_name)
                    elif call == "alt":
                        alt_records += 1
                        alt_qnames.add(read.query_name)
                    else:
                        other_or_uncalled += 1

                callable_depth = ref_records + alt_records
                rows.append(
                    {
                        "stage": stage,
                        "FN_site": site.label,
                        "target_covering_records": target_covering,
                        "callable_ref_records": ref_records,
                        "exact_alt_records": alt_records,
                        "other_or_uncalled_records": other_or_uncalled,
                        "callable_depth": callable_depth,
                        "exact_alt_AF_callable": (
                            f"{alt_records / callable_depth:.6f}"
                            if callable_depth
                            else "."
                        ),
                        "exact_alt_AF_target_covering": (
                            f"{alt_records / target_covering:.6f}"
                            if target_covering
                            else "."
                        ),
                        "ref_qnames": len(ref_qnames),
                        "exact_alt_qnames": len(alt_qnames),
                    }
                )

    fields = [
        "stage",
        "FN_site",
        "target_covering_records",
        "callable_ref_records",
        "exact_alt_records",
        "other_or_uncalled_records",
        "callable_depth",
        "exact_alt_AF_callable",
        "exact_alt_AF_target_covering",
        "ref_qnames",
        "exact_alt_qnames",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
