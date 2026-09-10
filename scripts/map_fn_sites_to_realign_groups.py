#!/usr/bin/env python3
"""Map FN allele spans to all overlapping SDrecall RG target BEDs."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

from build_fn_truth_haplotypes import Site, read_sites


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", required=True, type=Path)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--output-overlaps", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    return parser.parse_args()


def rg_sort_key(rg: str) -> tuple[int, str]:
    match = re.fullmatch(r"RG(\d+)", rg)
    return (int(match.group(1)), rg) if match else (10**9, rg)


def main() -> None:
    args = parse_args()
    sites = read_sites(args.sites)
    target_beds = sorted((args.run_dir / "realign_groups").glob("RG*/RG*.fc_target.bed"))
    if not target_beds:
        raise FileNotFoundError(f"No RG target BEDs under {args.run_dir}")

    sites_by_chrom: dict[str, list[tuple[Site, int, int]]] = defaultdict(list)
    for site in sites:
        start0 = site.pos - 1
        sites_by_chrom[site.chrom].append((site, start0, start0 + len(site.ref)))

    overlap_rows: list[dict[str, object]] = []
    rgs_by_site: dict[str, set[str]] = defaultdict(set)
    for bed_path in target_beds:
        rg = bed_path.parent.name
        with bed_path.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                chrom, start_text, end_text = fields[:3]
                start0, end0 = int(start_text), int(end_text)
                for site, site_start0, site_end0 in sites_by_chrom.get(chrom, []):
                    if site_start0 >= end0 or site_end0 <= start0:
                        continue
                    rgs_by_site[site.label].add(rg)
                    overlap_rows.append(
                        {
                            "FN_site": site.label,
                            "RG": rg,
                            "site_interval": f"{site.chrom}:{site_start0 + 1}-{site_end0}",
                            "target_interval": f"{chrom}:{start0 + 1}-{end0}",
                            "fc_target_bed": bed_path,
                        }
                    )

    missing = [site.label for site in sites if not rgs_by_site[site.label]]
    if missing:
        raise ValueError(f"FN sites outside every RG target BED: {missing}")

    args.output_overlaps.parent.mkdir(parents=True, exist_ok=True)
    args.output_summary.parent.mkdir(parents=True, exist_ok=True)
    overlap_fields = [
        "FN_site",
        "RG",
        "site_interval",
        "target_interval",
        "fc_target_bed",
    ]
    with args.output_overlaps.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=overlap_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(overlap_rows)

    with args.output_summary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["FN_site", "relevant_RGs"], delimiter="\t"
        )
        writer.writeheader()
        for site in sites:
            writer.writerow(
                {
                    "FN_site": site.label,
                    "relevant_RGs": ",".join(
                        sorted(rgs_by_site[site.label], key=rg_sort_key)
                    ),
                }
            )


if __name__ == "__main__":
    main()
