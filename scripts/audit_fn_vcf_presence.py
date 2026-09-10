#!/usr/bin/env python3
"""Report exact FN allele presence across raw, clean, and merged VCF stages."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pysam

from build_fn_truth_haplotypes import read_sites


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", required=True, type=Path)
    parser.add_argument(
        "--vcf",
        action="append",
        required=True,
        help="Named VCF as LABEL=PATH; repeat for multiple stages",
    )
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def parse_named_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise ValueError(f"Expected LABEL=PATH, observed {value!r}")
    label, path = value.split("=", 1)
    return label, Path(path)


def resolve_vcf_contig(vcf: pysam.VariantFile, chrom: str) -> str:
    contigs = vcf.header.contigs
    if chrom in contigs:
        return chrom
    aliases = [chrom.removeprefix("chr")]
    if not chrom.startswith("chr"):
        aliases.append(f"chr{chrom}")
    if chrom in {"chrM", "M"}:
        aliases.extend(["MT", "chrMT"])
    for alias in aliases:
        if alias in contigs:
            return alias
    raise ValueError(f"Contig {chrom!r} is absent from VCF")


def sample_gt(record: pysam.VariantRecord) -> str:
    samples = list(record.samples)
    if not samples:
        return "."
    sample = record.samples[samples[0]]
    gt = sample.get("GT") or ()
    separator = "|" if sample.phased else "/"
    return separator.join("." if allele is None else str(allele) for allele in gt)


def main() -> None:
    args = parse_args()
    sites = read_sites(args.sites)
    rows: list[dict[str, object]] = []

    for stage, vcf_path in map(parse_named_path, args.vcf):
        with pysam.VariantFile(str(vcf_path)) as vcf:
            for site in sites:
                vcf_chrom = resolve_vcf_contig(vcf, site.chrom)
                exact_records: list[pysam.VariantRecord] = []
                position_records: list[str] = []
                for record in vcf.fetch(vcf_chrom, site.pos - 1, site.pos):
                    if record.pos != site.pos:
                        continue
                    alts = tuple(record.alts or ())
                    position_records.append(
                        f"{record.ref}>{','.join(alts)}:{sample_gt(record)}"
                    )
                    if record.ref == site.ref and site.alt in alts:
                        exact_records.append(record)
                rows.append(
                    {
                        "stage": stage,
                        "FN_site": site.label,
                        "exact_allele_present": "yes" if exact_records else "no",
                        "exact_record_count": len(exact_records),
                        "exact_record_GT": (
                            ";".join(sample_gt(record) for record in exact_records)
                            if exact_records
                            else "."
                        ),
                        "position_records": ";".join(position_records) or ".",
                    }
                )

    fields = [
        "stage",
        "FN_site",
        "exact_allele_present",
        "exact_record_count",
        "exact_record_GT",
        "position_records",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
