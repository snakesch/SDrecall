#!/usr/bin/env python3
"""Exact, genotype-aware comparison of two normalized single-sample VCFs.

The benchmark wrapper normalizes, atomizes, splits multiallelic records, and
restricts both callsets to the same callable BED before invoking this script.
Consequently, an allele match is the exact tuple CHROM/POS/REF/ALT and a
genotype-aware match additionally requires the same alternate-allele dosage.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path

import vcfpy


@dataclass(frozen=True, order=True)
class AlleleKey:
    chrom: str
    pos: int
    ref: str
    alt: str


@dataclass(frozen=True, order=True)
class GenotypeKey:
    allele: AlleleKey
    alt_dosage: int


@dataclass
class Callset:
    sample: str
    alleles: set[AlleleKey]
    genotypes: set[GenotypeKey]
    skipped_symbolic: int
    skipped_non_variant_gt: int


def parse_alt_dosage(gt: object, source: Path, chrom: str, pos: int) -> int | None:
    if not isinstance(gt, str) or gt in {"", ".", "./.", ".|."}:
        return None

    allele_tokens = re.split(r"[/|]", gt)
    if not allele_tokens or any(token == "." for token in allele_tokens):
        return None

    try:
        alleles = [int(token) for token in allele_tokens]
    except ValueError as exc:
        raise ValueError(f"invalid GT={gt!r} at {chrom}:{pos} in {source}") from exc

    if any(allele < 0 or allele > 1 for allele in alleles):
        raise ValueError(
            f"multiallelic GT={gt!r} at {chrom}:{pos} in {source}; "
            "normalize with bcftools norm -m -both first"
        )
    return sum(allele == 1 for allele in alleles)


def load_callset(path: Path) -> Callset:
    reader = vcfpy.Reader.from_path(str(path))
    samples = reader.header.samples.names
    if len(samples) != 1:
        reader.close()
        raise ValueError(f"expected one sample in {path}, found {samples}")

    alleles: set[AlleleKey] = set()
    genotypes: set[GenotypeKey] = set()
    skipped_symbolic = 0
    skipped_non_variant_gt = 0

    try:
        for record in reader:
            if len(record.ALT) != 1:
                raise ValueError(
                    f"expected one ALT at {record.CHROM}:{record.POS} in {path}; "
                    "normalize with bcftools norm -m -both first"
                )

            alt = record.ALT[0]
            if not isinstance(alt, vcfpy.Substitution):
                skipped_symbolic += 1
                continue

            dosage = parse_alt_dosage(
                record.calls[0].data.get("GT"), path, record.CHROM, record.POS
            )
            if dosage is None or dosage == 0:
                skipped_non_variant_gt += 1
                continue

            allele = AlleleKey(
                chrom=str(record.CHROM),
                pos=int(record.POS),
                ref=str(record.REF).upper(),
                alt=str(alt.value).upper(),
            )
            alleles.add(allele)
            genotypes.add(GenotypeKey(allele=allele, alt_dosage=dosage))
    finally:
        reader.close()

    return Callset(
        sample=samples[0],
        alleles=alleles,
        genotypes=genotypes,
        skipped_symbolic=skipped_symbolic,
        skipped_non_variant_gt=skipped_non_variant_gt,
    )


def ratio(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else math.nan


def f1(precision: float, recall: float) -> float:
    if math.isnan(precision) or math.isnan(recall):
        return math.nan
    if precision + recall == 0:
        return 0.0
    return 2 * precision * recall / (precision + recall)


def metric_counts(called: set[object], truth: set[object]) -> dict[str, int | float]:
    true_positive = len(called & truth)
    false_positive = len(called - truth)
    false_negative = len(truth - called)
    precision = ratio(true_positive, true_positive + false_positive)
    recall = ratio(true_positive, true_positive + false_negative)
    return {
        "called": len(called),
        "truth": len(truth),
        "true_positive": true_positive,
        "false_positive": false_positive,
        "false_negative": false_negative,
        "precision": precision,
        "recall": recall,
        "f1": f1(precision, recall),
    }


def evaluate(called: Callset, truth: Callset) -> dict[str, int | float]:
    genotype = metric_counts(called.genotypes, truth.genotypes)
    allele = metric_counts(called.alleles, truth.alleles)
    called_by_allele: dict[AlleleKey, set[int]] = {}
    truth_by_allele: dict[AlleleKey, set[int]] = {}
    for key in called.genotypes:
        called_by_allele.setdefault(key.allele, set()).add(key.alt_dosage)
    for key in truth.genotypes:
        truth_by_allele.setdefault(key.allele, set()).add(key.alt_dosage)
    discordant = sum(
        not (called_by_allele[allele] & truth_by_allele[allele])
        for allele in called.alleles & truth.alleles
    )
    return {
        **{f"genotype_{key}": value for key, value in genotype.items()},
        **{f"allele_{key}": value for key, value in allele.items()},
        "genotype_discordant_shared_alleles": discordant,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare normalized, biallelic, single-sample VCFs using exact "
            "allele and alternate-dosage matching."
        )
    )
    parser.add_argument("called_vcf", type=Path)
    parser.add_argument("truth_vcf", type=Path)
    parser.add_argument("output_tsv", type=Path)
    parser.add_argument("--label", default="callset")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    called = load_callset(args.called_vcf)
    truth = load_callset(args.truth_vcf)
    if called.sample != truth.sample:
        raise ValueError(
            f"sample mismatch: called={called.sample!r}, truth={truth.sample!r}"
        )

    metrics = evaluate(called, truth)
    row: dict[str, object] = {
        "label": args.label,
        "sample": called.sample,
        "called_vcf": str(args.called_vcf.resolve()),
        "truth_vcf": str(args.truth_vcf.resolve()),
        **metrics,
        "called_skipped_symbolic": called.skipped_symbolic,
        "truth_skipped_symbolic": truth.skipped_symbolic,
        "called_skipped_non_variant_gt": called.skipped_non_variant_gt,
        "truth_skipped_non_variant_gt": truth.skipped_non_variant_gt,
    }

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=row.keys(), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

    print(
        f"{args.label}: genotype precision={metrics['genotype_precision']:.6g} "
        f"recall={metrics['genotype_recall']:.6g}; "
        f"allele precision={metrics['allele_precision']:.6g} "
        f"recall={metrics['allele_recall']:.6g}"
    )


if __name__ == "__main__":
    main()
