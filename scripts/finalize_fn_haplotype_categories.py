#!/usr/bin/env python3
"""Finalize exact-allele FN categories from correct-target haplotype traces."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Evidence:
    source: str
    marker_len: int
    summary: dict[str, str]
    per_rg: list[dict[str, str]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-table", required=True, type=Path)
    parser.add_argument("--per-rg-alt-summary", required=True, type=Path)
    parser.add_argument("--truth-71-summary", required=True, type=Path)
    parser.add_argument("--truth-71-per-rg", required=True, type=Path)
    parser.add_argument("--truth-51-summary", required=True, type=Path)
    parser.add_argument("--truth-51-per-rg", required=True, type=Path)
    parser.add_argument("--observed-71-summary", required=True, type=Path)
    parser.add_argument("--observed-71-per-rg", required=True, type=Path)
    parser.add_argument("--output-table", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    parser.add_argument("--all-genotype-table", type=Path)
    parser.add_argument("--output-all-genotype-table", type=Path)
    parser.add_argument("--output-all-genotype-summary", type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def evidence_map(
    source: str,
    marker_len: int,
    summary_path: Path,
    per_rg_path: Path,
) -> dict[str, Evidence]:
    summaries = {row["FN_site"]: row for row in read_tsv(summary_path)}
    per_rg_by_site: dict[str, list[dict[str, str]]] = {}
    for row in read_tsv(per_rg_path):
        per_rg_by_site.setdefault(row["FN_site"], []).append(row)
    return {
        site: Evidence(source, marker_len, summary, per_rg_by_site.get(site, []))
        for site, summary in summaries.items()
    }


def numeric(value: str) -> float:
    return float(value) if value not in {"", "."} else -1.0


def best_rg(evidence: Evidence) -> dict[str, str]:
    if not evidence.per_rg:
        return {}
    return max(
        evidence.per_rg,
        key=lambda row: (
            numeric(row["haplotype_fraction_of_full_window_spanning"]),
            numeric(row["haplotype_fraction_of_target_covering"]),
        ),
    )


def site_label(row: dict[str, str]) -> str:
    return f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}"


def main() -> None:
    args = parse_args()
    sources = [
        evidence_map(
            "GIAB_ALT_only_71mer",
            71,
            args.truth_71_summary,
            args.truth_71_per_rg,
        ),
        evidence_map(
            "GIAB_ALT_only_51mer",
            51,
            args.truth_51_summary,
            args.truth_51_per_rg,
        ),
        evidence_map(
            "HiFi_observed_ALT_haplotype_71mer",
            71,
            args.observed_71_summary,
            args.observed_71_per_rg,
        ),
    ]

    evidence_by_site: dict[str, list[Evidence]] = {}
    for source in sources:
        for site, evidence in source.items():
            if (
                evidence.summary["target_trace_status"]
                == "full_truth_haplotype_reaches_correct_target"
                and int(evidence.summary["input_haplotype_qnames_reaching_target"]) > 0
            ):
                evidence_by_site.setdefault(site, []).append(evidence)

    base_rows = read_tsv(args.base_table)
    per_rg_alt = {row["site"]: row for row in read_tsv(args.per_rg_alt_summary)}
    output_rows: list[dict[str, str]] = []
    verified_sites = 0
    for row in base_rows:
        site = site_label(row)
        exact_support = per_rg_alt.get(site)
        if exact_support is None:
            raise ValueError(f"Missing per-RG exact support for {site}")
        supporting_rgs = exact_support["supporting_rgs"].split(";", 1)[0]
        if supporting_rgs:
            rg_name, fraction, _exact_af = supporting_rgs.split(":", 2)
            exact_alt, exact_dp = fraction.split("/", 1)
        else:
            rg_name, exact_alt, exact_dp, exact_af = ".", "0", "0", "0.0000"
        row = dict(row)
        row["best_rg"] = rg_name
        row["best_rg_alt"] = exact_alt
        row["best_rg_dp"] = exact_dp
        row["max_rg_af"] = exact_support["max_per_rg_af"]
        candidates = evidence_by_site.get(site, [])
        selected = max(
            candidates,
            key=lambda item: (
                item.marker_len,
                item.source == "HiFi_observed_ALT_haplotype_71mer",
            ),
            default=None,
        )

        output = dict(row)
        output.update(
            {
                "haplotype_marker_source": ".",
                "haplotype_marker_len": ".",
                "haplotype_input_qnames": row.get("whole_haplotype_input_hits", "."),
                "haplotype_input_qnames_reaching_target": ".",
                "haplotype_best_rg": ".",
                "haplotype_best_rg_records": ".",
                "haplotype_best_rg_full_window_records": ".",
                "haplotype_best_rg_fraction": ".",
                "haplotype_target_verified": "no",
            }
        )

        if selected is not None:
            verified_sites += 1
            selected_rg = best_rg(selected)
            input_qnames = selected.summary["input_haplotype_qnames"]
            reaching_target = selected.summary["input_haplotype_qnames_reaching_target"]
            output["category"] = "Cat5_haplotype_verified_submerged"
            output["reason"] = (
                f"ALT-containing {selected.marker_len}-bp local haplotype is present in the "
                f"selected input; the same sequence/read names realign across the correct FN "
                f"target. Best exact-allele AD/DP is {row['best_rg_alt']}/{row['best_rg_dp']} "
                f"(AF={row['max_rg_af']}, below 0.20), and bcftools call -mv emits no allele; "
                f"the true haplotype is diluted by excess reference-like/paralog alignments."
            )
            output.update(
                {
                    "haplotype_marker_source": selected.source,
                    "haplotype_marker_len": str(selected.marker_len),
                    "haplotype_input_qnames": input_qnames,
                    "haplotype_input_qnames_reaching_target": reaching_target,
                    "haplotype_best_rg": selected_rg.get("RG", "."),
                    "haplotype_best_rg_records": selected_rg.get(
                        "haplotype_records_at_target", "."
                    ),
                    "haplotype_best_rg_full_window_records": selected_rg.get(
                        "full_window_spanning_records", "."
                    ),
                    "haplotype_best_rg_fraction": selected_rg.get(
                        "haplotype_fraction_of_full_window_spanning", "."
                    ),
                    "haplotype_target_verified": "yes",
                }
            )
        output_rows.append(output)

    if verified_sites != 15:
        raise ValueError(f"Expected 15 haplotype-verified sites, found {verified_sites}")

    args.output_table.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(output_rows[0])
    with args.output_table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    counts = Counter(row["category"] for row in output_rows)
    with args.output_summary.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["category", "count", "fraction_of_16_allele_FNs"])
        for category, count in counts.most_common():
            writer.writerow([category, count, f"{count / len(output_rows):.6f}"])

    all_genotype_args = [
        args.all_genotype_table,
        args.output_all_genotype_table,
        args.output_all_genotype_summary,
    ]
    if any(all_genotype_args) and not all(all_genotype_args):
        raise ValueError(
            "--all-genotype-table, --output-all-genotype-table, and "
            "--output-all-genotype-summary must be supplied together"
        )
    if args.all_genotype_table:
        exact_by_site = {site_label(row): row for row in output_rows}
        all_rows = read_tsv(args.all_genotype_table)
        for row in all_rows:
            site = site_label(row)
            exact = exact_by_site.get(site)
            row["allele_fn_category"] = exact["category"] if exact else "."
            row["allele_fn_reason"] = exact["reason"] if exact else "."
            if exact is None:
                continue
            if exact["category"] == "Cat5_haplotype_verified_submerged":
                row["root_cause_category"] = "haplotype_verified_submerged_before_calling"
                row["mechanistic_category"] = (
                    "true_haplotype_diluted_below_caller_threshold"
                )
            elif exact["category"] == "Cat1_raw_only_input_loss":
                row["root_cause_category"] = "selected_input_haplotype_loss"
                row["mechanistic_category"] = "recoverable_haplotype_absent_from_selected_input"

        args.output_all_genotype_table.parent.mkdir(parents=True, exist_ok=True)
        all_fields = list(all_rows[0])
        with args.output_all_genotype_table.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=all_fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(all_rows)

        all_counts = Counter(row["root_cause_category"] for row in all_rows)
        with args.output_all_genotype_summary.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["root_cause_category", "count", "fraction_of_63_genotype_FNs"])
            for category, count in all_counts.most_common():
                writer.writerow([category, count, f"{count / len(all_rows):.6f}"])


if __name__ == "__main__":
    main()
