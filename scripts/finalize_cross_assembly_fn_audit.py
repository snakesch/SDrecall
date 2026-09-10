#!/usr/bin/env python3
"""Join strict haplotype audit stages into the established FN categories."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

from build_fn_truth_haplotypes import read_sites


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--sites", required=True, type=Path)
    parser.add_argument("--marker-selection", required=True, type=Path)
    parser.add_argument("--pre-fp-summary", required=True, type=Path)
    parser.add_argument("--pre-fp-per-rg-summary", required=True, type=Path)
    parser.add_argument("--fp-control-summary", required=True, type=Path)
    parser.add_argument("--pooled-alleles", required=True, type=Path)
    parser.add_argument("--vcf-presence", required=True, type=Path)
    parser.add_argument("--rg-summary", required=True, type=Path)
    parser.add_argument("--per-rg-alleles", required=True, type=Path)
    parser.add_argument("--recruitment-subcauses", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    return parser.parse_args()


def read_dicts(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def int_value(value: str) -> int:
    return int(value) if value not in {"", "."} else 0


def classify(
    pre: dict[str, str],
    recruitment_subcause: str | None,
    raw_exact_allele_present: str,
) -> tuple[str, str]:
    if int_value(pre["input_marker_qnames"]) == 0:
        return (
            "pending",
            "selected input has no verified haplotype; resolve as Cat1, Cat1b, or Cat1c using the upstream input audit",
        )

    pre_fate = pre["pre_fp_qname_fate"]
    if pre_fate == "input_marker_qnames_absent_from_relevant_raw_bams":
        if recruitment_subcause == "source_locus_outside_relevant_fc_nfc":
            return (
                "Cat2",
                "input haplotype is absent from the relevant realigned RG BAM because its source locus is outside the relevant FC/NFC network",
            )
        return (
            "Cat3",
            "input haplotype source is covered by the relevant FC/NFC network but the qname is lost before or during relevant RG extraction/realignment",
        )
    if pre_fate == "recruited_but_marker_not_at_strict_target":
        overlaps = int_value(pre["raw_exact_marker_qnames_overlapping_target"])
        if overlaps:
            return (
                "pending",
                "realigned marker-bearing reads overlap the FN target but do not retain the strict truth-haplotype coordinate map",
            )
        return (
            "Cat4",
            "input haplotype qnames are recruited but realign to a different paralogous locus rather than the FN target",
        )
    if pre_fate != "recruited_and_correctly_placed":
        return ("pending", pre_fate)

    if raw_exact_allele_present != "no":
        return (
            "pending",
            "strict haplotype reaches the target and the raw VCF contains the exact allele; investigate benchmark merge or normalization",
        )
    return (
        "Cat5",
        "strict truth haplotype reaches the correct pre-FP target at low true-haplotype fraction, but the exact raw VCF allele is not emitted",
    )


def main() -> None:
    args = parse_args()
    sites = read_sites(args.sites)

    marker_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_dicts(args.marker_selection):
        marker_rows[row["FN_site"]].append(row)
    pre_rows = {row["FN_site"]: row for row in read_dicts(args.pre_fp_summary)}
    pre_rg_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_dicts(args.pre_fp_per_rg_summary):
        pre_rg_rows[row["FN_site"]].append(row)
    fp_rows = {row["FN_site"]: row for row in read_dicts(args.fp_control_summary)}
    pooled_rows = {
        (row["stage"], row["FN_site"]): row for row in read_dicts(args.pooled_alleles)
    }
    vcf_rows = {
        (row["stage"], row["FN_site"]): row for row in read_dicts(args.vcf_presence)
    }
    rgs_by_site = {
        row["FN_site"]: [rg for rg in row["relevant_RGs"].split(",") if rg]
        for row in read_dicts(args.rg_summary)
    }
    per_rg_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_dicts(args.per_rg_alleles):
        if row["stage"] in rgs_by_site.get(row["FN_site"], []):
            per_rg_rows[row["FN_site"]].append(row)
    recruitment_subcauses = {
        row["FN_site"]: row["recruitment_subcause"]
        for row in read_dicts(args.recruitment_subcauses)
    }

    output_rows: list[dict[str, object]] = []
    categories: Counter[str] = Counter()
    for site in sites:
        marker = marker_rows[site.label]
        pre = pre_rows[site.label]
        fp = fp_rows[site.label]
        raw = pooled_rows[("pooled_raw", site.label)]
        clean = pooled_rows[("pooled_clean", site.label)]
        raw_vcf = vcf_rows[("raw", site.label)]
        clean_vcf = vcf_rows[("clean", site.label)]
        merged_vcf = vcf_rows[("raw_clean_merged", site.label)]
        final_vcf = vcf_rows[("sdrecall_deepvariant_benchmark", site.label)]
        category, evidence = classify(
            pre,
            recruitment_subcauses.get(site.label),
            raw_vcf["exact_allele_present"],
        )
        categories[category] += 1

        site_pre_rg_rows = pre_rg_rows[site.label]
        max_true_haplotype_rg = max(
            site_pre_rg_rows,
            key=lambda row: float(row["haplotype_fraction_of_target_covering"]),
            default=None,
        )

        rg_rows = per_rg_rows[site.label]
        max_rg = max(
            rg_rows,
            key=lambda row: (
                float(row["exact_alt_AF_target_covering"])
                if row["exact_alt_AF_target_covering"] != "."
                else -1.0
            ),
            default=None,
        )
        output_rows.append(
            {
                "assembly": args.assembly,
                "FN_site": site.label,
                "truth_gt": next(
                    (
                        row.get("truth_gt", ".")
                        for row in read_dicts(args.sites)
                        if row.get("chrom") == site.chrom
                        and row.get("pos") == str(site.pos)
                        and row.get("ref") == site.ref
                        and row.get("alt") == site.alt
                    ),
                    ".",
                ),
                "selected_marker_lengths": ",".join(
                    sorted({row["selected_marker_len"] for row in marker}, key=int)
                ),
                "hifi_target_qnames": max(
                    int_value(row["target_haplotype_qnames"]) for row in marker
                ),
                "input_marker_qnames": int_value(pre["input_marker_qnames"]),
                "pre_fp_qnames_present": int_value(pre["raw_qnames_present"]),
                "pre_fp_marker_qnames_overlapping_target": int_value(
                    pre["raw_exact_marker_qnames_overlapping_target"]
                ),
                "pre_fp_strict_target_qnames": int_value(
                    pre["raw_exact_marker_qnames_at_strict_target"]
                ),
                "pre_fp_strict_target_records": int_value(
                    fp["pre_fp_input_haplotype_records_at_target"]
                ),
                "pre_fp_marker_mapping_span": pre["raw_exact_marker_mapping_span"],
                "fp_correct_qnames": int_value(fp["fp_correct_qnames"]),
                "fp_mismap_qnames": int_value(fp["fp_mismap_qnames"]),
                "fp_lowqual_qnames": int_value(fp["fp_lowqual_qnames"]),
                "clean_strict_input_haplotype_qnames": int_value(
                    fp["clean_input_haplotype_qnames_at_target"]
                ),
                "clean_strict_input_haplotype_records": int_value(
                    fp["clean_input_haplotype_records_at_target"]
                ),
                "pooled_raw_exact_AD": int_value(raw["exact_alt_records"]),
                "pooled_raw_target_DP": int_value(raw["target_covering_records"]),
                "pooled_clean_exact_AD": int_value(clean["exact_alt_records"]),
                "pooled_clean_target_DP": int_value(clean["target_covering_records"]),
                "pre_fp_true_haplotype_fraction": (
                    f"{int_value(fp['pre_fp_input_haplotype_records_at_target']) / int_value(raw['target_covering_records']):.6f}"
                    if int_value(raw["target_covering_records"])
                    else "."
                ),
                "clean_true_haplotype_fraction": (
                    f"{int_value(fp['clean_input_haplotype_records_at_target']) / int_value(clean['target_covering_records']):.6f}"
                    if int_value(clean["target_covering_records"])
                    else "."
                ),
                "max_true_haplotype_RG": (
                    max_true_haplotype_rg["RG"] if max_true_haplotype_rg else "."
                ),
                "max_per_RG_true_haplotype_records": (
                    int_value(max_true_haplotype_rg["haplotype_records_at_target"])
                    if max_true_haplotype_rg
                    else 0
                ),
                "max_per_RG_target_DP": (
                    int_value(max_true_haplotype_rg["target_covering_records"])
                    if max_true_haplotype_rg
                    else 0
                ),
                "max_per_RG_true_haplotype_fraction": (
                    max_true_haplotype_rg["haplotype_fraction_of_target_covering"]
                    if max_true_haplotype_rg
                    else "."
                ),
                "max_relevant_RG": max_rg["stage"] if max_rg else ".",
                "max_relevant_RG_exact_AD": (
                    int_value(max_rg["exact_alt_records"]) if max_rg else 0
                ),
                "max_relevant_RG_target_DP": (
                    int_value(max_rg["target_covering_records"]) if max_rg else 0
                ),
                "max_relevant_RG_exact_AF": (
                    max_rg["exact_alt_AF_target_covering"] if max_rg else "."
                ),
                "raw_exact_allele_in_vcf": raw_vcf["exact_allele_present"],
                "clean_exact_allele_in_vcf": clean_vcf["exact_allele_present"],
                "raw_clean_exact_allele_in_vcf": merged_vcf["exact_allele_present"],
                "benchmark_merged_exact_allele_in_vcf": final_vcf[
                    "exact_allele_present"
                ],
                "same_position_production_records": merged_vcf["position_records"],
                "downstream_fp_control_fate": fp["fp_control_fate"],
                "FN_category": category,
                "category_evidence": evidence,
            }
        )

    fields = list(output_rows[0])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    args.output_summary.parent.mkdir(parents=True, exist_ok=True)
    with args.output_summary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["assembly", "FN_category", "count"], delimiter="\t"
        )
        writer.writeheader()
        for category, count in sorted(categories.items()):
            writer.writerow(
                {"assembly": args.assembly, "FN_category": category, "count": count}
            )


if __name__ == "__main__":
    main()
