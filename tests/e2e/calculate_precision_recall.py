#!/usr/bin/env python3
"""
calculate_precision_recall.py — Standalone VCF precision/recall calculator

Compares a called VCF against a gold-standard VCF, reporting:
  - Recall  = TP / (TP + FN)
  - Precision = TP / (TP + FP)

Always runs BOTH comparison modes in a single invocation:
  - "site": match on CHROM + POS + REF + ALT only (ignore genotype)
  - "genotype": match on CHROM + POS + REF + ALT + GT (het vs hom-alt matters)

Output TSV contains one row per mode so downstream consumers see the
difference between site-level recovery and zygosity-accurate recovery.

Designed for the Rust SDrecall e2e benchmark. No dependency on
external annotation (gnomAD, CADD, etc.) — pure site-level comparison.
"""

import argparse
import sys
from collections import namedtuple
from pathlib import Path

try:
    import pysam
except ImportError:
    sys.exit("ERROR: pysam is required. Install with: pip install pysam")


VariantSite = namedtuple("VariantSite", ["chrom", "pos", "ref", "alt"])
VariantGT = namedtuple("VariantGT", ["chrom", "pos", "ref", "alt", "gt"])


def classify_gt(gt_tuple) -> int:
    """
    Classify genotype:
      0 = hom-ref or missing
      1 = het
      2 = hom-alt
    """
    if gt_tuple is None:
        return 0
    alleles = [a for a in gt_tuple if a is not None]
    if len(alleles) == 0:
        return 0
    if all(a == 0 for a in alleles):
        return 0
    if all(a > 0 for a in alleles) and len(set(alleles)) == 1:
        return 2
    if any(a > 0 for a in alleles):
        return 1
    return 0


def parse_vcf(vcf_path: str, exclude_filters: set = None) -> tuple:
    """
    Parse a VCF file and return TWO sets:
      - site_variants: set of VariantSite (CHROM, POS, REF, ALT)
      - gt_variants:   set of VariantGT   (CHROM, POS, REF, ALT, gt_type)

    Only non-ref genotype variants (het=1, hom-alt=2) are included.
    Variants matching any filter token in `exclude_filters` are skipped.
    """
    if exclude_filters is None:
        exclude_filters = set()

    site_variants = set()
    gt_variants = set()
    vcf = pysam.VariantFile(vcf_path)

    for rec in vcf:
        chrom = rec.chrom
        pos = rec.pos
        ref = rec.ref

        # Skip variants with excluded filter tokens
        if exclude_filters:
            rec_filters = set(rec.filter)
            if rec_filters & exclude_filters:
                continue

        # Determine genotype
        gt_type = 0
        if rec.samples:
            sample = rec.samples[0]
            gt = sample.get("GT", None)
            gt_type = classify_gt(gt)

            # Also check AD: if ref_AD > 0 and alt_AD == 0, it's effectively ref
            if "AD" in sample and gt_type > 0:
                ad = sample["AD"]
                if ad is not None and len(ad) >= 2:
                    if ad[0] > 0 and ad[1] == 0:
                        gt_type = 0

        if gt_type == 0:
            continue

        for alt in rec.alts or []:
            if alt == "*":
                continue
            site_variants.add(VariantSite(chrom, pos, ref, alt))
            gt_variants.add(VariantGT(chrom, pos, ref, alt, gt_type))

    vcf.close()
    return site_variants, gt_variants


def parse_vcf_gold(vcf_path: str) -> tuple:
    """
    Parse gold-standard VCF. Similar to parse_vcf but more permissive:
    includes all non-ref variants regardless of filter status.
    """
    site_variants = set()
    gt_variants = set()
    vcf = pysam.VariantFile(vcf_path)

    for rec in vcf:
        chrom = rec.chrom
        pos = rec.pos
        ref = rec.ref

        gt_type = 0
        if rec.samples:
            sample = rec.samples[0]
            gt = sample.get("GT", None)
            gt_type = classify_gt(gt)

        # Gold VCF: if there's no sample column, treat as het (genotype unknown)
        if gt_type == 0 and not rec.samples:
            gt_type = 1

        # Skip hom-ref
        if gt_type == 0:
            continue

        for alt in rec.alts or []:
            if alt == "*":
                continue
            site_variants.add(VariantSite(chrom, pos, ref, alt))
            gt_variants.add(VariantGT(chrom, pos, ref, alt, gt_type))

    vcf.close()
    return site_variants, gt_variants


def calculate_metrics(gold: set, called: set) -> dict:
    """Calculate precision, recall, F1 and variant counts."""
    tp = gold & called
    fn = gold - called
    fp = called - gold

    tp_count = len(tp)
    fn_count = len(fn)
    fp_count = len(fp)

    recall = tp_count / (tp_count + fn_count) if (tp_count + fn_count) > 0 else 0.0
    precision = tp_count / (tp_count + fp_count) if (tp_count + fp_count) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    return {
        "tp": tp_count,
        "fn": fn_count,
        "fp": fp_count,
        "recall": recall,
        "precision": precision,
        "f1": f1,
        "gold_total": len(gold),
        "called_total": len(called),
    }


def calculate_gt_concordance(gold_gt: set, called_gt: set, gold_site: set, called_site: set) -> dict:
    """
    Among site-level TPs, how many have correct zygosity?
    Returns het→hom and hom→het error counts.
    """
    site_tp = gold_site & called_site

    het_to_hom = 0
    hom_to_het = 0
    gt_concordant = 0

    gold_gt_map = {(v.chrom, v.pos, v.ref, v.alt): v.gt for v in gold_gt}
    called_gt_map = {(v.chrom, v.pos, v.ref, v.alt): v.gt for v in called_gt}

    for site in site_tp:
        key = (site.chrom, site.pos, site.ref, site.alt)
        g_gt = gold_gt_map.get(key, 0)
        c_gt = called_gt_map.get(key, 0)

        if g_gt == c_gt:
            gt_concordant += 1
        elif g_gt == 1 and c_gt == 2:
            het_to_hom += 1
        elif g_gt == 2 and c_gt == 1:
            hom_to_het += 1

    return {
        "gt_concordant": gt_concordant,
        "het_to_hom": het_to_hom,
        "hom_to_het": hom_to_het,
        "site_tp_total": len(site_tp),
        "gt_accuracy": gt_concordant / len(site_tp) if len(site_tp) > 0 else 0.0,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Calculate precision and recall between gold and called VCFs (both site-level and genotype-level)"
    )
    parser.add_argument("--gold", required=True, help="Gold-standard VCF path")
    parser.add_argument("--called", required=True, help="Called VCF path")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("--assembly", default="unknown", help="Assembly tag (hg19/hg38/t2t)")
    parser.add_argument("--caller", default="SDrecall", help="Caller tag for output")
    parser.add_argument(
        "--exclude-filters",
        nargs="*",
        default=["MISALIGNED"],
        help="FILTER tokens to exclude from called VCF (default: MISALIGNED). "
             "Variants with ANY of these tokens are skipped.",
    )
    args = parser.parse_args()

    exclude_set = set(args.exclude_filters) if args.exclude_filters else set()

    print(f"[INFO] Gold VCF:   {args.gold}", file=sys.stderr)
    print(f"[INFO] Called VCF: {args.called}", file=sys.stderr)
    print(f"[INFO] Assembly:   {args.assembly}", file=sys.stderr)
    print(f"[INFO] Exclude filters: {exclude_set or '(none)'}", file=sys.stderr)
    print(f"[INFO] Running BOTH site-level and genotype-level comparisons", file=sys.stderr)

    # Parse both VCFs
    gold_sites, gold_gts = parse_vcf_gold(args.gold)
    called_sites, called_gts = parse_vcf(args.called, exclude_filters=exclude_set)

    print(f"[INFO] Gold variants:   {len(gold_sites)} sites, {len(gold_gts)} genotyped", file=sys.stderr)
    print(f"[INFO] Called variants: {len(called_sites)} sites, {len(called_gts)} genotyped", file=sys.stderr)

    # Mode 1: Site-level (ignore zygosity)
    site_metrics = calculate_metrics(gold_sites, called_sites)
    print(f"[INFO] SITE-LEVEL:     TP={site_metrics['tp']}, FN={site_metrics['fn']}, FP={site_metrics['fp']}", file=sys.stderr)
    print(f"[INFO]   Recall={site_metrics['recall']:.4f}, Precision={site_metrics['precision']:.4f}, F1={site_metrics['f1']:.4f}", file=sys.stderr)

    # Mode 2: Genotype-level (zygosity must match)
    gt_metrics = calculate_metrics(gold_gts, called_gts)
    print(f"[INFO] GENOTYPE-LEVEL: TP={gt_metrics['tp']}, FN={gt_metrics['fn']}, FP={gt_metrics['fp']}", file=sys.stderr)
    print(f"[INFO]   Recall={gt_metrics['recall']:.4f}, Precision={gt_metrics['precision']:.4f}, F1={gt_metrics['f1']:.4f}", file=sys.stderr)

    # Zygosity concordance among site-level TPs
    gt_conc = calculate_gt_concordance(gold_gts, called_gts, gold_sites, called_sites)
    print(f"[INFO] GT CONCORDANCE: {gt_conc['gt_concordant']}/{gt_conc['site_tp_total']} = {gt_conc['gt_accuracy']:.4f}", file=sys.stderr)
    print(f"[INFO]   het→hom: {gt_conc['het_to_hom']}, hom→het: {gt_conc['hom_to_het']}", file=sys.stderr)

    # Write output TSV (two rows: one per mode)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sample_id = Path(args.called).stem.split(".")[0]

    headers = [
        "sample", "assembly", "caller", "mode",
        "recall", "precision", "f1",
        "tp", "fn", "fp",
        "gold_total", "called_total",
        "gt_accuracy", "het_to_hom", "hom_to_het",
    ]

    with open(output_path, "w") as f:
        f.write("\t".join(headers) + "\n")

        # Row 1: site-level (ignore zygosity)
        f.write("\t".join([
            sample_id, args.assembly, args.caller, "site",
            f"{site_metrics['recall']:.6f}", f"{site_metrics['precision']:.6f}", f"{site_metrics['f1']:.6f}",
            str(site_metrics["tp"]), str(site_metrics["fn"]), str(site_metrics["fp"]),
            str(site_metrics["gold_total"]), str(site_metrics["called_total"]),
            f"{gt_conc['gt_accuracy']:.6f}", str(gt_conc["het_to_hom"]), str(gt_conc["hom_to_het"]),
        ]) + "\n")

        # Row 2: genotype-level (zygosity matters)
        f.write("\t".join([
            sample_id, args.assembly, args.caller, "genotype",
            f"{gt_metrics['recall']:.6f}", f"{gt_metrics['precision']:.6f}", f"{gt_metrics['f1']:.6f}",
            str(gt_metrics["tp"]), str(gt_metrics["fn"]), str(gt_metrics["fp"]),
            str(gt_metrics["gold_total"]), str(gt_metrics["called_total"]),
            f"{gt_conc['gt_accuracy']:.6f}", str(gt_conc["het_to_hom"]), str(gt_conc["hom_to_het"]),
        ]) + "\n")

    print(f"[INFO] Results written to {args.output}", file=sys.stderr)

    # Summary pass/fail
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"  SITE-LEVEL:     Recall={site_metrics['recall']:.4f}  Precision={site_metrics['precision']:.4f}", file=sys.stderr)
    print(f"  GENOTYPE-LEVEL: Recall={gt_metrics['recall']:.4f}  Precision={gt_metrics['precision']:.4f}", file=sys.stderr)
    print(f"  GT ACCURACY:    {gt_conc['gt_accuracy']:.4f} ({gt_conc['gt_concordant']}/{gt_conc['site_tp_total']})", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)

    # Pass/fail on site-level thresholds
    if site_metrics["recall"] >= 0.93 and site_metrics["precision"] >= 0.45:
        print("[PASS] Site-level thresholds met (recall>=93%, precision>=45%).", file=sys.stderr)
        return 0
    else:
        print("[FAIL] Site-level thresholds NOT met.", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
