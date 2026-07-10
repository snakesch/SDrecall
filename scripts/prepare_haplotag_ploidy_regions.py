#!/usr/bin/env python3
"""Split a target BED into ploidy-aware region sets for WhatsHap haplotag."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path


def read_bed3(path: str) -> list[tuple[str, int, int, list[str]]]:
    intervals: list[tuple[str, int, int, list[str]]] = []
    with open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom, start, end = fields[:3]
            intervals.append((chrom, int(start), int(end), fields[3:]))
    return intervals


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def subtract(interval: tuple[int, int], masks: list[tuple[int, int]]) -> list[tuple[int, int]]:
    segments = [interval]
    for mask_start, mask_end in masks:
        next_segments: list[tuple[int, int]] = []
        for start, end in segments:
            if mask_end <= start or mask_start >= end:
                next_segments.append((start, end))
                continue
            if start < mask_start:
                next_segments.append((start, mask_start))
            if mask_end < end:
                next_segments.append((mask_end, end))
        segments = next_segments
        if not segments:
            break
    return segments


def intersections(
    interval: tuple[int, int], masks: list[tuple[int, int]]
) -> list[tuple[int, int]]:
    start, end = interval
    out = []
    for mask_start, mask_end in masks:
        left = max(start, mask_start)
        right = min(end, mask_end)
        if left < right:
            out.append((left, right))
    return out


def is_autosome(chrom: str) -> bool:
    if chrom.startswith("chr"):
        value = chrom[3:]
    else:
        value = chrom
    return value.isdigit() and 1 <= int(value) <= 22


def is_sex(chrom: str) -> bool:
    return chrom in {"X", "Y", "chrX", "chrY"}


def write_bed(path: Path, rows: list[tuple[str, int, int, str]]) -> None:
    with path.open("wt", encoding="utf-8") as handle:
        for chrom, start, end, label in rows:
            if start < end:
                handle.write(f"{chrom}\t{start}\t{end}\t{label}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-bed", required=True)
    parser.add_argument("--par-bed", required=True)
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()

    target = read_bed3(args.target_bed)
    par_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    par_labels: dict[tuple[str, int, int], str] = {}
    for chrom, start, end, extra in read_bed3(args.par_bed):
        par_by_chrom[chrom].append((start, end))
        par_labels[(chrom, start, end)] = extra[0] if extra else "PAR"
    par_by_chrom = {
        chrom: merge_intervals(intervals) for chrom, intervals in par_by_chrom.items()
    }

    autosomes: list[tuple[str, int, int, str]] = []
    sex_par: list[tuple[str, int, int, str]] = []
    sex_nonpar: list[tuple[str, int, int, str]] = []
    skipped: list[tuple[str, int, int, str]] = []

    for chrom, start, end, _extra in target:
        if is_autosome(chrom):
            autosomes.append((chrom, start, end, "ploidy2_autosome"))
        elif is_sex(chrom):
            masks = par_by_chrom.get(chrom, [])
            for left, right in intersections((start, end), masks):
                sex_par.append((chrom, left, right, "ploidy2_PAR"))
            for left, right in subtract((start, end), masks):
                sex_nonpar.append((chrom, left, right, "ploidy1_nonPAR"))
        else:
            skipped.append((chrom, start, end, "skip_non_primary_or_missing_vcf"))

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = {
        "autosomes_ploidy2": autosomes,
        "sex_PAR_ploidy2": sex_par,
        "sex_nonPAR_ploidy1": sex_nonpar,
        "skip_non_primary": skipped,
    }
    for suffix, rows in outputs.items():
        write_bed(prefix.with_suffix(f".{suffix}.bed"), rows)

    summary_path = prefix.with_suffix(".summary.tsv")
    with summary_path.open("wt", encoding="utf-8") as handle:
        handle.write("region_set\tintervals\tbp\n")
        for suffix, rows in outputs.items():
            bp = sum(end - start for _chrom, start, end, _label in rows)
            handle.write(f"{suffix}\t{len(rows)}\t{bp}\n")

    print(summary_path)


if __name__ == "__main__":
    main()
