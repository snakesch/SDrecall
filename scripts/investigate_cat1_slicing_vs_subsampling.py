#!/usr/bin/env python3
"""Split Cat1 FN sites into slicing-loss vs subsampling-loss evidence.

This uses the corrected FN category table. Cat1 means raw 300x has the
ALT-carrying read but the SDrecall input BAM does not. The pre-input BAM
preparation first slices the raw BAM by a merged SD BED, samples qnames from
that slice, and then extracts those qnames. Therefore:

* raw ALT loci covered by the slicing BED point to subsampling/qname loss;
* raw ALT loci outside the slicing BED point to slicing loss.

For slicing-loss loci, the script also checks whether the locus existed in
earlier WGAC/SD pair tables, which helps localize loss to CIGAR fragmentation,
trim/filter/highsim, or final merged BED generation.
"""

from __future__ import annotations

import bisect
import csv
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
FN_TABLE = ROOT / "docs/analysis/fn_sites/full_fn_categories_updated_after_cat1c_retry.tsv"
OUT_SITE = ROOT / "docs/analysis/fn_sites/cat1_slicing_vs_subsampling.tsv"
OUT_LOCUS = ROOT / "docs/analysis/fn_sites/cat1_raw_locus_slicing_detail.tsv"

PUBLIC_SD = Path("/paedyl01/disk1/yangyxt/public_data/SD_from_SEDEF/hg19/test_BISER")

ASSEMBLY_PATHS = {
    "hg19": {
        "raw_wgac_pair": PUBLIC_SD / "WGAC.hg19.bed",
        "cigar_pair": PUBLIC_SD / "WGAC.hg19.cigar.bed",
        "cigar_trim_pair": PUBLIC_SD / "WGAC.hg19.cigar.trim.bed",
        "cigar_trimmed_pair": PUBLIC_SD / "WGAC.hg19.cigar.trimmed.bed",
        "cigar_trim_homo_expanded_highsim_pair": PUBLIC_SD / "WGAC.hg19.cigar.trim.homo.expanded.highsim.bed",
        "cigar_trimmed_homo_expanded_highsim_pair": PUBLIC_SD / "WGAC.hg19.cigar.trimmed.homo.expanded.highsim.bed",
        # prepare_GIAB_bam_file.sh passes the chr-prefixed BED as the hg19
        # default, then strips chr before sambamba slice because the raw BAM is
        # hs37d5.300x.bam. This stripchr file is interval-identical after
        # chromosome-name normalization.
        "final_merged_slicing_bed": PUBLIC_SD / "WGAC.hg19.cigar.trimmed.homo.expanded.merged.stripchr.bed",
    },
    "hg38": {
        "raw_wgac_pair": PUBLIC_SD / "WGAC.hg38.bed",
        "cigar_pair": PUBLIC_SD / "WGAC.hg38.cigar.bed",
        "cigar_trim_pair": PUBLIC_SD / "WGAC.hg38.cigar.trim.bed",
        "cigar_trimmed_pair": PUBLIC_SD / "WGAC.hg38.cigar.trimmed.bed",
        "cigar_trim_homo_expanded_highsim_pair": PUBLIC_SD / "WGAC.hg38.cigar.trim.homo.expanded.highsim.bed",
        "cigar_trimmed_homo_expanded_highsim_pair": PUBLIC_SD / "WGAC.hg38.cigar.trimmed.homo.expanded.highsim.bed",
        "final_merged_slicing_bed": PUBLIC_SD / "WGAC.hg38.cigar.trimmed.homo.expanded.merged.bed",
    },
}

PAIR_STAGES = [
    "raw_wgac_pair",
    "cigar_pair",
    "cigar_trim_pair",
    "cigar_trimmed_pair",
    "cigar_trim_homo_expanded_highsim_pair",
    "cigar_trimmed_homo_expanded_highsim_pair",
]
ALL_STAGES = PAIR_STAGES + ["final_merged_slicing_bed"]


@dataclass(frozen=True)
class RawLocus:
    qname: str
    chrom: str
    pos_1based: int
    mapq: str

    @property
    def pos0(self) -> int:
        return self.pos_1based - 1

    @property
    def label(self) -> str:
        return f"{self.chrom}:{self.pos_1based}"


class IntervalIndex:
    def __init__(self) -> None:
        self._intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
        self._starts: dict[str, list[int]] = {}
        self._ends: dict[str, list[int]] = {}

    def add(self, chrom: str, start: int, end: int) -> None:
        if end <= start:
            return
        self._intervals[normalize_chrom(chrom)].append((start, end))

    def finalize(self) -> None:
        for chrom, intervals in self._intervals.items():
            intervals.sort()
            merged: list[tuple[int, int]] = []
            for start, end in intervals:
                if not merged or start > merged[-1][1]:
                    merged.append((start, end))
                elif end > merged[-1][1]:
                    merged[-1] = (merged[-1][0], end)
            self._intervals[chrom] = merged
            self._starts[chrom] = [x[0] for x in merged]
            self._ends[chrom] = [x[1] for x in merged]

    def contains_point(self, chrom: str, pos0: int) -> bool:
        chrom = normalize_chrom(chrom)
        starts = self._starts.get(chrom)
        if not starts:
            return False
        idx = bisect.bisect_right(starts, pos0) - 1
        return idx >= 0 and self._intervals[chrom][idx][0] <= pos0 < self._intervals[chrom][idx][1]

    def overlaps_window(self, chrom: str, start0: int, end0: int) -> bool:
        chrom = normalize_chrom(chrom)
        starts = self._starts.get(chrom)
        if not starts:
            return False
        idx = bisect.bisect_left(starts, end0)
        if idx > 0 and self._intervals[chrom][idx - 1][1] > start0:
            return True
        return idx < len(starts) and self._intervals[chrom][idx][0] < end0

    def nearest(self, chrom: str, pos0: int) -> tuple[str, int | str]:
        chrom = normalize_chrom(chrom)
        intervals = self._intervals.get(chrom)
        starts = self._starts.get(chrom)
        if not intervals or not starts:
            return "", "NA"
        idx = bisect.bisect_right(starts, pos0) - 1
        candidates: list[tuple[int, tuple[int, int]]] = []
        if idx >= 0:
            start, end = intervals[idx]
            dist = 0 if start <= pos0 < end else pos0 - end + 1
            candidates.append((dist, (start, end)))
        if idx + 1 < len(intervals):
            start, end = intervals[idx + 1]
            dist = 0 if start <= pos0 < end else start - pos0
            candidates.append((dist, (start, end)))
        dist, (start, end) = min(candidates, key=lambda item: item[0])
        return f"{chrom}:{start}-{end}", dist


def normalize_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.startswith("chr"):
        return chrom
    return f"chr{chrom}"


RAW_LOCUS_RE = re.compile(r"(?P<qname>.+)@(?P<chrom>[^:]+):(?P<pos>\d+):MAPQ(?P<mapq>[^;]+)$")


def parse_raw_loci(value: str) -> list[RawLocus]:
    if not value or value == ".":
        return []
    loci: list[RawLocus] = []
    for item in value.split(";"):
        item = item.strip()
        if not item:
            continue
        match = RAW_LOCUS_RE.match(item)
        if not match:
            raise ValueError(f"Cannot parse raw_loci entry: {item}")
        loci.append(
            RawLocus(
                qname=match.group("qname"),
                chrom=match.group("chrom"),
                pos_1based=int(match.group("pos")),
                mapq=match.group("mapq"),
            )
        )
    return loci


def load_bed_index(path: Path, both_pair_sides: bool) -> IntervalIndex:
    idx = IntervalIndex()
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            idx.add(fields[0], int(fields[1]), int(fields[2]))
            if both_pair_sides and len(fields) >= 6:
                idx.add(fields[3], int(fields[4]), int(fields[5]))
    idx.finalize()
    return idx


def cluster_loci(loci: Iterable[RawLocus], max_gap: int = 1000) -> str:
    grouped: dict[str, list[int]] = defaultdict(list)
    for locus in loci:
        grouped[normalize_chrom(locus.chrom)].append(locus.pos_1based)

    pieces: list[str] = []
    for chrom in sorted(grouped):
        positions = sorted(grouped[chrom])
        start = prev = positions[0]
        count = 1
        for pos in positions[1:]:
            if pos - prev <= max_gap:
                prev = pos
                count += 1
            else:
                pieces.append(f"{chrom}:{start}-{prev}:n={count}")
                start = prev = pos
                count = 1
        pieces.append(f"{chrom}:{start}-{prev}:n={count}")
    return ";".join(pieces)


def classify_cat1(point_covered: int, total: int) -> str:
    if total == 0:
        return "Cat1_no_raw_locus"
    if point_covered == total:
        return "Cat1_subsampling_loss_likely"
    if point_covered == 0:
        return "Cat1_slicing_loss_likely"
    return "Cat1_mixed_slicing_and_subsampling_likely"


def stage_inference(stage_counts: dict[str, int], total: int, slicing_covered: int) -> str:
    if total == 0:
        return "no_raw_locus"
    if slicing_covered == total:
        return "slicing_bed_contains_all_raw_alt_loci; input loss is consistent with subsampling/qname selection"
    if slicing_covered > 0:
        return "mixed: some raw ALT loci are in the slicing BED and some are outside it"

    if stage_counts["raw_wgac_pair"] == 0:
        return "slicing_loss; raw ALT loci are absent even from original WGAC pair intervals"
    if stage_counts["cigar_pair"] == 0:
        return "slicing_loss; original WGAC pair intervals contain the loci but CIGAR-fragmented intervals do not"

    previous = "cigar_pair"
    for stage in [
        "cigar_trim_pair",
        "cigar_trimmed_pair",
        "cigar_trim_homo_expanded_highsim_pair",
        "cigar_trimmed_homo_expanded_highsim_pair",
        "final_merged_slicing_bed",
    ]:
        if stage_counts[stage] == 0:
            return f"slicing_loss; loci present through {previous} but absent from {stage}"
        previous = stage
    return "slicing_loss; pair-level evidence exists but final slicing coverage check failed unexpectedly"


def main() -> None:
    stage_indexes: dict[str, dict[str, IntervalIndex]] = {}
    for assembly, paths in ASSEMBLY_PATHS.items():
        stage_indexes[assembly] = {}
        for stage, path in paths.items():
            both_pair_sides = stage in PAIR_STAGES
            stage_indexes[assembly][stage] = load_bed_index(path, both_pair_sides=both_pair_sides)

    site_rows: list[dict[str, str | int | float]] = []
    locus_rows: list[dict[str, str | int]] = []

    with FN_TABLE.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            assembly = row["assembly"]
            if assembly not in {"hg19", "hg38"}:
                continue
            if row["final_category"] != "Cat1_raw_only_input_loss":
                continue

            loci = parse_raw_loci(row["raw_loci"])
            total = len(loci)
            stage_counts: dict[str, int] = {}
            final_idx = stage_indexes[assembly]["final_merged_slicing_bed"]
            point_covered = 0
            padded_covered = 0

            for stage in ALL_STAGES:
                idx = stage_indexes[assembly][stage]
                stage_counts[stage] = sum(1 for locus in loci if idx.contains_point(locus.chrom, locus.pos0))

            for locus in loci:
                point_in_final = final_idx.contains_point(locus.chrom, locus.pos0)
                padded_in_final = final_idx.overlaps_window(locus.chrom, max(0, locus.pos0 - 250), locus.pos0 + 251)
                if point_in_final:
                    point_covered += 1
                if padded_in_final:
                    padded_covered += 1
                nearest_interval, nearest_distance = final_idx.nearest(locus.chrom, locus.pos0)
                locus_row: dict[str, str | int] = {
                    "assembly": assembly,
                    "FN_site": row["FN_site"],
                    "qname": locus.qname,
                    "raw_chrom": normalize_chrom(locus.chrom),
                    "raw_pos_1based": locus.pos_1based,
                    "raw_mapq": locus.mapq,
                    "point_in_final_slicing_bed": int(point_in_final),
                    "window_250bp_overlaps_final_slicing_bed": int(padded_in_final),
                    "nearest_final_slicing_interval": nearest_interval,
                    "nearest_final_slicing_distance_bp": nearest_distance,
                }
                for stage in ALL_STAGES:
                    locus_row[stage] = int(stage_indexes[assembly][stage].contains_point(locus.chrom, locus.pos0))
                locus_rows.append(locus_row)

            cat1_subcause = classify_cat1(point_covered, total)
            stage_summary = ";".join(f"{stage}={stage_counts[stage]}/{total}" for stage in ALL_STAGES)
            site_rows.append(
                {
                    "assembly": assembly,
                    "FN_site": row["FN_site"],
                    "raw_reads": row["raw_reads"],
                    "input_reads": row["input_reads"],
                    "raw_qnames": len({locus.qname for locus in loci}),
                    "raw_loci": total,
                    "raw_locus_clusters_1kb": cluster_loci(loci),
                    "final_slicing_point_covered_loci": point_covered,
                    "final_slicing_point_uncovered_loci": total - point_covered,
                    "final_slicing_point_coverage_fraction": f"{point_covered / total:.4f}" if total else "NA",
                    "final_slicing_250bp_window_covered_loci": padded_covered,
                    "final_slicing_250bp_window_coverage_fraction": f"{padded_covered / total:.4f}" if total else "NA",
                    "cat1_subcause": cat1_subcause,
                    "stage_covered_counts": stage_summary,
                    "stage_loss_inference": stage_inference(stage_counts, total, point_covered),
                    "fc_groups": row["fc_groups"],
                }
            )

    site_fields = [
        "assembly",
        "FN_site",
        "raw_reads",
        "input_reads",
        "raw_qnames",
        "raw_loci",
        "raw_locus_clusters_1kb",
        "final_slicing_point_covered_loci",
        "final_slicing_point_uncovered_loci",
        "final_slicing_point_coverage_fraction",
        "final_slicing_250bp_window_covered_loci",
        "final_slicing_250bp_window_coverage_fraction",
        "cat1_subcause",
        "stage_covered_counts",
        "stage_loss_inference",
        "fc_groups",
    ]
    locus_fields = [
        "assembly",
        "FN_site",
        "qname",
        "raw_chrom",
        "raw_pos_1based",
        "raw_mapq",
        "point_in_final_slicing_bed",
        "window_250bp_overlaps_final_slicing_bed",
        "nearest_final_slicing_interval",
        "nearest_final_slicing_distance_bp",
        *ALL_STAGES,
    ]

    with OUT_SITE.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=site_fields)
        writer.writeheader()
        writer.writerows(site_rows)

    with OUT_LOCUS.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=locus_fields)
        writer.writeheader()
        writer.writerows(locus_rows)

    print(f"Wrote {OUT_SITE} ({len(site_rows)} Cat1 sites)")
    print(f"Wrote {OUT_LOCUS} ({len(locus_rows)} raw ALT loci)")
    print("Cat1 subcause counts:")
    for (assembly, subcause), count in sorted(Counter((r["assembly"], r["cat1_subcause"]) for r in site_rows).items()):
        print(f"{assembly}\t{subcause}\t{count}")


if __name__ == "__main__":
    main()
