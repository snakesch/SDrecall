#!/usr/bin/env python3
"""Find qnames whose realigned records overlap two genomically remote FCs."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path


ROOT = Path("/paedyl01/disk1/yangyxt/SDrecall-rust-migration")
AUDIT_ROOT = ROOT / "test_tmp/nfc_fc_conflict_fix_20260723"
RUN_ROOTS = {
    "hg19": AUDIT_ROOT / "runs/hg19/HG002_hg19_nfcfix23_SDrecall",
    "hg38": AUDIT_ROOT / "runs/hg38/HG002_hg38_nfcfix23_SDrecall",
    "chm13": AUDIT_ROOT / "runs/t2t/HG002_chm13_nfcfix23_SDrecall",
}
TRACE_PATHS = {
    assembly: AUDIT_ROOT / "haplotypes" / assembly / "pre_fp_input_qname_trace.tsv"
    for assembly in RUN_ROOTS
}
DEFAULT_QNAMES = AUDIT_ROOT / "cat4_qname_audit/cat4_qname_summary.tsv"
RG_SUFFIX = re.compile(r":RG\d+$")


@dataclass(frozen=True)
class Interval:
    chrom: str
    start: int
    end: int
    strand: str = "."

    def overlaps(self, other: "Interval") -> bool:
        return (
            self.chrom == other.chrom
            and self.start < other.end
            and other.start < self.end
        )

    def remote_from(self, other: "Interval", threshold: float) -> bool:
        if self.chrom != other.chrom:
            return True
        if self.overlaps(other):
            return False
        gap = (
            other.start - self.end
            if self.end <= other.start
            else self.start - other.end
        )
        return gap > threshold

    def display(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}({self.strand})"


@dataclass(frozen=True)
class Placement:
    record_id: str
    interval: Interval
    mapq: int
    cigar: str
    marker: bool

    def display(self) -> str:
        marker = "marker" if self.marker else "mate"
        return (
            f"{self.record_id}:{self.interval.display()}:"
            f"MQ{self.mapq}:{self.cigar}:{marker}"
        )


@dataclass(frozen=True)
class PairKey:
    assembly: str
    allele: str
    rg: str
    fc_a: str
    interval_a: Interval
    fc_b: str
    interval_b: Interval


@dataclass
class PairEvidence:
    key: PairKey
    qnames: set[str] = field(default_factory=set)
    placements_a: set[str] = field(default_factory=set)
    placements_b: set[str] = field(default_factory=set)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qnames", type=Path, default=DEFAULT_QNAMES)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--threshold",
        action="append",
        default=[],
        metavar="ASSEMBLY=BP",
        help="remote-gap threshold; defaults to production median fragment size",
    )
    return parser.parse_args()


def canonical_qname(value: str) -> str:
    return RG_SUFFIX.sub("", value.strip())


def thresholds(assignments: list[str]) -> dict[str, float]:
    values = {"hg19": 572.35, "hg38": 572.35, "chm13": 563.30}
    for assignment in assignments:
        assembly, separator, number = assignment.partition("=")
        if not separator or assembly not in values:
            raise ValueError(f"invalid --threshold {assignment!r}")
        values[assembly] = float(number)
    return values


def read_fc_catalog(assembly: str, rgs: set[str]) -> dict[str, list[Interval]]:
    catalog: dict[str, list[Interval]] = defaultdict(list)
    for rg in sorted(rgs, key=lambda value: int(value[2:])):
        path = (
            RUN_ROOTS[assembly]
            / "realign_groups"
            / rg
            / f"{rg}_related_homo_regions.bed"
        )
        with path.open() as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 7 or not fields[6].startswith("FC:"):
                    continue
                catalog[fields[6]].append(
                    Interval(
                        fields[0],
                        int(fields[1]),
                        int(fields[2]),
                        fields[5],
                    )
                )
    return catalog


def read_needed_qnames(
    path: Path,
) -> tuple[set[tuple[str, str, str]], dict[str, set[str]]]:
    needed: set[tuple[str, str, str]] = set()
    rgs: dict[str, set[str]] = defaultdict(set)
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            assembly = row["Assembly"]
            allele = row["FN allele"]
            qname = canonical_qname(row["ALT-haplotype qname"])
            needed.add((assembly, allele, qname))
            for value in row["Realigned pair records"].split(";"):
                rg = value.split(":", 1)[0]
                if re.fullmatch(r"RG\d+", rg):
                    rgs[assembly].add(rg)
    return needed, rgs


def read_trace(
    assembly: str,
    needed: set[tuple[str, str, str]],
) -> dict[tuple[str, str, str, str], list[Placement]]:
    rows: dict[tuple[str, str, str, str], list[Placement]] = defaultdict(list)
    with TRACE_PATHS[assembly].open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key3 = (
                assembly,
                row["FN_site"],
                canonical_qname(row["read_qname"]),
            )
            if key3 not in needed:
                continue
            rows[(*key3, row["RG"])].append(
                Placement(
                    record_id=row["read_record_id"],
                    interval=Interval(
                        row["mapping_chrom"],
                        int(row["mapping_start"]) - 1,
                        int(row["mapping_end"]),
                    ),
                    mapq=int(row["MAPQ"]),
                    cigar=row["CIGAR"],
                    marker=row["exact_marker_in_record"] == "yes",
                )
            )
    return rows


def ordered_pair(
    left_label: str,
    left_interval: Interval,
    right_label: str,
    right_interval: Interval,
) -> tuple[str, Interval, str, Interval, bool]:
    left_key = (left_interval.chrom, left_interval.start, left_interval.end, left_label)
    right_key = (
        right_interval.chrom,
        right_interval.start,
        right_interval.end,
        right_label,
    )
    if left_key <= right_key:
        return left_label, left_interval, right_label, right_interval, False
    return right_label, right_interval, left_label, left_interval, True


def main() -> None:
    args = parse_args()
    remote_thresholds = thresholds(args.threshold)
    needed, rgs = read_needed_qnames(args.qnames)
    evidence: dict[PairKey, PairEvidence] = {}

    for assembly in sorted(RUN_ROOTS):
        assembly_needed = {key for key in needed if key[0] == assembly}
        catalog = read_fc_catalog(assembly, rgs[assembly])
        trace = read_trace(assembly, assembly_needed)
        for (observed_assembly, allele, qname, rg), placements in trace.items():
            assigned: list[tuple[Placement, str, Interval]] = []
            for placement in placements:
                for label, intervals in catalog.items():
                    if not label.startswith(f"FC:{rg}_"):
                        continue
                    for interval in intervals:
                        if placement.interval.overlaps(interval):
                            assigned.append((placement, label, interval))

            for left_index, (left, left_label, left_fc) in enumerate(assigned):
                for right, right_label, right_fc in assigned[left_index + 1 :]:
                    if left.record_id == right.record_id and left.interval == right.interval:
                        continue
                    if left_label == right_label:
                        continue
                    if not left_fc.remote_from(
                        right_fc, remote_thresholds[observed_assembly]
                    ):
                        continue
                    fc_a, interval_a, fc_b, interval_b, swapped = ordered_pair(
                        left_label, left_fc, right_label, right_fc
                    )
                    key = PairKey(
                        observed_assembly,
                        allele,
                        rg,
                        fc_a,
                        interval_a,
                        fc_b,
                        interval_b,
                    )
                    item = evidence.setdefault(key, PairEvidence(key))
                    item.qnames.add(qname)
                    if swapped:
                        item.placements_a.add(right.display())
                        item.placements_b.add(left.display())
                    else:
                        item.placements_a.add(left.display())
                        item.placements_b.add(right.display())

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "Assembly",
        "FN allele",
        "RG",
        "FC A",
        "FC A interval",
        "FC B",
        "FC B interval",
        "Supporting qname count",
        "Supporting qnames",
        "FC A record placements",
        "FC B record placements",
    ]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n"
        )
        writer.writeheader()
        for key in sorted(
            evidence,
            key=lambda item: (
                item.assembly,
                item.allele,
                item.rg,
                item.interval_a.chrom,
                item.interval_a.start,
                item.interval_b.chrom,
                item.interval_b.start,
            ),
        ):
            item = evidence[key]
            writer.writerow(
                {
                    "Assembly": key.assembly,
                    "FN allele": key.allele,
                    "RG": key.rg,
                    "FC A": key.fc_a,
                    "FC A interval": key.interval_a.display(),
                    "FC B": key.fc_b,
                    "FC B interval": key.interval_b.display(),
                    "Supporting qname count": len(item.qnames),
                    "Supporting qnames": ";".join(sorted(item.qnames)),
                    "FC A record placements": ";".join(sorted(item.placements_a)),
                    "FC B record placements": ";".join(sorted(item.placements_b)),
                }
            )

    by_assembly: dict[str, int] = defaultdict(int)
    for key in evidence:
        by_assembly[key.assembly] += 1
    for assembly in sorted(RUN_ROOTS):
        print(f"{assembly}_remote_fc_pairs={by_assembly[assembly]}")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
