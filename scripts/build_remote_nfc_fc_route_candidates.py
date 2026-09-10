#!/usr/bin/env python3
"""Resolve the actual remote NFC interval behind each Cat4 qname pathway.

The production FASTQ records only qnames.  ``cat4_qname_extraction_paths.tsv``
records the eligible NFC owner label but not which one of that owner's remote
NFC intervals overlapped the qname's original input-BAM alignment.  This script
joins those two pieces so graph auditing uses:

    exact remote NFC interval -> destination FC reached after realignment

rather than incorrectly substituting the owner's local FC interval.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path("/paedyl01/disk1/yangyxt/SDrecall-rust-migration")
sys.path.insert(0, str(ROOT / "scripts"))

from update_fn_qname_extraction_audit import (  # noqa: E402
    INPUT_SAMS,
    Interval,
    canonical_qname,
    parse_sam,
)


AUDIT_ROOT = ROOT / "test_tmp/nfc_fc_conflict_fix_20260723"
RUN_ROOTS = {
    "hg19": AUDIT_ROOT / "runs/hg19/HG002_hg19_nfcfix23_SDrecall",
    "hg38": AUDIT_ROOT / "runs/hg38/HG002_hg38_nfcfix23_SDrecall",
    "chm13": AUDIT_ROOT / "runs/t2t/HG002_chm13_nfcfix23_SDrecall",
}
DEFAULT_PATHS = AUDIT_ROOT / "cat4_qname_audit/cat4_qname_extraction_paths.tsv"


@dataclass(frozen=True)
class TaggedRegion:
    kind: str
    rg: str
    label: str
    interval: Interval
    strand: str


@dataclass(frozen=True)
class CandidateKey:
    assembly: str
    allele: str
    rg: str
    owner_fc: str
    owner_interval: TaggedRegion
    source_nfc: str
    source_interval: TaggedRegion
    destination_fc: str
    destination_interval: TaggedRegion


@dataclass
class Candidate:
    key: CandidateKey
    qnames: set[str] = field(default_factory=set)
    extraction_sources: set[str] = field(default_factory=set)
    path_provenance: set[str] = field(default_factory=set)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paths", type=Path, default=DEFAULT_PATHS)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def read_regions(assembly: str, rg: str) -> list[TaggedRegion]:
    path = RUN_ROOTS[assembly] / "realign_groups" / rg / f"{rg}_related_homo_regions.bed"
    if not path.is_file():
        raise FileNotFoundError(path)
    regions: list[TaggedRegion] = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7 or ":" not in fields[6]:
                continue
            kind, pair = fields[6].split(":", 1)
            if kind not in {"FC", "NFC"}:
                continue
            regions.append(
                TaggedRegion(
                    kind=kind,
                    rg=pair.split("_", 1)[0],
                    label=fields[6],
                    interval=Interval(fields[0], int(fields[1]), int(fields[2])),
                    strand=fields[5],
                )
            )
    return regions


def display(region: TaggedRegion) -> str:
    interval = region.interval
    return f"{interval.chrom}:{interval.start}-{interval.end}({region.strand})"


def distance_class(left_region: TaggedRegion, right_region: TaggedRegion) -> str:
    left = left_region.interval
    right = right_region.interval
    if left.chrom != right.chrom:
        return "different_chromosome"
    overlap = min(left.end, right.end) - max(left.start, right.start)
    if overlap > 0:
        return f"overlap:{overlap}"
    gap = right.start - left.end if left.end <= right.start else left.start - right.end
    return f"gap:{gap}"


def main() -> None:
    args = parse_args()
    with args.paths.open() as handle:
        path_rows = list(csv.DictReader(handle, delimiter="\t"))

    needed_assemblies = sorted(
        {
            row["Assembly"]
            for row in path_rows
            if row["Extraction source"].startswith("NFC:")
            and row["Path outcome"] == "other FC(s) only"
        }
    )
    input_records = {
        assembly: parse_sam(INPUT_SAMS[assembly]) for assembly in needed_assemblies
    }
    region_cache: dict[tuple[str, str], list[PairRegion]] = {}
    candidates: dict[CandidateKey, Candidate] = {}
    unresolved: list[tuple[str, str, str, str]] = []

    for row in path_rows:
        source_label = row["Extraction source"]
        if not source_label.startswith("NFC:") or row["Path outcome"] != "other FC(s) only":
            continue
        assembly = row["Assembly"]
        allele = row["FN allele"]
        rg = row["RG"]
        qname = canonical_qname(row["ALT-haplotype qname"])
        regions = region_cache.setdefault((assembly, rg), read_regions(assembly, rg))
        source_regions = [
            region
            for region in regions
            if region.kind == "NFC" and region.label == source_label
        ]
        owner_label = source_label.replace("NFC:", "FC:", 1)
        owner_regions = [
            region
            for region in regions
            if region.kind == "FC" and region.label == owner_label
        ]
        destination_labels = [
            value for value in row["Other FC labels"].split(";") if value
        ]
        destination_regions = {
            label: [
                region
                for region in regions
                if region.kind == "FC" and region.label == label
            ]
            for label in destination_labels
        }

        records = input_records[assembly].get(qname, [])
        record_intervals = [
            record.interval for record in records if record.interval is not None
        ]
        matched_sources = [
            region
            for region in source_regions
            if any(
                interval.overlaps(region.interval)
                for interval in record_intervals
            )
        ]
        if not matched_sources:
            unresolved.append((assembly, allele, qname, source_label))
            continue

        for owner in owner_regions:
            for source in matched_sources:
                for destination_label, destinations in destination_regions.items():
                    for destination in destinations:
                        key = CandidateKey(
                            assembly=assembly,
                            allele=allele,
                            rg=rg,
                            owner_fc=owner_label,
                            owner_interval=owner,
                            source_nfc=source_label,
                            source_interval=source,
                            destination_fc=destination_label,
                            destination_interval=destination,
                        )
                        candidate = candidates.setdefault(key, Candidate(key))
                        candidate.qnames.add(qname)
                        candidate.extraction_sources.add(source_label)
                        candidate.path_provenance.add(row["Path provenance"])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "Assembly",
        "FN allele",
        "RG",
        "Owner FC",
        "Owner FC interval",
        "Exact remote NFC",
        "Exact remote NFC interval",
        "Destination FC",
        "Destination FC interval",
        "Remote NFC to destination relationship",
        "Supporting qname count",
        "Supporting qnames",
        "Path provenance",
    ]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n"
        )
        writer.writeheader()
        for key in sorted(
            candidates,
            key=lambda item: (
                item.assembly,
                item.allele,
                item.rg,
                item.owner_fc,
                item.source_interval.interval.chrom,
                item.source_interval.interval.start,
                item.destination_fc,
            ),
        ):
            candidate = candidates[key]
            writer.writerow(
                {
                    "Assembly": key.assembly,
                    "FN allele": key.allele,
                    "RG": key.rg,
                    "Owner FC": key.owner_fc,
                    "Owner FC interval": display(key.owner_interval),
                    "Exact remote NFC": key.source_nfc,
                    "Exact remote NFC interval": display(key.source_interval),
                    "Destination FC": key.destination_fc,
                    "Destination FC interval": display(key.destination_interval),
                    "Remote NFC to destination relationship": distance_class(
                        key.source_interval, key.destination_interval
                    ),
                    "Supporting qname count": len(candidate.qnames),
                    "Supporting qnames": ";".join(sorted(candidate.qnames)),
                    "Path provenance": ";".join(sorted(candidate.path_provenance)),
                }
            )

    print(f"candidate_rows={len(candidates)}")
    print(f"unresolved_paths={len(unresolved)}")
    for item in unresolved[:20]:
        print("unresolved\t" + "\t".join(item))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
