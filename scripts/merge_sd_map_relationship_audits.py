#!/usr/bin/env python3
"""Merge the required hg19, hg38, and CHM13 relationship-audit tables."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_ROOT = (
    ROOT / "test_tmp/sd_map_relationship_audit_20260724/by_assembly"
)
ASSEMBLIES = ("hg19", "hg38", "chm13")
TABLES = (
    "assembly_summary.tsv",
    "relationship_candidates.tsv",
    "relationship_summary.tsv",
    "stage_audit.tsv",
    "stage_provenance.tsv",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--output-dir", type=Path)
    return parser.parse_args()


def merge_table(input_root: Path, output_dir: Path, filename: str) -> int:
    fieldnames: list[str] | None = None
    merged_rows: list[dict[str, str]] = []
    for assembly in ASSEMBLIES:
        path = input_root / assembly / filename
        if not path.is_file():
            raise FileNotFoundError(path)
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"missing header: {path}")
            if fieldnames is None:
                fieldnames = reader.fieldnames
            elif reader.fieldnames != fieldnames:
                raise ValueError(f"header mismatch: {path}")
            for row in reader:
                if row.get("Assembly") != assembly:
                    raise ValueError(
                        f"unexpected Assembly in {path}: {row.get('Assembly')!r}"
                    )
                merged_rows.append(row)

    if fieldnames is None:
        raise AssertionError(f"no header found for {filename}")
    with (output_dir / filename).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(merged_rows)
    return len(merged_rows)


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir or args.input_root.parent / "combined"
    output_dir.mkdir(parents=True, exist_ok=True)
    for filename in TABLES:
        rows = merge_table(args.input_root, output_dir, filename)
        print(f"{filename}\t{rows}\t{output_dir / filename}")


if __name__ == "__main__":
    main()
