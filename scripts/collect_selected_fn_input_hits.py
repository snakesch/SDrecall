#!/usr/bin/env python3
"""Collect scanner hits from the tier that supplied each selected FN marker."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected-definitions", required=True, type=Path)
    parser.add_argument("--definitions", action="append", required=True, type=Path)
    parser.add_argument("--input-hits", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or ()), list(reader)


def main() -> None:
    args = parse_args()
    if len(args.definitions) != len(args.input_hits):
        raise ValueError("Each --definitions file needs one --input-hits file")

    _, selected_rows = read_rows(args.selected_definitions)
    selected_sequences: dict[str, set[str]] = {}
    for row in selected_rows:
        selected_sequences.setdefault(row["FN_site"], set()).add(row["kmer_seq"])

    output_fields: list[str] | None = None
    output_rows: list[dict[str, str]] = []
    matched_sites: set[str] = set()
    seen_rows: set[tuple[str, ...]] = set()

    for definitions_path, hits_path in zip(
        args.definitions, args.input_hits, strict=True
    ):
        _, definition_rows = read_rows(definitions_path)
        tier_sites = {
            row["FN_site"]
            for row in definition_rows
            if row["kmer_seq"] in selected_sequences.get(row["FN_site"], set())
        }
        matched_sites.update(tier_sites)

        hit_fields, hit_rows = read_rows(hits_path)
        if output_fields is None:
            output_fields = hit_fields
        elif hit_fields != output_fields:
            raise ValueError("Scanner hit headers differ between marker tiers")
        for row in hit_rows:
            if row["FN_site"] not in tier_sites:
                continue
            row_key = tuple(row[field] for field in hit_fields)
            if row_key in seen_rows:
                continue
            seen_rows.add(row_key)
            output_rows.append(row)

    missing_tiers = sorted(set(selected_sequences) - matched_sites)
    if missing_tiers:
        raise ValueError(f"Selected markers not found in source tiers: {missing_tiers}")

    assert output_fields is not None
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)


if __name__ == "__main__":
    main()
