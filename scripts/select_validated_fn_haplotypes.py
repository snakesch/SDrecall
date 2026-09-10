#!/usr/bin/env python3
"""Select the longest target-validated marker for each FN truth haplotype."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--definitions", action="append", required=True, type=Path)
    parser.add_argument("--validation-summary", action="append", required=True, type=Path)
    parser.add_argument("--input-hits", action="append", type=Path)
    parser.add_argument("--output-definitions", required=True, type=Path)
    parser.add_argument("--output-kmers", required=True, type=Path)
    parser.add_argument("--output-selection", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if len(args.definitions) != len(args.validation_summary):
        raise ValueError("Each --definitions file needs one --validation-summary file")
    if args.input_hits and len(args.input_hits) != len(args.definitions):
        raise ValueError("Each --definitions file needs one --input-hits file")

    validated_counts: dict[tuple[str, str, int], tuple[int, int, str]] = {}
    input_counts: dict[tuple[str, int], int] = {}
    definition_rows: list[dict[str, str]] = []
    definition_fields: list[str] | None = None

    input_paths = args.input_hits or [None] * len(args.definitions)
    for definitions_path, summary_path, input_path in zip(
        args.definitions, args.validation_summary, input_paths, strict=True
    ):
        with summary_path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                key = (row["FN_site"], row["haplotype"], int(row["marker_len"]))
                validated_counts[key] = (
                    int(row["target_haplotype_qnames"]),
                    int(row["target_haplotype_records"]),
                    row["source"],
                )

        tier_input_qnames: dict[str, set[str]] = {}
        if input_path is not None:
            with input_path.open() as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    tier_input_qnames.setdefault(row["FN_site"], set()).add(
                        row["read_qname"]
                    )

        with definitions_path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if definition_fields is None:
                definition_fields = list(reader.fieldnames or ())
            elif list(reader.fieldnames or ()) != definition_fields:
                raise ValueError("Definition headers differ between marker lengths")
            for row in reader:
                definition_rows.append(row)
                input_counts[(row["FN_site"], int(row["kmer_len"]))] = len(
                    tier_input_qnames.get(row["FN_site"], set())
                )

    selected: dict[tuple[str, str], dict[str, str]] = {}
    for row in definition_rows:
        marker_key = (row["FN_site"], row["haplotype"], int(row["kmer_len"]))
        qnames, _, _ = validated_counts.get(marker_key, (0, 0, "."))
        if qnames == 0:
            continue
        haplotype_key = (row["FN_site"], row["haplotype"])
        previous = selected.get(haplotype_key)
        candidate_score = (
            input_counts.get((row["FN_site"], int(row["kmer_len"])), 0) > 0,
            int(row["kmer_len"]),
        )
        previous_score = (
            input_counts.get(
                (previous["FN_site"], int(previous["kmer_len"])), 0
            )
            > 0,
            int(previous["kmer_len"]),
        ) if previous is not None else None
        if previous_score is None or candidate_score > previous_score:
            selected[haplotype_key] = row

    expected_haplotypes = {
        (row["FN_site"], row["haplotype"]) for row in definition_rows
    }
    missing = sorted(expected_haplotypes - set(selected))
    if missing:
        raise ValueError(f"No target-validated marker for haplotypes: {missing}")

    selected_rows = [selected[key] for key in sorted(selected)]
    args.output_definitions.parent.mkdir(parents=True, exist_ok=True)
    args.output_kmers.parent.mkdir(parents=True, exist_ok=True)
    args.output_selection.parent.mkdir(parents=True, exist_ok=True)

    assert definition_fields is not None
    with args.output_definitions.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=definition_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(selected_rows)

    with args.output_kmers.open("w", newline="") as handle:
        fields = ["assembly", "FN_site", "haplotype", "kmer_seq", "kmer_len"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows({field: row[field] for field in fields} for row in selected_rows)

    selection_fields = [
        "FN_site",
        "haplotype",
        "selected_marker_len",
        "validation_source",
        "target_haplotype_qnames",
        "target_haplotype_records",
        "selected_input_haplotype_qnames",
        "selection_reason",
    ]
    with args.output_selection.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=selection_fields, delimiter="\t")
        writer.writeheader()
        for row in selected_rows:
            marker_key = (row["FN_site"], row["haplotype"], int(row["kmer_len"]))
            qnames, records, source = validated_counts[marker_key]
            input_qnames = input_counts.get(
                (row["FN_site"], int(row["kmer_len"])), 0
            )
            writer.writerow(
                {
                    "FN_site": row["FN_site"],
                    "haplotype": row["haplotype"],
                    "selected_marker_len": row["kmer_len"],
                    "validation_source": source,
                    "target_haplotype_qnames": qnames,
                    "target_haplotype_records": records,
                    "selected_input_haplotype_qnames": input_qnames,
                    "selection_reason": (
                        "longest_hifi_and_input_validated_marker"
                        if input_qnames > 0
                        else "longest_hifi_validated_marker_no_input_hit"
                    ),
                }
            )


if __name__ == "__main__":
    main()
