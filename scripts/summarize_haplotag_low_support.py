#!/usr/bin/env python3
"""Summarize low-support WhatsHap haplotag blocks from large TSV files.

The WhatsHap haplotag list has one row per primary alignment, so paired-end
mates usually appear as two rows with the same read name. This script reports
candidate (chromosome, phaseset, haplotype) groups whose row count could
correspond to <=N read templates, then performs an exact unique-readname count
only for those candidates.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl


KEYS = ["chromosome", "phaseset", "haplotype"]
READNAME = "#readname"
ALL_COLUMNS = [READNAME, "haplotype", "phaseset", "chromosome"]


def collect_streaming(lazy_frame: pl.LazyFrame) -> pl.DataFrame:
    try:
        return lazy_frame.collect(engine="streaming")
    except TypeError:
        return lazy_frame.collect(streaming=True)


def scan_haplotags(path: str) -> pl.LazyFrame:
    return pl.scan_csv(
        path,
        separator="\t",
        has_header=True,
        schema_overrides={column: pl.String for column in ALL_COLUMNS},
        infer_schema_length=0,
    )


def summarize(path: str, max_pairs: int, output_prefix: str | None) -> None:
    max_candidate_rows = max_pairs * 2
    basename = Path(path).name

    print(f"\n== {basename} ==", flush=True)
    print(f"candidate row cutoff: <= {max_candidate_rows} rows", flush=True)

    base = scan_haplotags(path).select(ALL_COLUMNS)
    tagged = base.filter(
        (pl.col("haplotype") != "none")
        & (pl.col("phaseset") != "none")
        & pl.col("haplotype").is_not_null()
        & pl.col("phaseset").is_not_null()
    )

    group_rows_lf = tagged.group_by(KEYS).agg(pl.len().alias("alignment_rows"))
    print("scanning row counts per chromosome/phaseset/haplotype ...", flush=True)
    group_rows = collect_streaming(group_rows_lf)

    support_hist = (
        group_rows.with_columns(
            pl.when(pl.col("alignment_rows") <= max_candidate_rows)
            .then(pl.lit(f"<= {max_candidate_rows} rows"))
            .when(pl.col("alignment_rows") <= 20)
            .then(pl.lit("7-20 rows"))
            .when(pl.col("alignment_rows") <= 100)
            .then(pl.lit("21-100 rows"))
            .otherwise(pl.lit(">100 rows"))
            .alias("row_count_bin")
        )
        .group_by("row_count_bin")
        .agg(pl.len().alias("haplotype_groups"))
        .sort("row_count_bin")
    )

    candidates = group_rows.filter(pl.col("alignment_rows") <= max_candidate_rows).sort(KEYS)

    print("row-count distribution across tagged haplotype groups:")
    print(support_hist)
    print(
        f"candidate groups needing exact unique-readname count: {candidates.height}",
        flush=True,
    )

    if candidates.is_empty():
        print(f"exact groups with <= {max_pairs} unique read pairs/templates: 0")
        return

    print("scanning exact unique read names for candidate groups ...", flush=True)
    candidate_keys = candidates.select(KEYS)
    exact_lf = (
        tagged.join(candidate_keys.lazy(), on=KEYS, how="inner")
        .select([READNAME, *KEYS])
        .unique()
        .group_by(KEYS)
        .agg(pl.len().alias("unique_read_pairs"))
        .join(candidates.lazy(), on=KEYS, how="left")
        .filter(pl.col("unique_read_pairs") <= max_pairs)
        .sort(["unique_read_pairs", "alignment_rows", *KEYS])
    )
    exact = collect_streaming(exact_lf)

    print(f"exact groups with <= {max_pairs} unique read pairs/templates: {exact.height}")
    if exact.height:
        print(exact.head(50))

    if output_prefix:
        out = f"{output_prefix}.{Path(path).stem}.low_support.tsv"
        exact.write_csv(out, separator="\t")
        print(f"wrote {out}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("haplotag_tsv", nargs="+")
    parser.add_argument("--max-pairs", type=int, default=3)
    parser.add_argument("--output-prefix")
    args = parser.parse_args()

    for path in args.haplotag_tsv:
        summarize(path, args.max_pairs, args.output_prefix)


if __name__ == "__main__":
    main()
