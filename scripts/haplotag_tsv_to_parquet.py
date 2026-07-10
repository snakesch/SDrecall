#!/usr/bin/env python3
"""Convert WhatsHap haplotag TSV output to partitioned ZSTD Parquet.

The default mode is row preserving: every TSV row becomes one Parquet row with a
compact schema. The writer uses Polars batches and writes per-chromosome part
files incrementally, avoiding the high memory use of a single partitioned sink
on multi-billion-row haplotag lists.

Exact template/read-pair compaction requires grouping by read name across
billions of rows. That should be done from BAM HP/PS tags during selection or
with an external-sort workflow, not as a naive global in-memory group-by.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import polars as pl


TSV_COLUMNS = ["#readname", "haplotype", "phaseset", "chromosome"]


def haplotag_scan(path: str, assembly: str, tagged_only: bool) -> pl.LazyFrame:
    hp = (
        pl.col("haplotype")
        .replace({"none": "0"})
        .str.strip_prefix("H")
        .cast(pl.UInt8, strict=False)
        .fill_null(0)
    )
    phase_set = (
        pl.when(pl.col("phaseset") == "none")
        .then(None)
        .otherwise(pl.col("phaseset"))
        .cast(pl.UInt64, strict=False)
    )

    frame = (
        pl.scan_csv(
            path,
            separator="\t",
            has_header=True,
            schema_overrides={column: pl.String for column in TSV_COLUMNS},
            infer_schema_length=0,
        )
        .select(
            [
                pl.lit(assembly).alias("assembly"),
                pl.col("chromosome"),
                pl.col("#readname").alias("read_name"),
                hp.alias("hp"),
                phase_set.alias("phase_set"),
                pl.lit(1).cast(pl.UInt16).alias("alignment_rows"),
            ]
        )
    )
    if tagged_only:
        frame = frame.filter(pl.col("hp") > 0)
    return frame


def safe_partition_value(value: object) -> str:
    text = str(value)
    return re.sub(r"[^A-Za-z0-9_.-]", "_", text)


def transform_batch(batch: pl.DataFrame, assembly: str, tagged_only: bool) -> pl.DataFrame:
    hp = (
        pl.col("haplotype")
        .replace({"none": "0"})
        .str.strip_prefix("H")
        .cast(pl.UInt8, strict=False)
        .fill_null(0)
    )
    phase_set = (
        pl.when(pl.col("phaseset") == "none")
        .then(None)
        .otherwise(pl.col("phaseset"))
        .cast(pl.UInt64, strict=False)
    )
    out = batch.select(
        [
            pl.lit(assembly).alias("assembly"),
            pl.col("chromosome"),
            pl.col("#readname").alias("read_name"),
            hp.alias("hp"),
            phase_set.alias("phase_set"),
            pl.lit(1).cast(pl.UInt16).alias("alignment_rows"),
        ]
    )
    if tagged_only:
        out = out.filter(pl.col("hp") > 0)
    return out


def write_rows(
    path: str,
    assembly: str,
    out_dir: str,
    tagged_only: bool,
    max_rows_per_file: int,
    compression_level: int | None,
    batch_rows: int,
    progress_batches: int,
) -> None:
    del max_rows_per_file  # Batches are the part-file boundary in this writer.
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    counters: dict[str, int] = {}
    manifest_rows: list[tuple[str, str, int, str]] = []

    reader = pl.read_csv_batched(
        path,
        separator="\t",
        has_header=True,
        columns=TSV_COLUMNS,
        schema_overrides={column: pl.String for column in TSV_COLUMNS},
        infer_schema_length=0,
        batch_size=batch_rows,
    )

    total_rows = 0
    batch_index = 0
    while True:
        raw_batches = reader.next_batches(1)
        if not raw_batches:
            break
        for raw_batch in raw_batches:
            batch_index += 1
            batch = transform_batch(raw_batch, assembly, tagged_only)
            total_rows += batch.height
            if not batch.is_empty():
                for key, part in batch.partition_by("chromosome", as_dict=True).items():
                    chrom = key[0] if isinstance(key, tuple) else key
                    chrom_text = safe_partition_value(chrom)
                    partition_dir = (
                        out_path / f"assembly={assembly}" / f"chromosome={chrom_text}"
                    )
                    partition_dir.mkdir(parents=True, exist_ok=True)
                    index = counters.get(chrom_text, 0)
                    counters[chrom_text] = index + 1
                    file_path = partition_dir / f"{index:08d}.parquet"
                    part.write_parquet(
                        file_path,
                        compression="zstd",
                        compression_level=compression_level,
                        statistics=True,
                        row_group_size=min(1_000_000, max(1, part.height)),
                    )
                    manifest_rows.append(
                        (assembly, str(chrom), part.height, str(file_path))
                    )
            if progress_batches and batch_index % progress_batches == 0:
                print(
                    f"{assembly}: batches={batch_index} rows={total_rows}",
                    flush=True,
                )

    manifest = out_path / f"_manifest.{assembly}.tsv"
    with manifest.open("wt", encoding="utf-8") as handle:
        handle.write("assembly\tchromosome\trows\tpath\n")
        for row in manifest_rows:
            handle.write("\t".join(map(str, row)) + "\n")
    (out_path / f"_SUCCESS.{assembly}").write_text("ok\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert WhatsHap haplotag TSV to partitioned ZSTD Parquet"
    )
    parser.add_argument("--input", required=True, help="WhatsHap haplotag TSV")
    parser.add_argument("--assembly", required=True, help="Assembly partition label")
    parser.add_argument(
        "--out-dir",
        required=True,
        help="Output directory, e.g. haplotags.parquet",
    )
    parser.add_argument(
        "--tagged-only",
        action="store_true",
        help="Only keep H1/H2/etc rows; drop haplotype=none rows",
    )
    parser.add_argument(
        "--max-rows-per-file",
        type=int,
        default=25_000_000,
        help=(
            "Deprecated compatibility option. The batch writer uses "
            "--batch-rows as the part-file boundary."
        ),
    )
    parser.add_argument(
        "--compression-level",
        type=int,
        default=6,
        help="ZSTD compression level. Default: 6",
    )
    parser.add_argument(
        "--batch-rows",
        type=int,
        default=5_000_000,
        help="Rows per Polars batch/part-file before chromosome splitting. Default: 5,000,000",
    )
    parser.add_argument(
        "--progress-batches",
        type=int,
        default=20,
        help="Print progress every N batches. Default: 20",
    )
    parser.add_argument(
        "--allow-existing",
        action="store_true",
        help="Allow writing into an existing dataset directory.",
    )
    args = parser.parse_args()

    out_path = Path(args.out_dir)
    if out_path.exists() and any(out_path.iterdir()) and not args.allow_existing:
        raise SystemExit(f"Output directory is not empty: {out_path}")

    print(f"polars={pl.__version__} index={pl.get_index_type()}", flush=True)
    print(f"input={args.input}", flush=True)
    print(f"assembly={args.assembly}", flush=True)
    print(f"out_dir={args.out_dir}", flush=True)
    print(f"tagged_only={args.tagged_only}", flush=True)
    write_rows(
        args.input,
        args.assembly,
        args.out_dir,
        args.tagged_only,
        args.max_rows_per_file,
        args.compression_level,
        args.batch_rows,
        args.progress_batches,
    )
    print("done", flush=True)


if __name__ == "__main__":
    main()
