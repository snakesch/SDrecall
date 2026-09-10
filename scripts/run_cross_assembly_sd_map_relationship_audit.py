#!/usr/bin/env python3
"""Run the SD-map relationship audit for all three assemblies in parallel."""

from __future__ import annotations

import argparse
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WORKER = ROOT / "scripts/audit_sd_map_relationship_provenance.py"
MERGER = ROOT / "scripts/merge_sd_map_relationship_audits.py"
DEFAULT_OUTPUT_ROOT = ROOT / "test_tmp/sd_map_relationship_audit_20260724"
ASSEMBLIES = ("hg19", "hg38", "chm13")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--qname-paths", type=Path)
    return parser.parse_args()


def run_assembly(
    assembly: str, output_root: Path, qname_paths: Path | None
) -> tuple[str, str]:
    output_dir = output_root / "by_assembly" / assembly
    command = [
        sys.executable,
        str(WORKER),
        "--assembly",
        assembly,
        "--output-dir",
        str(output_dir),
    ]
    if qname_paths is not None:
        command.extend(("--qname-paths", str(qname_paths)))
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"{assembly} audit failed ({result.returncode})\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return assembly, result.stdout.strip()


def main() -> None:
    args = parse_args()
    failures: list[str] = []
    with ThreadPoolExecutor(max_workers=len(ASSEMBLIES)) as executor:
        futures = {
            executor.submit(
                run_assembly, assembly, args.output_root, args.qname_paths
            ): assembly
            for assembly in ASSEMBLIES
        }
        for future in as_completed(futures):
            assembly = futures[future]
            try:
                _, output = future.result()
            except Exception as error:
                failures.append(f"{assembly}: {error}")
            else:
                print(output)

    if failures:
        raise RuntimeError("\n\n".join(failures))

    subprocess.run(
        [
            sys.executable,
            str(MERGER),
            "--input-root",
            str(args.output_root / "by_assembly"),
            "--output-dir",
            str(args.output_root / "combined"),
        ],
        cwd=ROOT,
        check=True,
    )


if __name__ == "__main__":
    main()
