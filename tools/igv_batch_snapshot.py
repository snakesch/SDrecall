#!/usr/bin/env python3
"""Generate IGV batch scripts for haplotype-grouped alignment snapshots.

Produces headless PNGs via xvfb-run + IGV, with reads grouped by HP tag
and labeled separators between haplotype groups.  Designed for validating
SDrecall's phasing/misalignment calls against the raw alignment.

Prerequisites (run once per session):
    source /etc/profile.d/easybuild.sh
    module load Xvfb/1.20.13-GCCcore-11.2.0 IGV/2.16.0-Java-11

    The --ref must be a LOCAL .fa/.fasta file (with .fai index).
    IGV genome IDs like "hg38" require internet access and will hang
    on the firewalled HPC.

Usage:
    # Single region
    python igv_batch_snapshot.py \\
        --bam realigned.bam --ref /path/to/hg38.fa \\
        --regions chr1:1633000-1635000 \\
        --outdir /paedyl01/disk1/yangyxt/test_tmp/igv_snapshots

    # Regions from BED file
    python igv_batch_snapshot.py \\
        --bam realigned.bam --ref /path/to/hg38.fa \\
        --bed islands.bed --pad 200 \\
        --outdir /paedyl01/disk1/yangyxt/test_tmp/igv_snapshots

    # Generate batch script only (for manual IGV run)
    python igv_batch_snapshot.py \\
        --bam realigned.bam --ref /path/to/hg38.fa \\
        --regions chr1:1633000-1635000 \\
        --outdir /tmp/igv_out --batch-only
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


# ── IGV preferences for alignment validation ──────────────────────────
# Snapshot shows: HP-grouped reads with labeled separators, base-level
# mismatches colored by nucleotide (A=green T=red C=blue G=orange),
# insertions marked, BQ shading on mismatches, MQ flagging, soft clips.
IGV_PREFERENCES = {
    # Group by HP tag — labeled separators between haplotype groups
    "SAM.GROUP_OPTION": "TAG",
    "SAM.GROUP_BY_TAG": "HP",
    # Expanded: every read gets its own row
    "SAM.DISPLAY_MODE": "EXPANDED",
    # Mismatches colored by nucleotide (always on)
    "SAM.SHOW_MISMATCHES": "TRUE",
    # BQ shading: mismatches with BQ<5 nearly invisible, BQ>=20 fully opaque
    "SAM.SHADE_BASE_QUALITY": "TRUE",
    "SAM.BASE_QUALITY_MIN": "5",
    "SAM.BASE_QUALITY_MAX": "20",
    # Insertions: purple I-markers + size labels
    "SAM.SHOW_INSERTION_MARKERS": "TRUE",
    "SAM.FLAG_LARGE_INDELS": "TRUE",
    "SAM.LARGE_INSERTIONS_THRESOLD": "1",
    # Soft clips visible (important for realigned reads)
    "SAM.SHOW_SOFT_CLIPPED": "TRUE",
    # MQ=0 reads flagged with transparency + solid outline
    "SAM.FLAG_ZERO_QUALITY": "TRUE",
    # MQ shading: low-MQ reads faded
    "SAM.SHADE_ALIGNMENT_BY": "MAPPING_QUALITY_LOW",
    "SAM.SHADE_QUALITY_LOW": "10",
    "SAM.SHADE_QUALITY_HIGH": "60",
    # Show ALL reads (no downsampling)
    "SAM.DOWNSAMPLE_READS": "FALSE",
    # Keep supplementary + secondary visible
    "SAM.FILTER_SUPPLEMENTARY_ALIGNMENTS": "FALSE",
    "SAM.FILTER_SECONDARY_ALIGNMENTS": "FALSE",
    # Coverage track on, junction track off
    "SAM.SHOW_COV_TRACK": "TRUE",
    "SAM.SHOW_JUNCTION_TRACK": "FALSE",
    # Center line for orientation
    "SAM.SHOW_CENTER_LINE": "TRUE",
}

# HPC module + Xvfb configuration
MODULES_CMD = (
    "source /etc/profile.d/easybuild.sh 2>/dev/null; "
    "module load Xvfb/1.20.13-GCCcore-11.2.0 IGV/2.16.0-Java-11 2>/dev/null"
)
HPC_ENV_FIX = (
    'export LD_LIBRARY_PATH="/usr/lib64:$LD_LIBRARY_PATH"; '
    'export GSETTINGS_BACKEND=memory; '
    'export GIO_EXTRA_MODULES=/usr/lib64/gio/modules; '
    'export XDG_DATA_DIRS="/usr/share:${XDG_DATA_DIRS:-/usr/local/share:/usr/share}"'
)


def parse_region(region_str: str) -> tuple[str, int, int]:
    """Parse 'chr:start-end' into (chrom, start, end)."""
    chrom, coords = region_str.split(":")
    start, end = coords.replace(",", "").split("-")
    return chrom, int(start), int(end)


def format_region(chrom: str, start: int, end: int) -> str:
    return f"{chrom}:{start}-{end}"


def regions_from_bed(bed_path: str, pad: int = 0) -> list[str]:
    """Read BED file, return list of 'chr:start-end' regions."""
    regions = []
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            chrom = parts[0]
            start = max(0, int(parts[1]) - pad)
            end = int(parts[2]) + pad
            regions.append(format_region(chrom, start, end))
    return regions


def generate_batch_script(
    bam_path: str,
    genome: str,
    regions: list[str],
    outdir: str,
    max_panel_height: int = 2000,
) -> str:
    """Generate an IGV batch script string."""
    lines = ["new"]
    lines.append(f"genome {genome}")

    for key, val in IGV_PREFERENCES.items():
        lines.append(f"preference {key} {val}")

    lines.append(f"maxPanelHeight {max_panel_height}")
    lines.append(f"load {bam_path}")
    lines.append(f"snapshotDirectory {outdir}")

    for region in regions:
        chrom, start, end = parse_region(region)
        safe_name = f"{chrom}_{start}_{end}.png"
        lines.append(f"goto {region}")
        lines.append("sort position")
        lines.append(f"snapshot {safe_name}")

    lines.append("exit")
    return "\n".join(lines) + "\n"


def run_igv_batch(
    batch_script_path: str,
    memory: str = "4g",
    display_size: str = "2560x1440x24",
    timeout: int = 300,
) -> int:
    """Run IGV batch script headlessly via xvfb-run on HPC."""
    shell_script = (
        f"{MODULES_CMD}; "
        f"{HPC_ENV_FIX}; "
        f'xvfb-run --auto-servernum --server-args="-screen 0 {display_size}" '
        f"igv.sh -b {batch_script_path}"
    )
    print(f"Running IGV (timeout {timeout}s)...", file=sys.stderr)
    try:
        result = subprocess.run(
            ["bash", "-c", shell_script],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        if result.stdout:
            print(result.stdout, file=sys.stderr)
        if result.returncode != 0 and result.stderr:
            print(f"IGV stderr:\n{result.stderr}", file=sys.stderr)
        return result.returncode
    except subprocess.TimeoutExpired:
        print(
            f"IGV timed out after {timeout}s. If using a genome ID (hg38), "
            "switch to a local FASTA path — the HPC firewall blocks genome downloads.",
            file=sys.stderr,
        )
        return 1


def main():
    parser = argparse.ArgumentParser(
        description="Generate IGV batch snapshots with HP-tag haplotype grouping",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--bam", required=True, help="BAM file path")
    parser.add_argument(
        "--ref",
        required=True,
        help="Reference genome: local FASTA path (must have .fai index)",
    )
    parser.add_argument(
        "--regions", nargs="+", help="Regions as chr:start-end (space-separated)"
    )
    parser.add_argument("--bed", help="BED file with regions")
    parser.add_argument(
        "--pad", type=int, default=0, help="Padding around each region (bp)"
    )
    parser.add_argument("--outdir", required=True, help="Output directory for PNGs")
    parser.add_argument(
        "--max-panel-height",
        type=int,
        default=2000,
        help="Max panel height in pixels (default: 2000)",
    )
    parser.add_argument(
        "--memory", default="4g", help="IGV JVM memory (default: 4g)"
    )
    parser.add_argument(
        "--display-size",
        default="2560x1440x24",
        help="Virtual display size (default: 2560x1440x24)",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=300,
        help="IGV process timeout in seconds (default: 300)",
    )
    parser.add_argument(
        "--batch-only",
        action="store_true",
        help="Only generate batch script, don't run IGV",
    )

    args = parser.parse_args()

    # Collect regions
    regions = []
    if args.regions:
        for r in args.regions:
            if args.pad:
                chrom, start, end = parse_region(r)
                start = max(0, start - args.pad)
                end += args.pad
                regions.append(format_region(chrom, start, end))
            else:
                regions.append(r)
    if args.bed:
        regions.extend(regions_from_bed(args.bed, pad=args.pad))
    if not regions:
        parser.error("Provide --regions or --bed")

    # Resolve paths
    bam_path = os.path.abspath(args.bam)
    ref = os.path.abspath(args.ref)
    if not os.path.isfile(ref):
        parser.error(f"Reference not found: {ref}")
    fai = ref + ".fai"
    if not os.path.isfile(fai):
        parser.error(f"Reference index not found: {fai}  (run: samtools faidx {ref})")

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    # Generate batch script
    script = generate_batch_script(
        bam_path, ref, regions, outdir, args.max_panel_height
    )
    batch_path = os.path.join(outdir, "igv_batch.txt")
    with open(batch_path, "w") as f:
        f.write(script)
    print(f"Batch script: {batch_path}", file=sys.stderr)

    if args.batch_only:
        print(script)
        return

    # Run IGV
    rc = run_igv_batch(
        batch_path,
        memory=args.memory,
        display_size=args.display_size,
        timeout=args.timeout,
    )
    if rc == 0:
        pngs = sorted(Path(outdir).glob("*.png"))
        print(f"Generated {len(pngs)} snapshots in {outdir}", file=sys.stderr)
        for p in pngs:
            print(p)
    else:
        print(f"IGV exited with code {rc}", file=sys.stderr)
        sys.exit(rc)


if __name__ == "__main__":
    main()
