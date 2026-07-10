#!/usr/bin/env python3
"""
FN Allele Haplotype Read Scanner — SOP for raw VCF false negative investigation.

For each FN site:
1. Build ALT 51mer from reference (with adjacent golden VCF variants if present)
2. If adjacent variants exist, refine 51mer using PacBio HiFi reads
3. Scan input BAM (forward-strand SEQ only — BAM SEQ is reference-forward)
4. Scan raw WGS BAM (if available)
5. Output: FN_site, read_qname, mapping_chrom, mapping_pos, MAPQ, source_BAM

Usage:
  python fn_haplotype_scanner.py \
    --fn-list FN_SITES.tsv \
    --ref REFERENCE.fasta \
    --gold GOLDEN.vcf.gz \
    --pacbio PACBIO.bam \
    --input-bam INPUT.bam \
    --raw-bam RAW_WGS.bam \
    --output OUTPUT.tsv

FN_SITES.tsv format (tab-separated, one per line):
  assembly  chrom  pos  ref  alt
"""

import argparse
import subprocess
import sys
import os
import collections
from typing import List, Tuple, Optional

def revcomp(s: str) -> str:
    c = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    return ''.join(c.get(b, 'N') for b in reversed(s))

def faidx_single(ref: str, chrom: str, start: int, end: int) -> str:
    """Fetch reference sequence (1-based inclusive)."""
    result = subprocess.run(
        ["samtools", "faidx", ref, f"{chrom}:{start}-{end}"],
        capture_output=True, text=True
    )
    lines = [l for l in result.stdout.split("\n") if l and not l.startswith(">")]
    return "".join(lines).upper()

def get_adjacent_variants(gold_vcf: str, chrom: str, pos: int, window: int = 25) -> list:
    """Find golden VCF variants within ±window bp of pos."""
    start = max(1, pos - window)
    end = pos + window
    result = subprocess.run(
        ["bcftools", "view", "-r", f"{chrom}:{start}-{end}", gold_vcf],
        capture_output=True, text=True
    )
    variants = []
    for line in result.stdout.split("\n"):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        v_pos = int(f[1])
        if v_pos == pos:
            continue  # skip the FN variant itself
        variants.append((f[0], v_pos, f[3], f[4]))
    return variants

def build_51mer_from_ref(ref: str, chrom: str, pos: int, ref_base: str, alt_base: str) -> str:
    """Build ALT 51mer centered on variant position."""
    start = max(1, pos - 25)
    seq = faidx_single(ref, chrom, start, pos + 25)
    if len(seq) < 51:
        return ""
    off = pos - start
    if seq[off] != ref_base:
        return ""
    return seq[:off] + alt_base + seq[off+1:]

def refine_51mer_from_pacbio(
    pacbio_bam: str, ref: str, chrom: str, pos: int,
    ref_base: str, alt_base: str, adjacent_variants: list
) -> List[str]:
    """If adjacent variants exist, extract actual ALT haplotype from PacBio reads."""
    if not adjacent_variants:
        return []

    # Get PacBio reads spanning the position
    start = max(1, pos - 25)
    result = subprocess.run(
        ["samtools", "view", pacbio_bam, f"{chrom}:{start}-{pos+25}"],
        capture_output=True, text=True
    )

    haplotypes = set()
    for line in result.stdout.split("\n"):
        if not line.strip():
            continue
        f = line.split("\t")
        seq = f[9].upper()  # BAM SEQ is reference-forward
        flag = int(f[1])
        read_pos = int(f[3])

        # Find offset of variant in this read
        cigar = f[5]
        # Simple: if read spans pos, extract 51bp window centered on pos
        if flag & 16:  # reverse strand — SEQ already RC'd to ref-forward by aligner
            pass

        # Check if read carries ALT at variant position
        # Parse CIGAR to find the base at genomic pos
        rel_off = pos - read_pos
        if rel_off < 0 or rel_off >= len(seq):
            continue

        if seq[rel_off] == alt_base:
            # Extract 51mer centered on variant
            mer_start = max(0, rel_off - 25)
            mer_end = min(len(seq), rel_off + 26)
            if mer_end - mer_start >= 51:
                haplotypes.add(seq[mer_start:mer_start+51])

    return list(haplotypes)

def scan_bam_for_kmers(
    bam_path: str, kmers: dict,  # kmer -> FN_site
    output_file: str,
    bam_label: str
):
    """Scan BAM for reads carrying any ALT 51mer.
    BAM SEQ is reference-forward, so only check forward kmers.

    Output per hit: FN_site, read_qname, mapping_chrom, mapping_pos, MAPQ, source_bam
    """
    proc = subprocess.Popen(
        ["samtools", "view", bam_path],
        stdout=subprocess.PIPE, text=True, bufsize=1024*1024
    )

    hits = []
    for line in proc.stdout:
        f = line.strip().split("\t")
        if len(f) < 10:
            continue
        qname = f[0]
        chrom = f[2]
        pos = int(f[3])
        mapq = int(f[4])
        seq = f[9].upper()

        for kmer, site in kmers.items():
            if kmer in seq:
                hits.append((site, qname, chrom, pos, mapq, bam_label))
                break  # one match per read

    proc.wait()

    with open(output_file, "a") as out:
        for site, qname, chrom, pos, mapq, label in hits:
            out.write(f"{site}\t{qname}\t{chrom}\t{pos}\t{mapq}\t{label}\n")

    return len(hits)

def main():
    parser = argparse.ArgumentParser(description="FN allele haplotype read scanner")
    parser.add_argument("--fn-list", required=True, help="TSV: assembly, chrom, pos, ref, alt")
    parser.add_argument("--ref", required=True, help="Reference FASTA")
    parser.add_argument("--gold", required=True, help="Golden VCF (gz)")
    parser.add_argument("--pacbio", required=True, help="PacBio HiFi BAM")
    parser.add_argument("--input-bam", required=True, help="SDrecall input BAM (sliced)")
    parser.add_argument("--raw-bam", default="", help="Raw WGS BAM (optional)")
    parser.add_argument("--output", required=True, help="Output TSV")
    parser.add_argument("--kmer-len", type=int, default=51, help="Kmer length")
    args = parser.parse_args()

    half = args.kmer_len // 2

    # Read FN sites
    fn_sites = []
    with open(args.fn_list) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                fn_sites.append((parts[0], parts[1], int(parts[2]), parts[3], parts[4]))

    print(f"# Loaded {len(fn_sites)} FN sites", file=sys.stderr)

    # Build 51mers
    all_kmers = {}  # forward_kmer -> FN_site
    pacbio_refined = 0
    for asm, chrom, pos, ref_b, alt_b in fn_sites:
        # Step 1: Build basic 51mer from reference
        kmer = build_51mer_from_ref(args.ref, chrom, pos, ref_b, alt_b)
        if not kmer:
            print(f"# SKIP {chrom}:{pos} (ref mismatch or too short)", file=sys.stderr)
            continue

        # Step 2: Check for adjacent variants
        adjacent = get_adjacent_variants(args.gold, chrom, pos, half)

        if adjacent:
            # Refine using PacBio
            refined = refine_51mer_from_pacbio(
                args.pacbio, args.ref, chrom, pos, ref_b, alt_b, adjacent
            )
            if refined:
                for h in refined:
                    all_kmers[h] = f"{chrom}:{pos}:{ref_b}:{alt_b}"
                pacbio_refined += 1
                print(f"# {chrom}:{pos} refined from PacBio ({len(refined)} haplotypes)", file=sys.stderr)
                continue

        # Use basic 51mer (no adjacent variants or PacBio refinement failed)
        all_kmers[kmer] = f"{chrom}:{pos}:{ref_b}:{alt_b}"

    print(f"# Built {len(all_kmers)} kmers ({pacbio_refined} PacBio-refined)", file=sys.stderr)

    # Write header
    with open(args.output, "w") as out:
        out.write("FN_site\tread_qname\tmapping_chrom\tmapping_pos\tMAPQ\tsource_BAM\n")

    # Step 3: Scan input BAM
    print(f"# Scanning input BAM: {args.input_bam}", file=sys.stderr)
    n = scan_bam_for_kmers(args.input_bam, all_kmers, args.output, "input_bam")
    print(f"# Input BAM: {n} hits", file=sys.stderr)

    # Step 4: Scan raw WGS BAM (if provided)
    if args.raw_bam and os.path.exists(args.raw_bam):
        print(f"# Scanning raw BAM: {args.raw_bam}", file=sys.stderr)
        n = scan_bam_for_kmers(args.raw_bam, all_kmers, args.output, "raw_bam")
        print(f"# Raw BAM: {n} hits", file=sys.stderr)

    # Step 5: Scan PacBio BAM (as last resort for Category 1)
    print(f"# Scanning PacBio BAM: {args.pacbio}", file=sys.stderr)
    n = scan_bam_for_kmers(args.pacbio, all_kmers, args.output, "pacbio")
    print(f"# PacBio BAM: {n} hits", file=sys.stderr)

    print(f"# DONE. Results: {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()
