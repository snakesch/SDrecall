#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    printf 'usage: %s BASELINE_RECALL_RESULTS CANDIDATE_RECALL_RESULTS\n' "$0" >&2
    exit 2
fi

baseline=$1
candidate=$2
artifacts=(
    HG002.sdrecall.clean.vcf.gz
    HG002.sdrecall.clean.vcf.gz.csi
    HG002.sdrecall.merged.vcf.gz
    HG002.sdrecall.merged.vcf.gz.csi
    HG002.sdrecall.vcf.gz
    HG002.sdrecall.vcf.gz.csi
)

checked=0
for artifact in "${artifacts[@]}"; do
    baseline_file=${baseline}/${artifact}
    candidate_file=${candidate}/${artifact}
    if [[ ! -e "$baseline_file" && ! -e "$candidate_file" ]]; then
        continue
    fi
    if [[ ! -f "$baseline_file" || ! -f "$candidate_file" ]]; then
        printf 'missing parity artifact: %s\n' "$artifact" >&2
        exit 1
    fi
    if ! cmp -s "$baseline_file" "$candidate_file"; then
        printf 'byte mismatch: %s\n' "$artifact" >&2
        sha256sum "$baseline_file" "$candidate_file" >&2
        exit 1
    fi
    sha256sum "$candidate_file"
    checked=$((checked + 1))
done

if [[ $checked -eq 0 ]]; then
    printf 'no parity artifacts found\n' >&2
    exit 1
fi

printf 'byte parity verified for %d artifacts\n' "$checked"

# BAM headers contain run-specific process IDs in @PG command lines. Compare
# all alignment records byte-for-byte as SAM text and compare headers after
# removing only those provenance lines.
baseline_bam=${baseline}/HG002.pooled.clean.bam
candidate_bam=${candidate}/HG002.pooled.clean.bam
for bam in "$baseline_bam" "$candidate_bam"; do
    if [[ ! -f "$bam" ]]; then
        printf 'missing parity artifact: %s\n' "$bam" >&2
        exit 1
    fi
done

if ! cmp -s \
    <(samtools view -H "$baseline_bam" | awk '$1 != "@PG"') \
    <(samtools view -H "$candidate_bam" | awk '$1 != "@PG"'); then
    printf 'BAM non-PG header mismatch\n' >&2
    exit 1
fi
if ! cmp -s <(samtools view "$baseline_bam") <(samtools view "$candidate_bam"); then
    printf 'BAM alignment-record mismatch\n' >&2
    exit 1
fi
printf 'canonical BAM parity verified\n'
