#!/usr/bin/env bash
# FN Investigation SOP — Reproducible script for all assemblies
#
# Scans for ALT allele haplotype-carrying reads in input BAM, raw WGS BAM,
# and PacBio HiFi BAM. Records read qname + mapping position for every hit.
#
# Usage: bash run_fn_investigation_sop.sh
#
# Output: /paedyl01/disk1/yangyxt/test_tmp/fn_sop_results/

set -eo pipefail
set +u
eval "$(conda shell.bash hook 2>/dev/null)" && conda activate SDrecall
set -eo pipefail

OUTROOT=/paedyl01/disk1/yangyxt/test_tmp/fn_sop_results
SCANNER=/paedyl01/disk1/yangyxt/SDrecall-rust-migration/scripts/fn_haplotype_scanner.py
mkdir -p "$OUTROOT"

# ── hg19 ──────────────────────────────────────────────────────────────────
echo "=== hg19 FN investigation ==="
mkdir -p "$OUTROOT/hg19"

# FN sites (extracted from raw VCF benchmark comparison)
# Format: assembly  chrom  pos  ref  alt
awk -F: '{print "hg19\t"$1"\t"$2"\t"$3"\t"$4}' /tmp/hg19_all_fns.txt > "$OUTROOT/hg19/fn_sites.tsv"

python3 "$SCANNER" \
  --fn-list "$OUTROOT/hg19/fn_sites.tsv" \
  --ref /paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta \
  --gold /paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg19/golden_vcfs/HG002_hg19_v4.2.1_benchmark.vcf.gz \
  --pacbio /paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/golden_vcfs/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh37.bam \
  --input-bam /paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/hg19/HG002.SD.deduped.bam \
  --raw-bam /paedyl01/disk1/yangyxt/wgs/GIAB_samples/raw_data/download_data/HG002.hs37d5.300x.sdrelated.bam \
  --output "$OUTROOT/hg19/alt_reads.tsv" \
  2>&1 | tee "$OUTROOT/hg19/scanner.log"

echo "hg19 done: $(wc -l < "$OUTROOT/hg19/alt_reads.tsv") lines"

# ── hg38 ──────────────────────────────────────────────────────────────────
echo "=== hg38 FN investigation ==="
mkdir -p "$OUTROOT/hg38"

# hg38 FN sites from earlier investigation
# (extract from the hg38 raw VCF comparison — adjust path as needed)
awk -F: '{print "hg38\t"$1"\t"$2"\t"$3"\t"$4}' \
  /paedyl01/disk1/yangyxt/test_tmp/fn_investigation_hg38_20260625_170718/isec_no_n/0000.vcf 2>/dev/null | \
  awk -F'\t' '{if(NF>=5) print}' > "$OUTROOT/hg38/fn_sites.tsv"

python3 "$SCANNER" \
  --fn-list "$OUTROOT/hg38/fn_sites.tsv" \
  --ref /paedyl01/disk1/yangyxt/indexed_genome/hg38/ucsc.hg38.fasta \
  --gold /paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_hg38_v4.2.1_benchmark.vcf.gz \
  --pacbio /paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/golden_vcfs/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam \
  --input-bam /paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/hg38/HG002.SD.deduped.bam \
  --raw-bam /paedyl01/disk1/yangyxt/wgs/GIAB_samples/raw_data/download_data/HG002.GRCh38.300x.bam \
  --output "$OUTROOT/hg38/alt_reads.tsv" \
  2>&1 | tee "$OUTROOT/hg38/scanner.log"

echo "hg38 done: $(wc -l < "$OUTROOT/hg38/alt_reads.tsv") lines"

# ── t2t/chm13 ─────────────────────────────────────────────────────────────
echo "=== t2t FN investigation ==="
mkdir -p "$OUTROOT/t2t"

awk -F: '{print "t2t\t"$1"\t"$2"\t"$3"\t"$4}' /tmp/t2t_all_fns.txt > "$OUTROOT/t2t/fn_sites.tsv"

python3 "$SCANNER" \
  --fn-list "$OUTROOT/t2t/fn_sites.tsv" \
  --ref /paedyl01/disk1/yangyxt/indexed_genome/chm13/chm13.draft_v1.1.fasta \
  --gold /paedyl01/disk1/yangyxt/wgs/GIAB_samples/test_recall_precision/t2t_test/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.cmrg.vcf.gz \
  --pacbio /paedyl01/disk1/yangyxt/wgs/GIAB_samples/test_recall_precision/t2t_test/HG002_PacBio-HiFi-Revio_20231031_48x_CHM13v2.0.bam \
  --input-bam /paedyl01/disk1/yangyxt/wgs/GIAB_samples/test_recall_precision/t2t_test/HG002.t2t.chm13.SD.pairs.bam \
  --raw-bam /paedyl01/disk1/yangyxt/wgs/GIAB_samples/test_recall_precision/t2t_test/HG002.t2t.chm13.v2.bam \
  --output "$OUTROOT/t2t/alt_reads.tsv" \
  2>&1 | tee "$OUTROOT/t2t/scanner.log"

echo "t2t done: $(wc -l < "$OUTROOT/t2t/alt_reads.tsv") lines"

echo "=== ALL DONE ==="
echo "Results: $OUTROOT/{hg19,hg38,t2t}/alt_reads.tsv"
echo "Each row: FN_site | read_qname | mapping_chrom | mapping_pos | MAPQ | source_BAM"
