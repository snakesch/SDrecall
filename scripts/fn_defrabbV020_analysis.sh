#!/usr/bin/env bash
# Correct FN analysis: restrict gold VCF to SDrecall FC target regions (the callable set)
# before computing FN. This matches what the pipeline actually attempts to call.
set -o pipefail
eval "$(conda shell.bash hook 2>/dev/null)" && conda activate SDrecall

DB=/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/defrabbV0.020_phased
OUT=/paedyl01/disk1/yangyxt/test_tmp/fn_defrabbV020
mkdir -p "$OUT"

run_asm() {
  local asm=$1 gold=$2 bed=$3 rawnorm=$4 fc_beds_glob=$5
  echo ""
  echo "########## $asm ##########"

  # 1. Build SDrecall FC target union (the callable region set)
  local fc=$(mktemp)
  cat $fc_beds_glob 2>/dev/null | sort -k1,1 -k2,2n | bedtools merge -i - > "$fc"
  local fc_kb=$(awk '{sum+=$3-$2}END{print sum/1000}' "$fc")
  echo "SDrecall FC target span: $fc_kb kb"

  # 2. Rename gold BED/VCF to chr-prefix if hg19
  local gbed="$bed"
  local gvcf="$gold"
  if [ "$asm" = "hg19" ]; then
    local gbed2=$(mktemp); sed 's/^\([0-9]\|X\|Y\|MT\)\t/chr\1\t/' "$bed" > "$gbed2"; gbed="$gbed2"
    local gvcf2=$(mktemp); bcftools annotate --rename-chrs <(for c in $(seq 1 22) X Y MT; do printf "%s\tchr%s\n" "$c" "$c"; done) "$gold" 2>/dev/null | bgzip > "$gvcf2"; tabix "$gvcf2"; gvcf="$gvcf2"
  fi

  # 3. Intersect FC target with gold benchmark bed -> callable SD subset
  local fc_callable=$(mktemp)
  bedtools intersect -a "$fc" -b "$gbed" > "$fc_callable"
  local call_kb=$(awk '{sum+=$3-$2}END{print sum/1000}' "$fc_callable")
  echo "FC within gold benchmark (callable) span: $call_kb kb"

  # 4. Restrict gold VCF to callable SD regions
  zcat "$gvcf" | grep -v '^#' | wc -l | xargs echo "gold total variants:"
  local gold_in_fc=$(mktemp)
  bedtools intersect -a "$gvcf" -b "$fc_callable" -header 2>/dev/null | bgzip > "$gold_in_fc"
  tabix "$gold_in_fc"
  local g_fc=$(zcat "$gold_in_fc" | grep -v '^#' | wc -l)
  echo "gold variants in callable SD regions: $g_fc"

  # 5. Normalize and split gold + raw.norm for key comparison
  local goldn="$OUT/${asm}.gold.norm.vcf.gz"
  local rawn="$OUT/${asm}.raw.norm.split.vcf.gz"
  zcat "$gold_in_fc" | bcftools norm -m -any 2>/dev/null | bgzip > "$goldn"; tabix "$goldn"
  bcftools norm -m -any "$rawnorm" 2>/dev/null | bgzip > "$rawn"; tabix "$rawn"
  local gnorm=$(zcat "$goldn" | grep -v '^#' | wc -l)
  local rnorm=$(zcat "$rawn" | grep -v '^#' | wc -l)
  echo "gold (normalized, in FC): $gnorm   raw.norm (split): $rnorm"

  # 6. Allele-aware key subtraction: FN = gold keys NOT in raw keys
  zcat "$rawn" | grep -v '^#' | awk -F'\t' '{print $1":"$2":"$4":"$5}' | sort -u > "$OUT/${asm}.raw.keys"
  zcat "$goldn" | grep -v '^#' | awk -F'\t' 'BEGIN{while((getline k<"'$OUT/${asm}.raw.keys'")>0) seen[k]=1}
    {key=$1":"$2":"$4":"$5; if(!(key in seen)) print}' > "$OUT/${asm}.FN.vcf"
  local fn=$(grep -vc '^#' "$OUT/${asm}.FN.vcf")
  # TP = gold keys IN raw
  local tp=$(zcat "$goldn" | grep -v '^#' | awk -F'\t' 'BEGIN{while((getline k<"'$OUT/${asm}.raw.keys'")>0) seen[k]=1}
    {key=$1":"$2":"$4":"$5; if(key in seen) c++} END{print c+0}')
  echo "FN (gold in FC not in raw.norm): $fn"
  echo "TP (shared): $tp"
  echo "recall = TP / gold_in_FC = $tp / $gnorm"
  rm -f "$fc" "$fc_callable" "$gold_in_fc"
}

run_asm hg19 \
  "$DB/GRCh37_HG2-T2TQ100-V1.1_smvar.vcf.gz" \
  "$DB/GRCh37_HG2-T2TQ100-V1.1_smvar.benchmark.bed" \
  /paedyl01/disk1/yangyxt/test_tmp/r9fix_hg19_eval/raw.norm.vcf.gz \
  "/paedyl01/disk1/yangyxt/test_tmp/rust_e2e_hg19/HG002_hg19_exome_SDrecall/realign_groups/RG*/RG*.fc_target.bed"

run_asm hg38 \
  "$DB/GRCh38_HG2-T2TQ100-V1.1_smvar.vcf.gz" \
  "$DB/GRCh38_HG2-T2TQ100-V1.1_smvar.benchmark.bed" \
  /paedyl01/disk1/yangyxt/test_tmp/hg38_alt_eval/raw.norm.vcf.gz \
  "/paedyl01/disk1/yangyxt/test_tmp/rust_e2e_hg38/HG002_hg38_exome_SDrecall/realign_groups/RG*/RG*.fc_target.bed"

run_asm t2t \
  "$DB/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.vcf.gz" \
  "$DB/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed" \
  /paedyl01/disk1/yangyxt/test_tmp/r9fix_t2t_eval/raw.norm.vcf.gz \
  "/paedyl01/disk1/yangyxt/test_tmp/rust_e2e_t2t/HG002_chm13_exome_SDrecall/realign_groups/RG*/RG*.fc_target.bed"

echo ""
echo "=== DONE — FN lists: $OUT/{hg19,hg38,t2t}.FN.vcf ==="
