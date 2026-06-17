#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# submit_e2e_benchmark_pbs.sh — Submit 3 parallel PBS jobs for E2E benchmarking
#
# Submits one job per assembly (hg19, hg38, T2T/CHM13) to the medium queue.
# Each job: 10 CPUs, 200GB RAM, 72h walltime.
#
# Usage:
#   bash tests/e2e/submit_e2e_benchmark_pbs.sh
#
# Each job runs independently and writes results + logs to:
#   /paedyl01/disk1/yangyxt/test_tmp/rust_e2e_{assembly}/
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
LOG_DIR="/paedyl01/disk1/yangyxt/test_tmp/rust_e2e_logs"
mkdir -p "${LOG_DIR}"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# --- Assembly-specific configuration -----------------------------------------

# hg19
HG19_INPUT_BAM="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/hg19/HG002.SD.deduped.bam"
HG19_REF="/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta"
HG19_GOLD_VCF="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg19/golden_vcfs/HG002_GRCh37_1_22_v4.2.1_benchmark.addchr.bench.no_lost_alt.vcf.gz"
HG19_BENCH_BED="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg19/golden_vcfs/HG002_GRCh37_1_22_v4.2.1_benchmark.addchr.bed"
HG19_TARGET_SD="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg19/HG002_hg19_exome_SDrecall/realign_groups/all_target_recall_SD_regions.bed"
HG19_OUTPUT="/paedyl01/disk1/yangyxt/test_tmp/rust_e2e_hg19"

# hg38
HG38_INPUT_BAM="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/hg38/HG002.SD.deduped.bam"
HG38_REF="/paedyl01/disk1/yangyxt/indexed_genome/hg38/ucsc.hg38.fasta"
HG38_GOLD_VCF="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_GRCh38_1_22_v4.2.1_benchmark.norm.bench.no_lost_alt.vcf.gz"
HG38_BENCH_BED="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_GRCh38_1_22_v4.2.1_benchmark.bed"
HG38_TARGET_SD="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002_hg38_exome_SDrecall/realign_groups/all_target_recall_SD_regions.bed"
HG38_OUTPUT="/paedyl01/disk1/yangyxt/test_tmp/rust_e2e_hg38"

# T2T (CHM13 v2)
T2T_INPUT_BAM="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/test_recall_precision/t2t_test/HG002.t2t.chm13.v2.bam"
T2T_REF="/paedyl01/disk1/yangyxt/indexed_genome/chm13/chm13.draft_v1.1.fasta"
T2T_GOLD_VCF="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/test_recall_precision/t2t_test/HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.cmrg.no_lost_alt.vcf.gz"
T2T_BENCH_BED="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/test_recall_precision/t2t_test/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed"
T2T_TARGET_SD="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/test_recall_precision/t2t_test/HG002_chm13_cmrg_SDrecall/realign_groups/all_target_recall_SD_regions.bed"
T2T_OUTPUT="/paedyl01/disk1/yangyxt/test_tmp/rust_e2e_t2t"

# --- PBS parameters ----------------------------------------------------------
QUEUE="medium"
NCPUS=10
MEM="200gb"
WALLTIME="72:00:00"
CONDA_ENV="SDrecall"

# --- Generate per-assembly job script ----------------------------------------

generate_job_script() {
    local assembly="$1"
    local input_bam="$2"
    local ref_genome="$3"
    local gold_vcf="$4"
    local bench_bed="$5"
    local target_sd="$6"
    local output_dir="$7"

    local job_script="${LOG_DIR}/e2e_${assembly}_${TIMESTAMP}.pbs"
    local log_file="${LOG_DIR}/e2e_${assembly}_${TIMESTAMP}.log"

    mkdir -p "${output_dir}"

    cat > "${job_script}" << 'HEREDOC_HEADER'
#!/bin/bash
#PBS -V
HEREDOC_HEADER

    cat >> "${job_script}" << EOF
#PBS -N e2e_${assembly}
#PBS -q ${QUEUE}
#PBS -l select=1:ncpus=${NCPUS}:mem=${MEM}
#PBS -l walltime=${WALLTIME}
#PBS -o ${log_file}
#PBS -e ${log_file}
#PBS -j oe

# ============================================================
# E2E Benchmark: ${assembly} — HG002
# Generated: ${TIMESTAMP}
# ============================================================

set -euo pipefail

echo "[\$(date)] Job started on \$(hostname)"
echo "[\$(date)] Assembly: ${assembly}"
echo "[\$(date)] CPUs: ${NCPUS}, Memory: ${MEM}"
echo ""

# Activate environment
source ~/.bashrc
eval "\$(conda shell.bash hook 2>/dev/null)"
conda activate ${CONDA_ENV}
export LIBCLANG_PATH=\$CONDA_PREFIX/lib
export OPENSSL_NO_VENDOR=1
export PKG_CONFIG_PATH=\$CONDA_PREFIX/lib/pkgconfig:\$PKG_CONFIG_PATH

# Input paths
INPUT_BAM="${input_bam}"
REF_GENOME="${ref_genome}"
GOLD_VCF="${gold_vcf}"
BENCH_BED="${bench_bed}"
TARGET_SD="${target_sd}"
OUTPUT_DIR="${output_dir}"
NCPUS_LOCAL=${NCPUS}
SCRIPT_DIR="${SCRIPT_DIR}"

echo "[\$(date)] Input BAM:    \${INPUT_BAM}"
echo "[\$(date)] Reference:    \${REF_GENOME}"
echo "[\$(date)] Gold VCF:     \${GOLD_VCF}"
echo "[\$(date)] Bench BED:    \${BENCH_BED}"
echo "[\$(date)] Target SD:    \${TARGET_SD}"
echo "[\$(date)] Output dir:   \${OUTPUT_DIR}"
echo ""

# ---- Step 1: Build Rust binary (if needed) ----
echo "[\$(date)] Step 1: Checking/building Rust binary..."
RUST_BIN="${PROJECT_ROOT}/rust_modules/target/release/sdrecall"
if [[ ! -x "\${RUST_BIN}" ]]; then
    echo "[\$(date)]   Building from source..."
    cd ${PROJECT_ROOT}/rust_modules
    cargo build --release --bin sdrecall 2>&1
    echo "[\$(date)]   Build complete."
else
    echo "[\$(date)]   Binary already exists: \${RUST_BIN}"
fi

# ---- Step 2: Run Rust SDrecall ----
echo "[\$(date)] Step 2: Running Rust SDrecall pipeline..."
RUST_LOG_FILE="\${OUTPUT_DIR}/sdrecall_run_${assembly}.log"

/usr/bin/time -v "\${RUST_BIN}" \\
    --bam "\${INPUT_BAM}" \\
    --ref "\${REF_GENOME}" \\
    --output-dir "\${OUTPUT_DIR}" \\
    --threads \${NCPUS_LOCAL} \\
    > "\${RUST_LOG_FILE}" 2>&1 || {
    echo "[\$(date)] WARNING: Rust SDrecall exited with error. Check \${RUST_LOG_FILE}"
}

echo "[\$(date)] Step 2 complete. Output:"
ls -lh "\${OUTPUT_DIR}/"*.bam 2>/dev/null || echo "  (no BAM files found)"

# ---- Step 3: Find the output clean BAM ----
RUST_CLEAN_BAM=\$(find "\${OUTPUT_DIR}" -name "*clean*bam" -o -name "*filtered*bam" | head -1)
if [[ -z "\${RUST_CLEAN_BAM}" ]]; then
    echo "[\$(date)] ERROR: No clean/filtered BAM found in \${OUTPUT_DIR}"
    echo "[\$(date)] Attempting to use pooled merged BAM as fallback..."
    RUST_CLEAN_BAM=\$(find "\${OUTPUT_DIR}" -name "*merged*bam" | head -1)
fi
if [[ -z "\${RUST_CLEAN_BAM}" ]]; then
    echo "[\$(date)] FATAL: No output BAM available. Aborting."
    exit 1
fi
echo "[\$(date)] Using clean BAM: \${RUST_CLEAN_BAM}"

# ---- Step 4: Variant calling (GATK HaplotypeCaller) ----
echo "[\$(date)] Step 4: Calling variants..."
CALLED_VCF="\${OUTPUT_DIR}/HG002.rust.${assembly}.called.vcf.gz"

gatk --java-options "-Xmx160G" HaplotypeCaller \\
    -R "\${REF_GENOME}" \\
    -I "\${RUST_CLEAN_BAM}" \\
    -L "\${TARGET_SD}" \\
    -O "\${CALLED_VCF}" \\
    --native-pair-hmm-threads \${NCPUS_LOCAL} \\
    2>&1 | tee "\${OUTPUT_DIR}/gatk_hc_${assembly}.log"

echo "[\$(date)] Step 4 complete: \$(ls -lh \${CALLED_VCF})"

# ---- Step 5: Intersect with benchmark region ----
echo "[\$(date)] Step 5: Normalizing and intersecting VCFs..."
BENCH_SD_BED="\${OUTPUT_DIR}/benchmark_SD_region.bed"
GOLD_TARGETED="\${OUTPUT_DIR}/gold_targeted.vcf.gz"
CALLED_TARGETED="\${OUTPUT_DIR}/called_targeted.vcf.gz"

bedtools intersect -a "\${TARGET_SD}" -b "\${BENCH_BED}" > "\${BENCH_SD_BED}"
echo "[\$(date)]   Benchmark SD region intervals: \$(wc -l < \${BENCH_SD_BED})"

bcftools view -R "\${BENCH_SD_BED}" "\${GOLD_VCF}" | \\
bcftools norm -m -both -f "\${REF_GENOME}" --multi-overlaps 0 -a -Ou - | \\
bcftools norm -d exact -Ou - | \\
bcftools view -i 'ALT!="*"' -Oz -o "\${GOLD_TARGETED}"
tabix -f -p vcf "\${GOLD_TARGETED}"

bcftools view -R "\${BENCH_SD_BED}" "\${CALLED_VCF}" | \\
bcftools norm -m -both -f "\${REF_GENOME}" --multi-overlaps 0 -a -Ou - | \\
bcftools norm -d exact -Ou - | \\
bcftools view -i 'ALT!="*"' -Oz -o "\${CALLED_TARGETED}"
tabix -f -p vcf "\${CALLED_TARGETED}"

echo "[\$(date)]   Gold targeted:   \$(bcftools stats \${GOLD_TARGETED} | grep 'number of records' | head -1)"
echo "[\$(date)]   Called targeted: \$(bcftools stats \${CALLED_TARGETED} | grep 'number of records' | head -1)"

# ---- Step 6: Calculate precision & recall ----
echo "[\$(date)] Step 6: Calculating precision & recall (site + genotype modes)..."
RESULTS_FILE="\${OUTPUT_DIR}/benchmark_results_${assembly}.tsv"

python3 "\${SCRIPT_DIR}/calculate_precision_recall.py" \\
    --gold "\${GOLD_TARGETED}" \\
    --called "\${CALLED_TARGETED}" \\
    --output "\${RESULTS_FILE}" \\
    --assembly "${assembly}" \\
    --caller "SDrecall_rust"

echo ""
echo "[\$(date)] ============================================"
echo "[\$(date)] RESULTS (${assembly}):"
echo "[\$(date)] ============================================"
column -t -s \$'\\t' "\${RESULTS_FILE}"
echo ""

# ---- Step 7: Summary ----
echo "[\$(date)] Job complete for assembly=${assembly}"
echo "[\$(date)] Results file: \${RESULTS_FILE}"
echo "[\$(date)] Full log: ${log_file}"
EOF

    echo "${job_script}"
}

# --- Submit jobs --------------------------------------------------------------

echo "=== Submitting E2E Benchmark PBS Jobs ==="
echo "Queue: ${QUEUE} | CPUs: ${NCPUS} | RAM: ${MEM} | Walltime: ${WALLTIME}"
echo ""

# Generate job scripts
JOB_HG19=$(generate_job_script "hg19" "${HG19_INPUT_BAM}" "${HG19_REF}" "${HG19_GOLD_VCF}" "${HG19_BENCH_BED}" "${HG19_TARGET_SD}" "${HG19_OUTPUT}")
JOB_HG38=$(generate_job_script "hg38" "${HG38_INPUT_BAM}" "${HG38_REF}" "${HG38_GOLD_VCF}" "${HG38_BENCH_BED}" "${HG38_TARGET_SD}" "${HG38_OUTPUT}")
JOB_T2T=$(generate_job_script "t2t" "${T2T_INPUT_BAM}" "${T2T_REF}" "${T2T_GOLD_VCF}" "${T2T_BENCH_BED}" "${T2T_TARGET_SD}" "${T2T_OUTPUT}")

echo "Generated PBS scripts:"
echo "  hg19: ${JOB_HG19}"
echo "  hg38: ${JOB_HG38}"
echo "  T2T:  ${JOB_T2T}"
echo ""

# Submit all 3
echo "Submitting hg19..."
JOB_ID_HG19=$(qsub "${JOB_HG19}")
echo "  Job ID: ${JOB_ID_HG19}"

echo "Submitting hg38..."
JOB_ID_HG38=$(qsub "${JOB_HG38}")
echo "  Job ID: ${JOB_ID_HG38}"

echo "Submitting T2T..."
JOB_ID_T2T=$(qsub "${JOB_T2T}")
echo "  Job ID: ${JOB_ID_T2T}"

echo ""
echo "=== All 3 jobs submitted ==="
echo ""
echo "Monitor with:"
echo "  qstat -u \$(whoami)"
echo ""
echo "Logs:"
echo "  tail -f ${LOG_DIR}/e2e_hg19_${TIMESTAMP}.log"
echo "  tail -f ${LOG_DIR}/e2e_hg38_${TIMESTAMP}.log"
echo "  tail -f ${LOG_DIR}/e2e_t2t_${TIMESTAMP}.log"
echo ""
echo "Results (after completion):"
echo "  cat ${HG19_OUTPUT}/benchmark_results_hg19.tsv"
echo "  cat ${HG38_OUTPUT}/benchmark_results_hg38.tsv"
echo "  cat ${T2T_OUTPUT}/benchmark_results_t2t.tsv"
echo ""
echo "Aggregate all results:"
echo "  head -1 ${HG38_OUTPUT}/benchmark_results_hg38.tsv > all_results.tsv"
echo "  tail -n +2 -q ${HG19_OUTPUT}/benchmark_results_hg19.tsv ${HG38_OUTPUT}/benchmark_results_hg38.tsv ${T2T_OUTPUT}/benchmark_results_t2t.tsv >> all_results.tsv"
echo "  column -t -s \$'\t' all_results.tsv"
