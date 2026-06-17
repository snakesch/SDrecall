#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# run_benchmark_hg002.sh — End-to-end benchmark for Rust SDrecall on HG002
#
# Usage:
#   bash tests/e2e/run_benchmark_hg002.sh [--rust-binary <path>] [--threads <N>]
#
# This script:
#   1. Runs the Rust SDrecall pipeline on HG002 SD BAM
#   2. Calls variants on the output clean BAM (via GATK or DeepVariant)
#   3. Intersects called VCF with benchmark region
#   4. Calculates precision/recall against GIAB gold VCF
#   5. Reports PASS/FAIL based on thresholds (recall≥93%, precision≥45%)
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# --- Configuration -----------------------------------------------------------

# Paths (override via environment variables if needed)
INPUT_BAM="${E2E_INPUT_BAM:-/paedyl01/disk1/yangyxt/wgs/GIAB_samples/aligned_results/hg38/HG002.SD.deduped.bam}"
REF_GENOME="${E2E_REF_GENOME:-/paedyl01/disk1/yangyxt/indexed_genome/hg38/ucsc.hg38.fasta}"
GOLD_VCF="${E2E_GOLD_VCF:-/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_GRCh38_1_22_v4.2.1_benchmark.norm.bench.no_lost_alt.vcf.gz}"
BENCH_REGION="${E2E_BENCH_REGION:-/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_GRCh38_1_22_v4.2.1_benchmark.bed}"
TARGET_SD_REGIONS="${E2E_TARGET_SD:-/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002_hg38_exome_SDrecall/realign_groups/all_target_recall_SD_regions.bed}"

# Output directory
OUTPUT_DIR="${E2E_OUTPUT_DIR:-/paedyl01/disk1/yangyxt/test_tmp/rust_e2e_HG002}"
THREADS="${E2E_THREADS:-8}"

# Thresholds
MIN_RECALL="0.93"
MIN_PRECISION="0.45"

# Rust binary (default: build from source)
DEFAULT_RUST_BINARY="${PROJECT_ROOT}/rust_modules/target/release/sdrecall"
RUST_BINARY="${E2E_RUST_BINARY:-${DEFAULT_RUST_BINARY}}"

# --- Argument parsing --------------------------------------------------------

while [[ $# -gt 0 ]]; do
    case $1 in
        --rust-binary) RUST_BINARY="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --help|-h)
            echo "Usage: $0 [--rust-binary <path>] [--threads <N>] [--output-dir <dir>]"
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# --- Utility functions -------------------------------------------------------

timestamp() { date +"%Y-%m-%d %H:%M:%S"; }

log() {
    echo "[$(timestamp)] [E2E] $*" >&2
}

check_file() {
    if [[ ! -f "$1" ]]; then
        log "ERROR: Required file not found: $1"
        exit 1
    fi
}

# --- Pre-flight checks -------------------------------------------------------

log "=== SDrecall Rust E2E Benchmark ==="
log "Input BAM:        ${INPUT_BAM}"
log "Reference:        ${REF_GENOME}"
log "Gold VCF:         ${GOLD_VCF}"
log "Benchmark region: ${BENCH_REGION}"
log "Target SD:        ${TARGET_SD_REGIONS}"
log "Output dir:       ${OUTPUT_DIR}"
log "Threads:          ${THREADS}"
log "Rust binary:      ${RUST_BINARY}"
log ""

check_file "${INPUT_BAM}"
check_file "${REF_GENOME}"
check_file "${GOLD_VCF}"
check_file "${BENCH_REGION}"
check_file "${TARGET_SD_REGIONS}"

mkdir -p "${OUTPUT_DIR}"

# --- Step 1: Build Rust binary -----------------------------------------------
# Always (re)build the DEFAULT binary from source so a stale prebuilt binary can
# never be silently reused (this exact trap once looked like a recurring code
# bug). cargo's incremental build is a fast no-op when the sources are unchanged.
# An explicit --rust-binary / E2E_RUST_BINARY override is treated as a
# user-provided binary and used as-is (no rebuild).

if [[ "${RUST_BINARY}" == "${DEFAULT_RUST_BINARY}" ]]; then
    log "Step 1: (Re)building Rust SDrecall binary from source (HEAD) to avoid stale-binary reuse..."
    (
        cd "${PROJECT_ROOT}/rust_modules"
        eval "$(conda shell.bash hook 2>/dev/null)" && conda activate SDrecall
        export LIBCLANG_PATH=$CONDA_PREFIX/lib
        export OPENSSL_NO_VENDOR=1
        export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH
        cargo build --release --bin sdrecall
    )
    check_file "${RUST_BINARY}"
else
    log "Step 1: Using user-provided Rust binary at ${RUST_BINARY} (no rebuild)"
    check_file "${RUST_BINARY}"
fi

# --- Step 2: Run Rust SDrecall pipeline --------------------------------------

RUST_CLEAN_BAM="${OUTPUT_DIR}/HG002.rust.clean.bam"
RUST_LOG_FILE="${OUTPUT_DIR}/rust_sdrecall_run.log"

log "Step 2: Running Rust SDrecall pipeline..."
/usr/bin/time -v "${RUST_BINARY}" \
    --bam "${INPUT_BAM}" \
    --ref "${REF_GENOME}" \
    --output-dir "${OUTPUT_DIR}" \
    --threads "${THREADS}" \
    2>&1 | tee "${RUST_LOG_FILE}"

# Verify output
if [[ ! -f "${RUST_CLEAN_BAM}" ]]; then
    log "WARNING: Expected output ${RUST_CLEAN_BAM} not found."
    log "Check logs at ${RUST_LOG_FILE}"
    log "Listing output directory:"
    ls -lh "${OUTPUT_DIR}/"
    log "Attempting to find any BAM output..."
    RUST_CLEAN_BAM=$(find "${OUTPUT_DIR}" -name "*.clean.bam" -o -name "*.filtered.bam" | head -1)
    if [[ -z "${RUST_CLEAN_BAM}" ]]; then
        log "ERROR: No clean BAM produced. Pipeline failed."
        exit 1
    fi
    log "Found: ${RUST_CLEAN_BAM}"
fi

# --- Step 3: Variant calling -------------------------------------------------

RUST_CALLED_VCF="${OUTPUT_DIR}/HG002.rust.called.vcf.gz"

log "Step 3: Calling variants with GATK HaplotypeCaller..."
gatk --java-options "-Xmx16G" HaplotypeCaller \
    -R "${REF_GENOME}" \
    -I "${RUST_CLEAN_BAM}" \
    -L "${TARGET_SD_REGIONS}" \
    -O "${RUST_CALLED_VCF}" \
    --native-pair-hmm-threads 4 \
    2>&1 | tee "${OUTPUT_DIR}/gatk_hc.log"

check_file "${RUST_CALLED_VCF}"

# --- Step 4: Prepare benchmark region and normalize VCFs ---------------------

BENCH_SD_BED="${OUTPUT_DIR}/benchmark_SD_region.bed"
GOLD_TARGETED="${OUTPUT_DIR}/gold_targeted.vcf.gz"
CALLED_TARGETED="${OUTPUT_DIR}/called_targeted.vcf.gz"

log "Step 4: Intersecting target SD with benchmark confidence region..."
bedtools intersect \
    -a "${TARGET_SD_REGIONS}" \
    -b "${BENCH_REGION}" \
    > "${BENCH_SD_BED}"

REGION_COUNT=$(wc -l < "${BENCH_SD_BED}")
log "  Benchmark SD region: ${REGION_COUNT} intervals"

log "Step 4b: Normalizing gold VCF..."
bcftools view -R "${BENCH_SD_BED}" "${GOLD_VCF}" | \
bcftools norm -m -both -f "${REF_GENOME}" --multi-overlaps 0 -a -Ou - | \
bcftools norm -d exact -Ou - | \
bcftools view -i 'ALT!="*"' -Oz -o "${GOLD_TARGETED}"
tabix -f -p vcf "${GOLD_TARGETED}"

log "Step 4c: Normalizing called VCF..."
bcftools view -R "${BENCH_SD_BED}" "${RUST_CALLED_VCF}" | \
bcftools norm -m -both -f "${REF_GENOME}" --multi-overlaps 0 -a -Ou - | \
bcftools norm -d exact -Ou - | \
bcftools view -i 'ALT!="*"' -Oz -o "${CALLED_TARGETED}"
tabix -f -p vcf "${CALLED_TARGETED}"

# --- Step 5: Calculate precision and recall ----------------------------------

RESULTS_FILE="${OUTPUT_DIR}/benchmark_results.tsv"

log "Step 5: Calculating precision and recall..."
python3 "${SCRIPT_DIR}/calculate_precision_recall.py" \
    --gold "${GOLD_TARGETED}" \
    --called "${CALLED_TARGETED}" \
    --output "${RESULTS_FILE}" \
    --mode site

# --- Step 6: Evaluate pass/fail ---------------------------------------------

log ""
log "=== BENCHMARK RESULTS ==="
cat "${RESULTS_FILE}"
log ""

RECALL=$(awk -F'\t' 'NR==2 {print $3}' "${RESULTS_FILE}")
PRECISION=$(awk -F'\t' 'NR==2 {print $4}' "${RESULTS_FILE}")

log "Recall:    ${RECALL} (threshold: >= ${MIN_RECALL})"
log "Precision: ${PRECISION} (threshold: >= ${MIN_PRECISION})"

PASS=true
if (( $(echo "${RECALL} < ${MIN_RECALL}" | bc -l) )); then
    log "FAIL: Recall ${RECALL} < ${MIN_RECALL}"
    PASS=false
fi
if (( $(echo "${PRECISION} < ${MIN_PRECISION}" | bc -l) )); then
    log "FAIL: Precision ${PRECISION} < ${MIN_PRECISION}"
    PASS=false
fi

if [[ "${PASS}" == "true" ]]; then
    log "========================================="
    log "  PASS: All thresholds met!"
    log "========================================="
    exit 0
else
    log "========================================="
    log "  FAIL: Thresholds not met."
    log "========================================="
    exit 1
fi
