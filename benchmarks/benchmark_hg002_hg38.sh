#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKSPACE_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)

SDRECALL_VCF=/paedyl01/disk1/yangyxt/SDrecall-test/results/sdrecall_rust_external_helper_limits_20260717/hg38/HG002_hg38_exome_avg50x_SDrecall/recall_results/HG002.sdrecall.vcf.gz
DEEPVARIANT_VCF=/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/HG002.deepvariant.sorted.vcf.gz
TRUTH_VCF=/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_hg38_v4.2.1_benchmark.vcf.gz
TRUTH_BED=/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_hg38_v4.2.1_benchmark.bed
REFERENCE=/paedyl01/disk1/yangyxt/indexed_genome/hg38/ucsc.hg38.fasta
RECALL_BED=/paedyl01/disk1/yangyxt/SDrecall-test/results/sdrecall_rust_external_helper_limits_20260717/hg38/HG002_hg38_exome_avg50x_SDrecall/realign_groups/all_target_recall_SD_regions.bed
VCF_OPS=${WORKSPACE_ROOT}/rust_modules/target/release/vcf-ops
EVALUATOR=${SCRIPT_DIR}/cal_prec_recall_vcf.py
PYTHON=python3
OUTPUT_DIR=/paedyl01/disk1/yangyxt/SDrecall-test/benchmarks/hg002_hg38_sdrecall_deepvariant_20260717
TMP_DIR=
THREADS=4
DRY_RUN=0
FORCE=0
ASSEMBLY_LABEL=hg38
TRUTH_LABEL=giab_v4.2.1
INCLUDE_MISALIGNED=0

usage() {
    cat <<'USAGE'
Benchmark a merged Rust SDrecall + DeepVariant HG002 callset against a truth
VCF within recall-SD regions that overlap the truth confident BED.

Activate the SDrecall environment before running:
  source ~/.bashrc
  mamba activate SDrecall

Usage:
  benchmark_hg002_hg38.sh [options]

Inputs (defaults point to the validated 2026-07-17 HG002 hg38 run):
  --sdrecall-vcf PATH       Rust SDrecall query VCF
  --deepvariant-vcf PATH    conventional caller/reference VCF
  --truth-vcf PATH          GIAB truth VCF
  --truth-bed PATH          GIAB confident-region BED
  --reference PATH          hg38 FASTA (requires PATH.fai)
  --recall-bed PATH         all_target_recall_SD_regions.bed
  --assembly-label NAME     output label such as hg19, hg38, or chm13
  --truth-label NAME        output label such as giab_v4.2.1 or defrabb_v0.020

Tools and outputs:
  --vcf-ops PATH            release vcf-ops executable
  --evaluator PATH          exact genotype-aware Python evaluator
  --python PATH             Python with vcfpy installed (default: python3)
  --output-dir PATH         durable benchmark output directory
  --tmp-dir PATH            temporary workspace (default: OUTPUT_DIR/tmp)
  --threads N               bcftools/vcf-ops threads, 1-255 (default: 4)
  --include-misaligned      retain SDrecall records tagged MISALIGNED
  --force                   permit overwriting this workflow's output files
  --dry-run                 validate inputs and print commands without writing
  -h, --help                show this help

The priority merge keeps DeepVariant records at matching alleles and adds
FILTER provenance tags DeepVariant and SDrecall. Query-only SDrecall records
are retained with the SDrecall tag. No heuristic genotype correction is
applied; later normalization can remap split/atomized genotypes.

By default, the evaluated merged callset excludes records carrying the
SDrecall MISALIGNED filter. Pass --include-misaligned for the sensitivity
audit that retains those records. Other FILTER values remain as provenance,
so neither mode is a generic PASS-only filter.
USAGE
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

need_value() {
    [[ $# -ge 2 && -n ${2:-} ]] || die "$1 requires a value"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sdrecall-vcf)
            need_value "$@"
            SDRECALL_VCF=$2
            shift 2
            ;;
        --deepvariant-vcf)
            need_value "$@"
            DEEPVARIANT_VCF=$2
            shift 2
            ;;
        --truth-vcf)
            need_value "$@"
            TRUTH_VCF=$2
            shift 2
            ;;
        --truth-bed)
            need_value "$@"
            TRUTH_BED=$2
            shift 2
            ;;
        --reference)
            need_value "$@"
            REFERENCE=$2
            shift 2
            ;;
        --recall-bed)
            need_value "$@"
            RECALL_BED=$2
            shift 2
            ;;
        --assembly-label)
            need_value "$@"
            ASSEMBLY_LABEL=$2
            shift 2
            ;;
        --truth-label)
            need_value "$@"
            TRUTH_LABEL=$2
            shift 2
            ;;
        --vcf-ops)
            need_value "$@"
            VCF_OPS=$2
            shift 2
            ;;
        --evaluator)
            need_value "$@"
            EVALUATOR=$2
            shift 2
            ;;
        --python)
            need_value "$@"
            PYTHON=$2
            shift 2
            ;;
        --output-dir)
            need_value "$@"
            OUTPUT_DIR=$2
            shift 2
            ;;
        --tmp-dir)
            need_value "$@"
            TMP_DIR=$2
            shift 2
            ;;
        --threads)
            need_value "$@"
            THREADS=$2
            shift 2
            ;;
        --include-misaligned)
            INCLUDE_MISALIGNED=1
            shift
            ;;
        --force)
            FORCE=1
            shift
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "unknown option: $1"
            ;;
    esac
done

[[ ${THREADS} =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer"
(( THREADS <= 255 )) || die "--threads must be <= 255 for vcf-ops"
[[ ${ASSEMBLY_LABEL} =~ ^[A-Za-z0-9._-]+$ ]] || die "--assembly-label contains unsafe characters"
[[ ${TRUTH_LABEL} =~ ^[A-Za-z0-9._-]+$ ]] || die "--truth-label contains unsafe characters"

OUTPUT_DIR=$(realpath -m -- "${OUTPUT_DIR}")
if [[ -z ${TMP_DIR} ]]; then
    TMP_DIR=${OUTPUT_DIR}/tmp
else
    TMP_DIR=$(realpath -m -- "${TMP_DIR}")
fi

LOG_DIR=${OUTPUT_DIR}/logs
VCF_DIR=${OUTPUT_DIR}/vcfs
METRICS_DIR=${OUTPUT_DIR}/metrics
REGION_DIR=${OUTPUT_DIR}/regions
RUN_LOG=${LOG_DIR}/benchmark.log
COMMAND_LOG=${LOG_DIR}/commands.log

OVERLAP_RAW=${REGION_DIR}/HG002.${ASSEMBLY_LABEL}.recall_x_truth_confident.raw.bed
OVERLAP_SORTED=${REGION_DIR}/HG002.${ASSEMBLY_LABEL}.recall_x_truth_confident.sorted.bed
CALLABLE_BED=${REGION_DIR}/HG002.${ASSEMBLY_LABEL}.recall_x_truth_confident.merged.bed
MERGED_RAW=${VCF_DIR}/HG002.sdrecall.deepvariant.merged.raw.vcf.gz
MERGED_NORM=${VCF_DIR}/HG002.sdrecall.deepvariant.callable.norm.vcf.gz
DEEPVARIANT_NORM=${VCF_DIR}/HG002.deepvariant.callable.norm.vcf.gz
TRUTH_NORM=${VCF_DIR}/HG002.${TRUTH_LABEL}.callable.norm.vcf.gz
MERGED_METRICS=${METRICS_DIR}/HG002.sdrecall_deepvariant_vs_${TRUTH_LABEL}.tsv
BASELINE_METRICS=${METRICS_DIR}/HG002.deepvariant_vs_${TRUTH_LABEL}.tsv
COMBINED_METRICS=${METRICS_DIR}/HG002.callset_comparison.tsv

for path in \
    "${SDRECALL_VCF}" \
    "${DEEPVARIANT_VCF}" \
    "${TRUTH_VCF}" \
    "${TRUTH_BED}" \
    "${REFERENCE}" \
    "${REFERENCE}.fai" \
    "${RECALL_BED}" \
    "${EVALUATOR}"; do
    [[ -f ${path} ]] || die "required file not found: ${path}"
done

if [[ ${VCF_OPS} == */* ]]; then
    [[ -x ${VCF_OPS} ]] || die "vcf-ops is not executable: ${VCF_OPS}"
else
    VCF_OPS=$(command -v -- "${VCF_OPS}") || die "vcf-ops not found: ${VCF_OPS}"
fi

if [[ ${PYTHON} == */* ]]; then
    [[ -x ${PYTHON} ]] || die "Python is not executable: ${PYTHON}"
else
    PYTHON=$(command -v -- "${PYTHON}") || die "Python not found: ${PYTHON}"
fi

for tool in bcftools bedtools awk cat find mkdir realpath sha256sum tee wc; do
    command -v -- "${tool}" >/dev/null || die "required command not found: ${tool}"
done

"${PYTHON}" -c 'import vcfpy' >/dev/null || die "${PYTHON} cannot import vcfpy"
for input_vcf in "${SDRECALL_VCF}" "${DEEPVARIANT_VCF}" "${TRUTH_VCF}"; do
    bcftools view -h "${input_vcf}" >/dev/null || die "invalid VCF: ${input_vcf}"
done

for indexed_vcf in "${DEEPVARIANT_VCF}" "${TRUTH_VCF}"; do
    [[ -f ${indexed_vcf}.csi || -f ${indexed_vcf}.tbi ]] || \
        die "indexed region access requires ${indexed_vcf}.csi or ${indexed_vcf}.tbi"
done

mapfile -t sdrecall_samples < <(bcftools query -l "${SDRECALL_VCF}")
mapfile -t deepvariant_samples < <(bcftools query -l "${DEEPVARIANT_VCF}")
mapfile -t truth_samples < <(bcftools query -l "${TRUTH_VCF}")
[[ ${#sdrecall_samples[@]} -eq 1 ]] || die "SDrecall VCF must contain one sample"
[[ ${#deepvariant_samples[@]} -eq 1 ]] || die "DeepVariant VCF must contain one sample"
[[ ${#truth_samples[@]} -eq 1 ]] || die "truth VCF must contain one sample"
[[ ${sdrecall_samples[0]} == "${deepvariant_samples[0]}" && \
   ${sdrecall_samples[0]} == "${truth_samples[0]}" ]] || \
    die "sample mismatch across input VCFs"

if (( ! DRY_RUN )); then
    if [[ -d ${OUTPUT_DIR} && -n $(find "${OUTPUT_DIR}" -mindepth 1 -maxdepth 1 -print -quit) && ${FORCE} -eq 0 ]]; then
        die "output directory is not empty; pass --force to reuse it: ${OUTPUT_DIR}"
    fi
    mkdir -p "${LOG_DIR}" "${VCF_DIR}" "${METRICS_DIR}" "${REGION_DIR}" "${TMP_DIR}"
    : > "${COMMAND_LOG}"
    : > "${RUN_LOG}"
    exec > >(tee -a "${RUN_LOG}") 2>&1
fi

print_command() {
    local token
    printf '$'
    for token in "$@"; do
        printf ' %q' "${token}"
    done
    printf '\n'
}

record_command() {
    print_command "$@"
    if (( ! DRY_RUN )); then
        print_command "$@" >> "${COMMAND_LOG}"
    fi
}

print_command_to_file() {
    local output=$1
    shift
    local token
    printf '$'
    for token in "$@"; do
        printf ' %q' "${token}"
    done
    printf ' > %q\n' "${output}"
}

run() {
    record_command "$@"
    if (( ! DRY_RUN )); then
        "$@"
    fi
}

run_to_file() {
    local output=$1
    shift
    print_command_to_file "${output}" "$@"
    if (( ! DRY_RUN )); then
        print_command_to_file "${output}" "$@" >> "${COMMAND_LOG}"
    fi
    if (( ! DRY_RUN )); then
        "$@" > "${output}"
    fi
}

normalize_and_subset() {
    local label=$1
    local input_vcf=$2
    local output_vcf=$3
    local include_expression=$4
    local work_dir=${TMP_DIR}/${label}
    local normalized_bcf=${work_dir}/normalized.bcf
    local deduplicated_bcf=${work_dir}/deduplicated.bcf
    local variant_bcf=${work_dir}/variant_only.bcf
    local sorted_vcf=${work_dir}/sorted.vcf.gz

    run mkdir -p "${work_dir}/sort"
    run bcftools norm \
        --threads "${THREADS}" \
        -m -both \
        -f "${REFERENCE}" \
        --multi-overlaps 0 \
        -a \
        -Ob \
        -o "${normalized_bcf}" \
        "${input_vcf}"
    run bcftools norm \
        --threads "${THREADS}" \
        -d exact \
        -Ob \
        -o "${deduplicated_bcf}" \
        "${normalized_bcf}"
    run bcftools filter \
        --threads "${THREADS}" \
        -i "${include_expression}" \
        -Ob \
        -o "${variant_bcf}" \
        "${deduplicated_bcf}"
    run bcftools sort \
        --temp-dir "${work_dir}/sort" \
        -Oz \
        -o "${sorted_vcf}" \
        "${variant_bcf}"
    run bcftools index -f "${sorted_vcf}"
    run bcftools view \
        --threads "${THREADS}" \
        -R "${CALLABLE_BED}" \
        --regions-overlap 1 \
        -Oz \
        -o "${output_vcf}" \
        "${sorted_vcf}"
    run bcftools index -f "${output_vcf}"
}

printf 'HG002 %s Rust benchmark workflow\n' "${ASSEMBLY_LABEL}"
printf 'sample=%s\ntruth_label=%s\ninclude_misaligned=%s\noutput_dir=%s\ntmp_dir=%s\nthreads=%s\n' \
    "${sdrecall_samples[0]}" "${TRUTH_LABEL}" "${INCLUDE_MISALIGNED}" \
    "${OUTPUT_DIR}" "${TMP_DIR}" "${THREADS}"
run bcftools --version
run bedtools --version
run "${PYTHON}" --version
run sha256sum "${VCF_OPS}" "${EVALUATOR}"

run mkdir -p "${LOG_DIR}" "${VCF_DIR}" "${METRICS_DIR}" "${REGION_DIR}" "${TMP_DIR}"
run_to_file "${OVERLAP_RAW}" bedtools intersect -a "${RECALL_BED}" -b "${TRUTH_BED}"
run_to_file "${OVERLAP_SORTED}" bedtools sort -i "${OVERLAP_RAW}"
run_to_file "${CALLABLE_BED}" bedtools merge -i "${OVERLAP_SORTED}"

if (( ! DRY_RUN )); then
    [[ -s ${CALLABLE_BED} ]] || die "recall/confident BED intersection is empty"
    callable_intervals=$(wc -l < "${CALLABLE_BED}")
    callable_bases=$(awk '{ total += $3 - $2 } END { print total + 0 }' "${CALLABLE_BED}")
    printf 'callable_intervals=%s\ncallable_bases=%s\n' \
        "${callable_intervals}" "${callable_bases}"
fi

run "${VCF_OPS}" \
    --log-level info \
    merge \
    --query-vcf "${SDRECALL_VCF}" \
    --reference-vcf "${DEEPVARIANT_VCF}" \
    --output-vcf "${MERGED_RAW}" \
    --ref-genome "${REFERENCE}" \
    --qv-tag SDrecall \
    --rv-tag DeepVariant \
    --threads "${THREADS}" \
    --tmp-dir "${TMP_DIR}"

BASE_VARIANT_EXPRESSION='ALT[0] != "*" && COUNT(GT="alt") > 0'
if (( INCLUDE_MISALIGNED )); then
    MERGED_VARIANT_EXPRESSION=${BASE_VARIANT_EXPRESSION}
else
    MERGED_VARIANT_EXPRESSION="${BASE_VARIANT_EXPRESSION} && FILTER!~\"MISALIGNED\""
fi

normalize_and_subset merged "${MERGED_RAW}" "${MERGED_NORM}" "${MERGED_VARIANT_EXPRESSION}"
normalize_and_subset deepvariant \
    "${DEEPVARIANT_VCF}" \
    "${DEEPVARIANT_NORM}" \
    "${BASE_VARIANT_EXPRESSION}"
normalize_and_subset truth "${TRUTH_VCF}" "${TRUTH_NORM}" "${BASE_VARIANT_EXPRESSION}"

run "${PYTHON}" "${EVALUATOR}" \
    "${MERGED_NORM}" \
    "${TRUTH_NORM}" \
    "${MERGED_METRICS}" \
    --label SDrecall+DeepVariant
run "${PYTHON}" "${EVALUATOR}" \
    "${DEEPVARIANT_NORM}" \
    "${TRUTH_NORM}" \
    "${BASELINE_METRICS}" \
    --label DeepVariant
run_to_file "${COMBINED_METRICS}" awk \
    'FNR == 1 && NR != 1 { next } { print }' \
    "${MERGED_METRICS}" \
    "${BASELINE_METRICS}"
run cat "${COMBINED_METRICS}"

if (( DRY_RUN )); then
    printf 'Dry run complete: no benchmark outputs were written.\n'
else
    printf 'Benchmark complete. Metrics: %s\n' "${COMBINED_METRICS}"
    printf 'Merged command/output log: %s\n' "${RUN_LOG}"
    printf 'Command manifest: %s\n' "${COMMAND_LOG}"
fi
