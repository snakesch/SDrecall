#!/usr/bin/env bash
set -euo pipefail

# Evaluate already completed Rust SDrecall HG002 runs against GIAB truth sets.
#
# This script intentionally does not run SDrecall and does not call variants
# from BAM. It compares the final Rust VCFs that already exist in the e2e
# output folders.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

OUT_DIR="${E2E_EVAL_OUT_DIR:-/paedyl01/disk1/yangyxt/test_tmp/rust_e2e_eval}"
ASSEMBLIES="hg38,hg19,t2t"
CALLER="SDrecall_rust"
EXCLUDE_FILTERS=("MISALIGNED")
VCF_STAGE="final"

usage() {
    cat <<EOF
Usage:
  bash tests/e2e/evaluate_existing_rust_outputs.sh [options]

Options:
  --assemblies LIST       Comma-separated subset: hg38,hg19,t2t (default: ${ASSEMBLIES})
  --output-dir DIR        Evaluation output directory (default: ${OUT_DIR})
  --caller NAME           Caller label in result TSV (default: ${CALLER})
  --vcf-stage STAGE       VCF to evaluate: raw, clean, merged, final
                          (default: ${VCF_STAGE})
  --exclude-filters LIST  Comma-separated FILTER tokens to skip in called VCF
                          (default: MISALIGNED; use 'none' to disable)
  --help                  Show this help

Inputs are the completed Rust output folders under /paedyl01/disk1/yangyxt/test_tmp/rust_e2e_*.
For each assembly the script:
  1. Intersects Rust all_target_recall_SD_regions.bed with the GIAB benchmark BED.
  2. Restricts gold and called VCFs to that intersection.
  3. Atomizes/splits/left-aligns both VCFs with bcftools norm.
  4. Runs calculate_precision_recall.py for site-level and genotype-level metrics.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --assemblies)
            ASSEMBLIES="$2"
            shift 2
            ;;
        --output-dir)
            OUT_DIR="$2"
            shift 2
            ;;
        --caller)
            CALLER="$2"
            shift 2
            ;;
        --vcf-stage)
            VCF_STAGE="$2"
            shift 2
            ;;
        --exclude-filters)
            IFS=',' read -r -a EXCLUDE_FILTERS <<< "$2"
            if [[ "$2" == "none" ]]; then
                EXCLUDE_FILTERS=()
            fi
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

timestamp() {
    date +"%Y-%m-%d %H:%M:%S"
}

log() {
    echo "[$(timestamp)] [eval] $*" >&2
}

require_file() {
    local path="$1"
    if [[ ! -f "$path" ]]; then
        log "ERROR: required file not found: $path"
        exit 1
    fi
}

require_tool() {
    local tool="$1"
    if ! command -v "$tool" >/dev/null 2>&1; then
        log "ERROR: required tool not found on PATH: $tool"
        exit 1
    fi
}

assembly_config() {
    local assembly="$1"

    case "$assembly" in
        hg38)
            RUN_DIR="/paedyl01/disk1/yangyxt/test_tmp/rust_e2e_hg38/HG002_hg38_exome_SDrecall"
            REF_GENOME="/paedyl01/disk1/yangyxt/indexed_genome/hg38/ucsc.hg38.fasta"
            GOLD_VCF="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_hg38_v4.2.1_benchmark.vcf.gz"
            BENCH_BED="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg38/golden_vcfs/HG002_hg38_v4.2.1_benchmark.bed"
            NORM_CHECK_MODE=()
            ;;
        hg19)
            RUN_DIR="/paedyl01/disk1/yangyxt/test_tmp/rust_e2e_hg19/HG002_hg19_exome_SDrecall"
            REF_GENOME="/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta"
            GOLD_VCF="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg19/golden_vcfs/HG002_hg19_v4.2.1_benchmark.vcf.gz"
            BENCH_BED="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/vcfs/hg19/golden_vcfs/HG002_hg19_v4.2.1_benchmark.bed"
            NORM_CHECK_MODE=()
            ;;
        t2t|chm13)
            RUN_DIR="/paedyl01/disk1/yangyxt/test_tmp/rust_e2e_t2t/HG002_chm13_exome_SDrecall"
            REF_GENOME="/paedyl01/disk1/yangyxt/indexed_genome/chm13/chm13.draft_v1.1.fasta"
            GOLD_VCF="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/test_recall_precision/t2t_test/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.norm.vcf.gz"
            BENCH_BED="/paedyl01/disk1/yangyxt/wgs/GIAB_samples/test_recall_precision/t2t_test/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed"
            # The historical CHM13 benchmark script used -c s for this truth set.
            NORM_CHECK_MODE=(-c s)
            ;;
        *)
            log "ERROR: unknown assembly: $assembly"
            exit 2
            ;;
    esac

    case "$VCF_STAGE" in
        raw)
            CALLED_VCF="${RUN_DIR}/recall_results/HG002.sdrecall.raw.vcf.gz"
            ;;
        clean)
            CALLED_VCF="${RUN_DIR}/recall_results/HG002.sdrecall.clean.vcf.gz"
            ;;
        merged)
            CALLED_VCF="${RUN_DIR}/recall_results/HG002.sdrecall.merged.vcf.gz"
            ;;
        final)
            CALLED_VCF="${RUN_DIR}/recall_results/HG002.sdrecall.vcf.gz"
            ;;
        *)
            log "ERROR: unknown VCF stage: $VCF_STAGE"
            exit 2
            ;;
    esac
    RECALL_BED="${RUN_DIR}/realign_groups/all_target_recall_SD_regions.bed"
}

normalize_vcf() {
    local input_vcf="$1"
    local region_bed="$2"
    local ref_genome="$3"
    local output_vcf="$4"
    shift 4
    local -a norm_check_mode=("$@")

    bcftools view -R "$region_bed" -Ou "$input_vcf" | \
        bcftools norm -m -both "${norm_check_mode[@]}" -f "$ref_genome" --multi-overlaps 0 -a -Ou - | \
        bcftools norm -d exact -Ou - | \
        bcftools view -i 'ALT!="*"' -Oz -o "$output_vcf"
    bcftools index -t -f "$output_vcf"
}

evaluate_one() {
    local assembly="$1"
    assembly_config "$assembly"

    local asm_label="$assembly"
    if [[ "$asm_label" == "chm13" ]]; then
        asm_label="t2t"
    fi

    require_file "$REF_GENOME"
    require_file "$GOLD_VCF"
    require_file "$BENCH_BED"
    require_file "$CALLED_VCF"
    require_file "$RECALL_BED"

    local asm_out="${OUT_DIR}/${VCF_STAGE}/${asm_label}"
    mkdir -p "$asm_out"

    local eval_bed="${asm_out}/rust_recall_x_benchmark.bed"
    local gold_targeted="${asm_out}/HG002.gold.targeted.norm.vcf.gz"
    local called_targeted="${asm_out}/HG002.sdrecall.${VCF_STAGE}.targeted.norm.vcf.gz"
    local results="${asm_out}/benchmark_results_${asm_label}_${VCF_STAGE}.tsv"

    log "Evaluating ${asm_label} (${VCF_STAGE} VCF)"
    log "  called VCF: $CALLED_VCF"
    log "  gold VCF:   $GOLD_VCF"
    log "  recall BED: $RECALL_BED"
    log "  bench BED:  $BENCH_BED"

    bedtools intersect -a "$RECALL_BED" -b "$BENCH_BED" | \
        LC_ALL=C sort -k1,1 -k2,2n -k3,3n > "$eval_bed"

    local interval_count
    interval_count="$(wc -l < "$eval_bed")"
    log "  evaluation intervals: ${interval_count}"
    if [[ "$interval_count" -eq 0 ]]; then
        log "ERROR: empty evaluation BED for ${asm_label}"
        return 1
    fi

    log "  normalizing gold VCF"
    normalize_vcf "$GOLD_VCF" "$eval_bed" "$REF_GENOME" "$gold_targeted" "${NORM_CHECK_MODE[@]}"

    log "  normalizing called VCF"
    normalize_vcf "$CALLED_VCF" "$eval_bed" "$REF_GENOME" "$called_targeted" "${NORM_CHECK_MODE[@]}"

    log "  gold records:   $(bcftools view -H "$gold_targeted" | wc -l)"
    log "  called records: $(bcftools view -H "$called_targeted" | wc -l)"

    local -a exclude_args=(--exclude-filters)
    if [[ "${#EXCLUDE_FILTERS[@]}" -gt 0 ]]; then
        exclude_args+=("${EXCLUDE_FILTERS[@]}")
    fi

    local calc_status=0
    set +e
    python3 "${SCRIPT_DIR}/calculate_precision_recall.py" \
        --gold "$gold_targeted" \
        --called "$called_targeted" \
        --output "$results" \
        --assembly "$asm_label" \
        --caller "$CALLER" \
        "${exclude_args[@]}"
    calc_status=$?
    set -e

    if [[ "$calc_status" -ne 0 ]]; then
        if [[ -s "$results" ]]; then
            log "  calculator returned ${calc_status}; metrics were written, so continuing"
            calc_status=0
        else
            log "  calculator failed before writing results"
        fi
    fi
    log "  results: $results"
    return "$calc_status"
}

main() {
    require_tool bedtools
    require_tool bcftools
    require_tool python3

    mkdir -p "$OUT_DIR"

    IFS=',' read -r -a assembly_list <<< "$ASSEMBLIES"

    local status=0
    local first_results=""
    local aggregate="${OUT_DIR}/all_results.tsv"
    rm -f "$aggregate"

    for assembly in "${assembly_list[@]}"; do
        assembly="${assembly//[[:space:]]/}"
        [[ -z "$assembly" ]] && continue

        if ! evaluate_one "$assembly"; then
            status=1
        fi

        local asm_label="$assembly"
        [[ "$asm_label" == "chm13" ]] && asm_label="t2t"
        local results="${OUT_DIR}/${VCF_STAGE}/${asm_label}/benchmark_results_${asm_label}_${VCF_STAGE}.tsv"
        if [[ -s "$results" ]]; then
            if [[ -z "$first_results" ]]; then
                first_results="$results"
                head -n 1 "$results" > "$aggregate"
            fi
            tail -n +2 "$results" >> "$aggregate"
        fi
    done

    if [[ -s "$aggregate" ]]; then
        log "Aggregated results: $aggregate"
        if command -v column >/dev/null 2>&1; then
            column -t -s $'\t' "$aggregate"
        else
            cat "$aggregate"
        fi
    fi

    return "$status"
}

main "$@"
