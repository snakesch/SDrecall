#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage:
  run_sex_gt_ploidy_haplotag.sh \
    --assembly hg19|hg38 \
    --ploidy-class haploid|diploid \
    --vcf defrabb.vcf.gz \
    --raw-bam HG002.300x.bam \
    --reference ref.fa \
    --out-dir out/hg19/sex_haploid \
    --parquet-out-dir out/parquet_by_ploidy/hg19.sex_haploid.haplotags.parquet

Runs WhatsHap haplotag on whole sex chromosomes using only the GT ploidy of
the VCF records to choose the WhatsHap ploidy:
  haploid records: GT has no '/' or '|' separator, run with --ploidy 1
  diploid records: GT has '/' or '|' separator, run with --ploidy 2

No PAR/non-PAR BED is used in this workflow.
USAGE
}

assembly=""
ploidy_class=""
vcf=""
raw_bam=""
reference=""
out_dir=""
parquet_out_dir=""
sample="HG002"
threads=25
ignore_read_groups=1
skip_existing=0

parquet_helper="/paedyl01/disk1/yangyxt/SDrecall-rust-migration/scripts/haplotag_tsv_to_parquet.py"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --assembly) assembly="$2"; shift 2 ;;
        --ploidy-class) ploidy_class="$2"; shift 2 ;;
        --vcf) vcf="$2"; shift 2 ;;
        --raw-bam) raw_bam="$2"; shift 2 ;;
        --reference) reference="$2"; shift 2 ;;
        --out-dir) out_dir="$2"; shift 2 ;;
        --parquet-out-dir) parquet_out_dir="$2"; shift 2 ;;
        --sample) sample="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        --no-ignore-read-groups) ignore_read_groups=0; shift ;;
        --skip-existing) skip_existing=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

require_file() {
    local path="$1"
    local label="$2"
    if [[ -z "${path}" ]]; then
        echo "Missing required option: ${label}" >&2
        exit 2
    fi
    if [[ ! -s "${path}" ]]; then
        echo "Required file is missing or empty (${label}): ${path}" >&2
        exit 2
    fi
}

require_cmd() {
    local cmd="$1"
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "Required command not found in PATH: ${cmd}" >&2
        exit 127
    fi
}

case "${assembly}" in
    hg19)
        sex_regions=(X Y)
        sex_region_csv="X,Y"
        ;;
    hg38)
        sex_regions=(chrX chrY)
        sex_region_csv="chrX,chrY"
        ;;
    *)
        echo "--assembly must be hg19 or hg38" >&2
        exit 2
        ;;
esac

case "${ploidy_class}" in
    haploid)
        ploidy=1
        gt_filter='GT!~"[|/]"'
        ;;
    diploid)
        ploidy=2
        gt_filter='GT~"[|/]"'
        ;;
    *)
        echo "--ploidy-class must be haploid or diploid" >&2
        exit 2
        ;;
esac

require_file "${vcf}" "--vcf"
require_file "${vcf}.tbi" "VCF tabix index"
require_file "${raw_bam}" "--raw-bam"
require_file "${raw_bam}.bai" "BAM index"
require_file "${reference}" "--reference"
require_file "${reference}.fai" "reference .fai"
require_file "${parquet_helper}" "haplotag_tsv_to_parquet.py"
if [[ -z "${out_dir}" || -z "${parquet_out_dir}" ]]; then
    echo "Missing one of --out-dir, --parquet-out-dir" >&2
    exit 2
fi

require_cmd bcftools
require_cmd whatshap
require_cmd samtools
require_cmd python3

mkdir -p "${out_dir}"/{logs,vcf,whatshap} "${parquet_out_dir}"

run_vcf="${out_dir}/vcf/${assembly}.sex_${ploidy_class}.vcf.gz"
tagged_bam="${out_dir}/whatshap/${assembly}.sex_${ploidy_class}.haplotagged.bam"
haplotags="${out_dir}/whatshap/${assembly}.sex_${ploidy_class}.haplotags.tsv"
regions_txt="${out_dir}/sex_regions.txt"

printf '%s\n' "${sex_regions[@]}" > "${regions_txt}"

if [[ "${skip_existing}" -eq 1 && -s "${run_vcf}" && -s "${run_vcf}.tbi" ]]; then
    echo "[skip] sex ${ploidy_class} VCF exists: ${run_vcf}" >&2
else
    echo "[run] bcftools sex ${ploidy_class} GT filter: ${gt_filter}" >&2
    bcftools view \
        -r "${sex_region_csv}" \
        -i "${gt_filter}" \
        -Oz \
        -o "${run_vcf}" \
        "${vcf}" \
        2>&1 | tee "${out_dir}/logs/00_filter_vcf.${ploidy_class}.log"
    bcftools index -f -t "${run_vcf}"
fi

variant_count="$(bcftools view -H "${run_vcf}" | wc -l | awk '{print $1}')"
echo "[info] ${assembly} sex_${ploidy_class} variants=${variant_count}" >&2
if [[ "${variant_count}" -eq 0 ]]; then
    echo "No ${ploidy_class} records found on sex chromosomes for ${assembly}" >&2
    exit 0
fi

haplotag_args=(
    haplotag
    --reference "${reference}"
    --output-haplotag-list "${haplotags}"
    -o "${tagged_bam}"
    --output-threads "${threads}"
    --sample "${sample}"
    --ploidy "${ploidy}"
    --skip-missing-contigs
)
for region in "${sex_regions[@]}"; do
    haplotag_args+=(--regions "${region}")
done
if [[ "${ignore_read_groups}" -eq 1 ]]; then
    haplotag_args+=(--ignore-read-groups)
fi
haplotag_args+=("${run_vcf}" "${raw_bam}")

if [[ "${skip_existing}" -eq 1 && -s "${tagged_bam}" && -s "${haplotags}" ]]; then
    echo "[skip] whatshap haplotag ${assembly}.sex_${ploidy_class}" >&2
else
    echo "[run] whatshap haplotag ${assembly}.sex_${ploidy_class} ploidy=${ploidy}" >&2
    echo "[regions] ${regions_txt}" >&2
    whatshap "${haplotag_args[@]}" 2>&1 | tee "${out_dir}/logs/01_whatshap_haplotag.log"
fi

if [[ -s "${tagged_bam}" && ! -e "${tagged_bam}.bai" ]]; then
    samtools index -@ "${threads}" "${tagged_bam}"
fi

if [[ "${skip_existing}" -eq 1 && -s "${parquet_out_dir}/_SUCCESS.${assembly}" ]]; then
    echo "[skip] parquet conversion ${assembly}.sex_${ploidy_class}" >&2
else
    python3 "${parquet_helper}" \
        --input "${haplotags}" \
        --assembly "${assembly}" \
        --out-dir "${parquet_out_dir}" \
        --allow-existing \
        --batch-rows 2000000 \
        --progress-batches 10 \
        2>&1 | tee "${out_dir}/logs/02_haplotags_to_parquet.log"
fi

cat > "${out_dir}/manifest.tsv" <<EOF
key	value
assembly	${assembly}
ploidy_class	${ploidy_class}
ploidy	${ploidy}
gt_filter	${gt_filter}
sex_regions	${sex_region_csv}
vcf	${vcf}
run_vcf	${run_vcf}
raw_bam	${raw_bam}
reference	${reference}
regions_txt	${regions_txt}
tagged_bam	${tagged_bam}
haplotags	${haplotags}
parquet_out_dir	${parquet_out_dir}
sample	${sample}
threads	${threads}
ignore_read_groups	${ignore_read_groups}
variant_count	${variant_count}
EOF

echo "Sex GT-ploidy haplotag complete: ${out_dir}/manifest.tsv" >&2
