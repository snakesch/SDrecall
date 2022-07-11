#!/bin/bash

# 6_makeIntrinsicVCF.sh
# Description: Call intrinsic variants for comparison in post-processing step.
# Author: Yang XT, She CH (2022)

source $(dirname $0)/src/miscellaneous.sh
source $(dirname $0)/src/errorHandling.sh
set -o errtrace
trap 'catch_exit_status $? $LINENO $0' ERR

[[ $# -eq 0 ]] && { $BASH_SOURCE --help; exit 2; }
while [[ $# -gt 0 ]]
do
    key="$1"
    case "$key" in
        --bed-dir|-b)
            BED_DIR="$2"
            shift
            shift
            ;;
        --ref-genome|-r)
            REF_GENOME="$2"
            shift
            shift
            ;;
        --ref-bed|-rb)
            TOTAL_BED="$2"
            shift
            shift
            ;;
        --thread|-t)
            NTHREADS="$2"
            shift
            shift
            ;;
        --read-length|-l)
            READ_LENGTH="$2"
            shift
            shift
            ;;
        --help|-h)
            echo "Extract multi-aligned reads from SD regions from BAM and prepare masked genomes"
            echo "Usage:    $BASH_SOURCE"
            echo "          --bed-dir | -b          | parent directory of homologous_regions/ and principal components/"
            echo "          --ref-genome | -r       | reference genome (hg19)"
            echo "          --ref-bed | -rb         | BED of reference genome"
            echo "          --read-length | -l      | window size for homologous regions"
            echo "          --thread | -t           | number of threads"
            echo "          --help | -h             | display this message and exit"
            exit 0
            ;;
        *)
            echo "Unknown option encountered - $1"
            shift
            ;;
    esac
done

# Create file lists
find "${BED_DIR}/homologous_regions" -name "*_related_homo_regions.bed" -type f > "${BED_DIR}/homologous_regions/all_sd_beds.tmp"
find "${BED_DIR}/principal_components" -name "*.bed" -type f > "${BED_DIR}/principal_components/all_pc_bed.tmp"

cat ${BED_DIR}/homologous_regions/all_sd_beds.tmp | xargs -n 1 -P${NTHREADS} -I{} $(dirname $0)/src/ref2fastq.sh {} ${REF_GENOME} ${READ_LENGTH}

_sd_regions=($(find ${BED_DIR}/homologous_regions -name "*.${READ_LENGTH}.fastq" -type f))
echo -e "$(timestamp) INFO: Extracted ${#_sd_regions[@]} FASTQs from reference genome."

[ -f "${BED_DIR}/principal_components/all_pc.reseq.${READ_LENGTH}.fastq" ] && rm -f "${BED_DIR}/principal_components/all_pc.reseq.${READ_LENGTH}.fastq"

for sd_region in "${_sd_regions[@]}"
do
    cat ${sd_region} >> "${BED_DIR}/principal_components/all_pc.reseq.${READ_LENGTH}.fastq"
done

_arr=($(cat "${BED_DIR}/principal_components/all_pc_bed.tmp"))

for pc_region in ${_arr[@]}
do
    cat ${pc_region} >> "${BED_DIR}/principal_components/all_pc.bed.tmp"
done
sort -V -k1,3 "${BED_DIR}/principal_components/all_pc.bed.tmp" | bedtools merge -i - > "${BED_DIR}/principal_components/all_pc.bed"

rm -f "${BED_DIR}/principal_components/all_pc.bed.tmp"

# Prepare masked genome
$(dirname $0)/src/buildMaskedGenome.sh "${REF_GENOME}" "${TOTAL_BED}" "${BED_DIR}/principal_components/all_pc.bed" all_pc "${BED_DIR}/masked_genome/"

# Masked alignment
[ -d "${BED_DIR}/masked_genome" ] || mkdir -p "${BED_DIR}/masked_genome"
time bwa mem -t "${NTHREADS}" -M -R "@RG\tID:all_pc\tLB:SureSelectXT Library Prep Kit\tPL:ILLUMINA\tPU:1064\tSM:all_pc" "${BED_DIR}/masked_genome/ucsc.hg19_all_pc_masked.fasta" "${BED_DIR}/principal_components/all_pc.reseq.${READ_LENGTH}.fastq" | samtools view -uSh -@ "${NTHREADS}" - | samtools sort -O bam -@ "${NTHREADS}" -o "${BED_DIR}/all_pc.realign.${READ_LENGTH}.bam" -
samtools index "${BED_DIR}/all_pc.realign.${READ_LENGTH}.bam"

# Variant calling
if [[ "${#_sd_regions[@]}" -gt 1 ]]
then
    bcftools mpileup -a FORMAT/AD,FORMAT/DP -f "${REF_GENOME}" "${BED_DIR}/all_pc.realign.${READ_LENGTH}.bam" | bcftools call -mv -Oz -o "${BED_DIR}/all_pc.realign.${READ_LENGTH}.vcf.gz"
    tabix -f -p vcf "${BED_DIR}/all_pc.realign.${READ_LENGTH}.vcf.gz"
    gatk LeftAlignAndTrimVariants \
    -V "${BED_DIR}/all_pc.realign.${READ_LENGTH}.vcf.gz" \
    -R "${REF_GENOME}" \
    --split-multi-allelics \
    -O "${BED_DIR}/all_pc.realign.${READ_LENGTH}.trim.vcf.gz"
else
    echo -e >&2 "$(timestamp) ERROR: Unable to extract FASTQs from reference genome (>1 is needed)."
    exit 1
fi

# Cleanup
rm -f "${BED_DIR}/homologous_regions/all_sd_beds.tmp"
rm -f "${BED_DIR}/principal_components/all_pc_bed.tmp"
rm -f "${BED_DIR}/principal_components/all_pc.bed"
rm -rf "${BED_DIR}/masked_genome"
rm -f "${BED_DIR}/all_pc.realign.${READ_LENGTH}.bam*"
rm -f "${BED_DIR}/all_pc.realign.${READ_LENGTH}.vcf.gz"
rm -f "${BED_DIR}/principal_components/all_pc.bed"
