#!/bin/bash

# 6_postProcessing.sh
# Description: This script compares gVCFs generated from 5_maskedAlignPolyVarCall.sh to an "intrinsic VCF" derived from the reference genome (hg19).
# Prerequisites:
# - bcftools
# - samtools
# - vcfpy
# - HTSlib

# Input handler
source $(dirname $(realpath -s $0))/src/miscellaneous.sh
source $(dirname $(realpath -s $0))/src/errorHandling.sh

set -o errtrace
trap 'catch_exit_status $? $LINENO $0' ERR
[[ $# -eq 0 ]] && { $BASH_SOURCE --help; exit 2; }
while [[ $# -gt 0 ]]
do
    key="$1"
    case "$key" in
        --vcfpath|-vp)
            VCFPATH="$2"
            shift
            shift
            ;;
        --regions|-r)
            REGIONS_PATH="$2"
            shift
            shift
            ;;
        --help|-h)
            echo "Postprocessing: Compare with intrinsic VCF"
            echo "Usage:    $BASH_SOURCE"
            echo "          --vcfpath | -p         | directory of VCFs in poor coverage regions"
            echo "          --regions | -r         | directory of all BED files of homologous regions"
            echo "          --help | -h            | display this message and exit"
            exit 0
            ;;
        *)
            echo "Unknown option encountered - $1"
            shift
            ;;
    esac
done

all_regions=($(find $REGIONS_PATH -name "*.bed" -type f -exec basename {} .bed \; | sed 's/_related_homo_regions//g'))
vcf_regions=($(find $VCFPATH -name "*.vcf.gz" -type f -exec basename {} .vcf.gz \; | grep -v g$ | cut -d_ -f2-))
valid_files=($(find $VCFPATH -name "*.vcf.gz" -type f | sort -uV | uniq))
samp_ID=$(basename ${valid_files[0]} | cut -d. -f1)

VALID=(${vcf_regions[@]})
INVALID=($(echo ${all_regions[@]} ${vcf_regions[@]} | tr ' ' '\n' | sort | uniq -u ))

echo -e "$(timestamp) INFO: Regions with valid VCF: ${VALID[@]}"
echo -e "$(timestamp) INFO: Regions with invalid VCF: ${INVALID[@]}"

find $VCFPATH -name "*.vcf.gz" -type f | grep -v "HC.g" | sort -uV | uniq > ${VCFPATH}/${samp_ID}.homo_region.tmp
bcftools concat -a --no-version -Ov -f ${VCFPATH}/${samp_ID}.homo_region.tmp -o ${VCFPATH}/${samp_ID}.homo_region.vcf
rm -f ${VCFPATH}/${samp_ID}.homo_region.tmp
bcftools sort -Oz ${VCFPATH}/${samp_ID}.homo_region.vcf -o ${VCFPATH}/${samp_ID}.homo_region.vcf.gz
if [ -f "${VCFPATH}/${samp_ID}.homo_region.vcf" ]; then rm -f "${VCFPATH}/${samp_ID}.homo_region.vcf"; fi
if ! isValidVCF ${VCFPATH}/${samp_ID}.homo_region.vcf.gz
then
    echo -e >&2 "$(timestamp) ERROR: Failed to concatenate VCFs of homologous regions."
    exit 1
fi

python3 $(dirname $(dirname $VCFPATH))/src/compare_with_intrinsic_vcf.py -qv ${VCFPATH}/${samp_ID}.homo_region.vcf.gz -ov ${VCFPATH}/${samp_ID}.homo_region.filtered.vcf.gz -iv "$(dirname $(dirname $VCFPATH))/ref/intrinsic_diff_calls.trimmed.bial.vcf.gz"
bcftools index ${VCFPATH}/${samp_ID}.homo_region.filtered.vcf.gz

# Cleanup (remove masked genomes & intermediate VCFs)
if [ -d $(dirname $VCFPATH)/masked_genome ]; then
    rm -rf $(dirname $VCFPATH)/masked_genome
fi

rm -f "${valid_files[@]}"
