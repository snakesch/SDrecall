#!/bin/bash

REF_GENOME=$1
TOTAL_BED=$2
REGION_BED=$3
PRIORITY_COMPONENT=$4
OUTPATH=$5

REF_GENOME="/home/louisshe/shortVariantVCF/ref/ucsc.hg19.fasta"
TOTAL_BED="/home/louisshe/shortVariantVCF/ref/ucsc.hg19.bed"
REGION_BED="/home/louisshe/shortVariantVCF/ref/homologous_regions/TERT_related_homo_regions.bed"
PRIORITY_COMPONENT="TERT"
OUTPATH="/home/louisshe/shortVariantVCF/out/masked_genome/"

# Prerequisites:
# - BEDTools
# - seqkit

source $(dirname $(realpath -s $0))/miscellaneous.sh
source $(dirname $(realpath -s $0))/errorHandling.sh
set -o errtrace
trap 'catch_exit_status $? $LINENO $0' ERR

[[ -f "${REF_GENOME}" ]] || exit 3;
REF_GENOME_NAME=$(basename ${REF_GENOME} .fasta)
if [[ -f "${OUTPATH}/${REF_GENOME_NAME}_${PRIORITY_COMPONENT}_masked.fasta" ]]
then
    # TODOs: Check file validity
    echo -e >&2 "$(timestamp) WARNING: Masked genome already exists for ${PRIORITY_COMPONENT}"
    exit 0
fi

# TODOs:
# - complementary BED
bedtools subtract -a ${TOTAL_BED} -b ${REGION_BED} > "${OUTPATH}/$(basename ${TOTAL_BED} .bed)_exclude_${PRIORITY_COMPONENT}.bed.tmp"

# - mask genome (main)
awk '{print $1;}' "${REGION_BED}" | sort -uV | uniq > "${OUTPATH}/$(basename ${REGION_BED} .bed).bed.tmp"
seqkit subseq -l 60 "${REF_GENOME}" "${OUTPATH}/$(basename ${REGION_BED} .bed).bed.tmp" > "${OUTPATH}/${REF_GENOME_NAME}_${PRIORITY_COMPONENT}_masked.fasta.tmp"
seqkit seq -l 60 -M "${OUTPATH}/$(basename ${TOTAL_BED} .bed)_exclude_${PRIORITY_COMPONENT}.bed.tmp" -n N "${OUTPATH}/${REF_GENOME_NAME}_${PRIORITY_COMPONENT}_masked.fasta.tmp" > "${OUTPATH}/${REF_GENOME_NAME}_${PRIORITY_COMPONENT}_masked.fasta"

rm -f ${OUTPATH}/$(basename ${TOTAL_BED} .bed)_exclude_${PRIORITY_COMPONENT}.bed.tmp ${OUTPATH}/$(basename ${REGION_BED} .bed).bed.tmp ${OUTPATH}/${REF_GENOME_NAME}_${PRIORITY_COMPONENT}_masked.fasta.tmp

# - check masked genome
# - build index for genome
# - validate masked genome
