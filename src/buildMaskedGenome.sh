#!/bin/bash

REF_GENOME=$1
TOTAL_BED=$2
REGION_BED=$3
PRIORITY_COMPONENT=$4
OUTPATH=$5

# Prerequisites:
# - BEDTools
# - seqkit
# - seqtk
# - GATK
# - bwa

source $(dirname $(realpath -s $0))/miscellaneous.sh
source $(dirname $(realpath -s $0))/errorHandling.sh
set -o errtrace
trap 'catch_exit_status $? $LINENO $0' ERR

REF_GENOME_NAME=$(basename ${REF_GENOME} .fasta)

if [[ -f "${OUTPATH}/${REF_GENOME_NAME}_${PRIORITY_COMPONENT}_masked.fasta" ]]
then
    validate_mgenome "${OUTPATH}/${REF_GENOME_NAME}_${PRIORITY_COMPONENT}_masked.fasta" "${REGION_BED}"
    echo -e >&2 "$(timestamp) WARNING: Masked genome already exists for ${PRIORITY_COMPONENT}"
    exit 0
fi

# - complementary BED
bedtools subtract -a ${TOTAL_BED} -b ${REGION_BED} > "${OUTPATH}/$(basename ${TOTAL_BED} .bed)_exclude_${PRIORITY_COMPONENT}.bed.tmp"

# - mask genome (main)
awk '{print $1;}' "${REGION_BED}" | sort -uV | uniq > "${OUTPATH}/$(basename ${REGION_BED} .bed).bed.tmp"
seqtk subseq -l 60 "${REF_GENOME}" "${OUTPATH}/$(basename ${REGION_BED} .bed).bed.tmp" > "${OUTPATH}/${REF_GENOME_NAME}_${PRIORITY_COMPONENT}_masked.fasta.tmp"
seqtk seq -l 60 -M "${OUTPATH}/$(basename ${TOTAL_BED} .bed)_exclude_${PRIORITY_COMPONENT}.bed.tmp" -n N "${OUTPATH}/${REF_GENOME_NAME}_${PRIORITY_COMPONENT}_masked.fasta.tmp" > "${OUTPATH}/${REF_GENOME_NAME}_${PRIORITY_COMPONENT}_masked.fasta"

rm -f ${OUTPATH}/$(basename ${TOTAL_BED} .bed)_exclude_${PRIORITY_COMPONENT}.bed.tmp ${OUTPATH}/$(basename ${REGION_BED} .bed).bed.tmp ${OUTPATH}/${REF_GENOME_NAME}_${PRIORITY_COMPONENT}_masked.fasta.tmp
index_genome "${OUTPATH}/${REF_GENOME_NAME}_${PRIORITY_COMPONENT}_masked.fasta" dict fai
validate_mgenome "${OUTPATH}/${REF_GENOME_NAME}_${PRIORITY_COMPONENT}_masked.fasta" "${REGION_BED}"
echo -e "$(timestamp) INFO: Masked genome created successfully. ("${OUTPATH}/${REF_GENOME_NAME}_${PRIORITY_COMPONENT}_masked.fasta")"
