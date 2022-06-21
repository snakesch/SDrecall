#!/bin/bash

# 4_buildSampleMap.sh
# Description: Extract XA tagged reads or reads with MQ < alpha; then create a masked genome for each region for downstream analysis.
# Note: This script is designed for multithread execution of mutliple regions.
# Author: Louis She (2022-06)
# Contact: louisshe@hku.hk

# Prerequisites:
# samtools v1.15
# BEDTools v2.27.1
# seqkit
# GATK
# bwa

source $(dirname $(realpath -s $0))/src/miscellaneous.sh
source $(dirname $(realpath -s $0))/src/errorHandling.sh

set -o errtrace
trap 'catch_exit_status $? $LINENO $0' ERR

[[ $# -eq 0 ]] && { $BASH_SOURCE --help; exit 2; }
while [[ $# -gt 0 ]]
do
    key="$1"
    case "$key" in
        --input-bam|-i)
            INPUT_BAM="$2"
            shift
            shift
            ;;
        --region-list|-l)
            REGION="$2"
            shift
            shift
            ;;
        --ref-genome|-r)
            REF_GENOME="$2"
            shift
            shift
            ;;
        --ref-bed|-rb)
            REF_BED="$2"
            shift
            shift
            ;;
        --mq-threshold|-mq)
            MQ_THRES="$2"
            shift
            shift
            ;;
        --thread|-t)
            NTHREADS="$2"
            shift
            shift
            ;;
        --out|-o)
            OUTPATH="$2"
            shift
            shift
            ;;
        --help|-h)
            echo "Extract multi-aligned readds from SD regions from BAM and prepare masked genomes"
            echo "Usage:    $BASH_SOURCE"
            echo "          --input-bam | -i        | path of indexed BAM file"
            echo "          --region-list | -l      | list of regions in BED format"
            echo "          --ref-genome | -r       | reference genome (hg19)"
            echo "          --ref-bed | -rb         | BED of reference genome"
            echo "          --mq-threshold | -mq    | threshold for mapping quality (MQ) (default: 30)"
            echo "          --out | -o              | output path"
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

# Input validation
[ -v INPUT_BAM ] || exit 2;
[ -v REGION ] || exit 2;
[ -v OUTPATH ] || exit 2;
[ -v REF_GENOME ] || exit 2;
[ -v REF_BED ] || exit 2;
[ -f "${INPUT_BAM}" ] || exit 3;
[ -f "${REGION}" ] || exit 3;
[ -f "${REF_GENOME}" ] || exit 3;
[ -f "${REF_BED}" ] || exit 3;

# Default thread num = 8
NTHREADS=${NTHREADS:-8}
MQ_THRES=${MQ_THRES:-30}
if [[ $( echo "$(samtools --version-only | cut -d. -f1-2) >= 1.15" | bc -l ) -eq 0 ]]; then
    echo -e >&2 "$(timestamp) ERROR: samtools version is not compatible. (Version>=1.15 is required)"
    exit 127
fi

if ! $(dirname $(realpath -s $0))/src/checkBAM.sh
then
    exit 5
fi

# Thread safety
REGION_PREFIX=$(basename $(basename $REGION .bed) _related_homo_regions)

# Sort and merge BED file by coordinates
bedtools sort -i "$REGION" | bedtools merge -i - | sort -k1,1 -V -s > ${OUTPATH}/${REGION_PREFIX}_cleaned.bed

# Extract SD regions from BAM
# Extract overlap (even if mate is not within the region)
samtools view -hb -@ ${NTHREADS} -P -L "$REGION" "$INPUT_BAM" > ${OUTPATH}/${REGION_PREFIX}_extracted.bam.tmp

# Extract XA reads OR MQ < MQ_THRES
samtools view -h -@ ${NTHREADS} ${OUTPATH}/${REGION_PREFIX}_extracted.bam.tmp | awk -F'\t' -v mq=${MQ_THRES} '($0 ~ /XA:Z/ || $5 < mq) {print $1}' | sort -uV | uniq > ${OUTPATH}/${REGION_PREFIX}_qnames.tmp
samtools view --qname-file ${OUTPATH}/${REGION_PREFIX}_qnames.tmp -h -b -@ ${NTHREADS} ${OUTPATH}/${REGION_PREFIX}_extracted.bam.tmp > ${OUTPATH}/${REGION_PREFIX}_extracted.bam
RECORDS=$(samtools view -c ${OUTPATH}/${REGION_PREFIX}_extracted.bam)
echo -e "$(timestamp) INFO: Total number of reads extracted : ${RECORDS}"
if [[ $RECORDS -eq 0 ]]
then
    echo -e >&2 "$(timestamp) ERROR: No reads extracted! (Criteria: XA / MQ < ${MQ_THRES})"
    echo -e >&2 "$(timestamp) ERROR: Hint: Check ${OUTPATH}/${REGION_PREFIX}_extracted.bam.tmp"
    exit 4
fi

# Sort BAM again by query names
samtools sort -n -o ${OUTPATH}/${REGION_PREFIX}_extracted_sorted.bam.tmp ${OUTPATH}/${REGION_PREFIX}_extracted.bam

# Generate paired-end FASTQ
[[ -d "${OUTPATH}/fastq/" ]] || mkdir -p "${OUTPATH}/fastq/"
samtools fastq -1 ${OUTPATH}/fastq/$(basename $INPUT_BAM | cut -d. -f1)_${REGION_PREFIX}-XA_1.fq -2 ${OUTPATH}/fastq/$(basename $INPUT_BAM | cut -d. -f1)_${REGION_PREFIX}-XA_2.fq -0 /dev/null -s /dev/null -n -@ 4 ${OUTPATH}/${REGION_PREFIX}_extracted_sorted.bam.tmp

gzip -f "${OUTPATH}/fastq/$(basename $INPUT_BAM | cut -d. -f1)_${REGION_PREFIX}-XA_1.fq"
gzip -f "${OUTPATH}/fastq/$(basename $INPUT_BAM | cut -d. -f1)_${REGION_PREFIX}-XA_2.fq"
echo -e "$(timestamp) INFO: Generated paired-end FASTQ file for ${REGION_PREFIX} for $(basename $INPUT_BAM | cut -d. -f1)"

if [[ $(zcat ${OUTPATH}/fastq/$(basename $INPUT_BAM | cut -d. -f1)_${REGION_PREFIX}-XA_1.fq.gz | wc -l) -lt 4 ]]
then
    echo -e >&2 "$(timestamp) WARNING: No multi-aligned reads for ${REGION_PREFIX} for $(basename $INPUT_BAM | cut -d. -f1)"
    rm -f ${OUTPATH}/fastq/$(basename $INPUT_BAM | cut -d. -f1)_${REGION_PREFIX}-XA_[1,2].fq
    rm -f ${OUTPATH}/${REGION_PREFIX}*.tmp
    exit 0
elif [[ $(zcat ${OUTPATH}/fastq/$(basename $INPUT_BAM | cut -d. -f1)_${REGION_PREFIX}-XA_2.fq.gz | wc -l) -lt 4 ]]
then
    echo -e >&2 "$(timestamp) WARNING: No multi-aligned reads for ${REGION_PREFIX} for $(basename $INPUT_BAM | cut -d. -f1)"
    rm -f ${OUTPATH}/fastq/$(basename $INPUT_BAM | cut -d. -f1)_${REGION_PREFIX}-XA_[1,2].fq
    rm -f ${OUTPATH}/${REGION_PREFIX}*.tmp
    exit 0
fi

# Check average depth of the region
depth=$(samtools depth -a -b "${REGION}" -g DUP,UNMAP,QCFAIL,SECONDARY "${INPUT_BAM}" -@ "${NTHREADS}" | awk -F'\t' '{sum+=$3} END{printf "%.2f", sum/NR; }')
if (( $(echo "$depth <= 0.0" | bc -l) ))
then
    rm -f ${OUTPATH}/fastq/$(basename $INPUT_BAM | cut -d. -f1)_${REGION_PREFIX}-XA_[1,2].fq
    rm -f ${OUTPATH}/${REGION_PREFIX}*.tmp
    echo -e >&2 "$(timestamp) WARNING: Zero coverage for ${REGION_PREFIX}."
    exit 0
else
    echo -e >&2 "$(timestamp) INFO: Coverage for ${REGION_PREFIX} : ${depth}"
fi

# prepare masked genome
[[ -d "${OUTPATH}/masked_genome/" ]] || mkdir -p "${OUTPATH}/masked_genome/"

$(dirname $(realpath -s $0))/src/buildMaskedGenome.sh "${REF_GENOME}" "${REF_BED}" "${REGION}" "${REGION_PREFIX}" "${OUTPATH}/masked_genome/"

# Cleanup
rm -f ${OUTPATH}/${REGION_PREFIX}_extracted.bam ${OUTPATH}/${REGION_PREFIX}_cleaned.bed
rm -f ${OUTPATH}/${REGION_PREFIX}*.tmp



