#!/bin/bash

# 4_buildSampleMap.sh
# Description:
# Note: This script is designed for multithread execution of mutliple regions.
# Author: Louis She (2022-06)
# Contact: louisshe@hku.hk

# Prerequisites:
# samtools v1.15
# BEDTools v2.27.1

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
            echo "Build a sample map (TBC)"
            echo "Usage:    $BASH_SOURCE"
            echo "          --input-bam | -i        | path of indexed BAM file"
            echo "          --region-list | -l      | list of regions in BED format"
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

# Default thread num = 8
NTHREADS=${NTHREADS:-8}
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
## 1. Extract overlap (even if mate is not within the region)
samtools view -hb -@ ${NTHREADS} -P -L "$REGION" "$INPUT_BAM" > ${OUTPATH}/${REGION_PREFIX}_extracted.bam.tmp

## 2. Extract XA reads OR MQ < 30
samtools view -h -@ ${NTHREADS} ${OUTPATH}/${REGION_PREFIX}_extracted.bam.tmp | awk -F'\t' '($0 ~ /XA:Z/ || $5 < 30) {print $0}' | samtools view -b > ${OUTPATH}/${REGION_PREFIX}_extracted.bam
RECORDS=$(samtools view -c ${OUTPATH}/${REGION_PREFIX}_extracted.bam)
echo -e "$(timestamp) INFO: Total number of reads extracted : ${RECORDS}"
if [[ $RECORDS -eq 0 ]]
then
    echo -e >&2 "$(timestamp) ERROR: No reads extracted! (Criteria: XA / MQ < 30)"
    echo -e >&2 "$(timestamp) ERROR: Hint: Check ${OUTPATH}/${REGION_PREFIX}_extracted.bam.tmp"
    exit 4
fi

## 3. Sort BAM again by query names
samtools sort -n -o ${OUTPATH}/${REGION_PREFIX}_extracted_sorted.bam.tmp ${OUTPATH}/${REGION_PREFIX}_extracted.bam

# Generate paired-end FASTQ
[[ -d "${OUTPATH}/fastq/" ]] || mkdir -p "${OUTPATH}/fastq/"
bedtools bamtofastq -i ${OUTPATH}/${REGION_PREFIX}_extracted_sorted.bam.tmp \
                    -fq ${OUTPATH}/fastq/$(basename $INPUT_BAM | cut -d. -f1)_${REGION_PREFIX}-XA_1.fq \
                    -fq2 ${OUTPATH}/fastq/$(basename $INPUT_BAM | cut -d. -f1)_${REGION_PREFIX}-XA_2.fq \
                    1>/dev/null
gzip -f ${OUTPATH}/fastq/$(basename $INPUT_BAM | cut -d. -f1)_${REGION_PREFIX}-XA_1.fq
gzip -f ${OUTPATH}/fastq/$(basename $INPUT_BAM | cut -d. -f1)_${REGION_PREFIX}-XA_2.fq
echo -e "$(timestamp) INFO: Generated paired-end FASTQ file for ${REGION_PREFIX} for $(basename $INPUT_BAM | cut -d. -f1)"

if [[ $(cat ${OUTPATH}/fastq/$(basename $INPUT_BAM | cut -d. -f1)_${REGION_PREFIX}-XA_1.fq | wc -l) -lt 4 ]]
then
    echo -e >&2 "$(timestamp) WARNING: No multi-aligned reads for ${REGION_PREFIX} for $(basename $INPUT_BAM | cut -d. -f1)"
    rm -f ${OUTPATH}/fastq/$(basename $INPUT_BAM | cut -d. -f1)_${REGION_PREFIX}-XA_[1,2].fq
    rm -f ${OUTPATH}/*.tmp
    exit 0
elif [[ $(cat ${OUTPATH}/fastq/$(basename $INPUT_BAM | cut -d. -f1)_${REGION_PREFIX}-XA_2.fq) -lt 4 ]]
then
    echo -e >&2 "$(timestamp) WARNING: No multi-aligned reads for ${REGION_PREFIX} for $(basename $INPUT_BAM | cut -d. -f1)"
    rm -f ${OUTPATH}/fastq/$(basename $INPUT_BAM | cut -d. -f1)_${REGION_PREFIX}-XA_[1,2].fq
    rm -f ${OUTPATH}/*.tmp
    exit 0
fi

# Check average depth of the region
depth=$(samtools depth -a -b "${REGION}" -g DUP,UNMAP,QCFAIL,SECONDARY "${INPUT_BAM}" -@ "${NTHREADS}" | awk -F'\t' '{sum+=$3} END{printf "%.2f", sum/NR; }')
if [[ "$depth" -le 0 ]]
then
    rm -f ${OUTPATH}/fastq/$(basename $INPUT_BAM | cut -d. -f1)_${REGION_PREFIX}-XA_[1,2].fq
    rm -f ${OUTPATH}/*.tmp
    echo -e >&2 "$(timestamp) WARNING: Zero coverage."
    exit 0
else
    echo -e >&2 "$(timestamp) INFO: Coverage for ${REGION_PREFIX} : ${depth}"
fi

# 2. prepare masked genome (in realign_masked_and_HC_multiploidy.sh prepare_masked_genome)
[[ -d "${OUTPATH}/masked_genome/" ]] || mkdir -p "${OUTPATH}/masked_genome/"


# 3. writes a result table showing file names of all intermediate outputs (writes header outside this script)
if [[ ! -f "${OUTPATH}/file_map.tsv" ]]
then
    echo -e "BED_file\tPriority_component\tgVCF_file\tfq_fwd\tfq_rev" > "${OUTPATH}/file_map.tsv"
fi

# Write a new record to file_map.tsv

# Cleanup
rm -f ${OUTPATH}/${REGION_PREFIX}_extracted.bam ${OUTPATH}/${REGION_PREFIX}_cleaned.bed
rm -f ${OUTPATH}/*.tmp



