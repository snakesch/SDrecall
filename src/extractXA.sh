#!/bin/bash

# This script extracts multi-mapped reads from BQSR BAM and writes extracted XA-tagged reads to FASTQ format.

# TODO:
# - check BAM validity
# - BAM to FASTQ conversion (extract XA tagged reads only)
# - Trim FASTQ files
# - Check for query name consistency

# Prerequisites:
# samtools
# bedtools
# Trim Galore!
# Java
# FastQC
# Cutadapt

source $(dirname $(realpath $0))/errorHandling.sh
source $(dirname $(realpath $0))/miscellaneous.sh

set -o errtrace
trap 'catch_exit_status $? $LINENO $0' ERR

SRC_PATH=$(dirname $(realpath $0))
TMP_DIR=$(dirname $(realpath $0))/tmp
[ -d "${TMP_DIR}" ] || mkdir -p "${TMP_DIR}"
LOG_PATH=${SRC_PATH}/../log/${0}.log
_THREADS=${_THREADS:-8} # TODO: Export _THREADS in alternative script

exec &> >(tee -a "${LOG_PATH}")
exec |& tee -a "${LOG_PATH}"

[ $# -eq 0 ] && { $BASH_SOURCE --help; exit 2; }
POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"

    case "$key" in
        --region-dir|-d)
            BED_REGION_DIR="$2"
            shift
            shift
            ;;
        --file|-f)
            FILE="$2"
            shift
            shift
            ;;
        --ref-genome|-r)
            REF_GENOME="$2"
            shift
            shift
            ;;
        --out|-o)
            OUTPATH="$2"
            shift
            shift
            ;;
        #--fwd|-1)
        #    FQ_FWD="$2"
        #    shift
        #    shift
        #    ;;
        #--rev|-2)
        #    FQ_REV="$2"
        #    shift
        #    shift
        #    ;;
        --help|-h)
            echo "Utility to extract XA tagged reads from BAM and writes to FASTQ format"
            echo "Usage:    $BASH_SOURCE"
            echo "          --region-dir | -d               | path of directory of BED files of homologous regions / regions to mask"
            echo "          --file | -f                     | BAM file to extract XA tagged reads"
            echo "          --ref-genome | -r               | path of reference genome (hg19)"
            echo "          --out | -o                      | output path"
            # echo "          --fwd | -1                      | path of forward FASTQ"
            # echo "          --rev | -2                      | path of reverse FASTQ"
            echo "          --help | -h                     | display this message"
            exit 0
            ;;
        *)
            echo "Unknown option encountered - $1"
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done
set -- "${POSITIONAL[@]}"

# Check if all required inputs are given
[[ -v REF_GENOME ]] || REF_GENOME="/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta"
[[ -v BED_REGION_DIR ]] || { echo -e >&2 "$(timestamp) ERROR: No mask region provided!" && exit 3; }
if [[ ! -v FILE ]] && [[ ! -v FQ_FWD ]]; then exit 3; fi

echo -e >&2 "$(timestamp) INFO: Executing on ${_THREADS} threads. More threads may be optimal for WGS data."

# Importing BED files of homologous regions
declare -a homolog_regions=( $(find ${BED_REGION_DIR} -mindepth 1 -type f -name "*_related_homo_regions.bed" -exec basename {} _related_homo_regions.bed \;) )
echo -e "$(timestamp) INFO: Loaded ${#homolog_regions[@]} homologous regions."

# prepare_map_file_per_region in GATK4.1 script.
# Original call:
# parallel --joblog ${bamf}/parallel_prepare_map_file.${sample_ID}.log --link --dry-run -j${threads} bash ${scrf/staging/paedyl01\/disk1\/yangyxt}/GATK4.             1_based_NGS_preprocessing_script.sh prepare_map_file_per_region {1} ${input_bam}      ${tmp_map_file} '>' ${bamf}/${sample_ID}.bqsr.{2}.genmap.log '2>&1' :::               "${homologous_regions[@]}" ::: "${hr_names[@]}"



# Case: < Input = BAM >
if [[ $FILE =~ *".bam" ]]; then

    # Check BAM validity
    ${SRC_PATH}/fileValidator.sh --file ${FILE} --type bam --ref-genome ${REF_GENOME}
    echo -e "$(timestamp) INFO: $(samtools view -c $FILE) reads detected in the input BAM file."
    SAMP_ID=$(samtools view -H $FILE | grep ^@RG | tr '\t' '\n' | grep ^SM | cut -d: -f2)

    # Merge BED files & extract exonic regions
    cat ${BED_REGION_DIR}/*.bed | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i stdin > ${TMP_DIR}/union.bed
    bedtools intersect -wo -a ${TMP_DIR}/union.bed -b ${EXONIC_REGIONS} > ${TMP_DIR}/union_exonic.bed

    # Extract exonic homologous regions & XA tagged reads
    samtools view -b -L ${TMP_DIR}/union_exonic.bed $FILE > ${TMP_DIR}/overlappedReads.bam

    # Cleanup to free storage
    rm -f ${TMP_DIR}/union.bed
    rm -f ${TMP_DIR}/union_exonic.bed

    if [ -z $(samtools view ${TMP_DIR}/overlappedReads.bam) ]; then
        echo -e >&2 "$(timestammp) ERROR: No reads found in homologous regions. Abort"
        rm -f ${TMP_DIR}/overlappedReads.bam
        exit 129;
    fi

    echo -e "$(timestamp) INFO: $(samtools view -c ${TMP_DIR}/overlappedReads.bam) reads fall within exonic homologous regions."
    samtools view -b -e '[XA]' ${TMP_DIR}/overlappedReads.bam | samtools sort - > ${TMP_DIR}/extractedReads.bam
    samtools index ${TMP_DIR}/extractedReads.bam

    # Cleanup to free storage
    rm -f ${TMP_DIR}/overlappedReads.bam
    if [ -z $(samtools view ${TMP_DIR}/extractedReads.bam) ]; then
        echo -e >&2 "$(timestamp) ERROR: No XA tags detected. Abort."
        rm -f ${TMP_DIR}/extractedReads.bam
        exit 129;
    fi

    echo -e "$(timestamp) INFO: $(samtools view -c ${TMP_DIR}/extractedReads.bam) XA tagged reads fall into exonic homologous regions."

    bedtools bamtofastq -i ${TMP_DIR}/extractedReads.bam -fq ${TMP_DIR}/${SAMP_ID}_1.fastq -fq2 ${TMP_DIR}/${SAMP_ID}_2.fastq
    gzip -c ${TMP_DIR}/${SAMP_ID}_1.fastq && gzip -c ${TMP_DIR}/${SAMP_ID}_2.fastq

    echo -e "$(timestamp) DEBUG: FASTQ files generated"

    trim_galore --paired --no_report_file --phred33 --fastqc_args "--outdir=${TMP_DIR}" -o ${TMP_DIR} ${TMP_DIR}/${SAMP_ID}_1.fastq.gz ${TMP_DIR}/${SAMP_ID}_2.fastq.gz

    # Only check query name consistency
    # Trim FASTQ and validation (DRAGMAP nomenclature; use seqtk)
    # < Not tested >
    zcat ${TMP_DIR}/${SAMP_ID}_1_val_1.*.gz | sed 's+/[1,2]++g' | gzip -c > ${OUTPATH}/${SAMP_ID}_1.fastq.gz
    zcat ${TMP_DIR}/${SAMP_ID}_2_val_2.*.gz | sed 's+/[1,2]++g' | gzip -c > ${OUTPATH}/${SAMP_ID}_2.fastq.gz

# Case: < Input = FASTQ > # Use FQ_FWD, FQ_REV
elif [ -v FQ_FWD ] && [ -v FQ_REV ]; then

    # Validate input FASTQ files
    # Extract XA reads & reads fall within exonic regions
    # Trim FASTQ
    echo -e >&2 "$(timestamp) ERROR: Not implemented.";
    exit 129;

# Case: < Input = Unknown file type >
else
    echo -e >&2 "$(timestamp) ERROR: Unknown file type - $FILE"
    exit 4;
fi

# Cleanup
set -- "${homolog_regions[@]}"

