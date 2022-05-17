#!/bin/bash

# This script extracts multi-mapped reads from BAM and writes extracted XA-tagged reads to FASTQ format.

# TODO:
# - check BAM validity
# - BAM to FASTQ conversion (extract XA tagged reads only)
# - Trim FASTQ files
# - Check for query name consistency

# Prerequisites:
#
#
#

source $(dirname $(realpath $0))/errorHandling.sh
source $(dirname $(realpath $0))/miscellaneous.sh

set -o errtrace
trap 'catch_exit_status $? $LINENO $0' ERR

SRC_PATH=$(dirname $(realpath $0))
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
[[ -v FILE ]] || exit 3;
[[ -v BED_REGION_DIR ]] || { echo -e >&2 "$(timestamp) ERROR: No mask region provided!" && exit 3; }

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

    # Prepare map file for each homologous region
    # Iterate over homolog_regions
    # Goal: Extract XA reads from BAM --> trimmed FASTQ (fwd + rev)
    # Func: prepare_fastq_files (only expect BAM input here)
    # Original call:
    # ./select_bam_by_regions.sh -i <bam> 0b
    # Initialize variables: - bam_ID = get_RG_SM
    # Merge BED files (intersectBed @ BedTools ) <-- writes to ~/tmp
    # Extract XA
    # Trim FASTQ
    # FASTQ validation - query name consistency

# Case: < Input = FASTQ >
elif [[ $FILE =~ *".fastq" ]] || [[ $FILE =~ *".fq" ]]; then
    echo -e >&2 "$(timestamp) ERROR: Not implemented.";
    exit 129;

# Case: < Input = Unknown file type >
else
    echo -e >&2 "$(timestamp) ERROR: Unknown file type - $FILE"
    exit 4;
fi

# Cleanup
set -- "${homolog_regions[@]}"

