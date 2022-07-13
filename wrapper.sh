#!/bin/bash

# wrapper.sh
# Description: This script is a wrapper script of SDrecall.
# Author: Yang XT, She CH (2022)

[[ $# -eq 0 ]] && { $BASH_SOURCE --help; exit 2; }
while [[ $# -gt 0 ]]
do
    key="$1"
    case "$key" in
        --input-bam|-i)
            BAM_PATH="$2"
            shift
            shift
            ;;
        --ref-bed|-rb)
            REF_BED="$2"
            shift
            shift
            ;;
        --ref-genome|-rg)
            REF_GENOME="$2"
            shift
            shift
            ;;
        --anno-ref|-a)
            ANNO_REF="$2"
            shift
            shift
            ;;
        --out|-o)
            OUTPATH="$2"
            shift
            shift
            ;;
        --fraglen|-f)
            FRAGLEN="$2"
            shift
            shift
            ;;
        --gaplen|-g)
            GAPLEN="$2"
            shift
            shift
            ;;
        --mq|-mq)
            MQ_THRES="$2"
            shift
            shift
            ;;
        --log)
            LOGLEVEL="$2"
            shift
            shift
            ;;
        --thread|-t)
            NTHREADS="$2"
            shift
            shift
            ;;
        --gene-list|-l)
            GENE_LIST="$2"
            shift
            shift
            ;;
        --read-length|-rl)
            READ_LENGTH="$2"
            shift
            shift
            ;;
        --help|-h)
            echo "A wrapper script of SDrecall"
            echo "Usage:    $BASH_SOURCE"
            echo -e "\nRequired:"
            echo "          --input-bam|-i          | input BQSR BAM file"
            echo "          --ref-bed|-rb           | BED file of reference genome"
            echo "          --ref-genome|-rg        | reference genome"
            echo "          --anno-ref|-a           | annotation table"
            echo "          --out|-o                | output directory (default: ./out)"
            echo -e "\nOptional:"
            echo "          --fraglen|-f            | fragment length FRAGLEN in CIGAR processing (default: 300)"
            echo "          --gaplen|-g             | small gap cutoff GAPLEN in CIGAR processing (default: 10)"
            echo "          --mq|-mq                | MQ threshold for extracting multi-aligned reads (default: 30)"
            echo "          --thread|-t             | number of threads (default: 8)"
            echo "          --log                   | log level (default: INFO)"
            echo "          --gene-list|-l          | list of genes of interest"
            echo "          --read-length | -rl     | window size for homologous regions"
            exit 0
            ;;
        *)
            echo "Unknown option encountered - $1"
            shift
            ;;
    esac
done

# Environment
[ -n "${BAM_PATH}" ] || exit 2;
[ -n "${REF_BED}" ] || exit 2;
[ -n "${REF_GENOME}" ] || exit 2;
[ -n "${ANNO_REF}" ] || exit 2;
ROOT=$(dirname $0)
NTHREADS=${NTHREADS:-8}
FRAGLEN=${FRAGLEN:-300}
GAPLEN=${GAPLEN:-10}
MQ_THRES=${MQ_THRES:-30}
OUTPATH=${OUTPATH:-${ROOT}/out}
LOGLEVEL=${LOGLEVEL:-INFO}
READ_LENGTH=${READ_LENGTH:-250}
export NTHREADS
source "${ROOT}/src/miscellaneous.sh"
source "${ROOT}/src/errorHandling.sh"
set -o errtrace
trap 'catch_exit_status $? $LINENO $0' ERR

# Helper function
function getThread () {
    local expected_thread=$1
    if [[ "$expected_thread" -gt "$NTHREADS" ]]
    then
        echo "$NTHREADS" && return 0
    else
        echo "$expected_thread" && return 0
    fi
}

# Create directories if needed
[[ -d "${ROOT}/out/fastq" ]] || mkdir -p "${ROOT}/out/fastq"
[[ -d "${ROOT}/ref/homologous_regions" ]] || mkdir -p "${ROOT}/ref/homologous_regions"
[[ -d "${ROOT}/ref/principal_components" ]] || mkdir -p "${ROOT}/ref/principal_components"

# Main
# $ROOT/1_biserFetch.sh --ref-genome "${REF_GENOME}" --out "${OUTPATH}" --thread $(getThread 8)

# if [ -f "${OUTPATH}/SD_hg19.trimmed.bed" ]
# then
#    rm -f "${OUTPATH}/SD_hg19.trimmed.bed"
# fi

# $ROOT/2_trimCIGAR.py -i "${OUTPATH}/SD_hg19.bed" -o "${OUTPATH}/SD_hg19.trimmed.bed" -v "${LOGLEVEL}" -f "${FRAGLEN}" -g "${GAPLEN}"

if [ -f "${GENE_LIST}" ]
then
    $ROOT/3_annotateExtract.py --input "${OUTPATH}/SD_hg19.trimmed.bed" --ref "${ANNO_REF}" --list "${GENE_LIST}" --verbose "${LOGLEVEL}" --out "${ROOT}/ref/"
else
    $ROOT/3_annotateExtract.py --input "${OUTPATH}/SD_hg19.trimmed.bed" --ref "${ANNO_REF}" --verbose "${LOGLEVEL}" --out "${ROOT}/ref"
    echo -e >&2 "$(timestamp) WARNING: Extracting homologous regions from ALL genes. This may take excessive time."
fi

find ${ROOT}/ref/homologous_regions/ -name "*.bed" | xargs -I{} -t -P3 bash ${ROOT}/4_preparation.sh --input-bam "${BAM_PATH}" --thread $(getThread 8) --region-list {} --out "${OUTPATH}" --ref-genome "${REF_GENOME}" --ref-bed "${REF_BED}" -mq "${MQ_THRES}"

for region in $(find ${ROOT}/ref/homologous_regions -name "*.bed" -type f | sort)
do
    $ROOT/5_maskedAlignPolyVarCall.sh --input-bam "${BAM_PATH}" --data "${OUTPATH}" --region-bed "${region}" --ref-genome "${REF_GENOME}" --thread $(getThread 14)
done

# Make intrinsic VCF
$ROOT/6_makeIntrinsicVCF.sh --bed-dir "${ROOT}/ref" --ref-genome "${REF_GENOME}" --ref-bed "${REF_BED}" --read-length "${READ_LENGTH}" --thread $(getThread 8)

$ROOT/7_postProcessing.sh --vcfpath "${OUTPATH}/vcf/" --regions "${ROOT}/ref/homologous_regions/" -iv "${ROOT}/ref/all_pc.realign.${READ_LENGTH}.trim.vcf.gz"

echo -e " ############# Part 2 done ############# "
# Deep cleaning ...
