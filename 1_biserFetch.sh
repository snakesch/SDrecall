#!/bin/bash

# 1_biserFetch.sh
# Description: This script extracts segmental duplication (SD) regions from given genome.
# Author: Louis She (2022-05)
# Contact: louisshe@hku.hk

# Prerequisites:
# samtools
# biser

source $(dirname $(realpath -s $0))/miscellaneous.sh
source $(dirname $(realpath -s $0))/errorHandling.sh

set -o errtrace
trap 'catch_exit_status $? $LINENO $0' ERR

[[ $# -eq 0 ]] && { $BASH_SOURCE --help; exit 2; }
while [[ $# -gt 0 ]]
do
    key="$1"
    case "$key" in
        --ref-genome|-r)
            REF_GENOME="$2"
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
            echo "Extract segmental duplication (SD) regions by BISER."
            echo "Usage:    $BASH_SOURCE"
            echo "          --ref-genome | -r          | reference genome (fasta)"
            echo "          --out | -o                 | output BED file path"
            echo "          --thread | -t              | number of threads to use"
            echo "          --help | -h                | display this message and exit"
            exit 0
            ;;
        *)
            echo "Unknown option encountered - $1"
            shift
            ;;
    esac
done

NTHREADS=${NTHREADS:-4}
samtools faidx ${REF_GENOME}
biser -o ${OUTPATH}/SD_hg19.bed -t ${NTHREADS} ${REF_GENOME}

