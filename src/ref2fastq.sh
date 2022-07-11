#!/bin/bash

# This script extracts subsequences from reference genome with a given BED file. Read length is optional.

BED_FILE=$1
REF_GENOME=$2
READ_LENGTH=$3

source $(dirname $BASH_SOURCE)/miscellaneous.sh
source $(dirname $BASH_SOURCE)/errorHandling.sh
set -o errtrace
trap 'catch_exit_status $? $LINENO $0' ERR

# Input checking
if [ ! -f "${BED_FILE}" ] || [ ! -f "${REF_GENOME}" ]
then
    exit 2;
fi

# BED preprocessing
sort -V -k1,3 "${BED_FILE}" | bedtools merge -i - | sort -k1,1n -k2,2n -k3,3n > "${BED_FILE}.merged"

if [[ ! ${READ_LENGTH} =~ ^[0-9]+$ ]] && [[ ! -z ${READ_LENGTH} ]]
then
    READ_LENGTH=250
elif [[ -z "${READ_LENGTH}" ]]
then
    READ_LENGTH=250
fi

if [[ ${READ_LENGTH} -gt 100 ]]
then
    bedtools makewindows -b "${BED_FILE}"  -w "${READ_LENGTH}" | bedtools merge -d -1 | \
    awk -F '\t' '{printf "%s\t%s\n", $0, $3-$2;}' | \
    awk -F '\t' '$NF >= 100 && $NF < '${READ_LENGTH}'{print;} \
                    $NF == '${READ_LENGTH}'{line = $0; chr = $1; start = $2; end = $3; \
                                        getline; \
                                        if (($NF < 100) && ($2 == end) && ($1 == chr)) printf "%s\t%s\t%s\t%s\n", $1, start, $3, $3 - start; \
                                        else printf "%s\n%s\n", line, $0; \
                                        next;} \
                    $NF < 100{print;}' > "${BED_FILE/.bed/.preproc.bed}"
else
    echo -e >&2 "$(timestamp) WARNING: Read length (${READ_LENGTH}) is too small (<100)."
    cp -f "${BED_FILE}" "${BED_FILE/.bed/.preproc.bed}"
fi

if [[ $(awk -F '\t' '{sum+=$3-$2;} END{print sum;}' "${BED_FILE/.bed/.preproc.bed}") -eq $(awk -F '\t' '{sum+=$3-$2;} END{print sum;}' "${BED_FILE}") ]]
then
    seqtk subseq "${REF_GENOME}" "${BED_FILE/.bed/.preproc.bed}" | \
    seqtk seq -F "F" - > "${BED_FILE/.bed/.refseq}.${READ_LENGTH}.fastq"
else
    echo -e >&2 "$(timestamp) ERROR: Realignment failed for BED file : ${BED_FILE}"
    echo -e >&2 "$(timestamp) ERROR: Does the BED file have reasonable read length (${READ_LENGTH})?"
fi


