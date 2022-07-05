#!/bin/bash

# Quick check for BAM file.

# Usage: ./checkBAM.sh [-t THREAD] FILE [FILES]
# Return value: 0 if OK, 1 otherwise
# if $( ./checkBAM.sh BAM ); then <case when the file(s) is OK>; fi

source $(dirname $(realpath -s $0))/miscellaneous.sh
source $(dirname $(realpath -s $0))/errorHandling.sh

if [[ "$1" == "-t" ]]
then
    NTHREADS="$2"
    shift 2
fi

NTHREADS=${NTHREADS:-4}

while [[ $# -gt 0 ]]
do
    FILE="$1"
    shift

    if [[ ! -f "$FILE" ]]; then exit 1; fi

    if [[ $(samtools view -c -@ ${NTHREADS} $FILE) -lt 10 ]]
    then
        echo -e >&2 "$(timestamp) ERROR: Insufficient records."
        exit 1
    fi

    if [[ $(samtools view -@ ${NTHREADS} "$FILE" | head -n5 | cut -f1 | grep -v ^@ | egrep -c '/1$|/2$') -gt 0 ]]
    then
        echo -e >&2 "$(timestamp) WARNING: BAM file contains unexpected query names."
     fi

    if [[ $(samtools stats -@ ${NTHREADS} "$FILE" | grep "is sorted" | cut -f3) -eq 0 ]]
    then
        echo -e >&2 "$(timestamp) ERROR: BAM is not sorted."
        exit 1
    fi

    if [[ "$FILE" -nt "${FILE}.bai" ]]
    then
        echo -e >&2 "$(timestamp) ERROR: BAM index file is out-of-date."
        exit 1
    fi
done
