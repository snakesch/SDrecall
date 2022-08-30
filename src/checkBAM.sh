#!/bin/bash

# Quick check for BAM file.

# Usage: ./checkBAM.sh [-t THREAD] FILE [FILES]
# Return value: 0 if OK, 1 otherwise
# if $( ./checkBAM.sh BAM ); then <case when the file(s) is OK>; fi

source $(dirname $(realpath -s $0))/miscellaneous.sh
source $(dirname $(realpath -s $0))/errorHandling.sh

set -x

while [[ $# -gt 0 ]]
do
    if [[ "$1" == "-t" ]]
    then
        NTHREADS="$2"
        shift 2
    fi
    NTHREADS=${NTHREADS:-4}

    FILE="$1"
    shift

    if [[ ! -f "$FILE" ]]; then exit 1; fi

    if [[ $(samtools view -c -@ ${NTHREADS} $FILE) -lt 10 ]]
    then
        echo -e >&2 "$(timestamp) ERROR: Insufficient records."
        exit 1
    fi

    if samtools head -n 5 "$FILE" | cut -f1 | egrep -c '/1$|/2$' > /dev/null
    then
        echo -e >&2 "$(timestamp) WARNING: BAM file contains unexpected query names."
    fi

    #if [[ $(samtools stats -@ ${NTHREADS} "$FILE" | grep "is sorted" | cut -f3) -eq 0 ]]
    order=$(samtools view -H "$FILE" | grep HD | cut -d: -f3)
    if [ ${order} == "unknown" ] || [ ${order} == "unsorted" ]
    then
        echo -e >&2 "$(timestamp) ERROR: BAM is not sorted."
        exit 1
    fi

    if [[ "$FILE" -nt "${FILE}.bai" ]]
    then
        echo -e >&2 "$(timestamp) ERROR: BAM index file is out-of-date."
        samtools index "${FILE}"
    fi
done
