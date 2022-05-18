#!/bin/bash

# Prerequisites:
# - GATK
# - SAMtools
# - HTSlib
# - BCFtools

source $(dirname $(realpath $0))/errorHandling.sh
source $(dirname $(realpath $0))/miscellaneous.sh
SRC_PATH=$(dirname $(realpath $0))

[ $# -eq 0 ] && { $BASH_SOURCE --help; exit 2; }
POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"

    case "$key" in
        --file)
            FILE="$2"
            shift
            shift
            ;;
        --type)
            TYPE="$2"
            shift
            shift
            ;;
        --ref-genome)
            REF_GENOME="$2"
            shift
            shift
            ;;
        -e|--expected-subjects)
            EXPECTED_SUB="$2"
            shift
            shift
            ;;
        -h|--help)
            echo "File validator for BAM and VCF"
            echo "Usage:    $BASH_SOURCE"
            echo "          --file                      | path of file to be validated"
            echo "          --type                      | type of file (bam/vcf/fastq)"
            echo "          --ref-genome                | path of reference genome file (hg19)"
            echo "          -e | --expected-subjects    | a file list of expected subject IDs (only for vcf)"
            echo "          -h | --help                 | display this message"
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

# Treat different file types separately
if [[ ${TYPE,,} == "bam" ]]; then
    if [ ! -f "$FILE" ]; then exit 5; fi

    # Default minimum line number set to 10
    line_cnt=$(samtools view -c $FILE)

    if [[ ${line_cnt} -lt 10 ]]; then
        echo -e >&2 "$(timestamp) ERROR: Insufficient records (record # = "${line_cnt}" in BAM file."
        exit 5;
    fi

    # Check query name disagreement from DRAGMAP
    if [[ $(samtools view "$FILE" | cut -f1 | grep -v ^@ | egrep -c '/1$|/2$') -gt 0 ]]; then
        echo -e >&2 "$(timestamp) WARNING: BAM file contains unexpected query names. Trying to format ..."
        samtools view "${FILE}" | cut -f1 | grep -v ^@ | egrep '/1$|/2$' | head -10
        samtools view "${FILE}" | sed 's+/1\t+\t+g;s+/2\t+\t+g' | \
        samtools sort -O bam -o ${FILE%%.bam}_sorted.bam -@ ${THREADS:-4} - | \
        samtools index ${FILE%%.bam}_sorted.bam
        echo -e "$(timestamp) INFO: Recoded BAM file as ${FILE%%.bam}_sorted.bam"

        gatk ValidateSamFile \
            -I ${FILE%%.bam}_sorted.bam \
            --INDEX_VALIDATION_STRINGENCY EXHAUSTIVE \
            -MO 1000000 > ${FILE%%.bam}.check
    fi

    # Check if BAM file has been sorted; check its index accordingly
    if [[ $(samtools stats "$FILE" | grep "is sorted" | cut -f3) -eq 1 ]]; then
        if [[ "$FILE" -nt "${FILE%%.bam}.bai" ]]; then samtools index "$FILE"; fi
    else
        samtools sort -O bam -o ${FILE%%.bam}_sorted.bam -@ ${_THREADS:-4} "$FILE"
        echo -e "$(timestamp) INFO: Recoded BAM file as ${FILE%%.bam}_sorted.bam"
        samtools index ${FILE%%.bam}_sorted.bam

        gatk ValidateSamFile \
            -I ${FILE%%.bam}_sorted.bam \
            --INDEX_VALIDATION_STRINGENCY EXHAUSTIVE \
            -MO 1000000 > ${FILE%%.bam}.check
    fi

    # Analyze .check file
    if [ -f "${FILE%%.bam}.check" ]; then
        while IFS=: read -r level empty message description
        do
            if [[ "$message" == "INVALID_INDEX_FILE_POINTER" ]]; then
                echo -e >&2 "$(timestamp) ERROR: "${FILE%%.bam}.bai" may be corrupted."
                exit 5; # <- expected manual re-indexing here as some subtle errors may present in source BAM file
            elif [[ "$message" == "MATE_NOT_FOUND" ]]; then
                echo -e >&2 "$(timestamp) WARNING: Singleton reads found in $FILE"
            elif [[ "$message" == "INVALID_PLATFORM_VALUE" ]]; then
                echo -e >&2 "$(timestamp) WARNING: $FILE may be a simulated golden BAM file without valid PL value."
                echo -e >&2 "$(timestamp) WARNING: Snippets of error messages: - "
                echo -e >&2 "$level"::"$message":"$description"
            else
                echo -e >&2 "$(timestamp) ERROR: $FILE may be corrupted; no simple fix available."
                exit 5;
            fi
        done < <(grep ^ERROR "${FILE%%.bam}.check");
        # <- .check file not removed for debugging purposes
    fi

elif [[ ${TYPE,,} == "vcf" ]]; then
    if [ ! -f "$FILE" ]; then exit 6; fi

    [ -v REF_GENOME ] || { echo -e >&2 "$(timestamp) ERROR: No reference genome provided." && exit 3; }

    # Check if the input file is in gzip format
    if [[ 0x$(head -c 2 $FILE | xxd -p) -eq 0x1f8b ]] && gunzip -t $FILE; then
        # Case: GZIP input (Note: BGZF follows the same header format as GZIP)
        echo -e "$(timestamp) INFO: GZIP file detected."

        # Check index file metadata
        if [[ $FILE -nt ${FILE}.tbi ]]; then
        tabix -f -p vcf $FILE || { echo -e >&2 "$(timestamp) ERROR: Input VCF may not be sorted by chr pos. Try bcftools sort."; exit 6; }
        fi
        echo -e "$(timestamp) INFO: Index file not up-to-date. A new index file has been created."

        # Check if the file is empty (i.e. no subject)
        if [[ $(bcftools query -l $FILE | wc -l) -eq 0 ]]; then
            echo -e >&2 "$(timestamp) ERROR: No subjects found in input VCF!"; exit 129;
        fi

        # GATK validation
        gatk ValidateVariants -R $REF_GENOME -V $FILE --validation-type-to-exclude ALL --verbosity DEBUG || { echo -e >&2 "$(timestamp) ERROR: Input file $FILE failed GATK validation." && exit 129; }

    elif [[ "$FILE" =~ \.[bgz,gz] ]]; then
        echo -e >&2 "$(timestamp) ERROR: Input file has corrupted."
        exit 6;

    else # Case: Uncompressed VCF
        if [[ $(stat -c %s $FILE) -le 30 ]]; then
            echo -e >&2 "$(timestamp) ERROR: Input is NULL or corrupted (bgzf)."
            exit 4;
        else
            # Check if the file is empty (i.e. no subject)
            if [[ $(bcftools query -l $FILE | wc -l) -eq 0 ]]; then
                echo -e >&2 "$(timestamp) ERROR: No subjects found in input VCF!"; exit 129;
            fi

            # GATK validation
            gatk ValidateVariants -R $REF_GENOME -V $FILE --validation-type-to-exclude ALL --verbosity DEBUG || { echo -e >&2 "$(timestamp) ERROR: Input file $FILE failed GATK validation." && exit 129; }
        fi
    fi
elif [[ ${TYPE,,} == "fastq" ]]; then
    if [[ $(${SRC_PATH}/fqValidator -i $FILE -t fastq) -eq 0 ]]; then
        echo -e >&2 "$(timestamp) ERROR: Input FASTQ file $FILE is invalid. Possibly due to invalid line number."
        exit 4;
    fi
elif [ ! -v TYPE ]; then
    exit 3;
else # Unrecognized file type
    echo -e >&2 "$(timestamp) ERROR: Unrecognized file type (BAM / VCF / FASTQ)."
    exit 4;
fi

# Compare subject IDs if the source VCF is valid
if [ ! -v EXPECTED_SUB ]; then
    echo -e >&2 "$(timestamp) WARNING: Expected list of subject IDs not provided. Skipping subject ID comparison."
    exit 0;
fi

if [[ ${TYPE,,} == "vcf" ]]; then
    dest="$EXPECTED_SUB"
    echo -e "$(timestamp) INFO: The following subjects are found both in the VCF and the expected list."
    grep -wf <(bcftools query -l "$FILE") "$dest"
    echo -e "$(timestamp) INFO: The following subjects are NOT found in the VCF but found in the expected list."
    grep -vwf <(bcftools query -l "$FILE") "$dest"
fi
