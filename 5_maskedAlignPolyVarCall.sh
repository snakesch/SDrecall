#!/bin/bash

# 5_maskedAlignPolyVarCall.sh
# Description: This script realigns BAM to masked genome, call for polyploidy variants if ploidy >= 2.
# Author: Yang XT, She CH (2022)

# Input handler
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
        --data|-d)
            DATAPATH="$2"
            shift
            shift
            ;;
        --region-bed|-b)
            REGION_BED="$2"
            shift
            shift
            ;;
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
        --help|-h)
            echo "Masked alignment and polyploidy variant calling"
            echo "Usage:    $BASH_SOURCE"
            echo "          --input-bam | -i        | input BQSR BAM file"
            echo "          --data | -d             | parent directory of masked_genome and fastq"
            echo "          --region-bed | -b       | BED file of selected SD regions"
            echo "          --ref-genome | -r       | reference genome"
            echo "          --thread | -t           | number of threads to use"
            echo "          --help | -h             | display this message and exit"
            exit 0
            ;;
        *)
            echo "Unknown option encountered - $1"
            shift
            ;;
    esac
done
NTHREADS=${NTHREADS:-10}
REGION=$(basename $(basename ${REGION_BED} .bed) _related_homo_regions)
declare -a fastq=($(find ${DATAPATH}/fastq/ -name "*${REGION}*.fq.gz" -type f | sort))
declare -a mgenome=($(find ${DATAPATH}/masked_genome/ -name "$(basename ${REF_GENOME} .fasta)_${REGION}_masked.fasta" -type f | sort))
if [[ "${#mgenome[@]}" -gt 1 ]]
then
    echo -e >&2 "$(timestamp) ERROR: More than 1 masked genome found for ${REGION}!"
    exit 129
elif [[ "${#mgenome[@]}" -eq 0 ]]
then
    echo -e >&2 "$(timestamp) WARNING: No masked genome found for ${REGION}! (Probably due to 0 coverage in input BAM file.)"
    exit 0
fi
[[ ${#fastq[@]} -eq 2 ]] || exit 3;

aligned_out=${INPUT_BAM%%.*}.only_${REGION}.bam
samp_ID=$(basename ${fastq[0]} | cut -d_ -f1)

[[ -f "${REF_GENOME%%.fasta}.dict" ]] || index_genome ${REF_GENOME} dict;

# Masked alignment
ori_depth=$(getDepth $INPUT_BAM $REGION_BED $NTHREADS)
ori_hq_depth=$(getDepth $INPUT_BAM $REGION_BED $NTHREADS HQ)
echo -e "$(timestamp) INFO: Average depth = ${ori_depth}"
echo -e "$(timestamp) INFO: Average depth (XA removed) = ${ori_hq_depth}"
echo -e "$(timestamp) INFO: Running BWA MEM on $NTHREADS threads for ${REGION}"
time bwa mem -t "${NTHREADS}" -M -R "@RG\tID:${samp_ID}\tLB:SureSelectXT Library Prep Kit\tPL:ILLUMINA\tPU:1064\tSM:${samp_ID}" "${mgenome}" "${fastq[0]}" "${fastq[1]}" | samtools view -uSh -@ "${NTHREADS}" - | samtools sort -O bam -@ "${NTHREADS}" -o "${aligned_out}"
samtools index "${aligned_out}"
echo -e "$(timestamp) INFO: BWA MEM completed for ${REGION}"
XA_depth=$(getDepth ${aligned_out} $REGION_BED $NTHREADS)
post_depth=$( echo "$XA_depth-$ori_hq_depth" | bc -l )
echo -e "$(timestamp) DEBUG: XA_depth = $XA_depth"
echo -e "$(timestamp) DEBUG: post_depth = $post_depth"
echo -e "$(timestamp) INFO: Read depth of ${REGION} region has increased by $( echo "$post_depth/$ori_depth" | bc -l ) fold."
echo -e "$(timestamp) INFO: Estimated ploidy: $(echo "scale=1;2*$post_depth/$ori_depth" | bc -l )"
PLOIDY=$(echo "scale=0;2*$post_depth/$ori_depth" | bc -l)

if [[ "${PLOIDY}" -ge 2 ]]
then
    # Call polyploid variants
    [[ -d "${DATAPATH}/vcf" ]] || mkdir -p "${DATAPATH}/vcf"
    echo -e "$(timestamp) INFO: Calling gVCFs from masked alignment ..."
    time gatk --java-options "-Xmx75G" HaplotypeCaller \
        --emit-ref-confidence GVCF \
        -R "${REF_GENOME}" \
        --ploidy "${PLOIDY}" \
        -I "${aligned_out}" \
        -O "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.HC.g.vcf.gz" \
        --bamout "${INPUT_BAM%%.*}.only_${REGION}_realigned.bam" \
        --min-base-quality-score 30 \
        --force-active \
        --allow-non-unique-kmers-in-ref \
        --debug-assembly \
        --kmer-size 21 \
        --max-num-haplotypes-in-population 256 \
        --bam-writer-type CALLED_HAPLOTYPES \
        --linked-de-bruijn-graph \
        -L "${REGION_BED}" \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -G StandardHCAnnotation \
        -imr OVERLAPPING_ONLY \
        -OBI true \
        --verbosity WARNING

    if [[ $(zgrep -v ^# "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.HC.g.vcf.gz" | wc -l) -eq 0 ]]
    then
        echo -e >&2 "$(timestamp) WARNING: No variants called by GATK HaplotypeCaller for ${REGION}"
        rm -f "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.HC.g.vcf.gz"
        exit 0
    fi

    gatk GenotypeGVCFs \
        -R "${REF_GENOME}" \
        -V "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.HC.g.vcf.gz" \
        -O "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz" \
        -G StandardAnnotation \
        -G AS_StandardAnnotation

    if ! isValidVCF "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz"
    then
        echo -e >&2 "$(timestamp) ERROR: GATK GenotypeGVCFs generated invalid VCF file ${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz"
        rm -f ${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz
        exit 129
    elif [[ $(zgrep -v ^# ${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz | wc -l) -eq 0 ]]
    then
        echo -e >&2 "$(timestamp) WARNING: No variant is called by GATK GenotypeGVCFs for ${REGION}."
        rm -f ${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz
        exit 0
    fi

    gatk LeftAlignAndTrimVariants \
        -R "${REF_GENOME}" \
        -V "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz" \
        -O "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.tmp.vcf.gz" \
        --max-leading-bases 2000 \
        --max-indel-length 2000 \
        --split-multi-allelics \
        --dont-trim-alleles

    rm -f "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz"
    zcat "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.tmp.vcf.gz" | awk -F'\t' '{if ($5 == "*") next; else print;}' | bgzip -c > "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz"
    rm -f "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.tmp.vcf.gz"
    tabix -f -p vcf "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz"

    [[ -f "${REF_GENOME}.fai" ]] || samtools faidx "${REF_GENOME}";
    awk -F'\t' '{printf "##contig=<ID=%s,length=%d>\n", $1, $2}' "${REF_GENOME}.fai" > "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.contig_header.tmp"
    zgrep "^##" "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz" | grep -v "^##contig=" > "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.other_header.tmp"
    zgrep "^#CHROM" "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz" > "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.chrom_header.tmp"
    
    cat "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.other_header.tmp" "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.contig_header.tmp" "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.chrom_header.tmp" > "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.new_header.tmp"
    rm -f "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.other_header.tmp" "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.contig_header.tmp" "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.chrom_header.tmp"
    
    bcftools reheader -h "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.new_header.tmp" -o "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.reheader.vcf.gz" "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz"
    rm -f "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz" "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.new_header.tmp" && mv "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.reheader.vcf.gz" "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz"

    python3 $(dirname $(realpath -s $0))/src/convert_vcf_to_diploidy.py -vp "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz" || echo -e >&2 "$(timestamp) ERROR: Failed to convert multiploid VCF to diploid VCF";

    pick_poorcov_region_bam -i "${INPUT_BAM}" -t "${REGION_BED}" -o "${INPUT_BAM/.bam}.${REGION}.txt" -d 15 && echo -e "$(timestamp) INFO: The poor coverage region in ${INPUT_BAM} overlapping with ${REGION_BED} is stored in "${INPUT_BAM/.bam}.${REGION}.txt
    bcftools view -R "${INPUT_BAM/.bam}.${REGION}.txt" -Oz -o "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.tmp.vcf.gz" "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz"
    bcftools sort -Oz "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.tmp.vcf.gz" -o "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz"
    rm -f "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.tmp.vcf.gz"
    if isValidVCF "${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz"; then
        echo -e "$(timestamp): INFO: In function ${FUNCNAME}: Picked out variants in poor coverage regions: ${DATAPATH}/vcf/${samp_ID}.only_${REGION}.vcf.gz"
    else
        >&2 echo -e "$(timestamp): ERROR: In function ${FUNCNAME}: Error in picking out variants in poor coverage regions: ${REGION}."
        exit 1
    fi
fi
