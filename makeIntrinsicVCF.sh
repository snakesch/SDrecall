#!/bin/bash

BED_DIR="$(dirname $0)/out"
REF_GENOME="$(dirname $0)/ref/ucsc.hg19.fasta"
TOTAL_BED="$(dirname $0)/ref/ucsc.hg19.bed"
NTHREADS=12
READ_LENGTH=250

source $(dirname $0)/src/miscellaneous.sh
source $(dirname $0)/src/errorHandling.sh
set -o errtrace
trap 'catch_exit_status $? $LINENO $0' ERR

# Create file lists
find "${BED_DIR}/homologous_regions" -name "*_related_homo_regions.bed" -type f > "${BED_DIR}/homologous_regions/all_sd_beds.tmp"
find "${BED_DIR}/principal_components" -name "*.bed" -type f > "${BED_DIR}/principal_components/all_pc_bed.tmp"

export ${REF_GENOME}
export ${READ_LENGTH}
cat ${BED_DIR}/homologous_regions/all_sd_beds.tmp | xargs -n 1 -P${NTHREADS} -I{} bash -c '$(dirname $0)/src/ref2fastq.sh {} ${REF_GENOME} ${READ_LENGTH}'

_num=$(find ${BED_DIR}/homologous_regions -name "*.${READ_LENGTH}.fastq" -type f | wc -l)
echo -e "$(timestamp) INFO: Extracted ${_num} from reference genome."
# Output = "${BED_FILE/.bed/.refseq}.${READ_LENGTH}.fastq"

cat "${BED_DIR}/principal_components/*.bed" | sort -V -k1,3 | bedtools merge -i - > "${BED_DIR}/principal_components/all_pc.bed"

# Prepare masked genome
$(dirname $0)/src/buildMaskedGenome.sh "${REF_GENOME}" "${TOTAL_BED}" "${BED_DIR}/principal_components/all_pc.bed" all_pc "${BED_DIR}/masked_genome/"

# Masked alignment
time bwa mem -t "${NTHREADS}" -M -R "@RG\tID:all_pc\tLB:SureSelectXT Library Prep Kit\tPL:ILLUMINA\tPU:1064\tSM:all_pc" "${BED_DIR}/masked_genome/ucsc.hg19_all_pc_masked.fasta" | samtools view -uSh -@ "${NTHREADS}" - | samtools sort -O bam -@ "${NTHREADS}" -o "${BED_DIR}/all_pc.realign.${READ_LENGTH}.bam" -
samtools index "${BED_DIR}/all_pc.realign.${READ_LENGTH}.bam"

echo -e "Done til bcftools mpileup!!!"
exit 0
# Variant calling
if [[ "${_num}" -gt 1 ]]
then
    bcftools mpileup -a FORMAT/AD,FORMAT/DP

fi

# Cleanup
rm -f "${BED_DIR}/homologous_regions/all_sd_beds.tmp" "${BED_DIR}/principal_components/all_pc_bed.tmp"
rm -f "${BED_DIR}/principal_components/all_pc.bed"