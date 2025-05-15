#! /usr/bin/env bash

SELF_PATH=$(realpath ${BASH_SOURCE[0]})
EXAMPLE_DIR=$(dirname $SELF_PATH)
PROJECT_DIR=$(dirname $EXAMPLE_DIR)
DATA_DIR="$PROJECT_DIR/data"

# Default values
OUTPUT_DIR="$(dirname $PROJECT_DIR)/test_output"
REF_GENOME="Null"
THREADS=10

# Parse command line arguments
usage() {
  echo "Usage: $0 [-o OUTPUT_DIR] [-r REF_GENOME] [-t THREADS]"
  echo "  -o  Output directory (default: $OUTPUT_DIR)"
  echo "  -r  Reference genome path (default: $REF_GENOME, the contig names must have the prefix 'chr')"
  echo "  -t  Number of threads (default: $THREADS)"
  echo "  -h  Show this help message"
  exit 1
}

while getopts "o:r:t:h" opt; do
  case $opt in
    o) OUTPUT_DIR="$OPTARG" ;;
    r) REF_GENOME="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h) usage ;;
    \?) usage ;;
  esac
done

[[ -d $OUTPUT_DIR ]] || { mkdir -p $OUTPUT_DIR || { echo "Cannot create output directory $OUTPUT_DIR"; exit 1 ; } }
[[ -f $REF_GENOME ]] || { echo "Reference genome $REF_GENOME does not exist"; exit 1 ; }

# Display the parameters being used
echo "Using parameters:"
echo "  Output directory: $OUTPUT_DIR"
echo "  Reference genome: $REF_GENOME"
echo "  Threads: $THREADS"
echo "Please wait around 15 min for the example run to complete, the real time progress log is $OUTPUT_DIR/test_SDrecall_CMRG_HG002.log"

# Run the command with the specified parameters
nohup $PROJECT_DIR/SDrecall \
--verbose DEBUG run \
-i $EXAMPLE_DIR/HG002.CMRG.hg38.test.bam \
-r $REF_GENOME \
-m $DATA_DIR/hg38/ref_SD/WGAC.hg38.cigar.trim.homo.expanded.highsim.bed.gz \
-o $OUTPUT_DIR \
-t $THREADS \
-b $EXAMPLE_DIR/HG002_hg38_CMRG_smallvar_v1.00.bed \
--ref_genome_tag hg38 \
--conventional_vcf $EXAMPLE_DIR/HG002.hg38.CMRG.deepvariant.vcf.gz \
--caller_name DeepVariant \
--target_tag CMRGTEST > $OUTPUT_DIR/test_SDrecall_CMRG_HG002.log 2>&1 && \
echo "Test completed successfully, now you can inspect the result VCF: $EXAMPLE_DIR/HG002.hg38.CMRG.deepvariant.sdrecall.CMRGTEST.vcf.gz, intermediate files in $OUTPUT_DIR/HG002_hg38_CMRGTEST_SDrecall/" && \
echo "Note that the added variants from SDrecall to the conventional callset are the variants only in the region specified by $OUTPUT_DIR/HG002_hg38_CMRGTEST_SDrecall/realign_groups/all_target_recall_SD_regions.bed" && \
echo "This region set is an intersection of the target region and the regions suffered from multi-alignments (covered by less than 10 high-MAPQ (MAPQ > 10) reads) according to the input BAM file" || \
{ echo "Test failed"; exit 1; }