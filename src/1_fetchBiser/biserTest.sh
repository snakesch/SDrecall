#!/bin/bash
#PBS -N biserTestRun
#PBS -l mem=110g,nodes=paedyl01:ppn=10,walltime=999:00:00
#PBS -q paedyl
#PBS -m abe
#PBS -M louisshe@hku.hk
#PBS -o biserTestRun.log
#PBS -j oe
#PBS -V

REF_GENOME="/home/louisshe/shortVariantVCF/ref/ucsc.hg19.fasta"
OUTPATH="/home/louisshe/ret_Xingtian"
TOOL="/home/louisshe/.local/bin"

module load samtools

samtools faidx ${REF_GENOME}
${TOOL}/biser -o ${OUTPATH}/SD_hg19.bed -t 8 ${REF_GENOME}

module unload samtools
