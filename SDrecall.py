#!/usr/bin/env python3

# SDrecall.py
# Description: A wrapper script that implements major steps in SDrecall.
# Author: Yang XT, She CH (2022)
# Contact: Xingtian Yang (u3005579@connect.hku.hk), Louis She (snakesch@connect.hku.hk)

import argparse
import subprocess
import os
import logging
import glob
import time
from src.utils import executeCmd
from src.geneAnnotation import *

def main():

    logger = logging.getLogger("root")
    start = time.time()
    
    # Check if we are in the correct directory
    if "1_getTrimmedSD.py" not in os.listdir():
        raise FileNotFoundError(f"Source scripts are not found. ")

    # Step 0: Get gene annotation table
#     if not os.path.exists(os.path.join(ROOT, 'ref', 'refgene_latest_anno.bed')):
#         report_coor(os.path.join(ROOT, 'ref'), genePanel_fp=args.list)
#         logger.info(f"****************** Step 0 completed in {(time.time() - start) / 60} minutes. ******************")
#     else:
#         logger.info(f"Found gene element coordinates BED in {os.path.join(ROOT, 'ref')}")

#     start = time.time()
#     # Step 1: Get homologous region and principal component coordinates
#     cmd = f"python3 1_getTrimmedSD.py --ref_genome {args.ref_genome} " \
#           f"--output {os.path.dirname(args.ref_genome)} --build {args.build} --list {args.list} " \
#           f"--table {os.path.join('ref', 'refgene_latest_anno.bed')} --thread {args.thread} --fraglen {args.fraglen} " \
#           f"--gaplen {args.gaplen} --verbose {args.verbose} "
#     if args.keep_trimmed:
#         cmd += "--keep_trimmed "
#     executeCmd(cmd)
#     logger.info(f"****************** Step 1 completed in {(time.time() - start) / 60} minutes. ******************")

#     start = time.time()
#     # Step 2: Preparation step for masked alignment
#     cmd = f"python3 2_preparation.py --input_bam {args.input_bam} --homo_dir {os.path.join(ROOT, 'ref/homologous_regions')} " \
#           f"--pc_dir {os.path.join(ROOT, 'ref/principal_components')} " \
#           f"--outpath {args.output} " \
#           f"--ref_genome {args.ref_genome} " \
#           f"--length {args.length} " \
#           f"--thread {args.thread} " \
#           f"--verbose {args.verbose} "
#     executeCmd(cmd)
#     logger.info(f"****************** Step 2 completed in {(time.time() - start) / 60} minutes. ******************")

    start = time.time()
    # Step 3: Masked re-alignment and variant calling of multiploid records
    cmd = f"python3 3_maskedAlignPolyVar.py --input_bam {args.input_bam} " \
          f"--bed_dir {os.path.join(ROOT, 'ref')} " \
          f"--intrinsic_vcf {os.path.join(args.output, 'vcf', 'all_pc.realign.' + str(args.length) + '.trim.vcf.gz')} " \
          f"--lower {args.lower} " \
          f"--masked_genomes {os.path.join(args.output, 'masked_genome')} " \
          f"--fastq_dir {os.path.join(args.output, 'fastq')} " \
          f"--ref_genome {args.ref_genome} " \
          f"--output_vcf {os.path.join(args.output, 'vcf')} " \
          f"--thread {args.thread} " \
          f"--verbose {args.verbose} "
    executeCmd(cmd)
    logger.info(f"****************** Step 3 completed in {(time.time() - start) / 60} minutes. ******************")

    start = time.time()
    # Step 4: Merge with prioritized VCF
    if not args.pvcf:
        return

    if not os.path.exists(args.pvcf):
        raise FileNotFoundError(f"Prioritized VCF not found {args.pvcf}. ")

    ovcf = glob.glob(os.path.join(args.output, 'vcf', '*filtered.vcf.gz'))
    if not ovcf:
        raise RuntimeError("Error in comparing with intrinsic VCF. ")

    cmd = f"python3 4_mergePVCF.py --pvcf {args.pvcf} " \
          f"--ovcf {ovcf[0]} " \
          f"--pv_tag GATK " \
          f"--ov_tag SDrecall " \
          f"--outpath {args.output} " \
          f"--verbose {args.verbose} "
    executeCmd(cmd)

    logger.info(f"****************** All processes completed in {time.time() - start // 60} minutes. ******************")

if __name__ == "__main__":

    # Argparse setup
    parser = argparse.ArgumentParser(description = "SDrecall wrapper.")
    parser._optionals.title = "Options"
    ROOT = os.path.dirname(__file__)
    parser.add_argument("-i", "--input_bam", type = str, required = True, help = "Input BAM file")
    parser.add_argument("-p", "--pvcf", type = str, default = None, help = "Path of prioritized VCF")
    parser.add_argument("-r", "--ref_genome", type = str, required = True, help = "reference genome to extract SD regions")
    parser.add_argument("-o", "--output", type = str, required = True, help = "output directory of resulting BED files")
    parser.add_argument("-b", "--build", type = str, required = True, help = "reference genome assembly")
    parser.add_argument("-l", "--list", type = str, required = True, help = "customized gene list")
    parser.add_argument("-t", "--thread", type = int, default = 8, help = "number of threads used with BISER (default = 8)")
    parser.add_argument("--length", type = int, default = 250, help = "BED window size for extracting reads from homologous regions (default: 250)")
    parser.add_argument("--lower", type = float, default = 1.1, help = "lower bound for unlikely intrinsic variants")
    parser.add_argument("-f", "--fraglen", type = int, default = 300, help = "expected fragment length (default: 300)")
    parser.add_argument("-g", "--gaplen", type = int, default = 10, help = "small gap cutoff value (default: 10)")
    parser.add_argument("--keep_trimmed", action = "store_true", help = "keep trimmed SD BED file for debugging")
    parser.add_argument("-v", "--verbose", type = str, default = "INFO", help = "verbosity level (default: INFO)")
    args = parser.parse_args()
    
    os.chdir(ROOT)
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P', level = args.verbose.upper())
    logger = logging.getLogger("root")
    logger.info(f"Working directory: {ROOT}")
    main()
