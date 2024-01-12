#!/usr/bin/env python3

# 2_preparation.py
# Description: Preparation for masked re-alignment.
# Author: Yang XT, She CH (2022)
# Contact: Xingtian Yang (u3005579@connect.hku.hk), Louis She (snakesch@connect.hku.hk)

import subprocess
import os
import logging
import pandas as pd
from multiprocessing import Pool
from itertools import repeat
import glob
import time
import argparse
from src.utils import executeCmd, get_gene_bed_fn
from src.maskGenome import *

def getHomoRelatedBam(infile, all_homo_bed, thread):
    """
    This function extracts reads mapped to any homologous regions and outputs a BAM.
    """

    logger = logging.getLogger("root")
    
    homorelated_out = infile + ".homorelated"
    _qname = infile + ".qname"
    _sorted_bam = infile + ".homorelated_sorted"
    
    if os.path.exists(_sorted_bam):
        return _sorted_bam
    
    if not os.path.exists(homorelated_out):
        cmd = f"samtools view -@ {thread} -L {all_homo_bed} {infile} | awk -F'\t' '{{print $1;}}' | sort -u > {_qname} "
        executeCmd(cmd)
        cmd = f"gatk FilterSamReads -I {infile} -O {homorelated_out} --FILTER includeReadList -RLF {_qname} -SO coordinate "
        executeCmd(cmd)
        os.remove(_qname)
        if subprocess.run(f"samtools view -c {homorelated_out}", shell=True, capture_output=True).stdout == '0':
            raise ValueError(f"No reads mapped to homologous regions in {infile}. ")
        logger.debug(f"Selected reads mapped to homologous regions in {infile}.")

    cmd = f"samtools sort -n -O bam -o {_sorted_bam} {homorelated_out} "
    executeCmd(cmd)
    os.remove(homorelated_out)

    return _sorted_bam

def extractFASTQ(infile, region_bed, outpath, mq):
    """
    This function extracts multi-aligned reads mapped to a given region and outputs reads in FASTQ format.

    infile should be in BAM format and only contain reads mapped to homologous regions and sorted by QNAME.

    It is designated for multithreading.

    """
    region_name = os.path.basename(region_bed).split("_r")[0]
    logger.debug(f"Detected region {region_name}. ")
    
    out_1 = os.path.join(outpath, os.path.basename(infile).split(".")[0] + "_" + region_name + "-XA_1.fq")
    out_2 = os.path.join(outpath, os.path.basename(infile).split(".")[0] + "_" + region_name + "-XA_2.fq")
    
    if os.path.exists(out_1) and os.path.exists(out_2):
        return

    tmp_bed = region_bed + ".tmp"
    cmd = f"bedtools sort -i {region_bed} | bedtools merge -i - > {tmp_bed}"
    executeCmd(cmd)

    tmp_qname = infile + "_" + region_name + ".qname"
    cmd = f"samtools view -L {tmp_bed} {infile} | sort -u > {tmp_qname} "
    executeCmd(cmd)

    os.remove(tmp_bed)

    if os.path.getsize(tmp_qname) == 0:
        os.remove(tmp_qname)
        logger.debug(f"No multi-aligned reads with MQ < 30 in {region_name}")
        return

    _extracted_bam = infile + "_only_" + region_name + ".bam"
    cmd = f"samtools view -N {tmp_qname} -o {_extracted_bam} {infile} "
    executeCmd(cmd)

    cmd = f"bedtools bamtofastq -i {_extracted_bam} -fq {out_1} -fq2 {out_2} >/dev/null 2>&1 && gzip -f {out_1} && gzip -f {out_2} "
    executeCmd(cmd)

    os.remove(_extracted_bam)
    os.remove(tmp_qname)

def mask_genome(bedf, ref_genome: str, outpath: str):

    rg = Genome(ref_genome)
    masked_genome = rg.mask(bedf)
    out = os.path.join(outpath, "masked_genome")
    cmd = f"mv {masked_genome} {out} "
    executeCmd(cmd)

def isMerged(df):
    """
    Equivalent to `bedtools merge -d 0`

    Returns:
    ---------
    Bool: True if the regions are merged by coordinates.

    """
    start = df.iloc[:, 1].tolist()
    end = df.iloc[:, 2].tolist()
    for i in range(len(start) - 1):
        if end[i] > start[i+1]:
            return False
    return True

def makeWindows(df, length):
    """
    This function takes a BED-like dataframe and makes windows of size LENGTH.

    Regions with length of trailing block < 100 are appended to the previous block.

    """
    start = df.iloc[:, 1].tolist()
    end = df.iloc[:, 2].tolist()
    ret_start, ret_end = [], []
    for i in range(len(start)):
        _start, _end = [], []
        cur_start, cur_end = start[i], end[i]
        short = True
        while cur_end - cur_start > length:
            _start.append(cur_start)
            _end.append(cur_start + length)
            cur_start += length
            short = False
        if not short and cur_end - cur_start < 100:
            _end[-1] = cur_end
        else:
            _start.append(cur_start)
            _end.append(cur_end)
        ret_start.append(_start)
        ret_end.append(_end)

    out_df = pd.DataFrame({
        "Chr": df.iloc[0,0],
        "Start": ret_start,
        "End": ret_end
    })

    return out_df.explode(["Start", "End"], ignore_index=True)

def getSubseq(bedf_path: str, fastq_path: str, length: int, ref_genome: str) -> None:

    if length <= 100:
        raise ValueError(f"Read length {length} is too small (<100). ")

    bedf = pd.read_csv(bedf_path, sep="\t", header=None)
    if not bedf.groupby(by=0).apply(isMerged).all():
        raise ValueError(f"Input BED file {bedf_path} is not sorted / merged. ")

    window_df = bedf.groupby(by=0, group_keys=False).apply(lambda x: makeWindows(x, 250)).reset_index(level=0, drop=True)
    window_df["Length"] = window_df.iloc[:, 2] - window_df.iloc[:, 1]
    bedf[3] = bedf[2] - bedf[1]
    if window_df.Length.sum() != bedf[3].sum():
        raise RuntimeError(f"Error in preprocessing {bedf_path}. Lengths of regions {window_df.Length.sum()} do not match. {bedf[3].sum()} ")
    window_out = bedf_path + ".preproc"
    window_df.to_csv(window_out, sep="\t", index=False, header=False)

    os.makedirs(fastq_path, exist_ok=True)
    fq_out = os.path.join(fastq_path, os.path.basename(bedf_path)[:-3] + str(length) + ".fastq")
    cmd = f"seqtk subseq {ref_genome} {window_out} | seqtk seq -F 'F' - > {fq_out} "
    executeCmd(cmd)
    logger.debug(f"Subsequence written to {fq_out}. ")
    os.remove(window_out)

    return

def sortBed(bedf):

    p = subprocess.run(f"cut -f1-3 {bedf} | sort -u -V -k1 -k2,2n | bedtools merge -i - > {bedf + '.tmp'} && mv {bedf + '.tmp'} {bedf} ", shell=True, capture_output=True)
    
    if p.returncode != 0:
        raise ValueError("BEDTools merge failed. Something wrong in the BED file?")

    return bedf

def getIntrinsicVcf(all_pc_fp, all_homo_regions_fp, fastq_dir, ref_genome, vcf_dir, length, nthreads):

    logger = logging.getLogger("root")
    
    # Genome masking
    all_pc_fp = sortBed(all_pc_fp)
    rg = Genome(ref_genome)
    pc_masked = rg.mask(all_pc_fp)

    if not os.path.exists(pc_masked + ".bwt"):
        cmd = f"bwa index {pc_masked} "
        executeCmd(cmd)
    else:
        logger.info(f"All priority component-masked genome index detected - {pc_masked + '.bwt'}")

    homo_sd = os.path.dirname(all_homo_regions_fp)
        
    # Get analytic regions BED
    pc_genes = (x.split("/")[-1].split(".")[0] for x in get_gene_bed_fn(os.path.dirname(all_pc_fp)))
    sd_genes = ( fn.split("_")[0] for fn in os.listdir(homo_sd) )
    analytic_genes = list(set(sd_genes) & set(pc_genes))
    fp = [ os.path.join(homo_sd, region + "_related_homo_regions.bed") for region in analytic_genes ]

    with Pool(nthreads) as p:
        p.starmap(getSubseq, zip(fp, repeat(fastq_dir), repeat(length), repeat(ref_genome)))

    out_seq = 0
    _out = os.path.join(homo_sd, "all_pc.reseq." + str(length) + ".fastq")
    with open(_out, "w") as ofs:
        for inf in glob.glob(os.path.join(fastq_dir, f"*.{length}.fastq")):
            with open(inf, "r") as ifs:
                for line in ifs:
                    ofs.write(line)
            os.remove(inf)
            out_seq += 1

    if out_seq <= 1:
        raise ValueError(f"Number of FASTQ files {out_seq} is fewer than expected. ")

    realigned_bam_out = os.path.join(vcf_dir, "all_pc.realign." + str(length) + ".bam")
    if not os.path.exists(realigned_bam_out):
        cmd = f"""bwa mem -t {nthreads} -M -R """ + r"""'@RG\tID:all_pc\tLB:SureSelectXT Library Prep Kit\tPL:ILLUMINA\tPU:1064\tSM:all_pc' """ \
              f"""{pc_masked} {_out} | samtools view -@ {nthreads} -uSh - | samtools sort -O bam -@ {nthreads} -o {realigned_bam_out} && """ \
              f"""samtools index {realigned_bam_out} """
        executeCmd(cmd)

    os.makedirs(vcf_dir, exist_ok=True)
    tmp_vcf = os.path.join(vcf_dir, os.path.basename(realigned_bam_out)[:-3] + "vcf.gz")
    cmd = f"bcftools mpileup -a FORMAT/AD,FORMAT/DP -f {ref_genome} {realigned_bam_out} |  " \
          f"bcftools call -mv -Oz -o {tmp_vcf} && tabix -p vcf {tmp_vcf} "
    executeCmd(cmd)
    vcf_out = os.path.join(vcf_dir, os.path.basename(realigned_bam_out)[:-3] + "trim.vcf.gz")
    cmd = f"gatk LeftAlignAndTrimVariants -V {tmp_vcf} -R {ref_genome} --split-multi-allelics -O {vcf_out} "
    executeCmd(cmd)
    
    logger.info(f"Writing intrinsic VCF to {vcf_out}")
    
    os.remove(realigned_bam_out)
    os.remove(tmp_vcf)
    os.remove(_out)

if __name__ == "__main__":

    # Argparse setup
    parser = argparse.ArgumentParser(description = "Preparation for masked re-alignment.")
    parser._optionals.title = "Options"
    ROOT = os.path.dirname(__file__)
    parser.add_argument("--input_bam", type = str, help = "input BAM file", required = True)
    parser.add_argument("--homo_dir", type = str, help = "directory of BED files of homologous regions", default = f"{ROOT}/ref/homologous_regions")
    parser.add_argument("--pc_dir", type = str, help = "directory of BED files of principal components", default = f"{ROOT}/ref/principal_components")
    parser.add_argument("-mq", "--mapping_quality", type = int, help = "MQ threshold", default = 30)
    parser.add_argument("--ref_genome", type = str, help = "reference genome", required = True)
    parser.add_argument("--outpath", type = str, help = "output directory for FASTQ and masked genomes", required = True)
    parser.add_argument("--length", type = int, help = "BED window size for extracting reads from homologous regions (default: 250)", default = 250)
    parser.add_argument("--thread", type = int, help = "number of threads (default: 8)", default = 8)
    parser.add_argument("-v", "--verbose", type = str, default = "INFO", help = "verbosity level (default: INFO)")
    args = parser.parse_args()
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P',
                        level = args.verbose.upper())
    
    logger = logging.getLogger("root")
    logger.debug(f"Working in {os.getcwd()}")

    fq_out = os.path.join(args.outpath, "fastq")
    vcf_out = os.path.join(args.outpath, "vcf")
    os.makedirs(fq_out, exist_ok=True)
    os.makedirs(vcf_out, exist_ok=True)
    os.makedirs(os.path.join(args.outpath, "masked_genome"), exist_ok=True)
    all_homo_bed = os.path.join(args.homo_dir, "all_homo_regions.bed")
    all_pc_bed = os.path.join(args.pc_dir, "all_pc.bed")
    if not os.path.exists(all_homo_bed):
        raise FileNotFoundError(f"{os.path.join(args.homo_dir, 'all_homo_regions.bed')} not found. ")
    if not os.path.exists(all_pc_bed):
        raise FileNotFoundError(f"{os.path.join(args.pc_dir, 'all_pc.bed')} not found. ")

    start = time.time()
    # Extract homologous region reads from BAM
    _sorted_bam = getHomoRelatedBam(args.input_bam, all_homo_bed, args.thread)
    logger.info(f"****************** BAM preprocessing completed in {time.time() - start:.2f} seconds ******************")

    # Extract reads from homologous regions
    start = time.time()
    region_beds = glob.glob(os.path.join(args.homo_dir, "*_related_homo_regions.bed"))

    with Pool(args.thread) as p:
        p.starmap(extractFASTQ, zip(repeat(_sorted_bam), region_beds, repeat(fq_out), repeat(args.mapping_quality)))
    logger.info(f"****************** FASTQ extraction completed in {time.time() - start:.2f} seconds ******************")

    start = time.time()
    overlap_bedfs = get_gene_bed_fn(args.pc_dir)
    with Pool(args.thread) as p:
        p.starmap(mask_genome, zip(overlap_bedfs, repeat(args.ref_genome), repeat(args.outpath)))
    logger.info(f"****************** Genome masking completed in {time.time() - start:.2f} seconds ******************")

    # Get intrinsic VCF for selecting unlikely intrinsic variants
    start = time.time()
    getIntrinsicVcf(all_pc_bed, all_homo_bed, fq_out, args.ref_genome, vcf_out, args.length, args.thread)
    logger.info(f"****************** Intrinsic VCF acquired in {time.time() - start:.2f} seconds ******************")

    os.remove(_sorted_bam)


