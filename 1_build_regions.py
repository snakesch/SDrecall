#!/usr/bin/env python3

# 1_getTrimmedSD.py
# Description: Extract homologous regions from trimmed SD regions from reference genome.
# Author: Yang XT, She CH (2022)
# Contact: Xingtian Yang (u3005579@connect.hku.hk), Louis She (snakesch@connect.hku.hk)

import os
import logging
import argparse
import shutil
import subprocess
import pandas as pd
from pandarallel import pandarallel as pa

from src.utils import executeCmd, timing
from src.trim_cigar import *

@timing
def biser(ref_genome: str, build: str, thread: int, output: str):
    """
    This function extracts SD regions from given reference genome using BISER.   
    """
    
    logger = logging.getLogger("root")
    logger.info(f"Extracting SD regions with {thread} threads.")
    logger.info(f"Reference genome asssembly : {build}")

    if shutil.which("conda"):
        _biser = subprocess.run("source activate SDrecall; conda list | grep biser", 
                      shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True).stdout
        _version = _biser.replace(" ", "").split("pypi")[0].lstrip("biser")
        logger.info(f"Running BISER version {_version}")
    
    if not os.path.exists(ref_genome + ".fai"):
        executeCmd(f"samtools faidx {ref_genome}")
    
    if shutil.which("biser"):
        executeCmd(f"biser -o {output} -t {thread} {ref_genome}")
    else:
        raise RuntimeError("BISER not found. ")

@timing       
def cigar_trim(infile: str, outfile: str, mismatch_rate=5, fraglen=300, gaplen=10, trace_trim=False):
    """
    This function implements trimming of CIGAR strings and creates a trimmed BED file.
    """
    import re
    
    logger = logging.getLogger("root")
    
    raw_df = pd.read_csv(infile, sep = "\t", header = None)
    
    # Filter SD regions with mismatch rate X < mismatch_rate
    raw_df = raw_df[raw_df.iloc[:, 13].apply(lambda x: float(re.split("\=|\;", x)[1]) < mismatch_rate)]
    raw_df = raw_df.reset_index(drop = True)
    logger.info(f"Selecting SD regions with mismatch rate X<{mismatch_rate}. ")
    
    logger.info("Total number of reads before trimming = {}".format(raw_df.shape[0]))
    ret = trim(raw_df, frag_len = fraglen, small_gap_cutoff = gaplen, trace=trace_trim)
    ret.to_csv(outfile, sep = "\t", index = False, header = None, mode = 'w')
    logger.info("CIGAR trimming completed")

@timing
def build_regions(infile, outpath, anno_table: str, keep=False):
    """
    This function annotates trimmed SD regions and extracts gene regions in gene list.
    
    Arguments:
    ----------
    infile: trimmed SD BED file
    outpath: output directory of BED files
    anno_table: path of annotation table (generated in previous step)
    keep: keep trimmed SD BED and annotated BED
    
    """
    
    import copy
    
    gene_idx = 16
    
    def merge_gene(df):
    
        if df.shape[0] == 1:
            return df

        out = df.head(1).copy()
        out[gene_idx] = ",".join(set(df.iloc[:, gene_idx].tolist()))

        return out 
    
    def build_per_gene(grp_df, gene):

        ### Output files (per gene)
        homo_out = os.path.join("homologous_regions", gene + "_related_homo_regions.bed")
        pc_out = os.path.join("principal_components", gene + ".bed")

        ### Homologous region BED
        homo_1, homo_2 = grp_df.iloc[:, [0,1,2,16]], grp_df.iloc[:, [3,4,5,16]]
        homo_2.columns = homo_1.columns
        homo_df = pd.concat([ homo_1, homo_2 ], axis=0, ignore_index=True) 
        homo_df.to_csv(homo_out + '.tmp', sep="\t", index=False, header=False)
        p = subprocess.run(f"bedtools sort -i {homo_out + '.tmp'} | bedtools merge -i - > {homo_out} ", shell=True, capture_output=False)
        if p.returncode != 0:
            raise RuntimeError("Error in building homologous regions. ")

        ### Principal component BED
        homo_1.to_csv(pc_out + '.tmp', sep="\t", index=False, header=False)
        p = subprocess.run(f"bedtools sort -i {pc_out + '.tmp'} | bedtools merge -i - > {pc_out} ", shell=True, capture_output=False)
        if p.returncode != 0:
            raise RuntimeError("Error in building principal components. ")

    logger = logging.getLogger("root")
    
    # Make two-way map for gene annotation
    trimmed_df = pd.read_csv(infile, sep = "\t", header = None, usecols = range(13))
    dup_df = trimmed_df.iloc[:, [3, 4, 5, 0, 1, 2, 6, 7, 9, 8, 10, 11, 12]]
    two_way_df = pd.concat([dup_df, trimmed_df], axis=0, ignore_index=True)
    two_way_df.to_csv(infile + ".expanded", sep="\t", index=False, header=False)
    
    # Annotate with refGene data
    anno_out = os.path.join(infile + ".expanded.anno")
    p = subprocess.run(f"tail -n+2 {anno_table} | bedtools intersect -a {infile + '.expanded'} -b - -wo | bedtools sort -i - > {anno_out} ", shell=True)
    if p.returncode != 0:
        raise RuntimeError("Failed to produce two-way map")
    ### Merge gene symbols
    anno_df = pd.read_csv(anno_out, sep="\t", header=None)
    cleaned_df = anno_df
#     cleaned_df = anno_df.groupby(by=list(range(7)), group_keys=False, sort=True).parallel_apply(merge_gene)
 
    # Create necessary directories
    os.makedirs(os.path.join(outpath, "principal_components"), exist_ok=True)
    os.makedirs(os.path.join(outpath, "homologous_regions"), exist_ok=True)
    os.chdir(outpath)
    
    # Output all homologous regions / principal components
    all_1, all_2 = cleaned_df.iloc[:, [0,1,2,16]], cleaned_df.iloc[:, [3,4,5,16]]
    all_2.columns = all_1.columns
    all_homo_df = pd.concat([all_1, all_2], axis=0, ignore_index=True)
    all_homo_df.to_csv("homologous_regions/all_homo_regions.bed", sep="\t", index=False, header=False)
    all_1.to_csv("principal_components/all_pc.bed", sep="\t", header=False, index=False)
    
    ### Output coordinates for each gene (pandarallel not working here)
    cleaned_df.groupby(by=gene_idx, group_keys=False, as_index=False).apply(lambda x: build_per_gene(x, x.name))
    
    subprocess.run("rm -f homologous_regions/*.tmp principal_components/*.tmp", shell=True, capture_output=False)

    if not keep:
        os.remove(anno_out)
        os.remove(infile)

if __name__ == "__main__":
    # Argparse setup
    parser = argparse.ArgumentParser(description = "Extract homologous trimmed SD regions from reference genome.")
    parser._optionals.title = "Options"
    parser.add_argument("-r", "--ref_genome", type = str, required = True, help = "reference genome to extract SD regions")
    parser.add_argument("-o", "--output", type = str, required = True, help = "output directory of resulting BED files")
    parser.add_argument("-b", "--build", type = str, required = True, help = "reference genome assembly")
    parser.add_argument("-a", "--table", type = str, required = True, help = "RefGene annotated table")
    parser.add_argument("-t", "--thread", type = int, default = 8, help = "number of threads used with BISER (default = 8)")
    parser.add_argument("-f", "--fraglen", type = int, default = 300, help = "expected fragment length (default: 300)")
    parser.add_argument("-m", "--mismatch_rate", type = int, default = 5, help = "mismatch rate threshold (default: 5)")
    parser.add_argument("-g", "--gaplen", type = int, default = 10, help = "small gap cutoff value (default: 10)")
    parser.add_argument("-v", "--verbose", type = str, default = "INFO", help = "verbosity level (default: INFO)")
    
    parser.add_argument("--trace_trim", action = "store_true", help = "print all trimming log")
    parser.add_argument("--keep_trimmed", action = "store_true", help = "keep trimmed SD BED file for debugging")
    args = parser.parse_args()
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P', level = args.verbose.upper())
    logger = logging.getLogger("root")
    
    pa.initialize(verbose=0, nb_workers=args.thread)
    
    # Extract SD regions    
    biser_out = os.path.join(args.output, "SD_" + args.build + ".bed")
    trimmed_out = os.path.join(args.output, "SD_" + args.build + "_trimmed.bed")
    if os.path.exists(trimmed_out):
        logger.info("Trimmed output detected. ")
    elif os.path.exists(biser_out): 
        logger.info(f"BISER output detected. Using {biser_out}")
        cigar_trim(biser_out, trimmed_out, mismatch_rate=args.mismatch_rate, fraglen=args.fraglen, gaplen=args.gaplen, trace_trim=args.trace_trim)
        os.remove(biser_out)
    else:
        biser(args.ref_genome, args.build, args.thread, biser_out)
        cigar_trim(biser_out, trimmed_out, mismatch_rate=args.mismatch_rate, fraglen=args.fraglen, gaplen=args.gaplen, trace_trim=args.trace_trim)
        os.remove(biser_out)
    logger.info("Building coordinates ... ")
    build_regions(trimmed_out, args.output, args.table, args.keep_trimmed)
    