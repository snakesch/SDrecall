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

from src.utils import executeCmd
from src.trim_cigar import *
from src.annotate_from_ref import *

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
        executeCmd(f"biser -o {output} -t {thread} --no-decomposition {ref_genome}")
    else:
        raise RuntimeError("BISER not found. ")
          
def main_trim(infile: str, outfile: str, fraglen=300, gaplen=10, logLevel="INFO"):
    """
    This function implements trimming of CIGAR strings and creates a trimmed BED file.
    """
    logger = logging.getLogger("root")
    
    raw_df = pd.read_csv(infile, sep = "\t", header = None)
    logger.info("Total number of reads before trimming = {}".format(raw_df.shape[0]))
    ret = trim(raw_df, frag_len = fraglen, small_gap_cutoff = gaplen, logLevel = logLevel)
    ret.to_csv(outfile, sep = "\t", index = False, header = None, mode = 'w')
    logger.info("****************** CIGAR trimming completed ******************")
    
def annotate_from_refgene(infile, outpath, anno_table: str, genelist: str, keep=False, genecol=16):
    """
    This function annotates trimmed SD regions and extracts gene regions in gene list.
    
    Arguments:
    ----------
    infile: trimmed SD BED file
    outpath: output directory of BED files
    anno_table: path of annotation table
    keep: keep trimmed SD BED file
    genecol: 0-based column index of gene names in annotated dataframe
    
    """
    
    def combine_pc(pc_dir):

        wd = os.getcwd()
        os.chdir(pc_dir)

        dfs = []
        for bedf in os.listdir():
            if not bedf.startswith("all"):
                cur_df = pd.read_csv(bedf, sep="\t", header=None)
                gene = bedf.split(".")[0]
                cur_df[3] = gene
                dfs.append(cur_df)
        df = pd.concat(dfs, ignore_index=True)
        df = df.sort_values(by=[0, 1, 2], ignore_index=True)
        df.to_csv("all_pc.bed", sep="\t", header=False, index=False)

        os.chdir(wd)
    
    logger = logging.getLogger("root")
    
    trimmed_df = pd.read_csv(infile, sep = "\t", header = None)
    dup_df = trimmed_df.iloc[:, [3, 4, 5, 0, 1, 2, 6, 7, 9, 8, 10, 11, 12, 13]]
    two_way_df = pd.concat([dup_df.T.reset_index(drop=True), trimmed_df.T.reset_index(drop=True)], axis=1, ignore_index=True).T
    two_way_out = os.path.join(infile + ".expanded")
    two_way_df.iloc[:, :13].to_csv(two_way_out, sep="\t", index=False, header=False)
    
    # Annotate with refGene data
    anno_out = os.path.join(infile + ".expanded.anno")
    annotate_from_ref(two_way_out, anno_table, anno_out)
    
    # Extract listed genes
    anno_df = pd.read_csv(anno_out, sep="\t", header=None)
    regions = set()
    with open(genelist, "r") as f:
        for line in f.readlines():
            gene = line.strip("\n").split(u"\xa0")[0]
            if " " not in gene and "Unknown" not in gene and "Large" not in gene and gene:
                regions.add(gene)
         
    cleaned_df = anno_df[anno_df[genecol].isin(regions)].drop_duplicates()
    cleaned_df = condense(cleaned_df)
    
    # Create necessary directories
    os.makedirs(outpath, exist_ok=True)
    os.makedirs(os.path.join(outpath, "principal_components"), exist_ok=True)
    os.makedirs(os.path.join(outpath, "homologous_regions"), exist_ok=True)
    by_gene = cleaned_df.groupby([genecol], as_index=False)
    all_homo_regions = by_gene.apply(generate_bed_per_gene, genecol, outpath)
    all_homo_regions = all_homo_regions.droplevel(level=0).reset_index(drop=True)
    logger.info("Writing all homo regions to " + os.path.join(outpath, "homologous_regions", "all_homo_regions.bed"))
    all_homo_regions.to_csv(os.path.join(outpath, "homologous_regions", "all_homo_regions.bed"), sep="\t", index=False, header=False)
    merge_bed_intervals(os.path.join(outpath, "homologous_regions", "all_homo_regions.bed"))
    
    # Re-annotate all_homo_regions.bed
    cmd = "tail -n+2 {ref} | bedtools intersect -a {bed} -b - -wao | cut -f 1,2,3,7 | sort -uV | uniq > {bed}.tmp; mv -f {bed}.tmp {bed}".format(
        bed=os.path.join(outpath, "homologous_regions", "all_homo_regions.bed"),
        ref=anno_table)
    code = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8").returncode
    if code:
        raise RuntimeError("Error in annotating " + os.path.join(outpath, "homologous_regions", "all_homo_regions.bed"))
    
    # Combine all principal component BED
    combine_pc(os.path.join(outpath, "principal_components"))
    
    # Remove intermediate files
    if not keep:
        os.remove(infile)
    os.remove(two_way_out)
    os.remove(anno_out)
    
if __name__ == "__main__":
    # Argparse setup
    parser = argparse.ArgumentParser(description = "Extract homologous trimmed SD regions from reference genome.")
    parser._optionals.title = "Options"
    parser.add_argument("-r", "--ref_genome", type = str, required = True, help = "reference genome to extract SD regions")
    parser.add_argument("-o", "--output", type = str, required = True, help = "output directory of resulting BED files")
    parser.add_argument("-b", "--build", type = str, required = True, help = "reference genome assembly")
    parser.add_argument("-l", "--list", type = str, required = True, help = "customized gene list")
    parser.add_argument("-a", "--table", type = str, required = True, help = "RefGene annotated table")
    parser.add_argument("-t", "--thread", type = int, default = 8, help = "number of threads used with BISER (default = 8)")
    parser.add_argument("-f", "--fraglen", type = int, default = 300, help = "expected fragment length (default: 300)")
    parser.add_argument("-g", "--gaplen", type = int, default = 10, help = "small gap cutoff value (default: 10)")
    parser.add_argument("--keep_trimmed", action = "store_true", help = "keep trimmed SD BED file for debugging")
    parser.add_argument("-v", "--verbose", type = str, default = "INFO", help = "verbosity level (default: INFO)")
    args = parser.parse_args()
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P', level = args.verbose.upper())
    logger = logging.getLogger("root")
    
    # Extract SD regions
    biser_out = os.path.join(args.output, "SD_" + args.build + ".bed")
    if not os.path.exists(biser_out):    
        biser(args.ref_genome, args.build, args.thread, biser_out)
    else:
        logger.info(f"BISER output detected. Using {biser_out}")
    
    # Trim SD regions
    trimmed_out = os.path.join(args.output, "SD_" + args.build + "_trimmed.bed")
    main_trim(biser_out, trimmed_out, fraglen=args.fraglen, gaplen=args.gaplen, logLevel=args.verbose.upper())
    os.remove(biser_out)
    
    annotate_from_refgene(trimmed_out, args.output, args.table, args.list, args.keep_trimmed, genecol=16)
