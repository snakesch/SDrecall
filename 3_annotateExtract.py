#!/usr/bin/env python3
import pandas as pd
import subprocess
import logging
import sys
import argparse
import re
import os

def makeTwoWay(INPATH: str, OUTPATH: str):

    trimmed_df = pd.read_csv(INPATH, sep = "\t", header = None)
    dup = trimmed_df[[3, 4, 5, 0, 1, 2, 6, 7, 9, 8, 10, 11, 12, 13]]
    out = pd.concat([trimmed_df, dup], axis = 0).reset_index(drop=True)
    out.drop(13, axis=1, inplace=True)
    out.to_csv(OUTPATH, sep="\t", index=False, header=None)

    return

def one_way(two_way_df, col_1 = [0, 1, 2, 16], col_2 = [3, 4, 5, 16]):
    '''
    We assume gene name is on the 17th column (0-based). Generate one-way map.
    '''
    df_1 = two_way_df.iloc[:, col_1]
    df_1.columns = ["chr", "start", "end", "gene"]
    df_2 = two_way_df.iloc[:, col_2]
    df_2.columns = ["chr", "start", "end", "gene"]

    return pd.concat([df_1, df_2], axis=0).drop_duplicates().reset_index(drop=True)

def intersectBED(sorted_df, to_path=None):
    """
    Note: This function works in a chromosome-wise and gene-wise manner. To_path should be a valid directory.
    """
    chrom = sorted_df["chr"].tolist()[0]
    start_lst = sorted_df["start"].tolist()
    end_lst = sorted_df["end"].tolist()
    gene = sorted_df["gene"].tolist()[0]
    if len(start_lst) == 1:
        return sorted_df
    prev, cur = 0, 1
    while cur < len(start_lst):
        if start_lst[prev] == start_lst[cur]:
            if end_lst[prev] > end_lst[cur]:
                end_lst.pop(cur)
                start_lst.pop(cur)
            else:
                end_lst.pop(prev)
                start_lst.pop(prev)
        elif end_lst[prev] == end_lst[cur]:
            if start_lst[prev] < start_lst[cur]:
                start_lst.pop(cur)
                end_lst.pop(cur)
            else:
                start_lst.pop(prev)
                end_lst.pop(prev)
        elif end_lst[prev] >= start_lst[cur]:
            end_lst.pop(prev)
            start_lst.pop(cur)
        else:
            prev += 1
            cur += 1
    out = pd.DataFrame({"chr": chrom, "start": start_lst, "end": end_lst, "gene": gene})
    if not to_path:
        return out
    else:
        os.makedirs(os.path.join(to_path), exist_ok=True)
        out.to_csv(os.path.join(to_path, gene + "_related_homo_regions.bed"), sep="\t", index=False, header=False)
        logging.info("Successfully extracted homologous regions to " + os.path.join(to_path, gene + "_related_homo_regions.bed") + ".")
        return

def annotate(INPATH: str, REF_REGIONS: str, OUTPATH: str):

    if subprocess.run("which bedtools", shell=True, capture_output=True).returncode != 0:
        logging.critical("BEDTools not detected. Abort.")
        raise ModuleNotFoundError("BEDTools not detected.")
    code = subprocess.run("sed 's/\t$//g' {} | bedtools intersect -a {} -b - -wo > {}".format(REF_REGIONS, INPATH, OUTPATH), shell=True).returncode
    if code != 0:
        logging.error("Error in bedtools intersect.")
        sys.exit(-1)
    return

def main():

    INPATH = args.input
    REF_REGIONS = args.ref
    REGION_PATH = args.list

    expanded_path = INPATH.rstrip("trimmed.bed") + ".homo.expanded.bed"
    anno_path = expanded_path.rstrip("bed") + "geneanno.bed"
    outpath = expanded_path.rstrip("bed") + "geneanno.region.bed"

    # Make two-way map
    two_way = makeTwoWay(INPATH, expanded_path)

    # Gene annotation
    annotate(expanded_path, REF_REGIONS, anno_path)
    logging.info("Writing annotated BED file to " + anno_path)
    anno_df = pd.read_csv(anno_path, sep="\t", header=None)
    if args.list:
        region_df = pd.read_csv(REGION_PATH, sep='\t', encoding="utf-8")
        region_df['Genetic defect'] = region_df['Genetic defect'].apply(lambda x: x[:x.index(u'\xa0')] if u'\xa0' in x else x)
        region_genes = set(region_df['Genetic defect'].tolist())
        cleaned_df = anno_df[anno_df[args.genecol].isin(region_genes)].drop_duplicates()
    else:
        cleaned_df = anno_df

    one_way_df = one_way(cleaned_df, col_1=[0, 1, 2, args.genecol], col_2=[3, 4, 5, args.genecol])
    sorted_df = one_way_df.sort_values(["chr", "start"]).reset_index(drop=True)
    grouped = sorted_df.groupby(["chr", "gene"], as_index=False).apply(intersectBED, to_path=args.out)
    total_bed = grouped.droplevel(level=0).reset_index(drop=True)
    os.makedirs(args.out, exist_ok=True)
    total_bed.to_csv(os.path.join(args.out, "all_homo_regions.bed"), sep="\t", index=False, header=False)
    logging.info("Successfully combined all homologous regions to " + os.path.join(args.out, "all_homo_regions.bed"))

    return

if __name__ == "__main__":
    # Argparse setup
    parser = argparse.ArgumentParser(description = "Make two-way map and extract genes of interest.")
    parser._optionals.title = "Options"
    parser.add_argument("-i", "--input", type = str, required = True, help = "trimmed BED file")
    parser.add_argument("-r", "--ref", type = str, required = True, help = "gene region list for annotation")
    parser.add_argument("-l", "--list", type = str, required = False, help = "list of genes/regions of interest")
    parser.add_argument("-c", "--genecol", type = int, default = 16, help = "0-based column index of \"Genetic defect\" in gene list (default: 16)")
    parser.add_argument("-o", "--out", type = str, required = True, help = "output path")
    parser.add_argument("-v", "--verbose", type = str, default = "INFO", help = "verbosity level (default: INFO)")
    args = parser.parse_args()
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P',
                        level = args.verbose.upper())
    main()
