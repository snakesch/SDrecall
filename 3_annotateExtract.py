#!/usr/bin/env python3
import pandas as pd
import subprocess
import logging
import sys
import argparse
import re
import os

# Shell utilities
def annotate(INPATH: str, REF_REGIONS: str, OUTPATH: str):

    if subprocess.run("which bedtools", shell=True, capture_output=True).returncode != 0:
        logging.critical("BEDTools not detected. Abort.")
        raise ModuleNotFoundError("BEDTools not detected.")
    code = subprocess.run("tail -n+2 {} | bedtools intersect -a {} -b - -wo > {}".format(REF_REGIONS, INPATH, OUTPATH), shell=True).returncode
    if code != 0:
        logging.error("Error in bedtools intersect.")
        sys.exit(-1)
    return

def update_bed_file(bed_df, bed_path):
    bed_df.drop_duplicates().to_csv(bed_path, sep="\t", index=False, header=False)
    cmd = "cat {tbed} | sort -V -k1,3 | cut -f 1-3 | bedtools merge -i - > {tbed}.tmp && mv {tbed}.tmp {tbed}".format(tbed=bed_path)
    code = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8").returncode
    if code:
        logging.error("Error in creating " + bed_path + ".")
        sys.exit(-1)
    return

def merge_bed_intervals(bed_path):
    cmd = "bedtools sort -i {bed} | bedtools merge -i - > {bed}.tmp && mv {bed}.tmp {bed}".format(bed=bed_path)
    code = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8").returncode
    if code:
        logging.error("Error in merging BED intervals.")
        sys.exit(-1)
    return   
    
# Python utilities
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

def condense(df, outpath="", gene_col=17):
    cols = df.columns.tolist()
    by_pairs = df.groupby(by=cols[:6])
    condensed_df = by_pairs.apply(condense_per_group, cols, gene_col)
    if len(outpath) > 0:
        condensed_df.to_csv(outpath, sep='\t', mode="a", header=False, index=False)
    else:
        return condensed_df
    
def condense_per_group(group_df, cols, gene_col):
    condense_col = cols[gene_col-1]
    group_df[condense_col] = ",".join(list(dict.fromkeys(group_df[condense_col].tolist())))
    return group_df.drop_duplicates()
   
def generate_bed_per_gene(group_df, gene_col, outpath):
    
    if group_df.shape[0] == 0:
        return
    gene = group_df.iloc[0, gene_col].replace(",", "_")
    
    # Homologous regions
    total_bed = one_way(group_df, col_1=[0, 1, 2, gene_col], col_2=[3, 4, 5, gene_col])
    total_bed.iloc[:,[1,2]] = total_bed.iloc[:,[1,2]].astype(int)
    update_bed_file(total_bed, os.path.join(outpath, "homologous_regions", gene + "_related_homo_regions.bed"))
    merge_bed_intervals(os.path.join(outpath, "homologous_regions", gene + "_related_homo_regions.bed"))
    
    # Principal components
    selected_bed = group_df.iloc[:, [0,1,2,gene_col]]
    selected_bed.iloc[:,[1,2]] = selected_bed.iloc[:,[1,2]].astype(int)
    update_bed_file(selected_bed, os.path.join(outpath, "principal_components", gene + ".bed"))
    merge_bed_intervals(os.path.join(outpath, "principal_components", gene + ".bed"))
    
    out = pd.read_csv(os.path.join(outpath, "homologous_regions", gene + "_related_homo_regions.bed"), sep="\t", header=None)

    return out

def main():

    GENE_COL=17

    expanded_path = args.input.rstrip("trimmed.bed") + ".homo.expanded.bed"
    anno_path = expanded_path.rstrip("bed") + "geneanno.bed"

    # Make two-way map
    two_way = makeTwoWay(args.input, expanded_path)

    # Gene annotation
    annotate(expanded_path, args.ref, anno_path)
    os.remove(expanded_path)
    logging.info("Writing annotated BED file to " + anno_path)
    anno_df = pd.read_csv(anno_path, sep="\t", header=None)
    if args.list:
        region_df = pd.read_csv(args.list, sep='\t', encoding="utf-8")
        region_df['Genetic defect'] = region_df['Genetic defect'].apply(lambda x: x[:x.index(u'\xa0')] if u'\xa0' in x else x)
        region_genes = set(region_df['Genetic defect'].tolist())
        cleaned_df = anno_df[anno_df[GENE_COL-1].isin(region_genes)].drop_duplicates()
        cleaned_df = condense(cleaned_df)
    else:
        cleaned_df = anno_df
        cleaned_df = condense(cleaned_df)

    # Condense twice and deploy
    os.remove(anno_path)
    
    # Create necessary directories
    os.makedirs(args.out, exist_ok=True)
    os.makedirs(os.path.join(args.out, "principal_components"), exist_ok=True)
    os.makedirs(os.path.join(args.out, "homologous_regions"), exist_ok=True)
    by_gene = cleaned_df.groupby([GENE_COL-1], as_index=False)
    all_homo_regions = by_gene.apply(generate_bed_per_gene, GENE_COL-1, args.out)
    all_homo_regions = all_homo_regions.droplevel(level=0).reset_index(drop=True)
    all_homo_regions.to_csv(os.path.join(args.out, "homologous_regions", "all_homo_regions.bed"), sep="\t", index=False, header=False)
    merge_bed_intervals(os.path.join(args.out, "homologous_regions", "all_homo_regions.bed"))
    
    # Re-annotate all_homo_regions.bed
    cmd = "bedtools intersect -a {bed} -b {ref} -wao | cut -f 1,2,3,7 | sort -uV | uniq > {bed}.tmp; mv -f {bed}.tmp {bed}".format(
        bed=os.path.join(args.out, "homologous_regions", "all_homo_regions.bed"),
        ref=args.ref)
    code = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8").returncode
    if code:
        logging.error("Error in annotating " + os.path.join(args.out, "homologous_regions", "all_homo_regions.bed"))
 
    return

if __name__ == "__main__":
    # Argparse setup
    parser = argparse.ArgumentParser(description = "Make two-way map and extract genes of interest.")
    parser._optionals.title = "Options"
    parser.add_argument("-i", "--input", type = str, required = True, help = "trimmed BED file")
    parser.add_argument("-r", "--ref", type = str, required = True, help = "gene region list for annotation")
    parser.add_argument("-l", "--list", type = str, required = False, help = "list of genes/regions of interest")
    parser.add_argument("-o", "--out", type = str, required = True, help = "output path")
    parser.add_argument("-v", "--verbose", type = str, default = "INFO", help = "verbosity level (default: INFO)")
    args = parser.parse_args()
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P',
                        level = args.verbose.upper())
    main()
