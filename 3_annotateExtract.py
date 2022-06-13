#!/usr/bin/env python3
import pandas as pd
import subprocess
import logging
import sys
import argparse

def makeTwoWay(INPATH: str, OUTPATH: str):

    trimmed_df = pd.read_csv(INPATH, sep = "\t", header = None)
    dup = trimmed_df[[3, 4, 5, 0, 1, 2, 6, 7, 9, 8, 10, 11, 12, 13]]
    out = pd.concat([trimmed_df, dup], axis = 0).reset_index(drop=True)
    out.drop(13, axis=1, inplace=True)
    out.to_csv(OUTPATH, sep="\t", index=False, header=None)

    return

def annotate(INPATH: str, REF_REGIONS: str, OUTPATH: str):

    # ! Change this when not working in HPC ! #
    if subprocess.run("which bedtools", shell=True, capture_output=True).returncode != 0:
        logging.critical("BEDTools not detected. Abort.")
        raise ModuleNotFoundError("BEDTools not detected.")
    code = subprocess.run("bedtools intersect -a {} -b {} -wo > {}".format(INPATH, REF_REGIONS, OUTPATH), shell=True).returncode
    if code != 0:
        logging.error("Error in bedtools intersect.")
        sys.exit(-1)
    return

def condense(df, gene_col = 17):

    cols = list(range(df.shape[1]))
    cols.pop(gene_col)
    out = df.groupby(cols)[gene_col].apply(",".join).reset_index()
    out = out[[*range(0, df.shape[1])]]
    logging.info("Original df shape : {}; Condensed df shape : {}".format(df.shape, out.shape))

    return out

def main():

    INPATH = args.input
    REF_REGIONS = args.ref
    PID_PATH = args.list

    expanded_path = INPATH.rstrip("trimmed.bed") + ".homo.expanded.bed"
    anno_path = expanded_path.rstrip("bed") + "geneanno.bed"
    raw_pid_path = expanded_path.rstrip("bed") + "geneanno.PID.bed"
    condensed_pid_path = expanded_path.rstrip("bed") + "geneanno.PID.condensed.bed"

    # Make two-way map
    two_way = makeTwoWay(INPATH, expanded_path)

    # Gene annotation
    annotate(expanded_path, REF_REGIONS, anno_path)

    if not args.list:
        return

    # Gene/Region filtering
    PID_df = pd.read_csv(PID_PATH, sep='\t', encoding="utf-8")
    PID_df['Genetic defect'] = PID_df['Genetic defect'].apply(lambda x: x[:x.index(u'\xa0')] if u'\xa0' in x else x)
    PID_genes = set(PID_df['Genetic defect'].tolist())

    anno_df = pd.read_csv(anno_path, sep="\t", header=None)
    raw_pid_df = anno_df[anno_df[args.genecol].isin(PID_genes)].drop_duplicates()
    raw_pid_df.to_csv(raw_pid_path, sep="\t", index=False, header=None)

    condensed_pid_df = condense(raw_pid_df, gene_col=args.genecol)
    condensed_pid_df.to_csv(condensed_pid_path, sep="\t", index=False, header=None)
    logging.info("Writing filtered BED to " + condensed_pid_path)

    return

if __name__ == "__main__":
    # Argparse setup
    parser = argparse.ArgumentParser(description = "Make two-way map and extract genes of interest.")
    parser._optionals.title = "Options"
    parser.add_argument("-i", "--input", type = str, required = True, help = "trimmed BED file")
    parser.add_argument("-r", "--ref", type = str, required = True, help = "gene region list for annotation")
    parser.add_argument("-l", "--list", type = str, required = False, help = "list of genes/regions of interest")
    parser.add_argument("-c", "--genecol", type = int, default = 16, help = "0-based column index of \"Genetic defect\" in gene list (default: 16)")
    parser.add_argument("-v", "--verbose", type = str, default = "INFO", help = "verbosity level (default: INFO)")
    args = parser.parse_args()
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P',
                        level = args.verbose.upper())
    main()
