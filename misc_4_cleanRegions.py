import pandas as pd
import logging
import argparse
import re
import os

# Argparse setup
parser = argparse.ArgumentParser(description = "Make two-way map and extract genes of interest.")
parser._optionals.title = "Options"
parser.add_argument("-i", "--input", type = str, required = True, help = "condensed & trimmed BED file")
parser.add_argument("-c", "--genecol", type = int, default = 16, help = "0-based column index of gene name in input file (default: 16)")
parser.add_argument("-o", "--out", type = str, required = True, help = "output path")
parser.add_argument("-v", "--verbose", type = str, default = "INFO", help = "verbosity level (default: INFO)")
args = parser.parse_args()
logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P',
                    level = args.verbose.upper())

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
    
def main():
    
    condensed_df = pd.read_csv(args.input, sep="\t", header=None)
    print(condensed_df)
    one_way_df = one_way(condensed_df, col_1=[0, 1, 2, args.genecol], col_2=[3, 4, 5, args.genecol])
    sorted_df = one_way_df.sort_values(["chr", "start"]).reset_index(drop=True)
    grouped = sorted_df.groupby(["chr", "gene"], as_index=False).apply(intersectBED)
    print(grouped)
    total_bed = grouped.droplevel(level=0).reset_index(drop=True)
    os.makedirs(args.out, exist_ok=True)
    total_bed.to_csv(os.path.join(args.out, "all_homo_regions.bed"), sep="\t", index=False, header=False)
    logging.info("Successfully combined all homologous regions to " + os.path.join(args.out, "all_homo_regions.bed"))
    
    return

if __name__ == "__main__":
    main()