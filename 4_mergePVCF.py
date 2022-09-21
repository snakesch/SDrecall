#!/usr/bin/env python3

# 4_mergePVCF.py
# Description: Merge prioritized VCF and original VCF.
# Author: Yang XT, She CH (2022)
# Contact: Xingtian Yang (u3005579@connect.hku.hk), Louis She (snakesch@connect.hku.hk)

import os
import pandas as pd
import logging
import argparse
import gzip
from src.utils import *
import time

class Variant:
    """
    This class creates a Variant object for each record in a graph.

    This is designated for VCFs with one sample.

    """

    import pandas as pd
    from difflib import SequenceMatcher

    def __init__(self, row):

        """
        Variant constructor based on a pandas dataframe

        Note:
        -------
        Rows should be taken from a dataframe sorted by chromosome (chrM, chr1, ..., chrX, chrY), or by the order in fai file

        Error:
        -------
        KeyError: No INSLEN in INS records
        ValueError: No sample calls

        """

        self.CHROM = row["#CHROM"]
        self.START = int(row["POS"])
        self.ID = row["ID"]
        self.FILTER = row["FILTER"]
        self.QUAL = row["QUAL"]
        self.REF = row["REF"]
        self.ALT = row["ALT"].split(",")[0] # Only consider the biallelic case
        self.idx = int(row["POS"])
        self._row = row

        # Data from INFO field
        self._INFO = {field.split("=")[0]: field.split("=")[1:] for field in row["INFO"].split(";")}
        self.END = int(self._INFO["END"][0]) if "END" in self._INFO.keys() else row["POS"]
        self.VARTYPE = self._INFO["SVTYPE"][0] if "SVTYPE" in self._INFO.keys() else "short_variant"
        if self.VARTYPE == "INS":
            try:
                self.INSEND = self.START + int(self._INFO["INSLEN"][0])
            except KeyError:
                raise KeyError(f"No INSLEN field in INS variant.\n{row} ")
            self.INSSEQ = self._INFO["INSSEQ"] if "INSSEQ" in self._INFO.keys() else "N"

        # Data from sample field
        if row.shape[0] - 9 > 1:
            logging.warning(f"Multiple samples found. Only consider the first sample. ")
        elif row.shape[0] == 9:
            raise ValueError(f"No sample found in the given VCF for the variants.\n{self.__repr__()}")

        _sample = [row["FORMAT"].split(":")] + [row[col].split(":") for col in range(9, row.shape[0], 1)]
        self.SAMPLE = {sample[0]: sample[1] for sample in zip(*_sample)}
        try:
            self.DP = float(self.SAMPLE["DP"])
        except (KeyError, TypeError, ValueError): # Case: DP not available or DP = "."
            self.DP = "."

        try:
            self.GT = self.SAMPLE["GT"] if self.DP != "." and self.DP != 0.0 else "."
        except:
            raise KeyError(f"GT not found in the variant:\n{self.__repr__()}")

        self.MISSING_GT = (self.GT == "./.") or (self.GT == ".|.") or (self.GT == ".")

        if "AD" in self.SAMPLE.keys():
            self.AD = list(map(int, self.SAMPLE["AD"].split(",")))
            if self.AD[0] > 0 and self.AD == 0:
                self.GT = "0"

    def __repr__(self):
        """
        String representation of the variant.
        """
        return f"{self.CHROM}:{self.START}-{self.END} {self.REF}>{self.ALT}\nType: {self.VARTYPE}\nGenotype: {self.GT}\nSample: {self.SAMPLE}"

    def getSeries(self):

        keys = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", self._row.index[9]]
        try:
            values = [self.CHROM, self.START, self.ID, self.REF, self.ALT, self.QUAL, self.FILTER,
                 ";".join([f"{key}={value[0]}" for key, value in self._INFO.items()]),
                 ":".join(self.SAMPLE.keys()),
                 ":".join(self.SAMPLE.values())]
        except:
            values = [self.CHROM, self.START, self.ID, self.REF, self.ALT, self.QUAL, self.FILTER,
                 ".",
                 ":".join(self.SAMPLE.keys()),
                 ":".join(self.SAMPLE.values())]

        return pd.Series(values, index=keys)

    def addFilter(self, *tag):

        tags = list(tag)

        of = self.FILTER.split(";") + tags
        self.FILTER = ";".join(filter(lambda x: x, of))

    def isAdjacent(self, other):

        s_thresh = 1

        if self == other:
            return True
        else:
            l, r = sorted([self, other])
            return (r.START - l.END <= s_thresh)

    def overlap_fraction(self, other):

        l, r = sorted([self, other])
        if r.START < l.END:
            return abs((r.START - l.END) / (r.END - l.START))
        else:
            return 0.0

    def __eq__(self, other):

        f_thresh = 0.95

        if (self.VARTYPE != other.VARTYPE) or (self.CHROM != other.CHROM):
            return False
        elif self.VARTYPE == other.VARTYPE == "short_variant":
            return (self.START == other.START) and (self.END == other.END) \
                    and (self.REF == other.REF) and (self.ALT == other.ALT)
        elif self.VARTYPE == other.VARTYPE == "BND":
            return self.isAdjacent(other)
        else:
            if self.VARTYPE != "INS":
                return (self.overlap_fraction(other) >= f_thresh)
            elif not self.isAdjacent(other):
                return False
            else:
                ofrac = self.overlap_fraction(other)
                if ofrac == 0.0:
                    return False
                elif self.INSSEQ != "N" and other.INSSEQ != "N":
                    seqsim = SequenceMatcher(None, self.INSSEQ, other.INSSEQ).ratio()
                    return (seqsim >= f_thresh)
                else:
                    return (ofrac >= f_thresh)

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        return self.idx < other.idx

    def __le__(self, other):
        return self.idx <= other.idx

    def __gt__(self, other):
        return self.idx > other.idx

    def __ge__(self, other):
        return self.idx >= other.idx

def bsearch(var, low, high, x) -> tuple:
    """
    A modified binary search function for finding SDrecall variants.

    Return:
    -------
    Found: Variant object, index of the object
    Not found: None

    """
    if high >= low:
        mid = (high + low) // 2
        if var[mid] == x:
            return var[mid], mid
        elif var[mid] > x:
            return bsearch(var, low, mid - 1, x)
        elif var[mid] < x:
            return bsearch(var, mid + 1, high, x)
        else: # Case: Multiple variants with same POS (Look for variants in proximity)
            if not var[mid+1] > var[mid]:
                return bsearch(var, mid + 1, high, x)
            elif not var[mid-1] < var[mid]:
                return bsearch(var, 0, mid - 1, x)
            else:
                return None, None
    else:
        return None, None

def inList(row, dup_ov_vars: list[tuple]):

    for tup in dup_ov_vars:
        if (row["POS"] == tup[0]) and (row["REF"] == tup[1]) and (row["ALT"] == tup[2]):
            return True
    return False

def pick_rep_rec(row, ov_vars, pv_tag, ov_tag, keep_both) -> pd.Series:
    """
    This function processes variants which are found in both VCFs.
    """
    pv_var = Variant(row)
    ov_var, idx = bsearch(ov_vars, 0, len(ov_vars) - 1, pv_var)
    if ov_var:
        _tup_key = (ov_var.START, ov_var.REF, ov_var.ALT)
        dup_ov_vars.append(_tup_key) # Keep indices of found items
        if keep_both: # Keep both records if specified
            if ov_var.GT[::2] != pv_var.GT[::2]: # Ignore phasing status 
                dup_ov_vars.remove(_tup_key) # Keep records: Do not drop
            elif pv_var.FILTER != "PASS":
                dup_ov_vars.remove(_tup_key) # Keep records: Do not drop
        if pv_var.GT == ov_var.GT:
            pv_var.addFilter(ov_tag, pv_tag)
            return pv_var.getSeries()
        elif pv_var.GT == "0/0" or pv_var.GT == "0|0" or "." in pv_var.GT:
            ov_var.addFilter(ov_tag, pv_tag)
            return ov_var.getSeries()
        elif ov_var.GT == "0/0" or ov_var.GT == "0|0"  or "." in ov_var.GT:
            pv_var.addFilter(ov_tag, pv_tag)
            return pv_var.getSeries()
        else:
            pv_var.addFilter(ov_tag, pv_tag)
            return pv_var.getSeries()
    else: # Not found in SDrecall VCF
        logging.debug(f"The row is not found : {pv_var.CHROM}:{pv_var.START} with GT {pv_var.GT}. ")
        if pv_var.GT == "0/0" or pv_var.GT == "0|0"  or "." in pv_var.GT:
            logging.debug(f"Removed row {pv_var.CHROM}:{pv_var.START} with GT {pv_var.GT}. ")
            row["#CHROM"] = None # Returning None will cause incompatibility in series/dataframe output
            row["POS"] = None # Returning None will cause incompatibility in series/dataframe output
            return row
        else:
            pv_var.addFilter(pv_tag)
            return pv_var.getSeries()

def concat_headers(header_1, header_2) -> list[str]:
    """
    This function takes two VCF headers and returns a merged header.
    """
    # Get original order
    fields = [ field.split("=")[0].strip("#") for field in header_1 ] + [ field.split("=")[0].strip("#") for field in header_2 ]
    _unique_fields = list(dict.fromkeys(fields))
    order = dict(zip(_unique_fields, range(len(_unique_fields))))

    # Merge headers
    merged_header = sorted(list(set(header_1 + header_2)), key=lambda x: order[x.split("=")[0].strip("#")])

    return merged_header

def merge(pv_vcf, ov_vcf, pv_tag, ov_tag, keep_both) -> pd.DataFrame:
    """
    This function merges prioritized VCF with original VCF.

    Arguments:
    ----------
    pv_vcf: path of prioritized VCF
    ov_vcf: path of original VCF
    pv_tag: FILTER tag for pv_vcf
    ov_tag: FILTER tag for ov_vcf
    keep_both: keep both records in merged output

    """
    # Load VCF data
    ov_header, ov_subjects, ov_df = loadVCF(ov_vcf)
    pv_header, pv_subjects, pv_df = loadVCF(pv_vcf)

    # Initialize an empty dataframe
    pv_ov = pd.DataFrame()

    pgroups = pv_df.groupby("#CHROM")
    for pgroup in pgroups:
        chrom = pgroup[0] # only process one chr at a time
        precs = pgroup[1].sort_values(by=["POS"])

        # Initialize an empty list for storing START coordinate, REF and ALT of found variants
        global dup_ov_vars
        dup_ov_vars = []

        logging.info(f"Processing chromosome {chrom} ...  ")
        ov_slice = ov_df[ov_df["#CHROM"] == chrom]
        ov_vars = sorted([ Variant(orec) for idx, orec in ov_slice.iterrows()])
        from_pv = precs.apply(pick_rep_rec, axis=1, args=(ov_vars, pv_tag, ov_tag, keep_both)).dropna()
        _filter = ov_slice.apply(inList, axis=1, args=(dup_ov_vars, ))
        from_ov = ov_slice[~_filter]
        if from_ov.shape[0] == 0:
            merged = from_pv
        else:
            from_ov["FILTER"] = from_ov["FILTER"].apply(lambda x: f"{x};{ov_tag}" if len(x) > 1 else ov_tag)
            merged = pd.concat([from_pv, from_ov], axis=0)
        if merged.shape[0] != 0:
            pv_ov = pd.concat([pv_ov, merged], axis=0).reset_index(drop=True)

    # Avoid decimals in POS
    pv_ov["POS"] = pv_ov["POS"].astype(int)

    return pv_ov.drop_duplicates()

def main(pvcf, ovcf, outpath, pv_tag, ov_tag, keep_both):
    """
    Main function for merging prioritized variants and original variants.

    Arguments:
    ----------
    pvcf: Prioritized VCF holding variants from all contigs
    ov_vcf: Original VCF (eg. VCF from SDrecall)
    ov_tag: tag used for variants called from ov_vcf
    pv_tag: tag used for variants called from pvcf
    outpath: absolute path of output VCF (gzipped)
    keep_both: keep both records in merged output (BOOL)

    """
    start = time.time()
    os.chdir(os.path.dirname(pvcf))

    ov_header, ov_subjects, ov_df = loadVCF(ovcf)
    pv_header, pv_subjects = loadVCF(pvcf, omit_record=True)

    # Write new VCF header
    ov_filter_head = f"##FILTER=<ID={ov_tag},Description='variants called from {ov_tag}'>"
    pv_filter_head = f"##FILTER=<ID={pv_tag},Description='variants called from {pv_tag}'>"
    ov_header.append(ov_filter_head)
    ov_header.append(pv_filter_head)
    merged_header = concat_headers(pv_header, ov_header)
    if os.path.exists(outpath):
        os.remove(outpath)
    with gzip.open(outpath, "wb") as f:
        f.write("\n".join(merged_header).encode())
        f.write("\n".encode())
        f.write("\t".join(ov_df.columns.tolist()).encode())
        f.write("\n".encode())

    # Process prioritized variants
    merged_df = merge(pvcf, ovcf, pv_tag, ov_tag, keep_both)

    try:
        if merged_df.shape[0] != 0:
            merged_df.to_csv(outpath, sep="\t", index=False, header=False, mode="a", compression="gzip")
    except AttributeError: # Case: merged_df == None
        raise ValueError("No variants in the VCFs! ")

    logging.info(f"****************** VCF merged in {time.time() - start:.2f} seconds ******************")

if __name__ == "__main__":

    # Argparse setup
    parser = argparse.ArgumentParser(description = "Merge prioritized VCF and original VCF.")
    parser._optionals.title = "Options"
    ROOT = os.path.dirname(__file__)
    parser.add_argument("--pvcf", type = str, help = "prioritized VCF (gz)", required = True)
    parser.add_argument("--ovcf", type = str, help = "original VCF (gz)", required = True)
    parser.add_argument("--pv_tag", type = str, help = "tag used for variants from prioritized VCF", required = True)
    parser.add_argument("--ov_tag", type = str, help = "tag used for variants from original VCF", required = True)
    parser.add_argument("--outpath", type = str, help = "absolute output path of merged VCF (gz)", required = True)
    parser.add_argument("--keep_both", action = "store_true", help = "keep both records in the merged output")
    parser.add_argument("-v", "--verbose", type = str, default = "INFO", help = "verbosity level (default: INFO)")
    args = parser.parse_args()
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P',
                        level = args.verbose.upper())
    logging.debug(f"Working in {ROOT}")

    main(args.pvcf, args.ovcf, args.outpath, args.pv_tag, args.ov_tag, args.keep_both)
