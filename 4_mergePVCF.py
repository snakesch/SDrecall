#!/usr/bin/env python3

import os
import pandas as pd
import logging
import argparse
import gzip
from src.utils import *
from pandarallel import pandarallel as pa
import time

class Variant:
    """
    This class creates a Variant object for each record in a graph.
    
    This is designated for VCFs with one sample.
    
    """
    
    import pandas as pd
    from difflib import SequenceMatcher
    
    def __init__(self, row, idx):
        
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
        except (KeyError, TypeError): # Case: DP not available or DP = "."
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

        keys = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", row.index[9]]
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
        if self.FILTER == "PASS":
            self.FILTER = ";".join(tags)
        else:
            of = self.FILTER.split(";") + tags
            self.FILTER = ";".join(of)     
    
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
    Binary search function for finding SDrecall variants.
    
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
        else:
            return bsearch(var, mid + 1, high, x)
    else:
        return None, None

def pick_rep_rec(row, ov_vars, pv_tag, ov_tag) -> pd.Series:
    """
    This function processes variants which are found in both VCFs.
    """  
    pv_var = Variant(row, row.name)
    ov_var, idx = bsearch(ov_vars, 0, len(ov_vars) - 1, pv_var)
    if ov_var:
        dup_ov_vars.append(idx) # Keep indices of found items
        if pv_var.GT == ov_var.GT:
            pv_var.addFilter(ov_tag, pv_tag)
            return pv_var.getSeries()
        elif "0" in pv_var.GT or "." in pv_var.GT:
            ov_var.addFilter(ov_tag, pv_tag)
            return ov_var.getSeries()
        elif "0" in ov_var.GT or "." in ov_var.GT:
            pv_var.addFilter(ov_tag, pv_tag)
            return pv_var.getSeries()
        else:
            pv_var.addFilter(ov_tag, pv_tag)
            return pv_var.getSeries()
    else: # Not found in SDrecall VCF
        if "0" in pv_var.GT or "." in pv_var.GT:
            return None
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

def merge(pv_vcf, ov_vcf, pv_tag, ov_tag):
    """
    This function merges prioritized VCF with original VCF on a chromosome-wise basis.
    """
    ov_header, ov_subjects, ov_df = loadVCF(ov_vcf)
    pv_header, pv_subjects = loadVCF(pv_vcf, omit_record=True)

    # Get chromosome number
    chrom = pv_vcf.split(".")[-1]
    logging.info(f"Processing chromosome {chrom} ... ")
    
    # Prepare ov variants
    ov_vars = []
    dup_ov_vars = []
    for idx, row in ov_df[ov_df["#CHROM"] == chrom].iterrows():
        ov_vars.append(Variant(row, idx))
    ov_vars = sorted(ov_vars)

    header = ov_df.columns
    processed = 0
    from_pv = pd.DataFrame()
    with pd.read_csv(pv_vcf, sep="\t", na_filter=False, engine="c", 
                     comment="#", header=None, names=header, compression="gzip", 
                     skiprows=32959110, chunksize=500000) as reader:
        for pv_chunk in reader:
            processed += 500000
            from_pv_tmp = pv_chunk.parallel_apply(pick_rep_rec, axis=1, args=(ov_vars, pv_tag, ov_tag,)).dropna()
            if from_pv_tmp.shape[0] != 0:
                from_pv = pd.concat([from_pv, from_pv_tmp], axis=0)
            logging.info(f"Processed {processed} variants. ")

    from_ov = ov_df[~ov_df.index.isin(dup_ov_vars)]
    pv_ov = pd.concat([from_pv, from_ov], axis=0).reset_index(drop=True)
    
    return pv_ov

def main(pvcf, ov_vcf, outpath, pv_tag, ov_tag, workers):
    """
    Main function for merging prioritized variants and original variants.
    
    Arguments:
    ----------
    pvcf: Prioritized VCF holding variants from all contigs
    ov_vcf: Original VCF (eg. VCF from SDrecall)
    ov_tag: tag used for variants called from ov_vcf
    pv_tag: tag used for variants called from pvcf
    outpath: absolute path of output VCF (gzipped)
    workers: number of cores to use
        
    """
    start = time.time()
    # Load prioritized VCF and original VCF (as file streams)
    os.chdir(os.path.dirname(pvcf))
    os.makedirs("tmp", exist_ok=True)

    ov_header, ov_subjects, ov_df = loadVCF(ov_vcf)
    pv_header, pv_subjects = loadVCF(pvcf, omit_record=True)
    all_chr = [ pcontig[13:].split(",")[0] for pcontig in pv_header if pcontig.startswith("##contig") ]

    for contig in all_chr:
        cmd = f"bcftools filter -r {contig} -Oz -o {os.path.join('tmp/', os.path.basename(pvcf) + '.') + str(contig)} {pvcf} "
        executeCmd(cmd)
    logging.info(f"****************** VCF splitting completed in {time.time() - start:.2f} seconds ******************")
    
    start = time.time()
    os.chdir("tmp/") # Descend into tmp/
    pv_vcfs = [pv_vcf for pv_vcf in os.listdir()]

    # Write new VCF header
    ov_filter_head = f"##FILTER=<ID={ov_tag},Description='variants called from {ov_tag}'>"
    pv_filter_head = f"##FILTER=<ID={pv_tag},Description='variants called from {pv_tag}'>"
    ov_header.append(ov_filter_head)
    ov_header.append(pv_filter_head)
    merged_header = concat_headers(pv_header, ov_header)
    with gzip.open(outpath) as f:
        f.write("\n".join(merged_header).encode())
        f.write("\n".encode())
        f.write("\t".join(ov_df.columns.tolist()).encode())
        f.write("\n".encode())

    # Process prioritized variants
    pa.initialize(progress_bar=False, verbose=0, nb_workers=workers)

    for pv_vcf in pv_vcfs:
        merged_df = merge(pv_vcf, ov_vcf, pv_tag, ov_tag)
        merged_df.to_csv(outpath, sep="\t", index=False, header=False, mode="a", compression="gzip")
    logging.info(f"****************** VCF merged in {time.time() - start:.2f} seconds ******************")
    # Cleanup
    os.chdir(os.path.dirname(pvcf))
    os.remove("tmp/")

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
    parser.add_argument("--thread", type = int, help = "number of threads (default: 8)", default = 8)
    parser.add_argument("-v", "--verbose", type = str, default = "INFO", help = "verbosity level (default: INFO)")
    args = parser.parse_args()
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P',
                        level = args.verbose.upper())
    logging.debug(f"Working in {ROOT}")  
    
    main(args.pvcf, args.ovcf, args.outpath, args.pv_tag, args.ov_tag, args.thread)