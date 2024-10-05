import pandas as pd
import numpy as np
import multiprocessing as mp
import logging
from scipy.stats import norm
from itertools import product, repeat

logger = logging.getLogger("SDrecall")

def extract_sd_coordinates_from_json(bam_json, avg_frag_size=400, std_frag_size=100, err_rate = 0.05, threads=1):
    
    df = pd.read_json(bam_json, lines=True)
    '''
    This function further refines paired SD coordinates, and identifies unpaired reads and discordant read pairs. 
    JSON follows typical format of SAM. 

    Returns:
    comma-separated SD pair coordinates e.g. chrX,1408045,1408193,chrY,1357307,1357878
    '''
    # Extract the align length on ref genome for the read 
    df["genome_consume_length"] = df['cigar'].str.extractall(r'(\d+)[MDN=X]')[0].astype(int).groupby(level=0).sum()
    
    df["read_length"] = df['cigar'].str.extractall(r'(\d+)[MIS=X]')[0].astype(int).groupby(level=0).sum()

    df["XA"] = df["tags"].apply(lambda d: d.get("XA", ""))
    
    df["pos_strand"] = df["flag"].astype(int) & 16 == 0
    
    df["first_in_pair"] = df["flag"].astype(int) & 64 == 0
    
    primary_alignment = df["flag"].astype(int) & 256 == 0
    
    # Only retain primary alignments for downstream analysis
    df = df.loc[primary_alignment, :]
    
    # Filter out the reads mapped to different chromosomes
    discordant_chr_pairs = set(df.loc[df["tlen"].astype(int) == 0, "qname"].drop_duplicates().to_list())
    df = df.loc[~df["qname"].isin(discordant_chr_pairs), :]
    
    # Drop unpaired reads
    group_counts = df.groupby('qname')['qname'].transform('count')
    df = df[group_counts == 2]

    # Calculate the coverage span of the fragment
    df["end_pos"] = 0
    pos_strand = df.loc[:, "pos_strand"]
    pnext_downstream = df.loc[:, "pos"].astype(int) <= df.loc[:, "pnext"].astype(int)
    pnext_upstream = np.logical_not(pnext_downstream)
    rev_strand = np.logical_not(pos_strand)
    df.loc[pos_strand & pnext_downstream, "end_pos"] = df.loc[pos_strand & pnext_downstream, "pos"].astype(int) + df.loc[pos_strand & pnext_downstream, "tlen"].astype(int).abs()
    df.loc[pos_strand & pnext_upstream, "end_pos"] = df.loc[pos_strand & pnext_upstream, "pnext"].astype(int) + df.loc[pos_strand & pnext_upstream, "tlen"].astype(int).abs()
    df.loc[rev_strand & pnext_upstream, "end_pos"] = df.loc[rev_strand & pnext_upstream, "pnext"].astype(int) + df.loc[rev_strand & pnext_upstream, "tlen"].astype(int).abs()
    df.loc[rev_strand & pnext_downstream, "end_pos"] = df.loc[rev_strand & pnext_downstream, "pos"].astype(int) + df.loc[rev_strand & pnext_downstream, "tlen"].astype(int).abs()
    null_fragment_interval = df["end_pos"].astype(int) <= df["pos"].astype(int)
    if df.loc[null_fragment_interval, :].shape[0] > 0:
        logger.warning("{} alignments have END larger than START (indicated by TLEN), possibly null genomic intervals:\n{}".format(df.loc[null_fragment_interval, :].shape[0],
                                                                                                                                   df.loc[null_fragment_interval, :].to_string(index=False)))
    df = df.loc[~null_fragment_interval, :]
    df.loc[:, "ref_cov_span"] = df.loc[:, "rname"] + ":" + df.loc[:, "pos"].astype(str) + "-" + df.loc[:, "end_pos"].astype(str)
    
    # Now we need to exclude the extreme tlen cause they might overlap with strcutrual variant event
    extreme_tlen_size = df["tlen"].astype(int).abs().apply(p_value, args=(avg_frag_size, std_frag_size,)) < 0.001
    if df.loc[extreme_tlen_size, :].shape[0] > 0:
        logger.warning("{} alignments have extreme tlen sizes:\n{}\n".format(df.loc[extreme_tlen_size, :].shape[0],
                                                                         df.loc[extreme_tlen_size, :].head(5).to_string(index=False)))
        df = df.loc[~extreme_tlen_size, :]
    
    # Drop the singleton reads
    group_counts = df.groupby('qname')['qname'].transform('count')
    df = df.loc[group_counts == 2, :]
    
    logger.info("After filtering, {} alignments remained.".format(df.shape[0]))
    
    sorted_df = df.sort_values("qname")
    chunks = split_dataframe(sorted_df, threads)

    with mp.Pool(min(threads, 4)) as pool:
        results = pool.starmap(extract_sd_pair_per_chunk, zip(chunks, repeat(avg_frag_size), repeat(std_frag_size), repeat(err_rate)))

    return "\n".join(results)
    
def extract_sd_pair_per_chunk(chunk_df, avg_frag_size=400, std_frag_size=100, err_rate=0.05):
    sd_pairs = chunk_df.groupby("qname")[["qname", "XA", "first_in_pair", "ref_cov_span", "tlen"]].apply(lambda rp: process_per_pair_return_tuple(rp, avg_frag_size=avg_frag_size, std_frag_size=std_frag_size, err_rate=err_rate))
    cleaned_sd_pairs = set(filter(lambda x: x is not None and x != '', sd_pairs.values.tolist()))
    return "\n".join(cleaned_sd_pairs)

def split_dataframe(df, num_chunks):
    
    num_rows = len(df)
    base_chunk_size = num_rows // num_chunks  # Base chunk size

    # If base_chunk_size is odd, decrement it by 1
    if base_chunk_size % 2 != 0:
        base_chunk_size -= 1

    # Create an array of chunk sizes
    chunk_sizes = np.full(num_chunks, base_chunk_size)

    # Distribute the remaining rows to the last chunk
    remainder = num_rows - base_chunk_size * num_chunks
    chunk_sizes[-1] += remainder

    # Split the DataFrame into chunks
    chunks = np.split(df, np.cumsum(chunk_sizes)[:-1])

    return chunks


def p_value(value, mean, std_dev):
    z_score = (value - mean)/std_dev
    p_value = 2 * (1- norm.cdf(abs(z_score)))
    return p_value


def calculate_tlen_per_pair(xa1, xa2):
    '''
    Takes a pair of XA reads as strings chr12,+125396989,148M,14 (chr, (strand) pos, cigar, edit distance)
    
    Returns:
    tlen: inferred template length (int)
    inferred_span: inferred genomic span (string; chr:start-end)
    '''
    import re

    xa1_list = xa1.strip().split(",")
    xa2_list = xa2.strip().split(",")
    chr1 = xa1_list[0]
    chr2 = xa2_list[0]

    # discordant pairs
    if chr1 != chr2:
        return 0, ""
    else:
        pos1 = int(xa1_list[1])
        pos2 = int(xa2_list[1])
        # improper pair with same orientation
        if pos1 * pos2 >= 0:
            return 0, ""
        else:
            inferred_span = chr1 + ":"
            if pos1 < pos2:
                rev_strand_leftmost = - pos1
                pos_strand_leftmost = pos2
                if rev_strand_leftmost >= pos_strand_leftmost:
                    tlen = rev_strand_leftmost - pos_strand_leftmost + sum(int(x) for x in re.findall(r'(\d+)[MDN=X]', xa1_list[2]))
                    inferred_span += str(pos_strand_leftmost) + "-" + str(pos_strand_leftmost + tlen)
                else:
                    tlen = pos_strand_leftmost - rev_strand_leftmost + sum(int(x) for x in re.findall(r'(\d+)[MDN=X]', xa2_list[2]))
                    inferred_span += str(rev_strand_leftmost) + "-" + str(rev_strand_leftmost + tlen)
            else:
                rev_strand_leftmost = - pos2
                pos_strand_leftmost = pos1
                if rev_strand_leftmost >= pos_strand_leftmost:
                    tlen = rev_strand_leftmost - pos_strand_leftmost + sum(int(x) for x in re.findall(r'(\d+)[MDN=X]', xa2_list[2]))
                    inferred_span += str(pos_strand_leftmost) + "-" + str(pos_strand_leftmost + tlen)
                else:
                    tlen = pos_strand_leftmost - rev_strand_leftmost + sum(int(x) for x in re.findall(r'(\d+)[MDN=X]', xa1_list[2]))
                    inferred_span += str(rev_strand_leftmost) + "-" + str(rev_strand_leftmost + tlen)
            return tlen, inferred_span
        

def process_per_pair_return_tuple(group_df, avg_frag_size = 400, std_frag_size = 100, err_rate = 0.05):

    assert len(group_df.loc[:, "qname"].drop_duplicates().values) == 1, f"Identified read pairs with different query names: \n{group_df.to_string()}\n"
    assert len(group_df) == 2, "This error should not be triggered. Previous qname filters are defective. "

    if (group_df.loc[:, "XA"].str.len() > 0).sum() < 2:
        return None
    else:
        group_df = group_df.assign(XA_list = group_df.loc[:, "XA"].str.strip(";").str.split(";"))
        
        fr_xas = group_df.loc[group_df.loc[:,"first_in_pair"], "XA_list"].values[0]
        rr_xas = group_df.loc[~group_df.loc[:, "first_in_pair"], "XA_list"].values[0]
        
        results = { rp: calculate_tlen_per_pair(*rp) for rp in product(fr_xas, rr_xas) } # read pair (rp) looks like: ("chr12,+125396989,148M,14", "chr12,+125396989,148M,14"), first from read1, second from read2
        # Valid read pairs should not have extreme template lengths, and should have one fwd and one rev
        valid_pairs = { k:v for k,v in results.items() if p_value(v[0], avg_frag_size, std_frag_size) >= err_rate and v[0] >= 100 and {k[0].split(",")[1][0], k[1].split(",")[1][0]} == {"+", "-"}}

        consolidated_coordinates = "\n".join([ ",".join(t).replace(":", ",").replace("-", ",") for t in product(group_df.loc[:, "ref_cov_span"].values, [v[1] for k, v in valid_pairs.items()]) ])
        return consolidated_coordinates

            