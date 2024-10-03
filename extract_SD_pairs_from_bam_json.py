import pandas as pd
import numpy as np
import re
import multiprocessing as mp
import logging
from scipy.stats import norm
from itertools import product, repeat

logger = logging.getLogger("SDrecall")

def main_parse_json_and_process(bam_json, avg_frag_size=400, std_frag_size=100, conf_level = 0.05, threads=1, read_SD_extract = True):
    df = pd.read_json(bam_json, lines=True) # This table is the json table 
    '''
    1. qname
    2. flag
    3. rname
    4. pos
    5. mapq
    6. cigar
    7. rnext
    8. pnext
    9. tlen
    10. seq
    11. qual: A list of integers, one integer stands for the phred-scale value of base quality
    12. tags: A dictionary per row
    '''
    
    # Extract the align length on ref genome for the read 
    df["genome_consume_length"] = df['cigar'].str.extractall('(\d+)[MDN=X]')[0].astype(int).groupby(level=0).sum()
    
    # Extract the align length on read
    df["read_length"] = df['cigar'].str.extractall('(\d+)[MIS=X]')[0].astype(int).groupby(level=0).sum()
    
    # Extract the XA infos to str
    df["XA"] = df["tags"].apply(lambda d: d.get("XA", ""))
    
    # Setup a column to tell whether the read is mapped to positive strand (python bitwise operator in a vectorized way)
    df["pos_strand"] = df["flag"].astype(int) & 16 == 0
    
    # Setup a column to tell whether the read is the first read in the pair
    df["first_in_pair"] = df["flag"].astype(int) & 64 == 0
    
    # Setup a column to tell whether this alignment row is primary alignment
    primary_alignment = df["flag"].astype(int) & 256 == 0
    
    # Only retain primary alignment records for downstream analysis
    df = df.loc[primary_alignment, :]
    
    # Filter out the read pairs mapped to different chromosomes
    inter_chr_pairs = set(df.loc[df["tlen"].astype(int) == 0, "qname"].drop_duplicates().to_list())
    df = df.loc[~df["qname"].isin(inter_chr_pairs), :]
    
    # Drop the singleton reads
    group_counts = df.groupby('qname')['qname'].transform('count')
    df = df[group_counts == 2]

    # Now that calculate the coverage span of the fragment
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
    logger.warning("These {} alignments suggest a null genomic interval after calculating the fragment coverage region in reference genome based on TLEN:\n{}\n".format(df.loc[null_fragment_interval, :].shape[0],
                                                                                                                                                                        df.loc[null_fragment_interval, :].to_string(index=False)))
    df = df.loc[~null_fragment_interval, :]
    df.loc[:, "ref_cov_span"] = df.loc[:, "rname"] + ":" + df.loc[:, "pos"].astype(str) + "-" + df.loc[:, "end_pos"].astype(str)
    
    # Now we need to exclude the extreme tlen cause they might overlap with strcutrual variant event
    extreme_tlen_size = df["tlen"].astype(int).abs().apply(p_value, args=(avg_frag_size, std_frag_size,)) < 0.001
    logger.warning("These {} alignments have an extreme tlen size:\n{}\n".format(df.loc[extreme_tlen_size, :].shape[0],
                                                                                 df.loc[extreme_tlen_size, :].iloc[:100, :].to_string(index=False)))
    df = df.loc[~extreme_tlen_size, :]
    
    # Drop the singleton reads
    group_counts = df.groupby('qname')['qname'].transform('count')
    df = df.loc[group_counts == 2, :]
    
    logger.info("After all the filtering, the remaining dataframe has {} rows and it looks like this:\n{}\n".format(df.shape[0],
                                                                                                                    df[:5].to_string(index=False)))
    
    sorted_df = df.sort_values("qname")

    chunks = split_dataframe(sorted_df, threads)
    pool = mp.Pool(threads)
    results = pool.starmap(extract_sd_pair_per_chunk, zip(chunks, repeat(avg_frag_size), repeat(std_frag_size), repeat(conf_level), repeat(read_SD_extract)))
    result_str = "\n".join(results)
    return result_str
    
        

def extract_sd_pair_per_chunk(chunk_df, avg_frag_size=400, std_frag_size=100, conf_level=0.05, read_SD_extract = True):
    result_sd_pairs = set()
    for i in range(0, len(chunk_df), 2):
        group_df = chunk_df.iloc[i:i+2, :]
        try:
            result_str = process_per_pair_return_tuple(group_df, avg_frag_size=avg_frag_size, std_frag_size=std_frag_size, conf_level=conf_level, output_read_SD=read_SD_extract)
        except AssertionError:
            logger.error(f"These few lines contain a singleton read: \n{chunk_df.iloc[i-2:i+2, :].to_string(index=False)}\n")
            raise AssertionError
        if result_str:
            result_sd_pairs.add(re.sub(r";|-|:", ",", result_str))
    return "\n".join(result_sd_pairs)


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
    # xa1 and xa2 format example: chr12,+125396989,148M,14 (chr, pos, cigar, edit_distance)
    xa1_list = xa1.strip().split(",")
    xa2_list = xa2.strip().split(",")
    chr1 = xa1_list[0]
    chr2 = xa2_list[0]
    if chr1 != chr2:
        return 0, ""
    else:
        pos1 = int(xa1_list[1])
        pos2 = int(xa2_list[1])
        if pos1 * pos2 >= 0:
            return 0, ""
        else:
            region_str = chr1 + ":"
            if pos1 < pos2:
                rev_strand_leftmost = - pos1
                pos_strand_leftmost = pos2
                if rev_strand_leftmost >= pos_strand_leftmost:
                    tlen = rev_strand_leftmost - pos_strand_leftmost + sum(int(x) for x in re.findall(r'(\d+)[MDN=X]', xa1_list[2]))
                    region_str += str(pos_strand_leftmost) + "-" + str(pos_strand_leftmost + tlen)
                else:
                    tlen = pos_strand_leftmost - rev_strand_leftmost + sum(int(x) for x in re.findall(r'(\d+)[MDN=X]', xa2_list[2]))
                    region_str += str(rev_strand_leftmost) + "-" + str(rev_strand_leftmost + tlen)
            else:
                rev_strand_leftmost = - pos2
                pos_strand_leftmost = pos1
                if rev_strand_leftmost >= pos_strand_leftmost:
                    tlen = rev_strand_leftmost - pos_strand_leftmost + sum(int(x) for x in re.findall(r'(\d+)[MDN=X]', xa2_list[2]))
                    region_str += str(pos_strand_leftmost) + "-" + str(pos_strand_leftmost + tlen)
                else:
                    tlen = pos_strand_leftmost - rev_strand_leftmost + sum(int(x) for x in re.findall(r'(\d+)[MDN=X]', xa1_list[2]))
                    region_str += str(rev_strand_leftmost) + "-" + str(rev_strand_leftmost + tlen)
            return tlen, region_str
        

def process_per_pair_return_tuple(group_df, avg_frag_size = 400, std_frag_size = 100, conf_level = 0.05, output_read_SD=True):
    try:
        assert len(group_df.loc[:, "qname"].drop_duplicates().values) == 1
    except AssertionError:
        logger.error(f"This two-line dataframe contains different query names: \n{group_df.to_string()}\n")
        raise AssertionError
    if (group_df.loc[:, "XA"].str.len() > 0).sum() < 2:
        return None
    else:
        group_df = group_df.assign(XA_list = group_df.loc[:, "XA"].str.strip(";").str.split(";"))
        tlen = abs(group_df.loc[:, "tlen"].values[0])
        assert len(group_df) == 2
        fr_xas = group_df.loc[group_df.loc[:,"first_in_pair"], "XA_list"].values[0]
        rr_xas = group_df.loc[~group_df.loc[:, "first_in_pair"], "XA_list"].values[0]
        
        results = { tuple: calculate_tlen_per_pair(*tuple) for tuple in product(fr_xas, rr_xas) } # Tuple looks like: ("chr12,+125396989,148M,14", "chr12,+125396989,148M,14"), first from read1, second from read2
        # Filter the pairs
        valid_pairs = { k:v for k,v in results.items() if p_value(v[0], avg_frag_size, std_frag_size) >= conf_level and v[0] >= 100 and {k[0].split(",")[1][0], k[1].split(",")[1][0]} == {"+", "-"}}
        
        if output_read_SD:
            valid_fr_xas = [ k[0] for k, v in valid_pairs.items() ]
            valid_rr_xas = [ k[1] for k, v in valid_pairs.items() ]
            
            # Now instead of report the fragment span corresponding ship, we report the read corresponding ship. Lets try this:
            # 1. Get the ref span of the read
            primary_read_end = group_df["pos"] + group_df["read_length"]
            group_df["primary_read_span"] = group_df["rname"] + ":" + group_df["pos"].astype(str) + "-" + primary_read_end.astype(float).astype(int).astype(str)
            
            # 2. Get the span of the XA alignments
            fr_xa_spans = [ fr_xa.split(",")[0] + ":" + fr_xa.split(",")[1].lstrip("+-") + "-" + str(int(float(fr_xa.split(",")[1].lstrip("+-")) + sum(int(i) for i in re.findall('(\d+)[MDN=X]', fr_xa.split(",")[-2])))) for fr_xa in valid_fr_xas ]
            rr_xa_spans = [ rr_xa.split(",")[0] + ":" + rr_xa.split(",")[1].lstrip("+-") + "-" + str(int(float(rr_xa.split(",")[1].lstrip("+-")) + sum(int(i) for i in re.findall('(\d+)[MDN=X]', rr_xa.split(",")[-2])))) for rr_xa in valid_rr_xas ]
            
            first_read_map = [";".join(t) for t in product(group_df.loc[group_df.loc[:,"first_in_pair"], "primary_read_span"].values, fr_xa_spans)]
            second_read_map = [";".join(t) for t in product(group_df.loc[~group_df.loc[:,"first_in_pair"], "primary_read_span"].values, rr_xa_spans)]
            
            result_str = "\n".join(first_read_map + second_read_map)
            return result_str
        else:
            result_str = "\n".join([ ";".join(t) for t in product(group_df.loc[:, "ref_cov_span"].values, [v[1] for k, v in valid_pairs.items()]) ])
            return result_str
            
            