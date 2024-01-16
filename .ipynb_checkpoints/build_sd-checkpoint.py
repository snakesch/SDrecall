#!/usr/bin/env python3

# This module consists of functions used in building SD maps. 

def build_SD_pair_from_pe_reads(pe_reads, avg_insert_size, std_insert_size, conf = 0.95, full = True):
    '''
    This function builds SD pairs from a pair of PE reads. 
    Specifying full = True will return a fully resolved SD map, with fully resolved SD regions; 
    specifying full = False will skip SD resolution, retaining the entire SD region.
    
    Expected return values:
    1) full = True: [('chr5:69348646-69348794', 'chr5:70224065-70224212'), ('chr5:69348646-69348794', 'chr5:70223709-70223808'), ('chr5:69348290-69348438', 'chr5:70224065-70224212'), ('chr5:69348290-69348438', 'chr5:70223709-70223808')]
    2) full = False: [('chr5:69348646-69348788', 'chr5:70223709-70224212'), ('chr5:69348290-69348788', 'chr5:70223709-70224212')]
    '''
    from itertools import product
    import re
    
    from multialign import xa_serialize
    from read_properties import calculate_tlen_per_pair
    from read_filter import p_value
    
    assert len(set(pe_reads["qname"].tolist())) == 1, f"Paired-end reads with different query names. {set(pe_reads['qname'].tolist())}"
    assert pe_reads.shape[0] == 2, f"Unpaired reads detected. ({pe_reads['qname']})"
    
    fr_xas, rr_xas = pe_reads[pe_reads["first_in_pair"]], pe_reads[~pe_reads["first_in_pair"]]
    fr_xa, rr_xa = fr_xas["XA"].apply(xa_serialize).values[0], rr_xas["XA"].apply(xa_serialize).values[0]
    
    ## Filter XA mappings
    valid_pairs = {}
    
    for fr, rr in product(fr_xa, rr_xa):
        tlen, genomic_span = calculate_tlen_per_pair(fr, rr)
        
        ## Check for extreme template lengths, direction of alternative alignments
        if p_value(tlen, avg_insert_size, std_insert_size) > 1 - conf and {fr[1][0], rr[1][0]} == {"-", "+"}:
            valid_pairs[(fr, rr)] = (tlen, genomic_span)
    
    ## Output genomic spans
    if full:
        valid_fr_xas, valid_rr_xas = [ k[0] for k in valid_pairs.keys() ], [ k[1] for k in valid_pairs.keys() ]
        ## Primary alignment statistics
        primary_alignment_end = pe_reads["pos"] + pe_reads["read_length"]
        primary_alignment_genomic_span = pe_reads["rname"] + ":" + pe_reads["pos"].astype(str) + "-" + primary_alignment_end.astype(str)
        
        ## Multialignment spans
        fr_xa_spans = [ fr_xa[0] + ":" + fr_xa[1].lstrip("+-") + "-" + str(int(fr_xa[1].lstrip("+-")) + sum(int(i) for i in re.findall('(\d+)[MDN=X]', fr_xa[-2]))) for fr_xa in valid_fr_xas ]
        rr_xa_spans = [ rr_xa[0] + ":" + rr_xa[1].lstrip("+-") + "-" + str(int(rr_xa[1].lstrip("+-")) + sum(int(i) for i in re.findall('(\d+)[MDN=X]', rr_xa[-2]))) for rr_xa in valid_rr_xas ]
        
        full_read_map = [ r for r in product(primary_alignment_genomic_span, fr_xa_spans + rr_xa_spans) ]
        return full_read_map
    else:        
        half_read_map = [ r for r in product(pe_reads["ref_cov_span"], [ v[1] for v in valid_pairs.values() ]) ] 
        return half_read_map
    
def groupby_qnames(large_chunk, avg_insert_size, std_insert_size, confidence, full):
           
    ret = large_chunk.groupby("qname").apply(build_SD_pair_from_pe_reads, avg_insert_size, std_insert_size, confidence, full)
    return ret

def validate_sd_map(sd_map):
    
    import pandas as pd
    
    ## Check contigs
    contig_pattern = r'(chr)*(X|Y|MT*|[0-9][0-9]*)$'
    sd_map = sd_map[sd_map["chr"].str.match(contig_pattern) & sd_map["chr_counterpart"].str.match(contig_pattern)]
    involved_contigs = set()
    involved_contigs.add(sd_map["chr"].unique())
    involved_contigs.add(sd_map["chr_counterpart"].unique())
    assert len(involved_contigs) <= 25, "More than 25 contigs present. Any alternative contigs?"
    
    

def build_sd_pairs(bam_file, confidence, threads, full=True):
    '''
    Expected output is a list containing mappings of primary alignments and multialignments.
    [('chr5:69348646-69348794', 'chr5:70224065-70224212'),
      ('chr5:69348646-69348794', 'chr5:70223709-70223808'),
      ('chr5:69348290-69348438', 'chr5:70224065-70224212'),
      ('chr5:69348290-69348438', 'chr5:70223709-70223808')]
    '''    
    from read_properties import get_fragment_dist
    from read_filter import p_value
    
    import multiprocessing as mp
    import pandas as pd
    from itertools import repeat
    
    # Step 1: Load JSON file
    multialignments = pd.read_json("selected_regions.json", lines=True)
    
    # Step 2: Filter out discordant pairs, unpaired reads, supplementary alignments and alignments with TLEN deviating significantly from expectation
    multialignments = multialignments[multialignments["flag"] & 1 == 1] # 0x1: paired reads
    multialignments = multialignments[multialignments["flag"] & 256 == 0] # 0x256: not primary alignments
    
    discordant_pairs = multialignments.loc[multialignments["tlen"] == 0, "qname"].drop_duplicates().to_list()
    multialignments = multialignments.loc[~multialignments["qname"].isin(discordant_pairs), :]
    
    # avg_insert_size, std_insert_size = get_fragment_dist(bam_file)
    avg_insert_size, std_insert_size = 563.0, 138.5
    extreme_reads = multialignments["tlen"].abs().apply(p_value, args = (avg_insert_size, std_insert_size,)) < 0.001
    multialignments = multialignments[~extreme_reads]
    
    # Step 3: Compute the length of bases aligned to in the reference genome and the read.
    multialignments["genome_consume_length"] = multialignments["cigar"].str.extractall('(\d+)[MDN=X]')[0].astype(int).groupby(level=0).sum()
    multialignments["read_length"] = multialignments["cigar"].str.extractall('(\d+)[MIS=X]')[0].astype(int).groupby(level=0).sum()
    
    # Step 4: Annotate XA, positive strand first in pair for easy manipulation
    multialignments["XA"] = multialignments["tags"].apply(lambda d: d.get("XA", ""))
    multialignments["pos_strand"] = multialignments["flag"] & 16 == 0
    multialignments["first_in_pair"] = multialignments["flag"] & 64 == 0
    
    # Step 5: Compute the region spanned by individual fragments
    multialignments["pnext_isDownstream"] = (multialignments["pos"] <= multialignments["pnext"])
    pnext_upstream, pnext_downstream = multialignments[~multialignments["pnext_isDownstream"]], multialignments[multialignments["pnext_isDownstream"]]
    
    pnext_upstream["end_pos"] = pnext_upstream["pnext"] + pnext_upstream["tlen"].abs()
    pnext_downstream["end_pos"] = pnext_downstream["pos"] + pnext_downstream["tlen"].abs()
    
    multialignments = pd.concat([pnext_upstream, pnext_downstream])
    multialignments["ref_cov_span"] = multialignments["rname"] + ":" + multialignments["pos"].astype(str) + "-" + multialignments["end_pos"].astype(str)
    
    ## Filter out null fragment intervals (END <= START)
    multialignments = multialignments[multialignments["end_pos"] > multialignments["pos"]]
    
    multialignments = multialignments.sort_values("qname")
    
    # Step 6: For each pair of PE reads, build a corresponding SD map
    n_read = multialignments.shape[0]
    n_chunk = n_read // 50
    chunks = [ multialignments.iloc[i:min(i+50, n_read), :] for i in range(0, n_read, 50) ]
    
    with mp.Pool(min(n_chunk + 1, threads)) as pool:
        ma = pool.starmap(groupby_qnames, zip(chunks, repeat(avg_insert_size), repeat(std_insert_size), repeat(confidence), repeat(full)))
    
    ## Filter the output regions
    ma_df = pd.DataFrame(ma[0])
    filtered_ma = [ m for m in ma_df[0].tolist() if m] 
    
    # Flatten the list of lists into a list of tuples
    flattened_ma = [item for sublist in filtered_ma for item in sublist]

    sd_map = pd.DataFrame(flattened_ma, columns=['col1', 'col2'])

    # Split the 'col1' and 'col2' values into separate columns
    sd_map[['chr', 'start', 'end']] = sd_map['col1'].str.split(':|-', expand=True)
    sd_map[['chr_counterpart', 'start_counterpart', 'end_counterpart']] = sd_map['col2'].str.split(':|-', expand=True)
    sd_map.drop(['col1', 'col2'], axis=1, inplace=True)

    sd_map = sd_map[['chr', 'start', 'end', 'chr_counterpart', 'start_counterpart', 'end_counterpart']]

    return sd_map