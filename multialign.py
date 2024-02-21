import logging
logger = logging.getLogger("root")

def get_multialign_regions(bam_file, min_mapq: int = 40, min_dp: int = 3, unambiguous_dp: int = 10, multialign_fraction = 0.7, 
                           target_region = None, use_inferred_coverage = True):
    '''
    For each base, we compute (A) raw sequencing coverage, (B) coverage of high MQ reads (MQ >= min_mapq), and (C) coverage of XA reads.
    We aims to extract bases satisfying:
    i.   A >= min_dp
    ii.  C/A >= multialign_fraction
    iii. B <= unambiguous_dp
    
    '''
    import multiprocess as mp
    import pandas as pd
    import pybedtools as pb
    
    if use_inferred_coverage:
        from coverages import calculate_inferred_coverage
        func = calculate_inferred_coverage
    else:
        from coverages import calculate_non_inferred_coverage
        func = calculate_non_inferred_coverage
    
    # Structure of coverage calculation parameters:
    #              | bam_file, min_mapq, filter_tags, target_region, use_inferred_coverage
    #              | ---------------------------------------------------------------------
    #     raw_depth| bam_file,        0,        None, target_region, use_inferred_coverage
    # high_mq_depth| bam_file, min_mapq,        None, target_region, use_inferred_coverage
    #  raw_xa_depth| bam_file,        0,      ["XA"], target_region, use_inferred_coverage
    #  raw_xs_depth| bam_file,        0,      ["XS"], target_region, use_inferred_coverage
    
    params = [(bam_file,        0,          [], target_region), 
              (bam_file, min_mapq,          [], target_region),
              (bam_file,        0,      ["XA"], target_region),
              (bam_file,        0,      ["XS"], target_region)]
    logger.info("Begin depth computation. ")
    with mp.Pool(4) as pool:
        raw_depth, high_mq_depth, raw_xa_depth, raw_xs_depth = pool.starmap(func, params)

    raw_depth = raw_depth.rename(columns={"depth": "raw_depth"})
    high_mq_depth = high_mq_depth.rename(columns={"depth": "high_MQ_depth"})
    raw_xa_depth = raw_xa_depth.rename(columns={"depth": "raw_XA_depth"})
    raw_xs_depth = raw_xs_depth.rename(columns={"depth": "raw_XS_depth"})
    logger.debug("Completed empirical depth calculations. ")
    merged_depth = raw_depth.merge(high_mq_depth, how="left", on=["chr", "position"]).merge(raw_xa_depth, how="left", on=["chr", "position"]).merge(raw_xs_depth, how="left", on=["chr", "position"]).drop_duplicates()
    merged_depth = merged_depth.fillna(0)
    logger.info("Depth calculation completed. ")
    
    # Now we filter the bases that do not fulfill the standards below: (AND logic)
    # 1. Has at least 5 reads covered. (no MQ considered)
    # 2. Has 60% of the overlapping reads with XA/XS tag
    # 3. Has no more than 10 MQ >= 50 reads covered

    min_depth = merged_depth["raw_depth"] >= min_dp
    most_xa = (merged_depth["raw_XA_depth"]/merged_depth["raw_depth"]) >= multialign_fraction
    most_xs = (merged_depth["raw_XS_depth"]/merged_depth["raw_depth"]) >= multialign_fraction
    not_enough_evidence = merged_depth["high_MQ_depth"] <= unambiguous_dp
        
    merged_depth = merged_depth.loc[min_depth & (most_xa | most_xs) & not_enough_evidence, :].drop_duplicates()

    merged_depth["start"] = merged_depth.loc[:, "position"] - 1
    merged_depth["end"] = merged_depth.loc[:, "position"]
    
    selected_regions = pb.BedTool.from_dataframe(merged_depth.loc[:, ["chr", "start", "end"]]).merge()
    selected_regions.saveas("selected_regions.bed")
    
    assert merged_depth.shape[0] > 0, "No valid multialigned regions. "
    
    return selected_regions

def bam_subset(bam_file, mq_cutoff = 20):
    '''
    This function extracts multialigned regions from per-sample BAM.
    '''
    import subprocess
    import os
    
    ## sambamba is significantly faster than pysam even without multithreading
    logger.debug("Begin BAM subsetting. ")
    
    # Step 1: Extract query name file
    cmd = f"module load sambamba; sambamba -q view -L selected_regions.bed -F '[SA] != null and mapping_quality <= {mq_cutoff + 10}' -f sam -t 8 {bam_file} | cut -f1 | sort -u > selected_regions.qname"
    
    proc = subprocess.run(cmd, shell = True, text = True, capture_output = True)
    
    assert proc.returncode == 0, "Sambamba error. "
    
    # Step 2: Extract relevant reads
    cmd = f"module load sambamba; samtools view -N selected_regions.qname -@ 7 -u {bam_file} | sambamba -q view -f json -o selected_regions.json -t 8 /dev/stdin"
    proc = subprocess.run(cmd, shell = True, text = True, capture_output = True)
    
    assert proc.returncode == 0, "Sambamba error. "
    
    if os.path.isfile("selected_regions.qname"):
        os.remove("selected_regions.qname")
    logger.debug("Completed BAM subsetting. ")
    
    return 

def xa_serialize(xa_string):
    '''
    This function takes a long XA string chrX,-1776046,148M,3;chrX,-1776046,148M,3; and serializes it to [('chrX', '-1776046', '148M', '3'), ('chrX', '-1776046', '148M', '3')].
    '''
    return [tuple(segment.split(",")) for segment in xa_string.split(";") if segment]