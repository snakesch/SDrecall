#!/usr/bin/env python3

# This module consists all different calculations of sequencing coverages

def calculate_inferred_coverage(bam_file, min_mapq = 10, 
                               filter_tags = [], target_region = "", filter_logic = "and",
                               genome_file = "ucsc.hg19.contigsize.genome"):
    import pybedtools as pb
    import pysam
    import os
    import pandas as pd
    
    from read_filter import read_evaluator
    
    '''
    This function is intended for identification of candidate realignment regions where template length >> 300bp.
    
    Some expected inputs: 
    filter_tags: a list of tags to search
    filter_logic: i. and: match all filter tags
                 ii.  or: match any one filter tags
                iii. not: not match any filter tags
    genome_file: provided, do not change
    target_region: BED file path of regions of interest
    '''
    
    # Load regions of interest if provided
    bamfp = pysam.AlignmentFile(bam_file, "rb")
    
    candidate_regions = []
    
    # Case: Valid regions provided
    if target_region and os.path.isfile(target_region):
        bed_obj = pb.BedTool(target_region).sort().merge().as_intervalfile()
        next_interval = next(bed_obj)
        try:
            while next_interval:
                for read in bamfp.fetch(next_interval.chrom, next_interval.start, next_interval.end):
                    candidate_region = read_evaluator(read, filter_tags, filter_logic, min_mapq)
                    if candidate_region != None:
                        candidate_regions.append(candidate_region)
                next_interval = next(bed_obj)
        except StopIteration:
            pass
    # Case: No region specified
    else:
        for read in bamfp.fetch():
            candidate_region = read_evaluator(read, filter_tags, filter_logic, min_mapq)
            if candidate_region != None:
                candidate_regions.append(candidate_region)
                
    candidate_regions = list(set(candidate_regions))
    bamfp.close()
    
    # Acquire base coverage statistics
    genome_coverage = pb.BedTool(candidate_regions).sort().genome_coverage(bg=True, g = genome_file)
    base_coverages = []
    for interval in genome_coverage:
        bases = [(interval.chrom, i, interval.name) for i in range(interval.start, interval.end+1, 1)]
        base_coverages.append(pd.DataFrame.from_records(bases))

    ret = pd.concat(base_coverages)
    ret.columns = ["chr", "position", "depth"]
    ret["depth"] = ret["depth"].astype(int)
    return ret

def calculate_non_inferred_coverage(bam_file, min_mapq: int = 10,
                                    filter_tags = [], target_region = "", thread: int = 8):
    '''
    Expected inputs:
    filter_tags: a list of tags to search
    target_region: BED file path of regions of interest
    '''
    import subprocess
    
    cmd = []
    if os.path.isfile(target_region):
        cmd.extend(["sambamba", "-q", "view", "-L", target_region])
        for tag in filter_tags:
            cmd.extend(["-F", f"'[{tag}] != null'"])
        cmd.extend(["|"])
    
    # Do we really need to include bases with zero coverage?
    cmd.extend(["samtools", "depth", "-J", "-q", "13", "-s", "-Q", str(min_mapq), bam_file])
    
    p = subprocess.run(cmd, text=True, capture_output=True)
    rows = p.stdout.strip().split("\n")
    data = [ row.split("\t") for row in rows ]
    
    ret = pd.DataFrame(data)
    ret.columns = ["chr", "position", "depth"]
    ret["depth"] = ret["depth"].astype(int)
    
    return ret
