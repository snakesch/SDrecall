#!/usr/bin/env python3

# This module consists all different calculations of sequencing coverages

def read_evaluator(read, filter_tags, filter_logic, min_mapq):
    '''
    This function retrieves query names of valid reads.
    '''
    assert isinstance(filter_tags, list), f"Filter tags must be a list of strings. ({filter_tags})"
    
    filter_pass = True
    
    # Retrieve alignments with filter e.g. XA
    if len(filter_tags) > 0:
        if filter_logic == "or":
            filter_pass = any([read.has_tag(tag) for tag in filter_tags])
        elif filter_logic == "and":
            filter_pass = all([read.has_tag(tag) for tag in filter_tags])
        elif filter_logic == "not":
            filter_pass = not any([read.has_tag(tag) for tag in filter_tags])
        else:
            logger.warning(f"Unsupported filter logic ({filter_logic}). Omitting filter ... ")
    
    # Assign filter_pass as True if no filter is imposed
    filter_pass = filter_pass and True

    # Select good reads
    candidate_reads = not (read.is_duplicate or read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail or read.mapping_quality >= min_mapq) and filter_pass
    
    if not candidate_reads:
        return
    
    chrom, start, end = read.reference_name, read.reference_start, read.reference_end
    
    if not read.mate_is_unmapped and read.is_proper_pair:
        # Process each read once only
        if read.is_read1:
            # mate is present on the same contig
            if read.next_reference_name == chrom:
                # Per pysam API, next_reference_end is not available
                return (chrom, min(start, read.next_reference_start), max(end, read.next_reference_start + 150), read.query_name)
            else:
                return (chrom, start, end, read.query_name)
    elif not read.is_proper_pair and read.is_paired:
        return (chrom, start, end, read.query_name)
    elif not read.is_paired:
        return (chrom, start, end, read.query_name)
    

def calculate_inferred_coverage(bam_file, min_mapq = 10, 
                               filter_tags = [], filter_logic = "and",
                               genome_file = "ucsc.hg19.contigsize.genome",
                               target_region = ""):
    import pybedtools as pb
    import pysam
    
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
    if os.path.isfile(target_region):
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
    return pd.concat(base_coverages)   

# Non-inferred depth
def calculate_non_inferred_coverage(bam_file, min_mapq = 10, 
                               filter_tags = [], thread: int = thread,
                               target_region = ""):
    import pysam
    import pybedtools as pb
    
    bamfp = pysam.AlignmentFile(bam_file, "rb")
    
    # Case: Valid regions provided
    if os.path.isfile(target_region):
        bed_obj = pb.BedTool(target_region).sort().merge().as_intervalfile()
        for read in bamfp.fetch()
    pass





## Initialization 
filter_tags = 
min_mapq = 10

target_region = ""
target_tag = target_tag or os.path.basename(target_region).split(".")[0]