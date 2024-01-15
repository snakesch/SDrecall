#!/usr/bin/env python3

# This module consists purpose-specific read filters.

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
    candidate_reads = not (read.is_duplicate or read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail or read.mapping_quality < min_mapq) and filter_pass
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