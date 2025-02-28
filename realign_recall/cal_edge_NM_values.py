import pysam
import numpy as np
import scipy.stats as stats

from src.log import logger

def calculate_NM_distribution_poisson(bam, conf_level=0.01, sample_size=3000000, logger=logger):
    '''
    The input vcf_rec is simply a tuple containing(chr, start, ref, alt), remember the coordinates from VCF are 1-indexed
    '''
    nm_array = np.empty(sample_size)
    bam_handle = pysam.AlignmentFile(bam, "rb")
    n = 0
    for read in bam_handle.fetch():
        if read.is_proper_pair and \
           not read.is_secondary and \
           not read.is_supplementary and \
           not read.is_duplicate and \
           not read.is_qcfail and \
           len([t for t in read.get_tags() if t[0] in ["XA"]]) == 0 and \
           read.mapping_quality == 60:
            gap_sizes = [t[1] for t in read.cigartuples if t[0] in [1,2] and t[1] > 1]
            max_gap_size = max(gap_sizes) if len(gap_sizes) > 0 else 0
            edit_dist = read.get_tag("NM")
            scatter_edit_dist = edit_dist - max_gap_size
            nm_array[n] = scatter_edit_dist
            n += 1
            if n >= sample_size:
                break

    bam_handle.close()

    # Deal with the edge case where the input BAM does not have enough reads
    if n < sample_size:
        logger.warning(f"Input BAM {bam} does not have enough reads, only {n} reads are used for NM distribution calculation")
        nm_array = nm_array[:n]
    
    nm_mean = np.mean(nm_array)
    cutoff = 0
    while stats.poisson.cdf(cutoff, nm_mean) < 1 - conf_level:
        cutoff += 1
    logger.info(f"NM distribution for {bam} has mean {nm_mean}, the value for {1-conf_level} percentile is {cutoff}")

    return cutoff, nm_mean