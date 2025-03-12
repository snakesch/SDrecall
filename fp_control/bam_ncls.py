import pysam
import numpy as np
from ncls import NCLS

from fp_control.numba_operators import fast_median, numba_sum
from src.log import logger

def overlapping_reads_iterator(ncls_dict, read_dict, chrom, start, end):
    """
    Generator function to lazily yield overlapping read objects.

    Parameters:
    -----------
    ncls_dict : dict
        Dictionary of NCLS objects, keyed by chromosome names.
    read_dict : dict
        Dictionary of read objects, keyed by query name indices.
    qname_dict : dict
        Dictionary mapping query name indices to query names.
    chrom : str
        Chromosome name.
    start : int
        Start position of the interval (0-based, inclusive).
    end : int
        End position of the interval (0-based, exclusive).

    Yields:
    -------
    pysam.AlignedSegment
        Overlapping read object.

    Notes:
    ------
    This function uses NCLS for efficient overlap queries and yields individual reads.
    It filters reads to ensure they actually overlap with the given interval.
    """

    # Perform an overlap query on the NCLS tree
    if chrom not in ncls_dict or ncls_dict[chrom] is None:
        yield from ()
    else:
        _, overlapping_read_qnames = ncls_dict[chrom].all_overlaps_both(np.array([start]), np.array([end]), np.array([0]))
        # overlapping_read_qnames = ncls_dict[chrom].find_overlap(start, end)
        for qname_idx in list(dict.fromkeys(overlapping_read_qnames)):
            for read in read_dict[qname_idx]:
                if read.reference_start < end and read.reference_end > start:
                    yield read



def overlap_qname_idx_iterator(ncls_dict, chrom, start, end):
    """
    Generator function to lazily yield query name indices of overlapping reads.

    Parameters:
    -----------
    ncls_dict : dict
        Dictionary of NCLS objects, keyed by chromosome names.
    chrom : str
        Chromosome name.
    start : int
        Start position of the interval (0-based, inclusive).
    end : int
        End position of the interval (0-based, exclusive).

    Yields:
    -------
    int
        Query name index of an overlapping read.

    Notes:
    ------
    This function uses NCLS for efficient overlap queries.
    The start position is inclusive, and the end position is exclusive.
    """

    # Perform an overlap query on the NCLS tree
    if chrom not in ncls_dict or ncls_dict[chrom] is None:
        yield from ()
    else:
        _, overlapping_read_qnames = ncls_dict[chrom].all_overlaps_both(np.array([start]), np.array([end]), np.array([0]))
        # overlapping_read_qnames = ncls_dict[chrom].find_overlap(start, end)
        for qname_idx in dict.fromkeys(overlapping_read_qnames):
            yield qname_idx



def calculate_mean_read_length(bam_file_path, sample_size=100000):
    bamfile = pysam.AlignmentFile(bam_file_path, "rb")

    total_length = 0
    read_count = 0

    for read in bamfile:
        total_length += read.query_length
        read_count += 1
        if read_count >= sample_size:
            break

    bamfile.close()

    if read_count == 0:
        return 0  # Handle the case where there are no reads

    mean_length = total_length / read_count
    return mean_length



def is_read_noisy(read, paired, mapq_filter, basequal_median_filter, filter_noisy = True):
    """Helper function to determine if a read is noisy based on various criteria."""
    if paired:
        # Use short circuit evaluation to speed up the function
        if read.is_secondary or \
           read.is_supplementary or \
           read.mapping_quality < mapq_filter or \
           read.is_qcfail or \
           read.is_unmapped or \
           read.reference_end is None or \
           read.reference_end - read.reference_start < 75 or \
           read.reference_name != read.next_reference_name or \
           not read.is_proper_pair or \
           read.query_sequence is None:
            return True
        
        if read.query_qualities is not None:
            return (fast_median(np.array(read.query_qualities, dtype=np.uint8)) <= basequal_median_filter or
                    numba_sum(np.array(read.query_qualities, dtype=np.uint8) < basequal_median_filter) >= 50) and filter_noisy
    else:
        if read.is_secondary or \
           read.is_supplementary or \
           read.is_duplicate or \
           read.is_unmapped or \
           read.mapping_quality < mapq_filter or \
           read.is_qcfail or \
           read.reference_end is None or \
           read.query_sequence is None:
            return True
        if read.query_qualities is not None:
            return (fast_median(np.array(read.query_qualities, dtype=np.uint8)) <= basequal_median_filter or
                    numba_sum(np.array(read.query_qualities, dtype=np.uint8) < basequal_median_filter) >= 40) and filter_noisy
    return False



def migrate_bam_to_ncls(bam_file,
                        mapq_filter = 25,
                        basequal_median_filter = 20,
                        paired = True,
                        filter_noisy = True,
                        logger = logger):
    """
    Migrate BAM file data to NCLS (Nested Containment List) format for efficient interval querying.

    This function filters reads based on quality criteria, and organizes the data
    into NCLS structures for fast overlap queries.

    It completely put all alignment info into memory, thus faster than pysam disk IO based query.

    Parameters:
    -----------
    bam_file : str
        Path to the input BAM file.
    mapq_filter : int, optional
        Minimum mapping quality threshold for reads (default is 25).
    basequal_median_filter : int, optional
        Minimum median base quality threshold for reads (default is 20).
    paired : bool, optional
        If True, process as paired-end data; if False, as single-end (default is True).
    logger : logging.Logger, optional
        Logger object for output messages.

    Returns:
    --------
    tuple
        A tuple containing:
        - ncls_dict : dict
            Dictionary of NCLS objects, keyed by chromosome names.
        - read_dict : dict
            Dictionary of read objects, keyed by query name indices.
        - qname_dict : dict
            Dictionary mapping query name indices to query names. (qname_idx -> qname) The qname_idx is determined by the growing variable n with the iteration of all reads
        - qname_idx_dict : dict
            Dictionary mapping query names to query name indices. (qname -> qname_idx)
        - noisy_qnames : set
            Set of query names identified as noisy and filtered out.

    Within an NCLS object:
    (start, end) ---ncls_dict[chrom]---> qname_idx(interval_idx) ---qname_dict---> qname
    qname_idx(interval_idx) ---read_dict---> read objects(pysam.AlignedSegment)
    qname ---qname_idx_dict---> qname_idx(interval_idx) (which is also the interval index)

    Notes:
    ------
    - The function uses pysam for BAM file processing and NCLS for interval data structure.
    - Reads are filtered based on mapping quality, median base quality, and other criteria.
    - For paired-end data, additional checks ensure proper pairing and orientation.
    - The NCLS structure allows for O(log n + m) time complexity for overlap queries, where n is the
      number of intervals and m is the number of overlaps found.

    Example:
    --------
    >>> bam_file = "path/to/aligned_reads.bam"
    >>> ncls_dict, read_dict, qname_dict, qname_idx_dict, noisy_qnames = migrate_bam_to_ncls(bam_file)
    >>> # Now you can use ncls_dict for efficient interval queries
    >>> chrom = "chr1"
    >>> start, end = 1000, 2000
    >>> overlapping_interval_indices = ncls_dict[chrom].find_overlap(start, end)
    >>> overlapping_reads = [read_dict[idx] for idx in overlapping_interval_indices]

    See Also:
    ---------
    pysam.AlignmentFile : For BAM file handling
    ncls.NCLS : Nested Containment List Structure for interval queries
    """

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        chroms = bam.references
        read_dict = {}
        qname_interval_dict = {chrom: {} for chrom in chroms}
        ncls_dict = {chrom: None for chrom in chroms}
        qname_idx_dict = {}
        qname_dict = {}

        n = 0
        noisy_qnames = set()

        for read in bam:
            qname = read.query_name
            if qname in noisy_qnames:
                continue

            if is_read_noisy(read, paired, mapq_filter, basequal_median_filter, filter_noisy):
                logger.debug(f"This qname {qname} is noisy. Skip it.")
                noisy_qnames.add(qname)
                continue

            chrom = read.reference_name
            start = read.reference_start
            end = read.reference_end

            if qname in qname_idx_dict:
                qname_idx = qname_idx_dict[qname]
            else:
                qname_idx = n
                qname_idx_dict[qname] = qname_idx
                qname_dict[qname_idx] = qname
                n += 1

            read_dict[qname_idx] = read_dict.get(qname_idx, []) + [read]

            if qname_idx in qname_interval_dict[chrom]:
                interval = qname_interval_dict[chrom][qname_idx]
                updated_interval = (min(start, int(interval[0])), max(end, int(interval[1])))
            else:
                updated_interval = (int(start), int(end))
            qname_interval_dict[chrom][qname_idx] = updated_interval

    for chrom in chroms:
        qname_intervals = qname_interval_dict[chrom]
        if len(qname_intervals) > 0:
            qname_indices = np.fromiter(qname_intervals.keys(), dtype=np.int32)
            starts = np.array([int(qname_intervals[qname_idx][0]) for qname_idx in qname_indices])
            ends = np.array([int(qname_intervals[qname_idx][1]) for qname_idx in qname_indices])
            ncls_dict[chrom] = NCLS(starts, ends, qname_indices)

    return ncls_dict, read_dict, qname_dict, qname_idx_dict, noisy_qnames