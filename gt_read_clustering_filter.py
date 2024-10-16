import os
# Below two env variables are validated to control the number of threads used by numpy and numba, without them parallel run this script might cause CPU leakage
# os.environ["OMP_NUM_THREADS"] = "24"
# os.environ["NUMBA_NUM_THREADS"] = "24"
# os.environ["OPENBLAS_NUM_THREADS"] = "20"
# os.environ["VECLIB_MAXIMUM_THREADS"] = "20"
# os.environ["NUMEXPR_NUM_THREADS"] = "20"
# os.environ["MKL_NUM_THREADS"] = "20"
# os.environ["TBB_NUM_THREADS"] = "24"

import gc
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
import graph_tool.all as gt
import numpy as np
import pandas as pd
import pybedtools as pb
import pysam
import logging
import argparse as ap
import numba
# numba.config.THREADING_LAYER = 'omp'
# numba.set_num_threads(4)
from numba import types, prange, get_num_threads
from io import StringIO
from scipy import sparse
from ncls import NCLS
from collections import defaultdict
from intervaltree import IntervalTree
from src.utils import executeCmd, prepare_tmp_file
from fp_control.clique_identify import clique_generator_per_component, find_components_inside_filtered_cliques

bash_utils_hub = "shell_utils.sh"


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(module)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)



class HashableAlignedSegment(pysam.AlignedSegment):
    def __init__(self, aligned_segment):
        # Copy attributes from the input AlignedSegment object
        self.read = aligned_segment

    def __getattr__(self, attr):
        return getattr(self.read, attr)

    def __eq__(self, other):
        if not isinstance(other, HashableAlignedSegment):
            return False
        return self.query_name == other.query_name and self.flag == other.flag

    def __hash__(self):
        return hash((self.query_name, self.flag))




def calculate_mean_read_length(bam_file_path, sample_size=1000000):
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



def migrate_bam_to_ncls(bam_file,
                        mapq_filter = 25,
                        basequal_median_filter = 20,
                        paired = True,
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

    bam = pysam.AlignmentFile(bam_file, "rb")
    chroms = bam.references
    # print(chroms)
    read_dict = {}
    qname_interval_dict = {chrom: dict({}) for chrom in chroms}
    ncls_dict = {chrom: None for chrom in chroms}
    qname_idx_dict = {}
    qname_dict = {}

    n = 0
    noisy_qnames = set([])

    for read in bam:
        if read.query_name in noisy_qnames:
            continue

        if paired:
            if read.is_secondary or \
               read.is_supplementary:
                logger.info(f"This alignment is not primary ")
                continue
            elif read.mapping_quality < mapq_filter or \
                 read.is_qcfail or \
                 read.reference_end is None or \
                 read.reference_end - read.reference_start < 80 or \
                 read.reference_name != read.next_reference_name or \
                 not read.is_proper_pair or \
                 fast_median(np.array(read.query_qualities, dtype=np.uint8)) <= basequal_median_filter or \
                 numba_sum(np.array(read.query_qualities, dtype=np.uint8) < basequal_median_filter) >= 50:
                logger.warning(f"This qname {read.query_name} is noisy. Skip it.")
                noisy_qnames.add(read.query_name)
                continue
        else:
            if read.is_secondary or \
               read.is_supplementary:
                logger.debug(f"This alignment is not primary")
                continue
            elif read.is_duplicate or \
                 read.mapping_quality < mapq_filter or \
                 read.is_qcfail or \
                 read.reference_end is None:
                logger.warning(f"This qname {read.query_name} is noisy. Skip it.")
                noisy_qnames.add(read.query_name)
                continue
            if read.query_qualities is not None:
                if fast_median(np.array(read.query_qualities, dtype=np.uint8)) <= basequal_median_filter or \
                   numba_sum(np.array(read.query_qualities, dtype=np.uint8) < basequal_median_filter) >= 40:
                    noisy_qnames.add(read.query_name)
                    continue

        chrom = read.reference_name
        start = read.reference_start
        end = read.reference_end
        qname = read.query_name

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
            qname_interval_dict[chrom][qname_idx] = updated_interval
        else:
            updated_interval = (int(start), int(end))
            qname_interval_dict[chrom][qname_idx] = updated_interval

    bam.close()

    for chrom in chroms:
        qname_intervals = qname_interval_dict[chrom]
        qname_indices = np.fromiter(qname_intervals.keys(), dtype=np.int32)
        if len(qname_intervals) > 0:
            # print(chrom, qnames.shape, type(qnames))
            starts = np.array([int(qname_intervals[qname_idx][0]) for qname_idx in qname_indices])
            ends = np.array([int(qname_intervals[qname_idx][1]) for qname_idx in qname_indices])
            ncls_dict[chrom] = NCLS(starts, ends, qname_indices)

    return ncls_dict, read_dict, qname_dict, qname_idx_dict, noisy_qnames



def get_overlapping_reads(ncls_dict, read_dict, qname_dict, chrom, start, end):
    """
    Retrieve overlapping reads for a given genomic interval using NCLS.

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

    Returns:
    --------
    list
        List of overlapping read objects (pysam.AlignedSegment).

    Notes:
    ------
    This function is deprecated and not used in the current implementation.
    Use overlapping_reads_generator instead for better memory efficiency.
    """

    if chrom not in ncls_dict:
        return []

    if ncls_dict[chrom] is None:
        return []

    # Perform an overlap query on the NCLS tree
    _, overlapping_read_qnames = ncls_dict[chrom].all_overlaps_both(np.array([start]), np.array([end]), np.array([0]))
    overlapping_reads = [read for qname_idx in list(dict.fromkeys(overlapping_read_qnames)) for read in read_dict[qname_idx] if read.reference_start < end and read.reference_end > start]

    return overlapping_reads


def lazy_get_overlapping_qname_idx(ncls_dict, read_dict, qname_dict, chrom, start, end):
    """
    Generator function to lazily yield query name indices of overlapping reads.

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
    int
        Query name index of an overlapping read.

    Notes:
    ------
    This function uses NCLS for efficient overlap queries.
    The start position is inclusive, and the end position is exclusive.
    """

    # Perform an overlap query on the NCLS tree
    if chrom not in ncls_dict:
        yield from ()
    elif ncls_dict[chrom] is None:
        yield from ()
    else:
        _, overlapping_read_qnames = ncls_dict[chrom].all_overlaps_both(np.array([start]), np.array([end]), np.array([0]))
        # overlapping_read_qnames = ncls_dict[chrom].find_overlap(start, end)
        for qname_idx in dict.fromkeys(overlapping_read_qnames):
            yield qname_idx




def check_edge(u, v, adj_set):
    '''
    1 means not sure
    2 means accept same haplotype
    -1 means reject same haplotype

    Cant use 0 here because python treats 0 as False
    '''
    return adj_set.get((u, v), False) or adj_set.get((v, u), False)







def get_hapvector_from_cigar(cigar_tuples,
                             query_sequence = None,
                             logger = logger):
    '''
    Convert CIGAR tuples to a haplotype vector with penalties for matches, mismatches, and gaps.
    Note that the length of the haplotype vector is strictly aligned with the covered reference genome span, only this way we can ensure that comparison between overlapping reads are strictly aligned to each other.
    Uses NumPy for efficient array operations.

    Example:
    --------
    >>> cigar_tuples = [(7, 10), (8, 1), (1, 2), (7, 5)]
    >>> query_sequence = "ACGTACGTACTATATCGA"
    >>> hapvector = get_hapvector_from_cigar(cigar_tuples, query_sequence)
    >>> print(hapvector)
    [1 1 1 1 1 1 1 1 1 1 99 1 1 1 1 1 ]

    Arguments:
    - cigar_tuples: List of tuples representing CIGAR operations. (generated from pysam.AlignedSegment.cigartuples)
    - query_sequence: The query sequence of the read. (generated from pysam.AlignedSegment.query_sequence)
    - logger: Logger for logging information or errors.

    Returns:
    - A NumPy array representing the haplotype vector.
    '''

    # Determine the length of the reference sequence consumption
    ref_length = sum([length for operation, length in cigar_tuples if operation in {0, 7, 8, 2, 3}])

    # Create a haplotype vector of zeros
    hapvector = np.empty(ref_length, dtype=np.int32)

    index = 0
    query_pos = 0
    insertion = False
    for operation, length in cigar_tuples:
        assert operation != 0, f"The CIGAR string require separate mismatch from match. But current one does not: {cigar_tuples}"
        if operation == 3:
            # N, skipping for ref seq
            if not insertion:
                hapvector[index:index + length] = 1
            else:
                hapvector[index] = insertion
                if length > 1:
                    hapvector[index+1:index + length] = 1
                insertion = False
            index += length
        if operation == 7:  # '=' Match (assuming operation code 7 for match)
            query_pos += length
            if not insertion:
                hapvector[index:index + length] = 1
            else:
                hapvector[index] = insertion
                if length > 1:
                    hapvector[index+1:index + length] = 1
                insertion = False
            index += length
        elif operation == 8:  # 'X' Mismatch
            true_mismatch = True
            if query_sequence is not None:
                base = query_sequence[query_pos:query_pos + length]
                if base == "N":
                    true_mismatch = False

            if not insertion:
                if true_mismatch:
                    hapvector[index:index + length] = -4
                else:
                    hapvector[index:index + length] = 1
            else:
                hapvector[index] = insertion
                if length > 1:
                    if true_mismatch:
                        hapvector[index+1:index + length] = -4
                    else:
                        hapvector[index+1:index + length] = 1
                insertion = False
            index += length
            query_pos += length
        elif operation == 1:  # 'I' Insertion
            query_pos += length
            if index > 0:  # Modify the base immediately before the insertion
                insertion = length * 4
        elif operation == 2:  # 'D' Deletion
            if not insertion:
                hapvector[index:index + length] = -6
            else:
                hapvector[index] = insertion
                if length > 1:
                    hapvector[index+1:index + length] = -6
                insertion = False
            index += length
        elif operation == 4:
            # Soft clipping
            query_pos += length

    return hapvector



def get_errorvector_from_cigar(read, cigar_tuples, logger=logger):
    '''
    read.reference_start position is 0-indexed and it is including the soft-clipped bases
    Generate an error vector from a read's CIGAR string and base qualities.

    This function creates an error vector representing the probability of errors
    at each reference position within the read's alignment. The original base qualities are extracted by the pysam.AlignedSegment.query_qualities method:

        read sequence base qualities, including soft clipped bases (None if not present).
        Quality scores are returned as a python array of unsigned chars. Note that this is not the ASCII-encoded value typically seen in FASTQ or SAM formatted files. Thus, no offset of 33 needs to be subtracted.
        Note that to set quality scores the sequence has to be set beforehand as this will determine the expected length of the quality score array.
        This method raises a ValueError if the length of the quality scores and the sequence are not the same.

    Parameters:
    -----------
    read : pysam.AlignedSegment
        The read object containing alignment information.
    cigar_tuples : list of tuples
        List of CIGAR operations, each tuple contains (operation, length).
    logger : logging.Logger, optional
        Logger object for output messages.

    Returns:
    --------
    numpy.ndarray
        An array of float values representing error probabilities for each base in the read's alignment that is aligned to the reference sequence.

    Notes:
    ------
    - The function assumes that the CIGAR string separates matches (7) from
      mismatches (8).
    - Insertions and deletions are initially marked with a placeholder value (99),
      which is later converted to 0 (zero err probability) in the error vector.
    - The final error vector is the raw error probability.

    Raises:
    -------
    AssertionError
        If the CIGAR string contains unexpected operations or if the resulting
        error vector length doesn't match the read's aligned length.
    AssertionError:
        If the returned numerical vector does not have the same length as the specified genomic interval size

    Example:
    --------
    >>> read = pysam.AlignedSegment()
    >>> read.cigartuples = [(7, 10), (8, 1), (1, 2), (7, 5)]
    >>> read.query_qualities = [30] * 18
    >>> error_vector = get_errorvector_from_cigar(read, read.cigartuples)
    >>> print(error_vector)
    [0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.0 0.001 0.001 0.001 0.001 0.001]
    '''
    errorvector = np.empty(read.reference_end - read.reference_start, dtype=float)
    # An array of phred-scaled integer quality scores
    base_qualities = np.array(read.query_qualities, dtype=np.int32)
    query_consume = 0
    ref_consume = 0

    for operation, length in cigar_tuples:
        assert operation != 0, f"The CIGAR string require separate mismatch from match. But current one does not: {cigar_tuples}"
        if operation == 3:
            #skip for ref
            errorvector[ref_consume:ref_consume + length] = 99
            ref_consume += length
        elif operation == 4:
            #Soft clipping
            query_consume += length
        elif operation == 7:
            # Match
            errorvector[ref_consume:ref_consume + length] = base_qualities[query_consume:query_consume + length]
            query_consume += length
            ref_consume += length
        elif operation == 8:
            # Mismatch
            errorvector[ref_consume:ref_consume + length] = base_qualities[query_consume:query_consume + length]
            query_consume += length
            ref_consume += length
        elif operation == 1:
            # Insertion, we assign 99 to the base at the reference position immediately followed by the insertion
            # The base quality does not affect the gaps appeared in the alignment, later 99 will be replaced by 0
            errorvector[ref_consume] = 99
            query_consume += length
        elif operation == 2:
            # Deletion, we assign 99 to the bases at the reference positions within the deletion
            # The base quality does not affect the gaps appeared in the alignment, later 99 will be replaced by 0
            errorvector[ref_consume:ref_consume + length] = np.array([99] * length)
            ref_consume += length

    assert errorvector.size == read.reference_end - read.reference_start, f"The error vector length is {len(errorvector)} while the read length is {read.reference_end - read.reference_start}. The cigar str is {read.cigarstring}"
    # Convert the phred_scaled base quality scores back to error probability, 99 will be replaced by 0
    errorvector = np.where(errorvector == 99, 0, 10**(-errorvector/10))

    return np.array(errorvector)



def get_read_id(read):
    """
    Generate a unique identifier for a read.

    This function creates a unique identifier for a read by combining its query name
    and flag. The flag is included because:

    1. It distinguishes between reads in a pair: For paired-end sequencing, each pair
       of reads will have the same query name but different flags.

    2. It provides information about the read's properties: The flag contains important
       information such as whether the read is paired, mapped, secondary alignment, etc.

    3. It ensures uniqueness: In cases where the same read might appear multiple times
       (e.g., secondary alignments), the flag helps to differentiate these instances.

    Parameters:
    -----------
    read : pysam.AlignedSegment
        The read object from which to generate the ID.

    Returns:
    --------
    str
        A string in the format "query_name:flag" that uniquely identifies the read.
    """
    return f"{read.query_name}:{read.flag}"




def prepare_ref_query_idx_map(qseq_ref_pos_arr):
    '''
    This function is used to create a mapping from reference genomic positions to query sequence positions
    input qseq_ref_pos_arr is an array of reference genomic positions corresponding to each base in the query sequence
    for example:
    query sequence: A,C,G,T,A,C,G,T (every base matched to the ref sequence)
    ref positions: 10, 11, 12, 13, 14, 15, 16, 17
    then the input array will be: [10, 11, 12, 13, 14, 15, 16, 17]
    The output dict will be: {10: 0, 11: 1, 12: 2, 13: 3, 14: 4, 15: 5, 16: 6, 17: 7}

    if there is a deletion in the cigar string (align query against the ref sequence):
    query sequence: A,C,G,T,A,C,G,T (Between 3rd and 4th base there contains a deletion)
    ref positions: 10, 11, 12, 15, 16, 17, 18, 19
    then the input array will be: [10, 11, 12, 15, 16, 17, 18, 19]
    The output dict will be: {10: 0, 11: 1, 12: 2, 15: 3, 16: 4, 17: 5, 18: 6, 19: 7}

    if there is an insertion in the cigar string (align ref against the query sequence):
    query sequence: A,C,G,T,A,C,G,T (Both the 3rd and the 4th base are inserted comparing to the ref sequence)
    ref positions: 10, 11, -1, -1, 12, 13, 14, 15
    then the input array will be: [10, 11, -1, -1, 12, 13, 14, 15]
    The output dict will be: {10: 0, 11: 1, 12: 4, 13: 5, 14: 6, 15: 7}

    ri stands for reference index, qi stands for query index
    '''
    numba_dict = {}
    for qi in range(qseq_ref_pos_arr.size):
        ri = qseq_ref_pos_arr[qi]
        if ri >= 0:
            numba_dict[ri] = qi
    return numba_dict



def get_interval_seq(read, interval_start, interval_end, read_ref_pos_dict = {}):
    """
    Extract the sequence of a read within a given genomic interval.

    Parameters:
    - read: pysam.AlignedSegment
        The read to extract sequence from.
    - start, end: int
        The start and end positions of the interval in genomic coordinates.
        both should be inclusive
    - read_ref_pos_dict: dict
        Dictionary mapping read positions to reference positions. read_id -> (ref_positions, qseq_ref_positions)
        where detailed explanation of the qseq_ref_positions and ref_positions can be found in the docstring of the prepare_ref_query_idx_map function

    Returns:
    - str: Extracted sequence.
    - dict: Updated read_ref_pos_dict.
    - numpy.ndarray: Array of reference genomic positions corresponding to each base in the extracted sequence. -1 means the base is not aligned to the reference.

    The detailed explanation of the returned qidx_ridx_arr can be found in the docstring of the prepare_ref_query_idx_map function
    """
    # Both interval_start and interval_end are inclusive
    # Input interval_start and interval_end are both 0-indexed
    # get_reference_positions() also returns 0-indexed positions (a list)
    read_id = get_read_id(read)

    if read_id in read_ref_pos_dict:
        ref_positions, qseq_ref_positions = read_ref_pos_dict[read_id]
    else:
        # Detailed explanation of the qseq_ref_positions can be found in the docstring of the prepare_ref_query_idx_map function
        qseq_ref_positions = np.array([ idx if idx is not None else -1 for idx in read.get_reference_positions(full_length=True) ], dtype=np.int32)
        ref_positions = prepare_ref_query_idx_map(qseq_ref_positions)
        read_ref_pos_dict[read_id] = (ref_positions, qseq_ref_positions)

    size = interval_end - interval_start + 1

    preseq = []
    interval_start_qidx = ref_positions.get(interval_start, None)
    if interval_start_qidx is not None:
        # logger.debug(f"Found the interval start {interval_start} is not in the reference positions: \n{ref_positions} of read {read.query_name}. Might locate in the middle of a deletion event.")
        while interval_start not in ref_positions and interval_start <= interval_end:
            interval_start += 1
        if interval_start > interval_end:
            return_seq = preseq
            logger.warning(f"The whole specified interval (size {size}) is in a deletion event for read {read_id}. Now the returned seq is: {return_seq}\n")
            return return_seq, read_ref_pos_dict, np.array([], dtype=np.int32)
        interval_start_qidx = ref_positions[interval_start]

    interval_end_qidx = ref_positions.get(interval_end, None)
    if not interval_end_qidx:
        # logger.debug(f"Found the interval end {interval_end} is not in the reference positions: \n{ref_positions} of read {read.query_name}. Might locate in the middle of a deletion event.")
        while interval_end not in ref_positions:
            # postseq = postseq + ["D"]
            interval_end -= 1
        interval_end_qidx = ref_positions[interval_end] + 1
    else:
        interval_end_qidx += 1

    interval_read_seq = read.query_sequence[interval_start_qidx:interval_end_qidx]
    # Detailed explanation of the returned qidx_ridx_arr can be found in the docstring of the prepare_ref_query_idx_map function
    qidx_ridx_arr = qseq_ref_positions[interval_start_qidx:interval_end_qidx]

    return interval_read_seq, read_ref_pos_dict, qidx_ridx_arr



def get_interval_seq_qual(read, interval_start, interval_end, read_ref_pos_dict = {}):
    '''
    Extract the sequence and base quality sequence of a read within a given genomic interval.

    Parameters:
    - read: pysam.AlignedSegment
        The read to extract sequence and base quality sequence from.
    - start, end: int
        The start and end positions of the interval in genomic coordinates.
        both should be inclusive and 0-indexed
    - read_ref_pos_dict: dict
        Dictionary mapping read positions to reference positions. read_id -> (ref_positions, qseq_ref_positions)
        where detailed explanation of the qseq_ref_positions and ref_positions can be found in the docstring of the prepare_ref_query_idx_map function
    '''
    # Both interval_start and interval_end are inclusive
    # Input interval_start and interval_end are both 0-indexed
    # get_reference_positions() also returns 0-indexed positions (a list)
    read_id = get_read_id(read)
    if read_id in read_ref_pos_dict:
        ref_positions, qseq_ref_positions = read_ref_pos_dict[read_id]
    else:
        qseq_ref_positions = np.array([ idx if idx is not None else -1 for idx in read.get_reference_positions(full_length=True) ], dtype=np.int32)
        ref_positions = prepare_ref_query_idx_map(qseq_ref_positions)
        read_ref_pos_dict[read_id] = (ref_positions, qseq_ref_positions)

    size = interval_end - interval_start + 1
    preseq = []
    qual_preseq = []

    interval_start_qidx = ref_positions.get(interval_start, None)
    if not interval_start_qidx:
        while interval_start not in ref_positions and interval_start <= interval_end:
            preseq = preseq + ["D"]
            qual_preseq = qual_preseq + [99]
            interval_start += 1
        if interval_start > interval_end:
            return preseq, qual_preseq, read_ref_pos_dict, np.array([], dtype=np.int32)
        interval_start_qidx = ref_positions[interval_start]

    interval_end_qidx = ref_positions.get(interval_end, None)
    if not interval_end_qidx:
        while interval_end not in ref_positions:
            interval_end -= 1
        interval_end_qidx = ref_positions[interval_end] + 1
    else:
        interval_end_qidx += 1

    interval_read_seq = read.query_sequence[interval_start_qidx:interval_end_qidx]
    interval_qual_seq = read.query_qualities[interval_start_qidx:interval_end_qidx]
    qidx_ridx_arr = qseq_ref_positions[interval_start_qidx:interval_end_qidx]

    return interval_read_seq, interval_qual_seq, read_ref_pos_dict, qidx_ridx_arr


def seq_err_det_stacked_bases(target_read,
                              position0,
                              nested_ad_dict,
                              read_ref_pos_dict = {},
                              logger = logger):
    '''
    This function is used to stat the stacked bases at position0 (meaning 0-indexed positions) to:
    1. Identify the base quality of the base in the target read
    2. Calculate the allele depth of the allele in target_read at position0 to deduce whether this lowqual base is sequencing error or not

    Detailed structure of the nested_ad_dict can be found in the docstring of the function
    As to read_ref_pos_dict, it is read_id -> (ref_positions, qseq_ref_positions)
    where detailed explanation of the qseq_ref_positions and ref_positions can be found in the docstring of the prepare_ref_query_idx_map function

    Returns:
    - bool: True if the mismatches can be explained by sequencing artifacts, False otherwise
    - dict: Updated read_ref_pos_dict
    '''

    # Get the base quality of the base at position0 in the target read, as well as the base at position0 (0-indexed) in the query read
    target_base, base_qual, read_ref_pos_dict, qr_idx_arr = get_interval_seq_qual(target_read, position0, position0, read_ref_pos_dict)
    base_qual = base_qual[0]
    target_base = target_base[0]
    if base_qual >= 10:
        return False, read_ref_pos_dict

    chrom = target_read.reference_name

    ad = nested_ad_dict.get(target_read.reference_name, {}).get(position0, {}).get(target_base, 0)
    dp = nested_ad_dict.get(target_read.reference_name, {}).get(position0, {}).get("DP", 0)

    if int(dp) == 0:
        logger.warning(f"The depth at position {position0} is 0. The AD is {ad}. The target base is {target_base}. The base quality is {base_qual}. The read is {target_read.query_name}.")
        return False, read_ref_pos_dict

    af = int(ad) / int(dp)

    if (af <= 0.02 or (int(ad) == 1 and int(dp) >= 10)) and base_qual < 10:
        return True, read_ref_pos_dict
    else:
        return False, read_ref_pos_dict



def tolerate_mismatches_two_seq(read1, read2,
                                abs_diff_indices,
                                nested_ad_dict,
                                read_ref_pos_dict,
                                logger = logger):
    '''
    This function is used to compare two sequences and identify the mismatch positions
    And to decide whether we are safe to determine the mismatches are originated from sequencing error or not.

    Use seq_err_det_stacked_bases() to determine if the mismatches can be explained by sequencing artifacts, details can be found in the docstring of the function
    '''

    tolerate_mismatches = []
    for diff_ind in abs_diff_indices:
        seq_err1, read_ref_pos_dict = seq_err_det_stacked_bases(read1,
                                                                diff_ind,
                                                                nested_ad_dict,
                                                                read_ref_pos_dict,
                                                                logger = logger)

        seq_err2, read_ref_pos_dict = seq_err_det_stacked_bases(read2,
                                                                diff_ind,
                                                                nested_ad_dict,
                                                                read_ref_pos_dict,
                                                                logger = logger)
        tolerate_mismatches.append(seq_err1 or seq_err2)

    if all(tolerate_mismatches):
        return True, read_ref_pos_dict, len(tolerate_mismatches)
    else:
        return False, read_ref_pos_dict, 0



def get_overlap_intervals(read_pair1, read_pair2):
    overlap_intervals = {}

    for r1 in read_pair1:
        for r2 in read_pair2:
            r1_start = r1.reference_start
            r2_end = r2.reference_end

            if r2_end - r1_start <= 0 or r2_end - r1_start >= 300:
                continue

            r1_end = r1.reference_end
            r2_start = r2.reference_start

            if r1_end - r2_start <= 0:
                continue

            overlap_start = max(r1_start, r2_start)
            overlap_end = min(r1_end, r2_end)

            overlap_intervals[(overlap_start, overlap_end)] = (r1, r2)

    return overlap_intervals




def determine_same_haplotype(read, other_read,
                             overlap_start, overlap_end,
                             read_hap_vectors = {},
                             nested_ad_dict = {},
                             read_ref_pos_dict = {},
                             mean_read_length = 148,
                             logger = logger):
    '''
    Determine if two reads belong to the same haplotype based on their overlapping region.

    This function compares two reads in their overlapping region to decide if they likely
    come from the same haplotype. It uses haplotype vectors, sequence comparison, and
    various heuristics to make this determination.

    The function uses a hierachical approach to determine if two reads belong to the same haplotype:
    1. Compare two reads' haplotype vectors within the overlapping region
        a. If two reads' haplotype vectors within the overlapping region are differed by more than 3 bases, we direct reject the possiblity of being the same haplotype and return False
        b. If two reads' haplotype vectors within the overlapping region are differed by at least one indel, we direct reject the possiblity of being the same haplotype and return False
    2. Extract the raw sequence of two reads within the overlapping region
        a. If two reads' raw sequence are identical, and the overlap span is greater than mean_read_length - 50, we consider the two reads belong to the same haplotype and return True
        b. If two reads' raw sequence are identical, and the overlap span is less than mean_read_length - 50, we consider not enough evidence to make a call and return np.nan
        c. If two reads' raw sequence are differed by at least one indel, we directly reject the possiblity of being the same haplotype and return False
        d. If two reads' raw sequence are differed by mismatches, we try to see if the mismatches can be explained by sequencing artifacts
           i. If the mismatches can be explained by sequencing artifacts and the overlap span is greater than mean_read_length - 50, we consider the two reads belong to the same haplotype and return True
           ii. If the mismatches can be explained by sequencing artifacts and the overlap span is less than mean_read_length - 50, we consider not enough evidence to make a call and return np.nan
           iii. If the mismatches cannot be explained by sequencing artifacts, we consider the two reads belong to different haplotypes and return False

    Parameters:
    - read, other_read: pysam.AlignedSegment
        The two reads to compare.
    - overlap_start, overlap_end: int
        The start and end positions of the overlapping region.
    - read_hap_vectors: dict
        Cache of haplotype vectors for reads. Structure: {read_id: haplotype_vector}
    - nested_ad_dict: dict
        Nested dictionary for allele depth information. Structure: {chrom: {pos: {base: depth}}}
    - read_ref_pos_dict: dict
        Dictionary mapping read positions to reference positions. Structure: {read_id: (ref_positions, qseq_ref_positions)}
        Detailed explanation of the qseq_ref_positions and ref_positions can be found in the docstring of the prepare_ref_query_idx_map function
    - mean_read_length: int
        Average read length, used for various calculations.
    - logger: logging.Logger
        Logger for output messages.

    Returns:
    - bool or float: True if same haplotype, False if different, np.nan if undetermined
    - dict: Updated read_ref_pos_dict
    - dict: Updated read_hap_vectors
    - float or None: Overlap span (weight) if determined to be same haplotype, else None
    '''

    read_id = get_read_id(read)
    start = read.reference_start
    end = read.reference_end

    # Below we extract the haplotype vector for both reads, which are an array of int values indicating the alignment status of each base within the reference genome span aligned by the read
    # Detailed explanation can be found in the docstring of the function get_hapvector_from_cigar
    if read_id in read_hap_vectors:
        read_hap_vector = read_hap_vectors[read_id]
    else:
        read_hap_vector = get_hapvector_from_cigar(read.cigartuples, read.query_sequence, logger = logger)
        read_hap_vectors[read_id] = read_hap_vector

    # We need a unique identifier for every read, which is a string composed by query name and the alignment flag.
    other_read_id = get_read_id(other_read)
    other_start = other_read.reference_start
    other_end = other_read.reference_end

    if other_read_id in read_hap_vectors:
        other_read_hap_vector = read_hap_vectors[other_read_id]
    else:
        other_read_hap_vector = get_hapvector_from_cigar(other_read.cigartuples, other_read.query_sequence, logger = logger)
        read_hap_vectors[other_read_id] = other_read_hap_vector

    # Extract the haplotype vector within the overlapping region
    interval_hap_vector = read_hap_vector[overlap_start - start:overlap_end - start]
    interval_other_hap_vector = other_read_hap_vector[overlap_start - other_start:overlap_end - other_start]

    # Calculate the overlap span
    overlap_span = overlap_end - overlap_start

    # Find the indices where two haplotype vectors differ
    diff_indices = numba_diff_indices(interval_hap_vector, interval_other_hap_vector)

    # First test whether there are too many mismatches between two reads that we wont tolerate
    if len(diff_indices) >= 3:
        # Cannot tolerate such mismatches
        # Short-circuit evaluation, if this condition is true, the following code will not be evaluated
        return False, read_ref_pos_dict, read_hap_vectors, None

    # Extract the positions in the read's haplotype vector where two haplotype vectors differ
    diff_pos_read = interval_hap_vector[diff_indices]
    # Extract the positions in the other read's haplotype vector where two haplotype vectors differ
    diff_pos_oread = interval_other_hap_vector[diff_indices]

    # Check if there are indels in the positions of the read's haplotype vector where two haplotype vectors differ
    read_indel_overlap = get_indel_bools(diff_pos_read)
    # Check if there are indels in the positions of the other read's haplotype vector where two haplotype vectors differ
    oread_indel_overlap = get_indel_bools(diff_pos_oread)

    # If there are indels in the positions of the read's haplotype vector where two haplotype vectors differ or the other read's haplotype vector where two haplotype vectors differ, then we cannot make a call on whether they belong to the same haplotype
    # Why we use hap vectors to compare ahead of using raw sequence?
    # The reason is that the haplotype vector comparison is much more efficient than string comparison
    # This is a short-circuit evaluation, if the first condition is true, the following condition will not be evaluated
    if numba_sum(read_indel_overlap) > 0 or numba_sum(oread_indel_overlap) > 0:
        return False, read_ref_pos_dict, read_hap_vectors, None

    # Now use overlap_start and overlap_end to extract the read sequence (a string series) of two reads within the overlap region
    # If two read sequences are not the same length, then it indicates they have difference in indels. A clear evidence to reject the possiblity of being the same haplotype
    read_seq, read_ref_pos_dict, qr_idx_arr = get_interval_seq(read, overlap_start, overlap_end - 1, read_ref_pos_dict)
    other_seq, read_ref_pos_dict, other_qr_idx_arr = get_interval_seq(other_read, overlap_start, overlap_end - 1, read_ref_pos_dict)

    # We need to ignore N during the comparison for read_seq and other_seq
    if "N" in read_seq or "N" in other_seq:
        base_dict = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
        read_seq_arr = np.array(list(map(base_dict.get, read_seq)), dtype=np.int8)
        other_seq_arr = np.array(list(map(base_dict.get, other_seq)), dtype=np.int8)
        total_match = compare_sequences(read_seq_arr, other_seq_arr, np.int8(base_dict["N"]))
    else:
        total_match = read_seq == other_seq
        read_seq_arr = None

    '''
    # code block for debugging
    # Allow us to review the extracted read sequence at specified overlap region
    print(read_seq, file = sys.stderr)
    print(other_seq, file = sys.stderr)
    print("\n", file = sys.stderr)
    logger.info(f"The read sequence is {read_seq} for {read_id} \nand the other read sequence is {other_seq} for other read_id {other_read_id} at region {overlap_span}")
    logger.info(f"The hap vector is {interval_hap_vector} for {read_id} \nand the other hap vector is {interval_other_hap_vector} for other read_id {other_read_id} at region {overlap_span}")
    assert interval_hap_vector.size == interval_other_hap_vector.size, f"The interval hap vector size {interval_hap_vector} is not the same as the other interval hap vector size {interval_other_hap_vector}"
    '''

    '''
    The following code block is used to calculate the edge weight between two reads and it is based on the overlapping span and the number of shared variants.
    1. overlap span is calculated based on the number of identical bases between two reads
    2. The number of shared variants is calculated by the function count_var
    3. The number of indels is calculated by the function count_continuous_indel_blocks

    We created the score array, this score array is used to assign weight to edges depending on the number of SNVs shared by two read pairs
    It will be later used by the function determine_same_haplotype to assign weight to edges
    '''

    identical_part = interval_hap_vector[interval_hap_vector == interval_other_hap_vector]
    overlap_span = identical_part.size
    var_count = count_var(identical_part)
    indel_num = count_continuous_indel_blocks(identical_part)
    # assert var_size is not None, f"The size of the variant should not be None, but the actual size is {var_size}, the input array is {interval_hap_vector}"

    snvcount_score_arr = np.array([mean_read_length * 1 - 50 + mean_read_length * i for i in range(50)])
    edge_weight = overlap_span + numba_sum(snvcount_score_arr[:var_count])
    edge_weight = overlap_span + mean_read_length * 3 * indel_num

    '''
    # Below is a code block for debugging
    # Allow us to review the read sequence of your interest query name at specified overlap region

    vis_qnames = ["HISEQ1:21:H9V1VADXX:2:1112:21127:38421:PC0",
                  "HISEQ1:26:HA2RRADXX:1:1113:11601:32503:PC0"]

    if read.query_name in vis_qnames and other_read.query_name in vis_qnames:
        logger.info(f"The read sequence is {interval_hap_vector.tolist()} for {read_id} \nand the other read sequence is {interval_other_hap_vector.tolist()} for other read_id {other_read_id} at region {overlap_start}, {overlap_end}. \nThe different indices are {diff_indices.tolist()} and the identical part is {identical_part.tolist()}")
        logger.info(f"The raw read sequence is {read_hap_vector.tolist()} for {read_id} with cigar {read.cigartuples} and query_sequence {read.query_sequence} \nand the other read sequence is {other_read_hap_vector.tolist()} for other read_id {other_read_id} with cigar {other_read.cigartuples} and query sequence {other_read.query_sequence}")
    '''

    if total_match:
        # Basically edge weight is overlap_span + indel_added_score + snv_added_score
        # If edge weight is greater than mean_read_length - 50 (minimum score), we consider the two reads belong to the same haplotype and return the edge weight.
        # This is just a heuristic cutoff and can be further tuned according to different sequencing profiles
        if edge_weight >= mean_read_length - 50:
            return True, read_ref_pos_dict, read_hap_vectors, edge_weight
        else:
            # If edge weight is smaller than mean_read_length - 50, we consider not enough evidence to make a call and return None edge weight
            return np.nan, read_ref_pos_dict, read_hap_vectors, None
    else:
        # Check if two reads within this overlapping region are of the same length
        # If not, we cannot make a conclusion on whether they belong to the same haplotype because they have a difference in indel
        if len(read_seq) != len(other_seq) or interval_hap_vector.size != interval_other_hap_vector.size:
            return False, read_ref_pos_dict, read_hap_vectors, None

        # Here, we know that though two reads are not identical, they do have the same length
        # Then we need to know whether the two mismatches are in-trans variants or in-cis variants, if its in-trans we cannot tolerate the two mismatches
        if read_seq_arr is None:
            base_dict = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
            read_seq_arr = np.array(list(map(base_dict.get, read_seq)), dtype=np.int8)
            other_seq_arr = np.array(list(map(base_dict.get, other_seq)), dtype=np.int8)

        # Use numba_diff_indices to find the location indices where two sequences differ
        q_diff_indices = numba_diff_indices(read_seq_arr, other_seq_arr)
        # Map the q_diff_indices to the reference positions
        r_diff_indices = qr_idx_arr[q_diff_indices]
        # If -1 in r_diff_indices, it means the two reads are not aligned properly and we cannot make a conclusion on whether they belong to the same haplotype
        if -1 in r_diff_indices:
            return False, read_ref_pos_dict, read_hap_vectors, None
        abs_diff_inds = r_diff_indices + overlap_start
        # See if the mismatches can be explained by sequencing artifacts
        tolerate, read_ref_pos_dict, tolerated_count = tolerate_mismatches_two_seq( read,
                                                                                    other_read,
                                                                                    abs_diff_inds,
                                                                                    nested_ad_dict,
                                                                                    read_ref_pos_dict,
                                                                                    logger = logger )

        if tolerate:
            edge_weight = edge_weight - tolerated_count * 20
            if edge_weight >= mean_read_length - 50:
                return True, read_ref_pos_dict, read_hap_vectors, edge_weight
            else:
                return np.nan, read_ref_pos_dict, read_hap_vectors, None
        else:
            return False, read_ref_pos_dict, read_hap_vectors, None



def pretty_print_matrix(matrix, precision=3):
    """
    Pretty prints a 2D NumPy array.

    Args:
        matrix (numpy.ndarray): The 2D array to be printed.
    """
    with np.printoptions(precision=precision, suppress=True):
        return "\n".join(["\t".join(map("{:.3f}".format, row)) for row in matrix])







@numba.njit(types.Tuple((types.float32[:], types.int32[:], types.int32[:], types.int32))(types.float32[:], types.int32[:], types.int32[:], types.boolean[:]), fastmath=True)
def sparse_mask_selection(data, indices, indptr, mask):
    new_size = numba_sum(mask)
    new_indptr = np.empty(new_size + 1, dtype = np.int32)
    new_data = np.empty(new_size, dtype = np.float32)
    new_indices = np.empty(new_size, dtype = np.int32)
    row_idx = 0
    new_indptr[row_idx] = 0
    for i in range(len(mask)):
        if mask[i]:
            start, end = indptr[i], indptr[i+1]
            for j in range(start, end):
                if mask[indices[j]]:
                    new_data[row_idx] = data[j]
                    new_indices[row_idx] = numba_sum(mask[:indices[j]])
            new_indptr[row_idx + 1] = len(new_data)
            row_idx += 1
    return new_data, new_indices, new_indptr, new_size





def graph_vertex_iter(vertex_indices, graph):
    for vid in vertex_indices:
        yield graph.vertex(vid)





def find_uncovered_regions(existing_intervals, new_interval):
    # Find overlapping intervals
    overlapping_intervals = existing_intervals.overlap(new_interval[0], new_interval[1])

    # Sort the overlapping intervals by start time
    sorted_intervals = sorted(overlapping_intervals, key=lambda x: x.begin)

    # Determine the uncovered regions
    uncovered_regions = []
    current_start = new_interval[0]

    if len(overlapping_intervals) == 0:
        return [new_interval]

    for interval in sorted_intervals:
        if interval.begin > current_start:
            uncovered_regions.append((current_start, interval.begin))
        current_start = max(current_start, interval.end)

    if current_start < new_interval[1]:
        uncovered_regions.append((current_start, new_interval[1]))

    return uncovered_regions



def stat_ad_dict(bam_file, logger = logger):
    bam_ad_file = f"{bam_file}.ad"
    cmd = f"""bcftools mpileup -Ou --no-reference -a FORMAT/AD --indels-2.0 -q 10 -Q 15 {bam_file} | \
              bcftools query -f '%CHROM\\t%POS\\t%ALT\\t[%AD]\\n' - > {bam_ad_file}"""
    executeCmd(cmd, logger = logger)
    ad_table = pd.read_table(bam_ad_file, header = None, sep = "\t", names = ["chrom", "pos", "alt", "ad"], dtype = {"chrom": str, "pos": int, "alt": str, "ad": str}, na_values=["", "<*>"])
    ad_expanded = ad_table["ad"].str.split(",", expand=True).replace({None: np.nan, "": np.nan, "0": np.nan}).dropna(axis = 1, how="all").astype(float)
    alt_expanded = ad_table["alt"].str.rstrip(",<*>").str.split(",", expand=True).replace({None: np.nan, "": np.nan}).dropna(axis = 1, how="all")

    if alt_expanded.shape[1] <= 1 or ad_expanded.shape[1] <= 1:
        logger.warning(f"No ALT allele found in this BAM file. Skip this entire script")
        logger.warning(f"Look at the two dataframes now: \n{alt_expanded.to_string(index=False)}\n{ad_expanded.to_string(index=False)}\n")
        return None

    logger.info("\n{}\n{}\n".format(ad_expanded.loc[~ad_expanded.iloc[:, -1].isna(), :][:10].to_string(index=False),
                                    alt_expanded.loc[~alt_expanded.iloc[:, -1].isna(), :][:10].to_string(index=False)))

    # After getting the table, we seek to organize the AD info by chromsome, pos, and alt into a nested dictionary
    # The dict structure is like:
    # nested_ad_dict[chrom][pos][alt] = ad
    # nested_ad_dict[chrom][pos]["DP"] = total_dp
    # Initialize the nested dictionary
    nested_ad_dict = {chrom: {} for chrom in ad_table["chrom"].unique()}
    column_width = alt_expanded.shape[1]

    # Iterate over the rows of ad_table to build the nested_ad_dict
    for i in range(len(ad_table)):
        chrom = ad_table.iloc[i, 0]
        outer_key = ad_table.iloc[i, 1]
        inner_key = alt_expanded.iloc[i, 0]
        value = ad_expanded.iloc[i, 0]

        if pd.isna(value):
            continue

        total_dp = value

        # Initialize the inner dictionary if the outer key is not present
        if outer_key not in nested_ad_dict[chrom]:
            nested_ad_dict[chrom][outer_key] = {}

        # Add the first pair of inner key-value
        nested_ad_dict[chrom][outer_key][inner_key] = value

        # Check for the second pair of inner key-value
        for c in range(1, column_width):
            if not pd.isna(ad_expanded.iloc[i, c]) and not pd.isna(alt_expanded.iloc[i, c]):
                nested_ad_dict[chrom][outer_key][alt_expanded.iloc[i, c]] = ad_expanded.iloc[i, c]
                total_dp += ad_expanded.iloc[i, c]

        nested_ad_dict[chrom][outer_key]["DP"] = total_dp

    return nested_ad_dict



def build_phasing_graph(bam_file,
                        ncls_dict,
                        ncls_read_dict,
                        ncls_qname_dict,
                        ncls_qname_idx_dict,
                        mean_read_length,
                        edge_weight_cutoff = 0.201,
                        logger = logger):
    '''
    Construct a phasing graph from BAM data for efficient haplotype identification.

    This function builds a graph where each vertex represents a read pair and each edge
    represents the confidence that two read pairs originated from the same haplotype.
    The graph is used for subsequent haplotype identification through clique finding.

    Parameters:
    -----------
    bam_file : str
        Path to the input BAM file.

    ncls_dict : dict
        Dictionary of NCLS objects, keyed by chromsome names.
        The ncls_dict looks like:
                                        { chromsome1: ncls_chr1,
                                          chromsome2: ncls_chr2,
                                          ...
                                          chromsomeN: ncls_chrN  }

    ncls_read_dict : dict
        Dictionary of read objects, keyed by query name indices (which are also the interval indices stored in the NCLS object)
        The ncls_read_dict looks like:
                                        { qname_idx1: [read1, read2],
                                          qname_idx2: [read1, read2]  }

    ncls_qname_dict : dict
        Dictionary mapping query name indices to query names. (qname_indices -> qnames)

    ncls_qname_idx_dict : dict
        Dictionary mapping query names to query name indices. (qnames -> qname_indices)

    mean_read_length : float
        Average read length in the BAM file.
    edge_weight_cutoff : float, optional
        Minimum edge weight threshold for including an edge in the graph (default is 0.201, a heuristic cutoff).
    logger : logging.Logger, optional
        Logger object for output messages.

    Returns:
    --------
    tuple
        A tuple containing:
        - phased_graph : graph_tool.Graph
            The constructed phasing graph.
        - weight_matrix : numpy.ndarray
            Matrix of edge weights between read pairs. (adjacency matrix)
        - qname_to_node : dict
            Dictionary mapping query names to graph vertices. (qnames -> graph_vertices_indices)
        - total_readhap_vector : dict
            Dictionary of read haplotype vectors. (read_ids -> haplotype_vectors), please refer read_ids to function get_read_id

    Notes:
    ------
    - The function uses the graph-tool library for efficient graph operations.
    - Edge weights are calculated based on two factors:
        1. The total overlapping span between two read pairs
        2. The number of variants (SNVs and small indels) shared by two read pairs
    - The resulting graph is used for subsequent haplotype identification through clique finding (Bron-Kerbosch algorithm).

    See Also:
    ---------
    migrate_bam_to_ncls : For creating the NCLS data structures.
    clique_generator_per_component : For identifying haplotypes in the constructed graph.
    '''

    logger.info(f"There are totally {len(ncls_read_dict)} pair of reads, mean read length is {mean_read_length}. with adequate mapping or base quality which can be used to build the graph")
    # Create an empty graph
    g = gt.Graph(directed = False)
    g.set_fast_edge_removal(fast = True)

    # Create a property map to store the query names for each node
    qname_prop = g.new_vertex_property("string")
    weight = g.new_edge_property("float")

    # When we try to identify a mismatch between two reads are due to a true variant or due to sequencing errors,
    # we need to know the allele depth at the mismatch site. This allows us to estimate the sequencing artifact likelihood by observing the allele's depth fraction at this site.
    # We would like an efficient way to extract AD info for all covered sites across this BAM. So we chose bcftools mpileup + bcftools query.
    # The output file is stored in the same directory as the input BAM file
    nested_ad_dict = stat_ad_dict(bam_file, logger = logger)
    if nested_ad_dict is None:
        return None, None, None, None

    '''
    Now we start to build an undirected graph of read pairs,
    The graph is constructed as follows:
    - Each read pair is represented by a vertex
    - Iterate through all the read pairs in the ncls_read_dict
       - For each read pair, extract its overlapping read pairs
       - Iterate through all the overlapping read pairs, skip the read pairs that have been inspected (undirected)

    '''

    # Create a dictionary to store tuples of vertex indices (v1, v2) as keys, and the value is a boolean indicating whether the pair of read pairs has been inspected
    # This is a temporary data hub to avoid repeated edge inspection for the same pair of reads
    qname_check_dict = {}

    # Initialize the weight matrix (adjacency matrix, which is diagonal symmetric) with np.eye
    total_qname_num = len(ncls_read_dict)
    weight_matrix = np.eye(total_qname_num, dtype=np.float32)

    # Create a dictionary to store the haplotype vectors corresponding to each read (read_ids -> haplotype_vectors), please refer read_ids to function get_read_id
    # This is a temporary data hub to avoid repeated hapvector extraction from the same read
    # The variable is primarily aggregated in the iterative calling of function determine_same_haplotype
    read_hap_vectors = {}

    # Create a dictionary to store the start and end positions for each read pair
    # This is a temporary data hub to avoid repeated read pair start and end position extraction from the same read
    # The variable is primarily aggregated in the iterative calling of function determine_same_haplotype
    read_ref_pos_dict = {}

    # Create a dictionary to map query names to their corresponding nodes (qnames -> graph_vertices_indices)
    # This is a temporary data hub to record whether a read pair has been added to the graph as a vertex(node)
    qname_to_node = {}

    # Create an IntervalTree to record the inspected overlapping intervals
    inspected_overlaps = IntervalTree()

    for qname_idx, paired_reads in ncls_read_dict.items():
        qname = ncls_qname_dict[qname_idx]
        assert len(dict.fromkeys(r.query_name for r in paired_reads)) == 1, f"The query names of the reads in the paired_reads are not the same: {paired_reads}"

        # Check if the query name already exists in the graph
        if qname not in qname_to_node:
            # Add a new node to the graph in case that the read pair has not been added to the graph
            qv = g.add_vertex()
            qname_prop[qv] = qname
            qname_to_node[qname] = int(qv)
        else:
            # Directly retrieve the existing vertex (graph-tool vertex object) from the graph
            qv = g.vertex(qname_to_node[qname])

        # Extract the overlapping span of the current iterating read pair
        # And get a generator of qname indices of the read pairs overlapping with the current read pair
        # Note that the span also include the inner gap between a pair of reads
        chrom = paired_reads[0].reference_name
        start = min(r.reference_start for r in paired_reads)
        end = max(r.reference_end for r in paired_reads)
        qidx_iter = lazy_get_overlapping_qname_idx(ncls_dict, ncls_read_dict, ncls_qname_dict, chrom, start, end)

        # Iterate through the overlapping read pairs
        for qidx in qidx_iter:
            other_reads = ncls_read_dict[qidx]
            # Get the query name
            other_qname = other_reads[0].query_name
            # When extracting overlapping reads from NCLS, the result qnames might contain the iterating qname itself
            if qname == other_qname:
                continue

            '''
            # Below is a code block for debugging, this allows you to check whether two read pairs overlap with each other during the inspection of them
            assert len(dict.fromkeys(r.query_name for r in other_reads)) == 1, f"The query names of the reads in the other_reads are not the same: {other_reads}"
            vis_qnames = ["HISEQ1:21:H9V1VADXX:2:1112:21127:38421:PC0",
                            "HISEQ1:26:HA2RRADXX:1:1113:11601:32503:PC0"]
            if qname in vis_qnames and other_qname in vis_qnames:
                logger.info(f"Found the qname {qname} and other_qname {other_qname} might overlap with each other.")
            '''

            # Check if the query name already exists in the graph
            if not other_qname in qname_to_node:
                # Add a new node to the graph
                oqv = g.add_vertex()
                qname_prop[oqv] = other_qname
                qname_to_node[other_qname] = int(oqv)
                # logger.debug(f"Added a new node {v} for qname {qname} to the graph. The current vertices are {g.get_vertices()}")
            else:
                oqv = g.vertex(qname_to_node[other_qname])

            # Check if the pair of read pairs has been inspected, if yes, skip the current iteration
            inspect_res = check_edge(int(qv), int(oqv), qname_check_dict)
            if inspect_res:
                continue

            # Proceeding to this step, we mark the pair of read pairs as inspected
            qname_check_dict[(int(qv), int(oqv))] = True

            # Inspect the overlap between the two pair of reads
            # There might be multiple overlapping intervals between two pair of reads
            overlap_intervals = get_overlap_intervals(paired_reads, other_reads)
            # logger.info(f"Found these overlap intervals for the read pair {qname} and {other_qname}: {overlap_intervals}")

            # If there is no overlap between two pair of reads, skip the current iteration
            # This might be some edge case when a read in pair A completely enclosed by the inner gap of pair B, while the other read in pair A exceeds the coverage span of pair B
            if len(overlap_intervals) == 0:
                continue

            '''
            There are four comparison results (C 2 2) to be recorded between two pair of reads
            - read1 in pair A and read1 in pair B comparison result
            - read1 in pair A and read2 in pair B comparison result
            - read2 in pair A and read1 in pair B comparison result
            - read2 in pair A and read2 in pair B comparison result

            There are also three possible results for the comparison of two reads
             1: meaning we found two reads from distinct pairs overlap with each other with enough size, (or share variants), so they should be in the same haplotype
             0: meaning we found two reads from distinct pairs do not overlap with each other or the overlapping span is too small (e.g < 100bp). So we do not have enough evidence to determine if they are in the same haplotype or not
            -1: meaning we found two reads from distinct pairs overlap with each other and we found clear evidence that they are not in the same haplotype (e.g. different variants)
            '''

            # Initialize the same_hap_bools with 4 zeros
            same_hap_bools = np.zeros(4, dtype=np.int32)

            n = 0
            pair_weight = None

            # Remember the overlap intervals we extracted for two pair of reads
            # Now we iterate through all the overlap intervals and use an IntervalTree to record the inspected intervals
            inspected_overlaps = IntervalTree()
            for (overlap_start, overlap_end), (read1, read2) in overlap_intervals.items():
                '''
                First we identify whehter the iterating overlap interval has some part of it that has been inspected
                If yes, we crop out the inspected part and only keep the uncovered overlapping span for downstream analysis
                Note that after cropping out the inspected part, the remaining span might not be continuous, leaving multiple uncovered intervals
                For example, the inspected intervals is [100, 200], [300, 400], [500, 600]
                If the current iterating overlap interval is [250, 550], then the remaining uncovered intervals are [250, 400] and [500, 550]
                The function IntervalTree.addi will return the following uncovered intervals: [250, 400], [500, 550]
                '''
                uncovered_overlaps = find_uncovered_regions(inspected_overlaps, (overlap_start, overlap_end))
                # logger.info(f"Found the uncovered regions {uncovered_overlaps} for the reads {read1.query_name} and {read2.query_name} in the overlap region ({overlap_start}-{overlap_end})")
                for uncovered_start, uncovered_end in uncovered_overlaps:
                    # Iterate over one uncovered interval at a time, and determine whether two read from distinct pairs overlapping at the current interval show evidence of being in the same haplotype or not
                    bool_res, read_ref_pos_dict, read_hap_vectors, read_weight = determine_same_haplotype(read1, read2,
                                                                                                          uncovered_start, uncovered_end,
                                                                                                          read_hap_vectors = read_hap_vectors,
                                                                                                          nested_ad_dict = nested_ad_dict,
                                                                                                          read_ref_pos_dict = read_ref_pos_dict,
                                                                                                          mean_read_length = mean_read_length,
                                                                                                          logger = logger)
                    # If the bool_res is NaN, we assign 0 to record the inspection result within this overlapping interval
                    # Otherwise, we assign 1 or -1 based on the bool_res
                    if np.isnan(bool_res):
                        same_hap_bools[n] = 0
                    elif bool_res:
                        same_hap_bools[n] = 1
                    else:
                        same_hap_bools[n] = -1

                    if read_weight is not None:
                        # Assign 1 (an aritificial very small number) to read_weight unless it a positive value
                        read_weight = read_weight if read_weight > 0 else 1
                        # Normalize the read_weight to adjust the scale of the weight, why denominate by mean_read_length * 10?
                        # mean_read_length * 10 is just a big enough number to make sure ur normal weight is at the scale between 0 and 1, u can understand it as a maximum normalization factor
                        norm_weight = read_weight/(mean_read_length * 10)
                        if pair_weight is None:
                            # If pair_weight is None, assign norm_weight to it
                            pair_weight = norm_weight
                        else:
                            # If pair_weight is not None, add norm_weight to it
                            pair_weight += norm_weight
                    else:
                        norm_weight = None

                # After inspecting the current uncovered interval, we add it to the inspected_overlaps
                # And adding 1 to n for next overlapping interval.
                inspected_overlaps.addi(overlap_start, overlap_end)
                n += 1

            # Now we iterated through all the overlapping intervals, we can trim the same_hap_bools to the actual number of overlapping intervals
            same_hap_bools = same_hap_bools[:n]

            # Because 1 in the adjacency matrix has a special meaning. The diagonal values are all 1 (representing self-self comparison)
            # So we add a very small number to the pair_weight to distinguish it from 1
            pair_weight = pair_weight if pair_weight != 1 else 1 + 1e-4

            # Below is just a debug code block
            # if qname in vis_qnames and other_qname in vis_qnames:
            #     logger.info(f"Found the qname {qname} (qv is {int(qv)}) and other_qname {other_qname} (oqv is {int(oqv)}) overlap with each other with pair_weight {pair_weight}. The same_hap_bools are {same_hap_bools}")

            # We need to consider all the overlapping intervals between two read pairs to decide the value we record inside the adjacency matrix
            if custom_all_numba(same_hap_bools):
                # logger.info(f"The same_hap_bools are {same_hap_bools}")
                if pair_weight > edge_weight_cutoff:
                    e = g.add_edge(qv, oqv)
                    weight[e] = pair_weight
                weight_matrix[int(qv), int(oqv)] = pair_weight
                weight_matrix[int(oqv), int(qv)] = pair_weight
            elif any_false_numba(same_hap_bools):
                # We found clear evidence that two read pairs are not in the same haplotype
                weight_matrix[int(qv), int(oqv)] = -1
                weight_matrix[int(oqv), int(qv)] = -1

    # Set the query name property for the graph
    g.vertex_properties["qname"] = qname_prop
    g.edge_properties["weight"] = weight

    logger.info(f"Now we finished building up the edges in the graph. There are currently {g.num_vertices()} vertices and {g.num_edges()} edges in the graph")
    return g, weight_matrix, qname_to_node, read_hap_vectors






def extract_var_pos(raw_vcf,
                    bam_region_bed,
                    logger = logger):
    cmd = f"bcftools view -R {bam_region_bed} -Ou {raw_vcf} | \
            bcftools sort -Ou - | \
            bcftools query -f '%CHROM\\t%POS\\t%END\\n' -"

    # ATTENTION! RETURNED coordinates are 1-indexed!
    var_pos_str = executeCmd(cmd, stdout_only = True, logger = logger)

    df = pd.read_table(StringIO(var_pos_str), header = None, sep = "\t", names = ["chrom", "start", "end"], dtype = {"chrom": str, "start": int, "end": int})
    logger.info(f"The variant positions are extracted from the VCF file. There are {df.shape[0]} variants in total. They looks like :\n{df[:5].to_string(index=False)}\n")
    gc.collect()

    return df



def record_haplotype_rank(haplotype_dict, mean_read_length = 150):
    starts = np.empty(len(haplotype_dict), dtype=np.int32)
    ends = np.empty(len(haplotype_dict), dtype=np.int32)
    hap_depths = np.empty(len(haplotype_dict), dtype=np.int32)
    hap_ids = np.empty(len(haplotype_dict), dtype=np.int32)
    var_counts = np.empty(len(haplotype_dict), dtype=np.int32)
    indel_counts = np.empty(len(haplotype_dict), dtype=np.int32)
    total_depth_col = np.empty(len(haplotype_dict), dtype=np.int32)

    i = 0
    total_depth = 0
    for hid in haplotype_dict:
        seq, reads, span, qnames = haplotype_dict[hid]
        starts[i] = span[0]
        ends[i] = span[1]
        depth = len(reads) * mean_read_length/len(seq)
        hap_depths[i] = depth
        hap_ids[i] = hid
        var_counts[i] = count_var(seq)
        indel_counts[i] = count_continuous_indel_blocks(seq)
        total_depth += depth
        i += 1

    total_depth_col.fill(total_depth)
    return np.column_stack((starts, ends, total_depth_col, hap_ids, hap_depths, var_counts, indel_counts))


def group_by_dict_optimized(vprop, vertices):
    grouped_keys = {}
    for v, rs in vertices.items():
        v_idx, q = v
        label = vprop[v_idx]
        if label not in grouped_keys:
            grouped_keys[label] = [[], [], []]
        grouped_keys[label][0].append(v_idx)
        grouped_keys[label][1].append(q)
        grouped_keys[label][2].append(rs)
    return grouped_keys




@numba.njit
def rank_unique_values(arr):
    # Step 1: Extract unique values and sort them
    unique_values = np.unique(arr)

    # Step 2: Create a mapping from each unique value to its rank
    value_to_rank = np.empty(unique_values.shape, dtype=np.int32)
    for i in range(unique_values.size):
        value_to_rank[i] = i + 1  # Ranks start from 1

    # Step 3: Apply the mapping to create a new array with the ranks
    ranks = np.empty(arr.shape, dtype=np.int32)
    for i in range(arr.size):
        for j in range(unique_values.size):
            if arr[i] == unique_values[j]:
                ranks[i] = value_to_rank[j]
                break
    return ranks



@numba.njit(types.float32[:](types.int32[:, :]), fastmath=True)
def calculate_coefficient(arr2d):
    # Span size is the weight for later mean value calculation
    rank = arr2d[:, 7]
    span = (arr2d[:, 1] - arr2d[:, 0])
    depth_frac = (1 - arr2d[:, 4]/arr2d[:, 2]).astype(types.float32)

    return rank * span * np.sqrt(depth_frac)




@numba.njit(types.int32[:,:](types.boolean[:]), fastmath=True)
def extract_true_stretches(bool_array):
    n = len(bool_array)

    # Pre-allocate maximum possible space (worst case: every element is True)
    stretches = np.zeros((n, 2), dtype=np.int32)
    count = 0

    if n == 0:
        return stretches[:count]

    in_stretch = False
    start_idx = 0

    for i in range(n):
        if bool_array[i]:
            if not in_stretch:
                in_stretch = True
                start_idx = i
        else:
            if in_stretch:
                stretches[count, 0] = start_idx
                stretches[count, 1] = i - 1
                count += 1
                in_stretch = False

    # Handle the case where the array ends with a True stretch
    if in_stretch:
        stretches[count, 0] = start_idx
        stretches[count, 1] = n - 1
        count += 1

    return stretches[:count]



@numba.njit(types.boolean(types.int32[:]), fastmath=True)
def judge_misalignment_by_extreme_vardensity(seq):
    five_vard = count_window_var_density(seq, padding_size = 42)
    six_vard = count_window_var_density(seq, padding_size = 65)
    read_vard = count_window_var_density(seq, padding_size = 74)
    # indel_count = count_continuous_indel_blocks(seq)
    if numba_sum(five_vard >= 5/85) > 0:
        select_bool = five_vard >= 5/85
        # pad the select_bool by 42 to both directions
        true_segments = extract_true_stretches(select_bool)
        max_indel_count = 0
        for i in range(true_segments.shape[0]):
            start = true_segments[i, 0] - 42
            end = true_segments[i, 1] + 42
            five_seq = seq[start:end]
            indel_count = count_continuous_indel_blocks(five_seq)
            max_indel_count = max(max_indel_count, indel_count)
        if max_indel_count > 1:
            return True
    elif numba_sum(six_vard >= 6/131) > 0:
        select_bool = six_vard >= 6/131
        true_segments = extract_true_stretches(select_bool)
        max_indel_count = 0
        for i in range(true_segments.shape[0]):
            start = true_segments[i, 0] - 65
            end = true_segments[i, 1] + 65
            five_seq = seq[start:end]
            indel_count = count_continuous_indel_blocks(five_seq)
            max_indel_count = max(max_indel_count, indel_count)
        if max_indel_count > 1:
            return True
    elif numba_sum(read_vard > 11/148) > 0:
        return True
    return False



@numba.njit
def find_indices(hap_ids, included_hapids):
    hapid_indices = np.empty(len(included_hapids), dtype=np.int32)
    for i in range(len(included_hapids)):
        hapid = included_hapids[i]
        indices = np.where(hap_ids == hapid)[0]
        if len(indices) > 0:
            hapid_indices[i] = indices[0]
    return hapid_indices







def calulate_coefficient_per_group(record_df, logger=logger):
    # Generate a rank for haplotypes
    record_df["rank"] = rank_unique_values(record_df["var_count"].to_numpy() * 50 + record_df["indel_count"].to_numpy())
    # logger.info(f"Before calculating the coefficient for this region, the dataframe looks like :\n{record_df[:10].to_string(index=False)}\n")
    record_df["coefficient"] = calculate_coefficient(record_df.loc[:, ["start",
                                                                       "end",
                                                                       "total_depth",
                                                                       "hap_id",
                                                                       "hap_depth",
                                                                       "var_count",
                                                                       "indel_count",
                                                                       "rank"]].to_numpy(dtype=np.int32))
    # logger.info(f"After calculating the coefficient for this region, the dataframe looks like :\n{record_df[:10].to_string(index=False)}\n")
    return record_df


def sweep_region_inspection(input_bam,
                            output_bed = None,
                            depth_cutoff = 5,
                            window_size = 120,
                            step_size = 30,
                            logger = logger):

    if output_bed is None:
        output_bed = prepare_tmp_file(suffix = ".bed").name

    cmd = f"samtools depth {input_bam} | \
            mawk -F '\\t' '$3 >= {depth_cutoff}{{printf \"%s\\t%s\\t%s\\n\", $1, $2-1, $2;}}' | \
            bedtools merge -i stdin | \
            bedtools makewindows -b stdin -w {window_size} -s {step_size} > {output_bed}"

    executeCmd(cmd, logger = logger)

    op_bed_obj = pb.BedTool(output_bed)
    return op_bed_obj






def main_function(bam,
                  output_bam = None,
                  filter_out_bam = None,
                  intrinsic_bam = None,
                  raw_vcf = None,
                  bam_region_bed = None,
                  max_varno = 5,
                  mapq_cutoff = 20,
                  basequal_median_cutoff = 15,
                  edge_weight_cutoff = 0.201,
                  logger=logger):
    '''
    Main function for processing and analyzing BAM files to identify and filter misaligned reads.

    This function performs the following main tasks:
    1. Migrates BAM files to NCLS format for efficient query of reads or read pairs overlapping a query genomic region
    2. Builds a read-pair graph from the BAM data for phasing, where each vertex represents a read pair and each edge represents the confidence that two read pairs are originated from the same haplotype
    3. Identifies haplotypes in the graph by finding non-overlapping maximal cliques iteratively (using an approximate but not exact algorithm for efficient clique search, basically we iteratively run Bron-Kerbosch algorithm)
    4. After read pair grouped into different haplotypes, we put them into a binary integer linear programming model to solve the haplotype-level misalignment with HiGHs solver.
    5. Generates output BAM files with filtered and annotate haplotype index assigned to each read pairs.

    Parameters:
    - bam (str): Path to the input BAM file
    - output_bam (str, optional): Path for the output filtered BAM file, if not specified, output_bam will be the input bam with .clean.bam suffix
    - filter_out_bam (str, optional): Path for the BAM file containing filtered-out reads, if not specified, we do not output the filtered BAM file
    - intrinsic_bam (str): Path to the intrinsic BAM file, generated earlier in the SDrecall workflow (basically it aligns the reference sequence of one SD to the other SD)
    - raw_vcf (str): Path to the raw VCF file, the variants detected in the raw pooled alignments.
    - bam_region_bed (str, optional): Path to the BED file defining covered regions of the processing bam file, if not specified, bam_region_bed will be generated with a name after the input bam with .coverage.bed suffix
    - max_varno (float): Maximum variant number allowed
    - mapq_cutoff (int): Mapping quality cutoff to be included in the analysis
    - basequal_median_cutoff (int): Base quality median cutoff to be included in the analysis (if the median BQ of a read is lower than this cutoff, we will discard this read because it is too noisy)
    - edge_weight_cutoff (float): Edge weight cutoff separating two rounds of BK clique searches
    - logger (Logger object): Logger for output messages

    Returns:
    - phased_graph (Graph object): The constructed phasing graph

    This function integrates various analysis steps including BAM processing,
    graph construction, haplotype identification, and read filtering to improve
    the quality of genomic alignments and identify potential misalignments.
    '''

    # Given the top 1% mismatch count per read (one structrual variant count as 1 mismatch)
    tmp_bam = bam.replace(".bam", ".tmp.bam")

    if output_bam is None:
        replace = True
        output_bam = bam.replace(".bam", ".clean.bam")
    else:
        replace = False

    max_varno = float(max_varno)

    if filter_out_bam is None:
        filter_out_bam = output_bam.replace(".bam", ".noise.bam")

    noisy_qnames = set([])
    mismap_qnames = set([])
    norm_qnames = {"":set([])}
    noisy_num = 0
    total_num = 0

    bam_ncls = migrate_bam_to_ncls(bam,
                                   mapq_filter = mapq_cutoff,
                                   basequal_median_filter = basequal_median_cutoff,
                                   logger=logger)

    # parse the results from the tuple returned by migrate_bam_to_ncls
    ncls_dict, read_dict, qname_dict, qname_idx_dict, total_lowqual_qnames = bam_ncls
    logger.info(f"Successfully migrated the BAM file {bam} to NCLS format\n\n")

    # Now migrate the intrinsic BAM file to NCLS format
    intrin_bam_ncls = migrate_bam_to_ncls(intrinsic_bam,
                                          mapq_filter = 0,
                                          basequal_median_filter = 0,
                                          paired = False,
                                          logger=logger)
    # Since intrinsic BAM reads are reference sequences, therefore there are no low quality reads
    intrin_bam_ncls = intrin_bam_ncls[:-1]
    logger.info(f"Successfully migrated the intrinsic BAM file {intrinsic_bam} to NCLS format\n")
    logger.info(f"Containing {len(intrin_bam_ncls[1])} reads in total.\n\n")
    bam_graph = bam.replace(".bam", ".phased.graphml")

    # Calculate the mean read length of the input bam file, which can be used for read pair similarity calculation
    mean_read_length = calculate_mean_read_length(bam)

    # Create the read-pair graph used for phasing
    # Detailed description of the graph construction can be found in the function docstring.
    phased_graph, weight_matrix, qname_to_node, total_readhap_vector = build_phasing_graph(bam,
                                                                                           ncls_dict,
                                                                                           read_dict,
                                                                                           qname_dict,
                                                                                           qname_idx_dict,
                                                                                           mean_read_length,
                                                                                           edge_weight_cutoff = edge_weight_cutoff,
                                                                                           logger = logger)
    if phased_graph is None:
        return None

    logger.info(f"Now succesfully built the phasing graph with {phased_graph.num_vertices()} vertices and {phased_graph.num_edges()} edges. Save it to {bam_graph}\n\n")
    # Now we need to extract the components in the phased graph
    phased_graph.save(bam_graph)

    # Now we need to do local phasing for each component in the graph. (Finding non-overlapping high edge weight cliques inside each component iteratively)
    logger.info(f"Now start finding haplotypes in the setup weight matrix, the numba parallel threads are set to {get_num_threads()}")
    total_cliques = clique_generator_per_component(phased_graph,
                                                    weight_matrix,
                                                    ew_cutoff = edge_weight_cutoff,
                                                    logger = logger)

    total_cliques = list(total_cliques)

    # clique_sep_component_idx = 0
    qname_hap_info = defaultdict(int)
    qname_hap_info = find_components_inside_filtered_cliques( total_cliques,
                                                              phased_graph,
                                                              qname_hap_info,
                                                              weight_matrix,
                                                              edge_weight_cutoff,
                                                              logger = logger )
    gc.collect()

    logger.info(f"The final components are {qname_hap_info}")
    hap_qname_info = defaultdict(set)
    for vid, hid in qname_hap_info.items():
        qname = phased_graph.vp.qname[vid]
        hap_qname_info[hid].add(qname)
    logger.info(f"The final haplotype clusters are {hap_qname_info}")

    # Inspect the raw BAM corresponding variants to get the high density regions
    # It's like active region identification for GATK HC
    if not bam_region_bed:
        bam_region_bed = bam.replace(".bam", ".coverage.bed")
        cmd = f"bash {bash_utils_hub} samtools_bam_coverage \
                -i {bam} \
                -d 0 \
                -o {bam_region_bed}"
        executeCmd(cmd, logger = logger)

    # tobe_inspected_regions, var_df = extract_var_pos(raw_vcf, bam_region_bed, padding_size = 30, density_cutoff = 1/50, logger = logger)
    var_df = extract_var_pos(raw_vcf, bam_region_bed, logger = logger)

    total_genomic_haps = {}
    total_readerr_vector = {}
    compare_haplotype_meta_tab = bam.replace(".bam", ".haplotype_meta.tsv")

    correct_qnames, mismap_qnames = inspect_by_haplotypes(bam,
                                                          hap_qname_info,
                                                          qname_hap_info,
                                                          bam_ncls,
                                                          intrin_bam_ncls,
                                                          phased_graph,
                                                          weight_matrix,
                                                          qname_to_node,
                                                          total_lowqual_qnames,
                                                          total_readhap_vector,
                                                          total_readerr_vector,
                                                          total_genomic_haps,
                                                          var_df,
                                                          compare_haplotype_meta_tab = compare_haplotype_meta_tab,
                                                          mean_read_length = mean_read_length,
                                                          logger = logger )

    assert len(correct_qnames & mismap_qnames) == 0, f"The correct_qnames and mismap_qnames have overlap: {correct_qnames & mismap_qnames}"

    logger.info(f"In total has found {len(hap_qname_info)} clique separated components (haplotypes) in the target inspected regions.")
    logger.info(f"We found {len(mismap_qnames)} read pairs that are likely to be misaligned in the target regions.\n {mismap_qnames}\n And {len(correct_qnames)} read pairs that are likely to be correctly aligned in the target regions.\n {correct_qnames}\n")

    with pysam.AlignmentFile(bam, "rb") as bam_handle:
        # Extract the sample name (SM) from the existing read groups in the header
        with pysam.AlignmentFile(tmp_bam, "wb", header = bam_handle.header) as tmp_handle:
            with pysam.AlignmentFile(output_bam, "wb", header = bam_handle.header) as output_handle:
                with pysam.AlignmentFile(filter_out_bam, "wb", header = bam_handle.header) as noisy_handle:
                    # viewed_regions = {}
                    for read in bam_handle:
                        if read.is_supplementary or \
                           read.is_secondary:
                           continue
                        total_num += 1
                        qname = read.query_name
                        hap_id = qname_hap_info.get(qname_to_node.get(qname, -1), "NA")
                        if qname in mismap_qnames:
                            hap_id = f"{hap_id}_HIGHVD"
                        if qname in total_lowqual_qnames:
                            hap_id = f"{hap_id}_LOWQUAL"
                        # Filter out reads with oddly high editing distance that breakthrough the cutoff
                        if read.is_mapped and read.mapping_quality > 10:
                            gap_sizes = [t[1] for t in read.cigartuples if t[0] in [1,2] and t[1] > 1]
                            max_gap_size = max(gap_sizes) if len(gap_sizes) > 0 else 0
                            edit_dist = read.get_tag("NM")
                            scatter_edit_dist = edit_dist - max_gap_size
                            if qname in total_lowqual_qnames:
                                pass
                            elif qname in noisy_qnames:
                                # logger.debug(f"Read {1 if read.is_read1 else 2} from pair {read.query_name} has an edit distance {scatter_edit_dist} which is too high to be likely correctly mapped")
                                noisy_handle.write(read)
                            elif qname in mismap_qnames:
                                # logger.debug(f"Read {1 if read.is_read1 else 2} from pair {read.query_name} has an edit distance {scatter_edit_dist} which is significantly higher than other reads so it is considered as noisy and misaligned")
                                noisy_handle.write(read)
                            elif scatter_edit_dist > max_varno and qname not in correct_qnames:
                                noisy_num += 1
                                logger.info(f"Read {1 if read.is_read1 else 2} from pair {read.query_name} has an edit distance {scatter_edit_dist} which is so high that it is considered as noisy and misaligned")
                                noisy_qnames.add(qname)
                                noisy_handle.write(read)
                            elif not read.is_secondary and \
                                not read.is_supplementary and \
                                not read.is_duplicate and \
                                not read.is_qcfail:
                                # logger.debug(f"Read {1 if read.is_read1 else 2} from pair {read.query_name} has an edit distance {read.get_tag('NM')} which is considered as normal")
                                qname_pair = norm_qnames.get(qname, set([]))
                                qname_pair.add(read)
                                norm_qnames[qname] = qname_pair
                        else:
                            logger.info(f"Read {1 if read.is_read1 else 2} from pair {read.query_name} is not mapped or has low mapping quality {read.mapping_quality}, skip this read")

                        if qname in noisy_qnames:
                            hap_id = f"{hap_id}_HIGHVD"
                        read.set_tag('HP', f'HAP_{hap_id}')
                        tmp_handle.write(read)

                # logger.warning(f"Check {check_odd_num} reads to find oddly high editing distance reads, {zero_odd_num} reads found no odd editing distance read pairs")
                for qname, pair in norm_qnames.items():
                    if (not qname in noisy_qnames) and (not qname in mismap_qnames):
                        for read in pair:
                            output_handle.write(read)

    logger.warning(f"Filtered out {len(noisy_qnames)} noisy read-pairs (Editing distance without the biggest gap > {max_varno}) and {len(mismap_qnames - noisy_qnames)} read-pairs with ODD high editing distance, remaining {len(set(norm_qnames.keys()) - noisy_qnames - mismap_qnames)} read-pairs from {bam} (with total {total_num} reads) and output to {output_bam}\n\n")

    # Replace the input BAM file with the tmp BAM file with modified RG tags for visualization of haplotype clusters
    executeCmd(f"samtools sort -O bam -o {bam} {tmp_bam} && samtools index {bam} && rm {tmp_bam}", logger=logger)

    if replace:
        logger.info(f"Replacing {bam} with {output_bam}")
        executeCmd(f"samtools sort -O bam -o {bam} {output_bam} && samtools index {bam} && rm {output_bam}", logger=logger)
    else:
        logger.info(f"Generated a new BAM file: {output_bam}")
        tmp_bam = output_bam.replace(".bam", ".tmp.bam")
        executeCmd(f"samtools sort -O bam -o {tmp_bam} {output_bam} && mv {tmp_bam} {output_bam} && samtools index {output_bam}", logger = logger)

    cmd = f"samtools sort -O bam -o {tmp_bam} {filter_out_bam} && mv {tmp_bam} {filter_out_bam} && samtools index {filter_out_bam}"
    executeCmd(cmd, logger = logger)
    logger.info(f"Generated two BAM files: {filter_out_bam} and {output_bam}.")
    return phased_graph



if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("-f", "--function", type=str, help="The function name", required=True)
    parser.add_argument("-a", "--arguments", type=str, help="The function's input arguments, delimited by semi-colon ;", required=False, default=None)
    parser.add_argument("-k", "--key_arguments", type=str, help="Keyword arguments for the function, delimited by semi-colon ;", required=False, default=None)

    args = parser.parse_args()
    try:
        fargs = [ convert_input_value(a) for a in args.arguments.split(";") ] if type(args.arguments) == str else []
        fkwargs = { t.split("=")[0]: convert_input_value(t.split("=")[1]) for t in args.key_arguments.split(";") } if type(args.key_arguments) == str else {}
        logger.info("Running function: {}, input args are {}, input kwargs are {}".format(args.function, fargs, fkwargs))
    except Exception as e:
        logger.error("Input argument does not meet the expected format, encounter Parsing error {}, Let's check the input:\n-f {}, -a {}, -k {}".format(
            e,
            args.function,
            args.arguments,
            args.key_arguments
        ))
        raise e

    globals()[args.function](*fargs, **fkwargs)




