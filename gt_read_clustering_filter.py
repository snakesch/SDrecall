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
from highspy import Highs
import numba
# numba.config.THREADING_LAYER = 'omp'
# numba.set_num_threads(4)
from numba import types, prange, get_num_threads
from io import StringIO
from scipy import sparse
from ncls import NCLS
from collections import defaultdict
from intervaltree import IntervalTree
from python_utils import convert_input_value
from src.utils import executeCmd, prepare_tmp_file

bash_utils_hub = "shell_utils.sh"

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(module)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)




@numba.njit
def fast_median(data):
    return np.median(data)


@numba.njit
def numba_sum(data):
    return np.sum(data)


@numba.njit(types.bool_[:](types.bool_[:], types.bool_[:]), fastmath = True)
def numba_and(arr1, arr2):
    return np.logical_and(arr1, arr2)

@numba.njit(types.bool_[:](types.bool_[:]), fastmath=True)
def numba_not(arr):
    return np.logical_not(arr)

@numba.njit(types.float32(types.float32[:], types.boolean[:]), fastmath=True)
def numba_max(data, index_mask=None):
    if index_mask is not None:
        return np.max(data[index_mask])
    else:
        return np.max(data)

@numba.njit(types.int32(types.float32[:,:], types.boolean[:], types.boolean[:]), fastmath=True)
def find_row_index_of_max(matrix, row_index_mask, column_index_mask):
    # Apply the masks to get the submatrix
    submatrix = matrix[row_index_mask, :][:, column_index_mask]

    # Find the index of the maximum value in the submatrix
    max_index = np.argmax(submatrix)

    # Convert the flat index to a row index
    row_index = max_index // submatrix.shape[1]

    # Return the actual row index in the original matrix
    return np.where(row_index_mask)[0][row_index]

@numba.njit
def custom_all_numba(array):
    has_true = False
    for value in array:
        if value == 0:
            continue
        elif value == 1:
            has_true = True
        elif value == -1:
            return False
    return has_true


@numba.njit
def any_false_numba(array):
    has_false = False
    for value in array:
        if value == -1:
            return True
    return has_false

@numba.njit(parallel=True)
def numba_isin(arr, uniq_vals):
    '''
    return the index (as booleans) of the array that their corresponding value is in another array
    '''
    mask = np.zeros(arr.shape, dtype=np.bool_)
    uniq_set = set(uniq_vals)
    for i in prange(arr.size):
        if arr[i] in uniq_set:
            mask[i] = True
    return mask

@numba.njit(types.boolean[:](types.int32, types.int32[:]), fastmath=True)
def boolean_mask(total_rows, row_indices):
    boolean_mask = np.zeros(total_rows, dtype=np.bool_)
    for index in row_indices:
        boolean_mask[index] = True
    return boolean_mask


@numba.njit(types.boolean[:](types.int32, types.int32[:]), fastmath=True)
def reverse_boolean_mask(total_rows, row_indices):
    boolean_mask = np.ones(total_rows, dtype=np.bool_)
    for index in row_indices:
        boolean_mask[index] = False
    return boolean_mask

@numba.njit(types.float32[:,:](types.float32[:,:], types.boolean[:]), fastmath=True)
def numba_ix(matrix, mask):
    return matrix[mask, :][:, mask]


@numba.njit(types.Tuple((types.int32, types.float32))(types.float32[:], types.boolean[:]), fastmath=True)
def numba_max_idx_mem(data, index_mask=None):
    max_val = -np.inf
    max_idx = -1

    if index_mask is None:
        for i in range(data.size):
            if data[i] >= max_val:
                max_idx = i
                max_val = data[i]
    else:
        for i in range(data.size):
            if index_mask[i]:
                if data[i] >= max_val:
                    max_idx = i
                    max_val = data[i]

    return max_idx, max_val



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

    This function processes a BAM file, filters reads based on quality criteria, and organizes the data
    into NCLS structures for fast overlap queries. It's particularly useful for handling genomic interval
    data in large-scale sequencing projects.

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
    >>> overlapping_reads = ncls_dict[chrom].find_overlap(start, end)

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
    Use lazy_get_overlapping_reads instead for better memory efficiency.
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



def lazy_get_overlapping_reads(ncls_dict, read_dict, qname_dict, chrom, start, end):
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
    if chrom not in ncls_dict:
        yield from ()
    elif ncls_dict[chrom] is None:
        yield from ()
    else:
        _, overlapping_read_qnames = ncls_dict[chrom].all_overlaps_both(np.array([start]), np.array([end]), np.array([0]))
        # overlapping_read_qnames = ncls_dict[chrom].find_overlap(start, end)
        for qname_idx in list(dict.fromkeys(overlapping_read_qnames)):
            for read in read_dict[qname_idx]:
                if read.reference_start < end and read.reference_end > start:
                    yield read


def check_edge(u, v, adj_set):
    '''
    1 means not sure
    2 means accept same haplotype
    -1 means reject same haplotype

    Cant use 0 here because python treats 0 as False
    '''
    return adj_set.get((u, v), False) or adj_set.get((v, u), False)


@numba.njit(types.int32(types.int32[:]), fastmath=True)
def count_continuous_indel_blocks(array):
    """
    This function counts the number of continuous (consecutive) blocks of negative numbers in a 1D array.

    Args:
    array (list or numpy array): Input array containing numbers.

    Returns:
    int: The number of continuous negative blocks.
    """
    is_var = (array == -6) | (array > 1)
    is_var_size = is_var.size
    if is_var_size == 0:
        return 0

    # Add False at the start and end to correctly count the transitions at the boundaries
    extended = np.empty(is_var_size + 2, dtype=np.bool_)
    extended[0] = False
    extended[-1] = False
    extended[1:-1] = is_var

    # Count drops from True to False which represent the end of a block of True values
    # For this method to work. Both ends cannot be True
    not_bools = numba_not(extended[1:])
    block_bools = numba_and(extended[:-1], not_bools)
    block_count = numba_sum(block_bools)

    return block_count



@numba.njit(types.int32(types.bool_[:]), fastmath=True)
def count_continuous_blocks(array):
    """
    This function counts the number of continuous (consecutive) blocks of True values in a 1D boolean array.

    Args:
    array (list or numpy array): Input array containing numbers.

    Returns:
    int: The number of continuous negative blocks.
    """
    is_var = array
    is_var_size = is_var.size
    if is_var_size == 0:
        return 0

    # Add False at the start and end to correctly count the transitions at the boundaries
    extended = np.empty(is_var_size + 2, dtype=np.bool_)
    extended[0] = False
    extended[-1] = False
    extended[1:-1] = is_var

    # Count drops from True to False which represent the end of a block of True values
    # For this method to work. Both ends cannot be True
    not_bools = numba_not(extended[1:])
    block_bools = numba_and(extended[:-1], not_bools)
    block_count = numba_sum(block_bools)

    return block_count



@numba.njit(types.int32(types.int32[:]), fastmath=True)
def count_snv(array):
    snv_bools = array == -4
    return count_continuous_blocks(snv_bools)


@numba.njit(types.int32(types.int32[:]), fastmath=True)
def count_var(array):
    return count_snv(array) + count_continuous_indel_blocks(array)



@numba.njit(types.boolean[:](types.int32), fastmath=True)
def create_default_true_mask(size):
    return np.ones(size, dtype=numba.bool_)


@numba.njit(types.boolean[:](types.int32), fastmath=True, parallel=True)
def para_create_default_true_mask(size):
    return np.ones(size, dtype=numba.bool_)


@numba.njit(types.boolean[:](types.int32), fastmath=True)
def create_default_false_mask(size):
    return np.zeros(size, dtype=numba.bool_)



@numba.njit(types.int32[:](types.int32, types.boolean[:]), fastmath=True, parallel=True)
def apply_index_mask(size, initial_index_mask):
    return np.arange(size, dtype=np.int32)[initial_index_mask]



@numba.njit(types.float32[:](types.int32[:], types.int32), fastmath=True)
def count_window_var_density(array, padding_size = 25):
    is_var = array != 1
    is_var_size = numba_sum(is_var)
    if is_var_size == 0:
        return np.zeros(array.size, dtype=np.float32)

    density_arr = np.empty(array.size, dtype = np.float32)
    for i in range(array.size):
        start = max(i-padding_size, 0)
        end = min(i+padding_size, array.size)
        iter_arr = array[start:end]
        var_count = count_var(iter_arr)
        density = var_count/(padding_size *2 + 1)
        density_arr[i] = density

    return density_arr



@numba.njit(types.Tuple((types.int32, types.int32))(types.int32[:], types.int32[:]), fastmath=True)
def ref_genome_similarity(query_read_vector,
                          genomic_hap_vector):
    '''
    This function is used to identify the haplotypes that the query read vector belongs to.
    The genomic_hap_vectors is a 2D numpy array
    The query_read_vector is a 1D numpy array
    The max_vector_distance is a float reflecting the max manhattan distance between a read and a genomic haplotype
    '''
    # assert len(genomic_hap_vectors.shape) == 2, f"The input genomic haplotype vectors should be a 2D numpy array, but the actual input is {genomic_hap_vectors}"
    # Remove the reference haplotype first
    if numba_sum(genomic_hap_vector != 1) == 0:
        return 0, 0

    # ref_read_vector = np.ones(query_read_vector.size, dtype=np.int32)
    # Calculate the Manhattan distance between the query read vector and each genomic haplotype
    alt_var_size = count_continuous_blocks(genomic_hap_vector != query_read_vector)
    var_size = count_var(query_read_vector)

    return var_size, alt_var_size




def get_hapvector_from_cigar(cigar_tuples,
                             query_sequence = None,
                             logger = logger):
    '''
    Convert CIGAR tuples to a haplotype vector with penalties for matches, mismatches, and gaps.
    Note that the length of the haplotype vector is strictly aligned to the covered reference genome span, only this way we can ensure that comparison between overlapping reads are strictly aligned to each other.
    Uses NumPy for efficient array operations.

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
    at each position in the read's alignment to the reference. It uses the CIGAR
    string to determine the alignment structure and the base qualities to estimate
    error probabilities.

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
        An array of float values representing error probabilities for each
        position in the read's alignment to the reference. Values range from
        0 to 1, where 0 indicates high confidence and 1 indicates low confidence
        or gaps in the alignment.

    Notes:
    ------
    - The function assumes that the CIGAR string separates matches (7) from
      mismatches (8).
    - Insertions and deletions are initially marked with a placeholder value (99),
      which is later converted to 0 (high confidence) in the error vector.
    - The final error vector is scaled based on the Phred quality scores,
      converting them to error probabilities.
    - The error probabilities are intended to represent the expected deviation
      from the original haplotype.

    Raises:
    -------
    AssertionError
        If the CIGAR string contains unexpected operations or if the resulting
        error vector length doesn't match the read's aligned length.

    Example:
    --------
    >>> read = pysam.AlignedSegment()
    >>> read.cigartuples = [(7, 10), (8, 1), (1, 2), (7, 5)]
    >>> read.query_qualities = [30] * 18
    >>> error_vector = get_errorvector_from_cigar(read, read.cigartuples)
    >>> print(error_vector)
    [0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.
     0.001 0.001 0.001 0.001 0.001]
    '''
    errorvector = np.empty(read.reference_end - read.reference_start, dtype=float)
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
            # Insertion, the base quality is the mean value of the inserted bases
            # The base quality does not affect the gaps appeared in the alignment, later 99 will be replaced by 0
            errorvector[ref_consume] = 99
            query_consume += length
        elif operation == 2:
            # Deletion, the base quality is the value of the base immediately upstream of the deletion
            # The base quality does not affect the gaps appeared in the alignment, later 99 will be replaced by 0
            errorvector[ref_consume:ref_consume + length] = np.array([99] * length)
            ref_consume += length

    assert errorvector.size == read.reference_end - read.reference_start, f"The error vector length is {len(errorvector)} while the read length is {read.reference_end - read.reference_start}. The cigar str is {read.cigarstring}"
    # Then we need to convert the phred_scaled base quality scores back to error probability
    # Since the error vector is actually the expected deviation from originality haplotype, we need to multiply with the score deviation between mismatch and match is 5 (match for 1, mismatch for -4)
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



# @numba.njit()
def prepare_ref_query_idx_map(qseq_ref_pos_arr):
    # numba_dict = Dict.empty(key_type=types.int32, value_type=types.int32)
    numba_dict = {}
    for qi in range(qseq_ref_pos_arr.size):
        ri = qseq_ref_pos_arr[qi]
        if ri >= 0:
            numba_dict[ri] = qi
    return numba_dict



def get_interval_seq(read, interval_start, interval_end, read_ref_pos_dict = {}):
    """
    Extract the sequence of a read for a given genomic interval.

    Parameters:
    - read: pysam.AlignedSegment
        The read to extract sequence from.
    - start, end: int
        The start and end positions of the interval in genomic coordinates.
        both should be inclusive
    - read_ref_pos_dict: dict
        Dictionary mapping read positions to reference positions.

    Returns:
    - str: Extracted sequence.
    - dict: Updated read_ref_pos_dict.
    - numpy.ndarray: Array of query positions corresponding to the extracted sequence.
    """
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

    # interval_ref_pos = ref_positions[interval_start_qidx: interval_end_qidx]
    # We need to also remove the inserted sequences
    interval_read_seq = read.query_sequence[interval_start_qidx:interval_end_qidx]
    qidx_ridx_arr = qseq_ref_positions[interval_start_qidx:interval_end_qidx]

    return interval_read_seq, read_ref_pos_dict, qidx_ridx_arr



def get_interval_seq_qual(read, interval_start, interval_end, read_ref_pos_dict = {}):
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
    2. Use the base at other reads to deduce whether this lowqual base is sequencing error or not
    '''

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



@numba.njit(types.bool_[:](types.int32[:]), fastmath=True)
def get_indel_bools(seq_arr):
    return (seq_arr > 1) | (seq_arr == -6)


@numba.njit
def numba_diff_indices(arr1, arr2):
    """
    Find the indices where two arrays differ.

    This function uses Numba for faster computation.

    Parameters:
    - arr1, arr2: numpy.ndarray
        The two arrays to compare.

    Returns:
    - numpy.ndarray: Indices where the arrays differ.
    """
    return np.where(arr1 != arr2)[0]



@numba.njit(types.bool_(types.int8[:], types.int8[:], types.int8), fastmath=True)
def compare_sequences(read_seq, other_seq, except_char):
    """
    Compare two sequences, handling 'N' bases.

    Parameters:
    - seq1, seq2: numpy.ndarray
        Arrays representing the sequences to compare.
    - N_value: int
        Integer representation of 'N' in the sequence arrays.

    Returns:
    - bool: True if sequences match (ignoring 'N's), False otherwise.
    """
    if read_seq.shape != other_seq.shape:
        return False

    total_match = True
    for i in range(read_seq.size):
        if read_seq[i] != other_seq[i]:
            if read_seq[i] != except_char and other_seq[i] != except_char:
                return False

    return total_match



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

    if read_id in read_hap_vectors:
        read_hap_vector = read_hap_vectors[read_id]
    else:
        read_hap_vector = get_hapvector_from_cigar(read.cigartuples, read.query_sequence, logger = logger)
        read_hap_vectors[read_id] = read_hap_vector

    other_read_id = get_read_id(other_read)
    other_start = other_read.reference_start
    other_end = other_read.reference_end

    if other_read_id in read_hap_vectors:
        other_read_hap_vector = read_hap_vectors[other_read_id]
    else:
        other_read_hap_vector = get_hapvector_from_cigar(other_read.cigartuples, other_read.query_sequence, logger = logger)
        read_hap_vectors[other_read_id] = other_read_hap_vector

    interval_hap_vector = read_hap_vector[overlap_start - start:overlap_end - start]
    interval_other_hap_vector = other_read_hap_vector[overlap_start - other_start:overlap_end - other_start]

    overlap_span = overlap_end - overlap_start
    # assert interval_hap_vector.size > 0, f"The interval hap vector is {interval_hap_vector} for read {read.query_name} ({start}-{read.reference_end}) and other read {other_read.query_name} ({other_start}-{other_read.reference_end}) and the overlap span is {overlap_start} to {overlap_end}"
    # assert interval_hap_vector.size == interval_other_hap_vector.size, f"The interval hap vector size {interval_hap_vector} is not the same as the other interval hap vector size {interval_other_hap_vector}"

    diff_indices = numba_diff_indices(interval_hap_vector, interval_other_hap_vector)

    # First test whether there are too many mismatches between two reads that we wont tolerate
    if len(diff_indices) >= 3:
        # Cannot tolerate such mismatches
        return False, read_ref_pos_dict, read_hap_vectors, None

    diff_pos_read = interval_hap_vector[diff_indices]
    diff_pos_oread = interval_other_hap_vector[diff_indices]

    read_indel_overlap = get_indel_bools(diff_pos_read)
    oread_indel_overlap = get_indel_bools(diff_pos_oread)

    if numba_sum(read_indel_overlap) > 0 or numba_sum(oread_indel_overlap) > 0:
        return False, read_ref_pos_dict, read_hap_vectors, None

    # Now use overlap_start and overlap_end to extract the read sequence of two reads within the overlap region
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
        if edge_weight >= mean_read_length - 50:
            return True, read_ref_pos_dict, read_hap_vectors, edge_weight
        else:
            return np.nan, read_ref_pos_dict, read_hap_vectors, None
    else:
        if len(read_seq) != len(other_seq) or interval_hap_vector.size != interval_other_hap_vector.size:
            return False, read_ref_pos_dict, read_hap_vectors, None

        # Here we know that there are only equal to or smaller than 2 mismatches between the two reads
        # Then we need to know whether the two mismatches are in-trans variants or in-cis variants, if its in-trans we cannot tolerate the two mismatches
        if read_seq_arr is None:
            base_dict = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
            read_seq_arr = np.array(list(map(base_dict.get, read_seq)), dtype=np.int8)
            other_seq_arr = np.array(list(map(base_dict.get, other_seq)), dtype=np.int8)

        q_diff_indices = numba_diff_indices(read_seq_arr, other_seq_arr)
        r_diff_indices = qr_idx_arr[q_diff_indices]
        if -1 in r_diff_indices:
            return False, read_ref_pos_dict, read_hap_vectors, None
        abs_diff_inds = r_diff_indices + overlap_start
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


@numba.njit(types.float32[:](types.float32[:, :], types.boolean[:], types.float32[:]), fastmath=True)
def row_wise_max_with_mask_nb(matrix, index_mask, mask_values):
    row_max_values = np.zeros(matrix.shape[0], dtype=np.float32)

    for i in range(matrix.shape[0]):
        if not index_mask[i]:
            continue

        max_val = -np.inf
        for j in range(matrix.shape[1]):
            if index_mask[j] and matrix[i, j] != mask_values[0] and i != j:
                max_val = max(max_val, matrix[i, j])

        row_max_values[i] = max_val if max_val != -np.inf else 0

    return row_max_values


@numba.njit(types.float32[:](types.float32[:, :], types.boolean[:], types.float32[:]), fastmath=True, parallel=True)
def para_row_wise_max_with_mask_nb(matrix, index_mask, mask_values):
    row_max_values = np.zeros(matrix.shape[0], dtype=np.float32)

    for i in numba.prange(matrix.shape[0]):
        if not index_mask[i]:
            continue

        max_val = -np.inf
        for j in range(matrix.shape[1]):
            if index_mask[j] and matrix[i, j] != mask_values[0] and i != j:
                max_val = max(max_val, matrix[i, j])

        row_max_values[i] = max_val if max_val != -np.inf else 0

    return row_max_values

def pretty_print_matrix(matrix, precision=3):
    """
    Pretty prints a 2D NumPy array.

    Args:
        matrix (numpy.ndarray): The 2D array to be printed.
    """
    with np.printoptions(precision=precision, suppress=True):
        return "\n".join(["\t".join(map("{:.3f}".format, row)) for row in matrix])



def iteratively_find_largest_cliques(selected_indices,
                                     weight_matrix,
                                     cutoff = 0.1,
                                     qname_dict = {},
                                     logger = logger):
    # Generate a index_mask array for the indices of rows in the weight_matrix that the indices are in the initial_vert_indices
    # weight_matrix = weight_matrix[np.ix_(initial_index_mask, initial_index_mask)]
    # weight_matrix[np.ix_(initial_index_mask, initial_index_mask)]
    if weight_matrix.shape[0] > 1000000:
        parallel = True
    else:
        parallel = False

    if parallel:
        index_mask = para_create_default_true_mask(weight_matrix.shape[0])
    else:
        index_mask = create_default_true_mask(weight_matrix.shape[0])

    sparse_matrix = sparse.csr_matrix(weight_matrix)
    sparse_matrix_data = sparse_matrix.data
    sparse_matrix_indices = sparse_matrix.indices
    sparse_matrix_indptr = sparse_matrix.indptr
    matrix_size = sparse_matrix.shape[0]

    assert selected_indices.size == matrix_size, f"The selected indices size is {selected_indices.size} and the matrix size is {matrix_size}"

    while True:
        logger.info(f"Start to find the largest clique inside a subgraph of {index_mask.size} nodes, left {numba_sum(index_mask)} nodes to be inspected.")
        if matrix_size == 0:
            logger.warning(f"Found that in this iteration, the subgraph size is already 0")
            break

        clique_indices, drop_mask = heuristic_find_largest_edge_weight_clique_sparse(sparse_matrix_data,
                                                                                     sparse_matrix_indices,
                                                                                     sparse_matrix_indptr,
                                                                                     matrix_size,
                                                                                     index_mask,
                                                                                     cutoff = cutoff)
        # This is utter important. It makes sure the returned vertex index are in the original total graph
        assert np.max(clique_indices) < selected_indices.size, f"The returned clique indices are out of the range of the selected indices. The max index is {np.max(clique_indices)} and the size of selected indices is {selected_indices.size}"
        raw_clique_indices = selected_indices[clique_indices]
        # logger.info(f"The returned qname idx are  {raw_clique_indices.tolist()}")
        # logger.info(f"Found a clique containing qnames {[(qname_dict[qname_idx], qname_idx) for qname_idx in raw_clique_indices]} composing the largest clique in the subgraph. Remaining {drop_mask.sum()} nodes in this iteration. ")
        remain_ind_count = numba_sum(drop_mask)
        logger.info(f"Found {clique_indices.size} composing the largest clique in the subgraph. Remaining {remain_ind_count} nodes after this iteration. ")
        if raw_clique_indices.size >= 1:
            yield set(raw_clique_indices)
        elif remain_ind_count > 0:
            logger.error(f"Found that in this iteration, the clique size is 0. Meaning highs decides to drop all the nodes to break all the pairs of -1 value row-col indices. Which is odd and intolerable")
            # logger.error(f"Lets check the input matrix for this iteration \n{pretty_print_matrix(weight_matrix[np.ix_(drop_mask, index_mask)])}\n")
            break

        if remain_ind_count == 0:
            logger.warning(f"Found that after this iteration, remain 0 nodes. All the nodes has been assigned to cliques")
            break

        # Prepare the index_mask for next_round
        # index_mask = numba_and(drop_mask, index_mask)
        selected_indices = selected_indices[drop_mask]
        # sparse_matrix = sparse_matrix[np.ix_(drop_mask, drop_mask)]
        # weight_matrix = numba_ix(weight_matrix, drop_mask)
        sparse_matrix = sparse_matrix[drop_mask, :][:, drop_mask]
        sparse_matrix_data = sparse_matrix.data
        sparse_matrix_indices = sparse_matrix.indices
        sparse_matrix_indptr = sparse_matrix.indptr
        matrix_size = sparse_matrix.shape[0]
        # sparse_matrix_data, sparse_matrix_indices, sparse_matrix_indptr, matrix_size = sparse_mask_selection( sparse_matrix_data,
        #                                                                                                       sparse_matrix_indices,
        #                                                                                                       sparse_matrix_indptr,
        #                                                                                                       drop_mask )

        if matrix_size < 1000000:
            parallel = False

        if parallel:
            index_mask = para_create_default_true_mask(matrix_size)
        else:
            index_mask = create_default_true_mask(matrix_size)



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


@numba.njit(types.Tuple((types.int32[:], types.boolean[:]))(types.float32[:, :], types.boolean[:], types.int32[:], types.float32), fastmath=True)
def heuristic_find_largest_edge_weight_clique(weight_matrix, initial_index_mask, raw_indices, cutoff = 0.1):
    # logger.info(f"\n{pretty_print_matrix(weight_matrix)}\n")
    if weight_matrix.shape[0] > 5000:
        parallel = True
    else:
        parallel = False

    if parallel:
        row_wise_maxs = para_row_wise_max_with_mask_nb(weight_matrix, initial_index_mask, mask_values = np.array([-1.0], dtype=np.float32))
    else:
        row_wise_maxs = row_wise_max_with_mask_nb(weight_matrix, initial_index_mask, mask_values = np.array([-1.0], dtype=np.float32))
    select_indices = np.empty(weight_matrix.shape[0], dtype=np.int32)

    if weight_matrix.shape[0] > 40000:
        parallel = True
    else:
        parallel = False

    if parallel:
        max_row_ind, max_value = numba_max_idx_mem(row_wise_maxs, initial_index_mask)
    else:
        max_row_ind, max_value = numba_max_idx_mem(row_wise_maxs, initial_index_mask)

    select_indices[0] = max_row_ind
    neg_1_mask = weight_matrix[max_row_ind] != -1
    neg_1_mask[max_row_ind] = False

    index_mask = numba_and(initial_index_mask, neg_1_mask)
    i = 1
    # print(f"Initially use {max_row_ind} (raw index is {raw_indices[max_row_ind]}) as the first selected index", file=sys.stderr)
    # selected_mask = create_default_false_mask(weight_matrix.shape[0])
    # selected_mask[max_row_ind] = True
    while numba_sum(index_mask) > 0:
        # logger.info(f"{max_row_ind}, {index_mask.tolist()}")
        # max_row_ind = find_row_index_of_max(weight_matrix, index_mask, selected_mask)
        if parallel:
            next_max_ind, next_max_value = numba_max_idx_mem(weight_matrix[max_row_ind], index_mask)
        else:
            next_max_ind, next_max_value = numba_max_idx_mem(weight_matrix[max_row_ind], index_mask)
        if next_max_value <= cutoff:
            # print(f"Next max value is {next_max_value} and its smaller than the cutoff {cutoff} while the corresponding index is {next_max_ind}, the raw index is {raw_indices[next_max_ind]}", file=sys.stderr)
            # Cannot be used to extend a phase group
            trial = 0
            while next_max_value <= cutoff and trial < i:
                # Try other selected_indice
                select_indice = select_indices[trial]
                if parallel:
                    next_max_ind, next_max_value = numba_max_idx_mem(weight_matrix[select_indice], index_mask)
                else:
                    next_max_ind, next_max_value = numba_max_idx_mem(weight_matrix[select_indice], index_mask)
                trial += 1

            if trial >= i and next_max_value <= cutoff:
                # Run out of paths
                # print(f"Run out of paths to extend the clique. The next max value is {next_max_value} and the cutoff is {cutoff}. Directly break", file=sys.stderr)
                break
            else:
                max_row_ind = next_max_ind

        else:
            # print(f"Found the next max value is {next_max_value} larger than the cutoff {cutoff} and the corresponding index is {next_max_ind}, raw index is {raw_indices[next_max_ind]}, last max raw index is {raw_indices[max_row_ind]}", file=sys.stderr)
            max_row_ind = next_max_ind


        select_indices[i] = max_row_ind
        # selected_mask[max_row_ind] = True
        i += 1

        if i >= weight_matrix.shape[0]:
            break

        neg_1_mask = weight_matrix[max_row_ind] != -1
        neg_1_mask[max_row_ind] = False
        # not_selected_mask = numba_not(selected_mask)
        index_mask = numba_and(index_mask, neg_1_mask)

    select_indices = select_indices[:i]
    remain_mask = reverse_boolean_mask(weight_matrix.shape[0], select_indices)
    remain_mask = numba_and(remain_mask, initial_index_mask)
    return select_indices, remain_mask


@numba.njit(types.float32[:](types.float32[:], types.int32[:], types.int32[:], types.int32, types.boolean[:], types.float32), fastmath=True)
def row_wise_max_with_mask_sparse(matrix_data,
                                  matrix_indices,
                                  matrix_indptr,
                                  matrix_size,
                                  index_mask,
                                  mask_value=-1):
    # assert sparse.isspmatrix_csr(matrix), "Input matrix must be in CSR format"
    row_max_values = np.zeros(matrix_size, dtype=np.float32)

    for i in range(matrix_size):
        if not index_mask[i]:
            continue

        row_start, row_end = matrix_indptr[i], matrix_indptr[i+1]
        row_data = matrix_data[row_start:row_end]
        row_cols = matrix_indices[row_start:row_end]

        index_mask_i = index_mask[i]
        index_mask[i] = False
        valid_mask = (row_data != mask_value) & index_mask[row_cols]
        index_mask[i] = index_mask_i
        valid_data = row_data[valid_mask]

        if valid_data.size > 0:
            row_max_values[i] = np.max(valid_data)

    return row_max_values



@numba.njit(types.float32[:](types.float32[:], types.int32[:], types.int32[:], types.int32, types.boolean[:], types.float32), fastmath=True, parallel=True)
def para_row_wise_max_with_mask_sparse(matrix_data,
                                       matrix_indices,
                                       matrix_indptr,
                                       matrix_size,
                                       index_mask,
                                       mask_value=-1):
    # assert sparse.isspmatrix_csr(matrix), "Input matrix must be in CSR format"
    row_max_values = np.zeros(matrix_size, dtype=np.float32)

    for i in prange(matrix_size):
        if not index_mask[i]:
            continue

        row_start, row_end = matrix_indptr[i], matrix_indptr[i+1]
        row_data = matrix_data[row_start:row_end]
        row_cols = matrix_indices[row_start:row_end]

        index_mask_i = index_mask[i]
        index_mask[i] = False
        valid_mask = (row_data != mask_value) & index_mask[row_cols]
        index_mask[i] = index_mask_i
        valid_data = row_data[valid_mask]

        if valid_data.size > 0:
            row_max_values[i] = np.max(valid_data)

    return row_max_values



@numba.njit(types.boolean[:](types.float32[:], types.int32[:], types.int32, types.float32), fastmath=True)
def efficient_mask(sarr_data,
                   sarr_indices,
                   sarr_size,
                   mask_value=-1.0):
    '''
    The input should be the data, indices, arr_size of a one row sparse_matrix instance
    sparse_matrix[row_index] operation:
    When you perform sparse_matrix[row_index] on a CSR matrix, it returns a csr_matrix representing that single row. This is a 1 x n sparse matrix.

    Attributes of the returned row:
    The returned row matrix has its own data, indices, and indptr attributes, but these are views of the original matrix's data. This means:
    data: Contains only the non-zero values for that specific row
    indices: Contains the column indices for the non-zero values in that row
    indptr: For a single row, this will always be [0, number_of_nonzeros_in_row]

    Efficiency implications:
    Because these are views, accessing sparse_matrix[row_index] is very efficient. It doesn't create a copy of the data, but rather provides a view into the original matrix's data structure.
    '''
    # Extract the row without converting the whole row to dense format
    neg_1_mask = np.ones(sarr_size, dtype=np.bool_)
    data_mask = sarr_data == mask_value
    neg_1_mask[sarr_indices[data_mask]] = False

    return neg_1_mask


@numba.njit(types.Tuple((types.int32, types.float32))(types.float32[:], types.int32[:], types.boolean[:]), fastmath=True)
def efficient_row_max(sarr_data,
                      sarr_indices,
                      index_mask):
    '''
    The input should be the data, indices, arr_size of a one row sparse_matrix instance
    sparse_matrix[row_index] operation:
    When you perform sparse_matrix[row_index] on a CSR matrix, it returns a csr_matrix representing that single row. This is a 1 x n sparse matrix.

    Attributes of the returned row:
    The returned row matrix has its own data, indices, and indptr attributes, but these are views of the original matrix's data. This means:
    data: Contains only the non-zero values for that specific row
    indices: Contains the column indices for the non-zero values in that row
    indptr: For a single row, this will always be [0, number_of_nonzeros_in_row]
    '''
    valid_mask = numba_and((sarr_data != -1), index_mask[sarr_indices])
    if np.any(valid_mask):
        # assert len(sarr_data) == len(valid_mask), f"The mask are not of the equal size of array data"
        max_ind_short, max_value = numba_max_idx_mem(sarr_data, valid_mask)
        # logger.info(f"max_ind_short is {max_ind_short}, max_value is {max_value}, valid_mask sum is {np.sum(valid_mask)}. Original sarr_data is {sarr_data.tolist()}")
        max_index = sarr_indices[max_ind_short]
        # logger.info(f"max_index is {max_index}, max_ind_short is {max_ind_short}")
        return max_index, max_value
    return -1, -2.0


@numba.njit(types.Tuple((types.int32[:], types.boolean[:]))(types.float32[:], types.int32[:], types.int32[:], types.int32, types.boolean[:], types.float32), fastmath=True)
def heuristic_find_largest_edge_weight_clique_sparse(matrix_data,
                                                     matrix_indices,
                                                     matrix_indptr,
                                                     matrix_size,
                                                     initial_index_mask,
                                                     cutoff=0.1):
    # assert sparse.isspmatrix_csr(weight_matrix), "Input matrix must be in CSR format"
    if matrix_size > 5000:
        parallel = True
    else:
        parallel = False

    if parallel:
        row_wise_maxs = para_row_wise_max_with_mask_sparse(matrix_data,
                                                           matrix_indices,
                                                           matrix_indptr,
                                                           matrix_size,
                                                           initial_index_mask,
                                                           -1.0)
    else:
        row_wise_maxs = row_wise_max_with_mask_sparse(matrix_data,
                                                      matrix_indices,
                                                      matrix_indptr,
                                                      matrix_size,
                                                      initial_index_mask,
                                                      -1.0)

    select_indices = np.empty(matrix_size, dtype=np.int32)
    # print(len(row_wise_maxs))

    max_row_ind, max_value = numba_max_idx_mem(row_wise_maxs, initial_index_mask)

    # print(max_row_ind, max_value)
    index_mask = initial_index_mask.copy()
    select_indices[0] = max_row_ind

    max_row_start = matrix_indptr[max_row_ind]
    max_row_end = matrix_indptr[max_row_ind+1]
    max_row_data = matrix_data[max_row_start: max_row_end]
    max_row_cols = matrix_indices[max_row_start: max_row_end]

    neg_1_mask = efficient_mask(max_row_data,
                                max_row_cols,
                                matrix_size,
                                -1.0)

    neg_1_mask[max_row_ind] = False  # Exclude the element itself
    index_mask = numba_and(index_mask, neg_1_mask)

    i = 1
    while np.any(index_mask):
        # print(f"max_row_ind is {max_row_ind}, while the total array size is {matrix_size}, and the index_mask sum is {np.sum(index_mask)}")
        next_max_ind, next_max_value = efficient_row_max(max_row_data,
                                                         max_row_cols,
                                                         index_mask)

        if next_max_value <= cutoff:
            trial = 0
            while next_max_value <= cutoff and trial < i:
                select_indice = select_indices[trial]
                select_row_start = matrix_indptr[select_indice]
                select_row_end = matrix_indptr[select_indice+1]
                select_row_data = matrix_data[select_row_start: select_row_end]
                select_row_cols = matrix_indices[select_row_start: select_row_end]
                next_max_ind, next_max_value = efficient_row_max(select_row_data,
                                                                 select_row_cols,
                                                                 index_mask)
                trial += 1
            if next_max_value <= cutoff and trial >= i:
                # The two condition cant be satisfied at the same time, which it will stay in while loop
                break
            else:
                max_row_ind = next_max_ind
        else:
            max_row_ind = next_max_ind

        select_indices[i] = max_row_ind
        i += 1

        if i >= matrix_size:
            break

        max_row_start = matrix_indptr[max_row_ind]
        max_row_end = matrix_indptr[max_row_ind+1]
        max_row_data = matrix_data[max_row_start: max_row_end]
        max_row_cols = matrix_indices[max_row_start: max_row_end]

        neg_1_mask = efficient_mask(max_row_data,
                                    max_row_cols,
                                    matrix_size,
                                    -1.0)

        neg_1_mask[max_row_ind] = False  # Exclude the element itself
        index_mask = numba_and(index_mask, neg_1_mask)


    select_indices = select_indices[:i]
    remain_mask = reverse_boolean_mask(matrix_size, select_indices)
    remain_mask = numba_and(remain_mask, initial_index_mask)
    return select_indices, remain_mask



def find_cliques_in_every_component(graph, weight_matrix, ew_cutoff = 0.101, logger = logger):
    # Find all components for this graph
    components, _ = gt.label_components(graph)
    component_dict = defaultdict(set)
    for v in graph.vertices():
        cid = components[v]
        component_dict[cid].add(np.int32(v))

    # Find row indices where the entire row is below or equal to 0.1
    small_row_mask = np.all(weight_matrix <= 0.1, axis=1)
    big_row_mask = np.logical_not(small_row_mask)
    small_row_indices = set(np.where(small_row_mask)[0])

    if not weight_matrix.flags['C_CONTIGUOUS']:
        weight_matrix = np.ascontiguousarray(weight_matrix)

    logger.info(f"Found {len(component_dict)} components in the graph. The big weight matrix has {numba_sum(big_row_mask)} rows and columns. The original weight_matrix has {weight_matrix.shape[0]} rows and columns and its contiguity in memory is {weight_matrix.flags['C_CONTIGUOUS']}. The small row indices are {small_row_indices}")

    # total_cliques = []
    for component_id, component_verts in component_dict.items():
        # logger.info(f"Before the iteration start, the 515, 216 cell value for weight matrix is {weight_matrix[515, 216]}")
        comp_index_mask = numba_isin(np.arange(weight_matrix.shape[0], dtype=np.int32), component_verts)
        comp_index_mask = numba_and(comp_index_mask, big_row_mask)
        selected_indices = apply_index_mask(weight_matrix.shape[0], comp_index_mask)
        big_weight_matrix = weight_matrix[np.ix_(comp_index_mask, comp_index_mask)]
        # big_weight_matrix = numba_ix(weight_matrix, comp_index_mask)
        logger.info(f"Start to find the largest clique in the component {component_id}, which contains {len(component_verts)} vertices. Contiguous? {big_weight_matrix.flags['C_CONTIGUOUS']}. Does the current matrix a view of the original one ? {big_weight_matrix.base is weight_matrix}")
        # logger.info(f"The selected indices are {selected_indices.tolist()}, here are the component vertices: {component_verts}")
        if len(component_verts) <= 3:
            small_row_indices.update(component_verts)
            logger.info(f"Adding the {len(component_verts)} vertices in the component {component_id} to the small row indices")
            continue
        else:
            cliques_iter = iteratively_find_largest_cliques(selected_indices, big_weight_matrix, cutoff = ew_cutoff, qname_dict = graph.vertex_properties['qname'], logger = logger)
            # logger.info(f"Found {len(cliques)} cliques in the component {component_id}\n")
            for clique in cliques_iter:
                logger.info(f"Receiving a clique containing {[graph.vertex_properties['qname'][qid] for qid in clique]} in the component {component_id}")
                if len(clique) <= 3:
                    logger.info(f"Found {[graph.vertex_properties['qname'][qid] for qid in clique]} read pairs in a very small clique.")
                    small_row_indices.update(clique)
                    logger.info(f"Adding the {len(clique)} vertices in the clique to the small row indices")
                    continue
                else:
                    logger.info(f"The clique is big enough to be directly yielded")
                    yield clique

    logger.info(f"Remaining {len(small_row_indices)} vertices that are not included in the cliques. Here we find cliques again among them:\n{small_row_indices}")
    # return total_cliques
    if len(small_row_indices) > 0:
        small_row_mask = numba_isin(np.arange(weight_matrix.shape[0], dtype=np.int32), small_row_indices)
        selected_indices = apply_index_mask(weight_matrix.shape[0], small_row_mask)
        small_weight_matrix = weight_matrix[np.ix_(small_row_mask, small_row_mask)]
        # small_weight_matrix = numba_ix(weight_matrix, small_row_mask)
        logger.info(f"Start to find the largest clique in the small weight matrix, which contains {numba_sum(small_row_mask)} rows and columns. Contiguous? {small_weight_matrix.flags['C_CONTIGUOUS']}")
        cliques_iter = iteratively_find_largest_cliques(selected_indices, small_weight_matrix, cutoff = ew_cutoff/2, qname_dict = graph.vertex_properties['qname'], logger = logger)
        for clique in cliques_iter:
            yield clique



# def refine_small_cliques_in_every_component(graph, weight_matrix, small_clique_verts, logger = logger):
#     # Find all components for this graph
#     components, _ = gt.label_components(graph)
#     component_dict = defaultdict(set)
#     for v in graph.vertices():
#         cid = components[v]
#         component_dict[cid].add(int(v))

#     for component_id, component_verts in component_dict.items():
#         logger.info(f"Start to find the largest clique in the component {component_id}, which contains {len(component_verts)} vertices")
#         if len(component_verts) <= 2:
#             yield set(component_verts)
#         else:
#             component_small_verts = set([v for v in component_verts if v in small_clique_verts])
#             cliques = iteratively_refine_small_cliques(component_small_verts, weight_matrix, logger = logger)
#             # logger.info(f"Found {len(cliques)} cliques in the component {component_id}\n")
#             for clique in cliques:
#                 yield clique


def graph_vertex_iter(vertex_indices, graph):
    for vid in vertex_indices:
        yield graph.vertex(vid)



def find_components_inside_filtered_cliques(final_cliques,
                                            graph,
                                            final_components,
                                            weight_matrix,
                                            ew_cutoff = 0.101,
                                            logger=logger):
    # Drop the zero weight edges for the graph before performing this function
    # For each clique, use vfilter to extract the subgraph, then use efilter to only view the weight > 0 edges
    # Then use gt.label_components to find the connected components
    # final_components = defaultdict(int)
    # last_c_idx = 0

    # vis_qnames = ["HISEQ1:59:HB66DADXX:1:1110:1670:77678:PC0",
    #               "HISEQ1:63:HB65FADXX:2:2214:15844:10511:PC0",
    #               "HISEQ1:61:HB66HADXX:1:1213:7009:61371:PC0",
    #               "HISEQ1:63:HB65FADXX:1:1114:2865:3570:PC357"]

    haplotype_idx = 0
    for clique in final_cliques:
        # Prepare a boolean v property map used for vfilt
        vertex_filter = graph.new_vertex_property("bool", val=False)
        # Set the filter to True for vertices in the list
        for vid in clique:
            vertex_filter[graph.vertex(vid)] = True

        subgraph = gt.GraphView(graph, vfilt = vertex_filter)
        clique = list(clique)
        # Supplement the small weight edges to the graph
        for i in range(len(clique)):
            for j in range(i, len(clique)):
                if i != j:
                    if weight_matrix[clique[i], clique[j]] > 0.05 and weight_matrix[clique[i], clique[j]] <= ew_cutoff:
                        logger.info(f"Adding an edge between {graph.vertex_properties['qname'][graph.vertex(clique[i])]}, {clique[i]} and {graph.vertex_properties['qname'][graph.vertex(clique[j])]}, {clique[j]} because their weight in matrix at {clique[i]} row and {clique[j]} column is {weight_matrix[clique[i], clique[j]]}")
                        v1 = subgraph.vertex(clique[i])
                        v2 = subgraph.vertex(clique[j])
                        subgraph.add_edge(v1, v2)

        component_labels, hist = gt.label_components(subgraph)
        components_dict = {}
        components_dict = defaultdict(set)
        for v in graph_vertex_iter(clique, subgraph):
            c_index = component_labels[v]
            components_dict[c_index].add(int(v))

        logger.info(f"Found {len(components_dict)} components in the clique {clique}")
        for c_index, vs in components_dict.items():
            for v in vs:
                final_components[v] = haplotype_idx
            qnames = "\n".join([graph.vertex_properties['qname'][graph.vertex(v)] for v in vs])
            logger.info(f"Assigning \n{qnames}\nto the component group {haplotype_idx}")
            haplotype_idx += 1

    # logger.info(f"This is the final components: {final_components}")
    return final_components



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
        Dictionary of read objects, keyed by query name indices.
        The ncls_read_dict looks like:
                                        { qname1: [read1, read2],
                                          qname2: [read1, read2]  }

    ncls_qname_dict : dict
        Dictionary mapping query name indices to query names. (qname_indices -> qnames)

    ncls_qname_idx_dict : dict
        Dictionary mapping query names to query name indices. (qnames -> qname_indices)

    mean_read_length : float
        Average read length in the BAM file.
    edge_weight_cutoff : float, optional
        Minimum edge weight threshold for including an edge in the graph (default is 0.201).
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
    find_cliques_in_every_component : For identifying haplotypes in the constructed graph.
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
        return None, None, None, None

    logger.info("\n{}\n{}\n".format(ad_expanded.loc[~ad_expanded.iloc[:, -1].isna(), :][:10].to_string(index=False),
                                    alt_expanded.loc[~alt_expanded.iloc[:, -1].isna(), :][:10].to_string(index=False)))

    # After getting the table, we seek to organize the AD info by chromsome, pos, and alt into a nested dictionary
    # The dict structure is like:
    # nested_ad_dict[chrom][pos][alt] = ad
    # nested_ad_dict[chrom][pos]["DP"] = total_dp
    # Initialize the nested dictionary
    nested_ad_dict = {chrom: {} for chrom in ad_table["chrom"].unique()}
    column_width = alt_expanded.shape[1]
    inspected_overlaps = IntervalTree()

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
                    if np.isnan(bool_res):
                        same_hap_bools[n] = 0
                    elif bool_res:
                        same_hap_bools[n] = 1
                    else:
                        same_hap_bools[n] = -1

                    if read_weight is not None:
                        read_weight = read_weight if read_weight > 0 else 1
                        norm_weight = read_weight/(mean_read_length * 10)
                        if pair_weight is None:
                            pair_weight = norm_weight
                        else:
                            pair_weight += norm_weight
                    else:
                        norm_weight = None

                inspected_overlaps.addi(overlap_start, overlap_end)
                n += 1
            same_hap_bools = same_hap_bools[:n]
            pair_weight = pair_weight if pair_weight != 1 else 1 + 1e-4
            # if qname in vis_qnames and other_qname in vis_qnames:
            #     logger.info(f"Found the qname {qname} (qv is {int(qv)}) and other_qname {other_qname} (oqv is {int(oqv)}) overlap with each other with pair_weight {pair_weight}. The same_hap_bools are {same_hap_bools}")
            if custom_all_numba(same_hap_bools):
                # logger.info(f"The same_hap_bools are {same_hap_bools}")
                if pair_weight > edge_weight_cutoff:
                    e = g.add_edge(qv, oqv)
                    weight[e] = pair_weight
                weight_matrix[int(qv), int(oqv)] = pair_weight
                weight_matrix[int(oqv), int(qv)] = pair_weight
            elif any_false_numba(same_hap_bools):
                weight_matrix[int(qv), int(oqv)] = -1
                weight_matrix[int(oqv), int(qv)] = -1

    # Set the query name property for the graph
    g.vertex_properties["qname"] = qname_prop
    g.edge_properties["weight"] = weight

    logger.info(f"Now we finished building up the edges in the graph. There are currently {g.num_vertices()} vertices and {g.num_edges()} edges in the graph")
    return g, weight_matrix, qname_to_node, read_hap_vectors


@numba.njit(types.int32[:](types.int32[:,:], types.float32[:,:], types.int32[:,:]), fastmath=True)
def assemble_consensus(seq_arrays, qual_arrays, reads):
    # Every read in the reads can have different length
    # Find the start and end position of the consensus sequence
    start_pos = reads[:, 0].min()
    end_pos = reads[:, 1].max()

    # logger.info(", ".join([f"{r.reference_start}-{r.reference_end}" for r in reads]))

    # Initialize the consensus sequence and quality scores with zeros
    consensus_seq = np.ones(end_pos - start_pos, dtype=np.int32)
    consensus_qual = np.full(end_pos - start_pos, 0.1, dtype=np.float32)

    for i in prange(len(seq_arrays)):
        '''
        For every iteration,
        seq and qual should have the same length

        Across iterations,
        seq and qual are not bound with the same length
        '''
        seq = seq_arrays[i]
        qual = qual_arrays[i]
        read = reads[i]
        start = read[0]
        # assert seq.size == qual.size, f"Found the hap_vector {list(seq)} ({seq.size}bp) and qual_vector {list(qual)} ({qual.size}bp) have different lengths for read {read} at position {start}"
        # Filter out NaN values
        non_na_values = numba_sum(seq >= -8)

        nona_seq = np.empty(non_na_values, dtype=np.int32)
        nona_qual = np.empty(non_na_values, dtype=np.float32)

        for j in prange(non_na_values):
            nona_seq[j] = seq[j]
            nona_qual[j] = qual[j]

        seq = nona_seq
        qual = nona_qual

        # Calculate the relative start and end positions of the current sequence
        rel_start = start - start_pos
        rel_end = rel_start + non_na_values

        # Create a mask for positions where the current sequence has higher quality
        # assert rel_end - rel_start == len(qual), f"Found the relative start and end positions (determined by {len(seq)}) of the current sequence are not equal to the quality vector size {len(qual)} for read {read} at position {start}"
        # qual here is actually err prob, so smaller the better
        mask = qual <= consensus_qual[rel_start:rel_end]
        mask = mask & (qual < 0.1)

        # Update the consensus sequence and quality scores where the current sequence has higher quality
        consensus_seq[rel_start:rel_end][mask] = seq[mask]
        consensus_qual[rel_start:rel_end][mask] = qual[mask]

    return consensus_seq


@numba.njit
def find_row_index(array, value):
    for i in prange(array.shape[0]):
        if array[i] == value:
            return i
    return -1


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

def extract_continuous_regions_dict(reads):
    # Sort reads by reference start position
    reads = sorted(reads, key=lambda read: read.reference_start)

    continuous_regions = {}
    current_region_reads = []
    current_start = None
    current_end = None

    for read in reads:
        read_start = read.reference_start
        read_end = read.reference_end

        if current_start is None:
            # Initialize the first region
            current_start = read_start
            current_end = read_end
            current_region_reads.append(read)
        else:
            if read_start <= current_end:
                # Extend the current region
                current_end = max(current_end, read_end)
                current_region_reads.append(read)
            else:
                # Save the current region and start a new one
                continuous_regions[(current_start, current_end)] = current_region_reads
                current_start = read_start
                current_end = read_end
                current_region_reads = [read]

    # Append the last region
    if current_region_reads:
        continuous_regions[(current_start, current_end)] = current_region_reads

    return continuous_regions


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



def identify_misalignment_per_region(region,
                                    bam_ncls,
                                    intrin_bam_ncls,
                                    phasing_graph,
                                    qname_hap_info,
                                    hap_qname_info,
                                    weight_matrix,
                                    qname_to_node,
                                    lowqual_qnames,
                                    clique_sep_component_idx,
                                    var_df,
                                    viewed_regions = {},
                                    total_hapvectors = {},
                                    total_errvectors = {},
                                    total_genomic_haps = {},
                                    mean_read_length = 148,
                                    logger = logger ):
    bam_ncls, read_pair_dict, qname_dict, qname_idx_dict = bam_ncls
    '''
    read pair dict use internal qname_idx (created when executing function migrate bam to ncls) as keys instead of qnames
    qname_dict is a dictionary that maps qname_idx to qname
    '''

    # qname_idx_dict = {qname: qname_idx for qname_idx, qname in qname_dict.items()}
    overlap_reads_iter = lazy_get_overlapping_reads(bam_ncls,
                                                    read_pair_dict,
                                                    qname_dict,
                                                    region[0],
                                                    region[1],
                                                    region[2])

    vertices = defaultdict(list)
    vert_inds = set()
    for read in overlap_reads_iter:
        qname = read.query_name
        if qname in lowqual_qnames or qname not in qname_to_node:
            continue
        vert_idx = qname_to_node[qname]
        vertices[(vert_idx, qname)].append(read)
        vert_inds.add(vert_idx)

    subgraphs = group_by_dict_optimized(qname_hap_info, vertices)

    if len(subgraphs) == 0:
        logger.error(f"No haplotype clusters are found for region {region}. These are the qnames found overlapping this region: {vertices}. Skip this region.\n")
        return None, total_hapvectors, total_errvectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx

    # For each subgraph, we need to generate the consensus sequence for the reads
    # Each subgraph is a connected component in the phasing graph and it represents a haplotype
    region_haplotype_info = {}
    # overlap_haplotype_info = {}
    inspect_results = []
    for subgraph_id, (subgraph_vertices, qnames, read_pair_lists) in subgraphs.items():
        reads = []
        for read_pair_list in read_pair_lists:
            reads += read_pair_list
        # Sort the reads by start positions
        continuous_covered_regions = extract_continuous_regions_dict(reads)

        for span, sreads in continuous_covered_regions.items():
            if span[0] <= region[1] and span[1] >= region[2]:
                region_haplotype_info[span] = (sreads, set([qname_to_node[r.query_name] for r in sreads]), subgraph_id, qnames)
            else:
                inspect_results.append((subgraph_id, qnames, span, sreads, (min(span[1], region[2]) - max(span[0], region[1])) / (region[2] - region[1])))

    if len(region_haplotype_info) <= 2:
        logger.info(f"At region {region}, only found {len(region_haplotype_info)} haplotype clusters enwrapping the whole region. Try recover some regions.\n")
        inspect_results = sorted(inspect_results, key = lambda x: x[4], reverse = True)
        if len(inspect_results) == 0:
            return None, total_hapvectors, total_errvectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx
        elif inspect_results[0][4] < 0.8:
            return None, total_hapvectors, total_errvectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx
        else:
            recover_results = [t for t in inspect_results if t[4] >= 0.8]
            recover_span = max(recover_results, key = lambda x: x[2][0])[2][0], min(recover_results, key = lambda x: x[2][1])[2][1]
            overlapping_span = (max(recover_span[0], region[1]), min(recover_span[1], region[2]))
            for t in recover_results:
                subgraph_id, qnames, span, sreads, _ = t
                region_haplotype_info[span] = (sreads, set([qname_to_node[r.query_name] for r in sreads]), subgraph_id, qnames)
            if len(region_haplotype_info) <= 2:
                logger.info(f"At region {region}, only found {len(region_haplotype_info)} haplotype clusters enwrapping the whole region even after the recover trial.\n")
                return None, total_hapvectors, total_errvectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx
    else:
        overlapping_span = (region[1], region[2])
    region_str = f"{read.reference_name}:{overlapping_span[0]}-{overlapping_span[1]}"

    final_clusters = {}
    for span, (reads, vert_inds, haplotype_idx, qnames) in region_haplotype_info.items():
        hap_vectors = []
        err_vectors = []
        assert all([qname_hap_info[vid] == haplotype_idx for vid in vert_inds]), f"The vertices {vert_inds} are not in the same connected component."
        if len(reads) == 0:
            logger.warning(f"No reads are found for the haplotype across {span}. Skip this haplotype cluster.")
            continue

        read_spans = np.empty((len(reads), 2), dtype=np.int32)
        hap_vectors = np.full((len(reads), 500), -10, dtype=np.int32)
        err_vectors = np.full((len(reads), 500), -10, dtype=np.float32)
        for i, r in enumerate(reads):
            read_spans[i, 0] = r.reference_start
            read_spans[i, 1] = r.reference_end
            rid = get_read_id(r)
            if rid in total_hapvectors:
                hap_vector = total_hapvectors[rid]
            else:
                hap_vector = get_hapvector_from_cigar(r.cigartuples, r.query_sequence, logger = logger)
                total_hapvectors[rid] = hap_vector

            hap_vectors[i, :hap_vector.size] = hap_vector

            if rid in total_errvectors:
                err_vector = total_errvectors[rid]
            else:
                err_vector = get_errorvector_from_cigar(r, r.cigartuples, logger = logger)
                total_errvectors[rid] = err_vector

            err_vectors[i, :err_vector.size] = err_vector
        # logger.info(f"Haplotype across {span} contains {len(reads)} reads.")
        consensus_sequence = assemble_consensus(hap_vectors, err_vectors, read_spans)
        # logger.info(f"The consensus sequence for haplotype {haplotype_idx} across {span} composed of {len(reads)} and the qnames are {[r.query_name for r in reads]}. The consensus sequence is \n{consensus_sequence}\n")
        overlapping_con_seq = consensus_sequence[overlapping_span[0] - span[0]:overlapping_span[1] - span[0] + 1]
        final_clusters[haplotype_idx] = (overlapping_con_seq, reads, overlapping_span, qnames)

    # logger.info(f"The Final haplotype dict looks like this \n{final_clusters}\n")
    if len(final_clusters) <= 2:
        logger.info(f"Only {len(final_clusters)} haplotype clusters are found for region {region_str}. Do not need to choose 2 haplotypes, Skip this region.\n")
        return None, total_hapvectors, total_errvectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx

    record_2d_arr = record_haplotype_rank( final_clusters, mean_read_length = mean_read_length )
    record_df = pd.DataFrame(record_2d_arr, columns = ["start", "end", "total_depth", "hap_id", "hap_depth", "var_count", "indel_count"])
    record_df["chrom"] = region[0]

    logger.info(f"Found {len(final_clusters)} haplotype clusters for region {region_str}. The dataframe recording the haplotypes in this region looks like :\n{record_df.to_string(index=False)}\n")
    return record_df, total_hapvectors, total_errvectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx


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




def lp_solve_remained_haplotypes(total_record_df,
                                 logger = logger):
    # total_record_df = total_record_df.loc[np.logical_not(total_record_df["extreme_vard"]), :]
    hap_id_coefficients = total_record_df.groupby(["hap_id"])["coefficient"].mean().reset_index()
    hap_id_coefficients = hap_id_coefficients.merge(total_record_df.loc[:, ["hap_id", "ref_genome_similarities"]], how="left", on="hap_id")
    hap_id_coefficients.loc[:, "coefficient"] = hap_id_coefficients["coefficient"] + hap_id_coefficients["ref_genome_similarities"]
    hap_id_coefficients.drop_duplicates(inplace=True)

    logger.info(f"The hap_id_coefficients looks like \n{hap_id_coefficients.to_string(index=False)}\n")
    # rank_3_count = total_record_df[total_record_df["rank"] <= 3].shape[0]
    # rank_1_count = total_record_df[total_record_df["rank"] <= 1].shape[0]
    hap_no = hap_id_coefficients.shape[0]
    '''
    Note that here the hap_id_coefficient table is a dataframe with columns ["hap_id", "coefficient"]
    The row index is used to represent the hap_id in the constraint matrix below, but the the value does not necessarily equal to hap_id
    Therefore, we need two dicts to map from row_index to hap_id and from hap_id to row_index
    '''

    assert len(hap_id_coefficients["hap_id"].unique()) == hap_no, f"The hap_id_coefficients table has {hap_no} rows, but the hap_id column has {len(hap_id_coefficients['hap_id'].unique())} unique values."
    index_to_hapid = dict(zip(range(hap_no), hap_id_coefficients["hap_id"].to_numpy()))
    hapid_to_index = dict(zip(hap_id_coefficients["hap_id"].to_numpy(), range(hap_no)))

    logger.info(f"Now the index_to_hapid dict looks like \n{index_to_hapid}\n")
    logger.info(f"Now the hapid_to_index dict looks like \n{hapid_to_index}\n")

    by_region = total_record_df.groupby(["chrom", "start", "end"])
    inspect_region_no = by_region.ngroups
    # Initialize the HiGHS model
    highs = Highs()

    # Set the output_flag to False to suppress output
    highs.setOptionValue('output_flag', False)
    # Create a constraint matrix
    logger.info(f"The constraint matrix will have a shape of {(inspect_region_no, hap_no)}")

    # Then we need to fill the constraint matrix A properly to correctly specify start indices, row_indices, and values in the future function calls.
    lower_bounds = np.zeros(hap_no, dtype=np.int32)
    upper_bounds = np.ones(hap_no, dtype=np.int32)

    # Set the output_flag to False to suppress output
    highs.setOptionValue('output_flag', True)

    col_costs = - hap_id_coefficients["coefficient"].to_numpy().astype(np.double)

    # Add Vars
    highs.addVars(hap_no, lower_bounds, upper_bounds)

    # Change column costs
    highs.changeColsCost(hap_no, np.arange(hap_no, dtype=np.int32), col_costs)

    # Set variable types to integer
    integrality = np.ones(hap_no, dtype=np.int32)  # 1 for integer variables
    highs.changeColsIntegrality(hap_no, np.arange(hap_no, dtype=np.int32), integrality)

    # Adding rows
    for name, group in by_region:
        included_hapids = group["hap_id"].unique()
        if included_hapids.size <= 2:
            continue
        rank_2_count = group[group["rank"] <= 2].shape[0]
        # Column index extraction
        # hapid_indices = [hap_id_coefficients[hap_id_coefficients["hap_id"] == hapid].index[0] for hapid in included_hapids]
        # hapid_indices = find_indices(hap_id_coefficients["hap_id"].to_numpy(), included_hapids)
        hapid_indices = [hapid_to_index[hapid] for hapid in included_hapids]
        lower_bound = included_hapids.size - rank_2_count
        upper_bound = max(included_hapids.size - 2, lower_bound)
        status = highs.addRow(0,
                              lower_bound,
                              included_hapids.size,
                              np.array(hapid_indices, dtype=np.int32),
                              np.ones(included_hapids.size, dtype=np.double))
        logger.info(f"The addrow status is {status}, the group included hap_ids are {included_hapids}, the corresponding hap_id indices are {hapid_indices} the lower bound is {lower_bound}, the upper bound is {upper_bound}.")

    # Run solver
    highs.run()
    solution = highs.getSolution()
    info = highs.getInfo()
    num_var = highs.getNumCol()
    model_status = highs.getModelStatus()
    basis = highs.getBasis()

    # Access and print the constraint matrix
    lp = highs.getLp()
    logger.info(f"column coefficients: {lp.col_cost_}")
    logger.info(f"Model status = {highs.modelStatusToString(model_status)}")
    logger.debug(f"Optimal objective = {info.objective_function_value}")
    logger.info(f"Iteration count = {info.simplex_iteration_count}")
    logger.info(f"Primal solution status = {highs.solutionStatusToString(info.primal_solution_status)}")

    # Get solution
    solution = np.array(highs.getSolution().col_value)
    # print(f"Solution is {solution}")
    select_hap_indices = np.where(solution == 0)[0]
    drop_hap_indices = np.where(solution == 1)[0]
    select_hap_ids = hap_id_coefficients.iloc[select_hap_indices, 0]
    drop_hap_ids = hap_id_coefficients.iloc[drop_hap_indices, 0]
    return set(select_hap_ids), set(drop_hap_ids), highs.modelStatusToString(model_status)


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


def inspect_by_haplotypes(input_bam,
                          hap_qname_info,
                          qname_hap_info,
                          bam_ncls,
                          intrin_bam_ncls,
                          phased_graph,
                          weight_matrix,
                          qname_to_node,
                          total_lowqual_qnames,
                          total_hap_vectors,
                          total_err_vectors,
                          total_genomic_haps,
                          var_df,
                          compare_haplotype_meta_tab = "",
                          mean_read_length = 150,
                          logger = logger):
    ncls_dict, read_dict, qname_dict, qname_idx_dict = bam_ncls
    record_dfs = []
    clique_sep_component_idx = 0
    ref_genome_similarities = defaultdict(dict)
    hid_extreme_vard = defaultdict(bool)
    hid_var_count = defaultdict(int)
    scatter_hid_dict = defaultdict(bool)
    hsid_hid_dict = defaultdict(int)
    hsid_seq_dict = {}
    phased_graph_qname_vprop = phased_graph.vertex_properties['qname']
    logger.info("All the haplotype IDs are :\n{}\n".format(list(hap_qname_info.keys())))

    for hid, qnames in hap_qname_info.items():
        logger.info(f"The haplotype {hid} contains {len(qnames)} read pairs in total.")

        qname_indices = [qname_idx_dict[qname] for qname in qnames]
        reads = [read for qname_idx in qname_indices for read in read_dict[qname_idx]]
        conregion_dict = extract_continuous_regions_dict(reads)

        if len(qnames) < 2:
            scatter_hid_dict[hid] = True

        high_vard = False
        mid_vard_dict = {}
        # Generate consensus sequence for each continuously covered region
        for span, reads in conregion_dict.items():
            # span end is exclusive
            logger.info(f"The haplotype {hid} contains {len(reads)} reads in the continuous region {span}.")
            chrom = reads[0].reference_name
            srids = set([])
            read_spans = np.empty((len(reads), 2), dtype=np.int32)
            hap_vectors = np.full((len(reads), 500), -10, dtype = np.int32)
            err_vectors = np.full((len(reads), 500), -10, dtype = np.float32)
            for i, r in enumerate(reads):
                read_spans[i, 0] = r.reference_start
                read_spans[i, 1] = r.reference_end
                rid = get_read_id(r)
                srids.add(rid)
                if rid in total_hap_vectors:
                    hap_vector = total_hap_vectors[rid]
                else:
                    hap_vector = get_hapvector_from_cigar(r.cigartuples, r.query_sequence, logger = logger)
                    total_hap_vectors[rid] = hap_vector

                hap_vectors[i, :hap_vector.size] = hap_vector

                if rid in total_err_vectors:
                    err_vector = total_err_vectors[rid]
                else:
                    err_vector = get_errorvector_from_cigar(r, r.cigartuples, logger = logger)
                    total_err_vectors[rid] = err_vector

                err_vectors[i, :err_vector.size] = err_vector

            consensus_sequence = assemble_consensus(hap_vectors, err_vectors, read_spans)
            extreme_vard = judge_misalignment_by_extreme_vardensity(consensus_sequence)
            var_count = count_var(consensus_sequence)
            hid_var_count[hid] += var_count
            logger.info(f"For haplotype {hid} at continuous region {span}, the consensus sequence is {consensus_sequence.tolist()}. The extreme variant density is {extreme_vard}.")
            hid_extreme_vard[hid] = extreme_vard or hid_extreme_vard[hid]
            # We need to extract the overlapping genomic haplotype vectors
            # Get the genomic haplotype vectors for the overlapping region
            max_similarity = 0
            for hread in lazy_get_overlapping_reads(*intrin_bam_ncls[:-1], chrom, span[0], span[1]):
                hread_start = hread.reference_start
                hread_end = hread.reference_end
                hread_qname = hread.query_name
                overlap_span = (max(hread_start, span[0]), min(hread_end, span[1]))
                # logger.info(f"The genomic sequence spanning {hread_start} to {hread_end}. And it overlaps with the haplotype at {overlap_span}")
                if overlap_span[1] - overlap_span[0] < 140:
                    continue
                hregion_str = f"{chrom}:{hread_start}-{hread_end}"
                hread_id = get_read_id(hread) + ":" + hregion_str
                # hread_cigars = tuple(hread.cigartuples)
                if hread_id in total_genomic_haps:
                    hread_hap_vector = total_genomic_haps[hread_id]
                else:
                    hread_hap_vector = get_hapvector_from_cigar(hread.cigartuples, logger = logger)
                    total_genomic_haps[hread_id] = hread_hap_vector

                try:
                    interval_genomic_hap = hread_hap_vector[overlap_span[0] - hread_start:overlap_span[1] - hread_start]
                except IndexError:
                    continue
                interval_con_seq = consensus_sequence[overlap_span[0] - span[0]:overlap_span[1] - span[0]]
                # logger.info(f"Inspect the conesensus sequence (length {len(consensus_sequence)})\n{consensus_sequence.tolist()}\nwith the genomic haplotype (length {len(interval_genomic_hap)})\n{interval_genomic_hap.tolist()}\n")
                varcount, alt_varcount = ref_genome_similarity(interval_con_seq, interval_genomic_hap)
                if hread_qname in ref_genome_similarities[hid]:
                    ref_genome_similarities[hid][hread_qname].append((varcount, alt_varcount))
                else:
                    ref_genome_similarities[hid][hread_qname] = [(varcount, alt_varcount)]

    ref_genome_str = '\n'.join([f"{k}: {v}" for k,v in ref_genome_similarities.items()])
    logger.info(f"ref genome similarities: {ref_genome_str}")
    logger.info(f"extreme variant density haplotypes: {hid_extreme_vard}")
    logger.info(f"scatter haplotypes: {scatter_hid_dict}")
    final_ref_seq_similarities = {}
    ref_seq_delta = {}
    for hid, gdict in ref_genome_similarities.items():
        max_similarity = 0
        for genomic_seq_qname, pairs in gdict.items():
            total_varcount = np.sum([t[0] for t in pairs])
            total_alt_varcount = np.sum([t[1] for t in pairs])
            similarity = total_varcount / total_alt_varcount if total_alt_varcount > 0 else total_varcount
            delta = total_varcount - total_alt_varcount
            # logger.info(f"For haplotype {hid}, comparing to mapping against the genomic sequence {genomic_seq_qname}. Variant count ratio is {similarity} and variant count delta is {delta}.")
            similarity_score = delta + similarity
            max_similarity = max(max_similarity, similarity_score)

        final_ref_seq_similarities[hid] = max_similarity if max_similarity > 0 else 0

    ref_genome_similarities = final_ref_seq_similarities
    logger.info(f"Final ref sequence similarities for all the haplotypes are {final_ref_seq_similarities}")

    sweep_region_bed = input_bam.replace(".bam", ".sweep.bed")
    sweep_regions = sweep_region_inspection(input_bam, sweep_region_bed, logger = logger)
    logger.info(f"Now the sweep regions are saved to {sweep_region_bed}.")
    for interval in sweep_regions:
        region = (interval.chrom, interval.start, interval.end)
        logger.info(f"Inspecting the region {region}.")
        record_df, total_hap_vectors, total_err_vectors, total_genomic_haps, qname_hap_info, clique_sep_component_idx = identify_misalignment_per_region(region,
                                                                                                                                                        bam_ncls,
                                                                                                                                                        intrin_bam_ncls,
                                                                                                                                                        phased_graph,
                                                                                                                                                        qname_hap_info,
                                                                                                                                                        hap_qname_info,
                                                                                                                                                        weight_matrix,
                                                                                                                                                        qname_to_node,
                                                                                                                                                        total_lowqual_qnames,
                                                                                                                                                        clique_sep_component_idx,
                                                                                                                                                        var_df,
                                                                                                                                                        total_hapvectors = total_hap_vectors,
                                                                                                                                                        total_errvectors = total_err_vectors,
                                                                                                                                                        total_genomic_haps = total_genomic_haps,
                                                                                                                                                        mean_read_length = mean_read_length,
                                                                                                                                                        logger = logger )
        if record_df is not None:
            record_dfs.append(record_df)

    failed_lp = False
    if len(record_dfs) > 0:
        total_record_df = pd.concat(record_dfs)
        record_hapids = set(total_record_df["hap_id"].unique())
        logger.info(f"The haplotypes that have been inspected are {record_hapids}.")
        # Is there any other haplotype that has not been inspected?
        remained_hapids = set(hap_qname_info.keys()) - record_hapids
        if len(remained_hapids) > 0:
            remained_hapid_dict = { k:v for k,v in hap_qname_info.items() if k in remained_hapids }
            logger.info(f"The haplotypes that have not been inspected are {remained_hapids}. Their qnames are: \n{remained_hapid_dict}\n")

        total_record_df["extreme_vard"] = total_record_df["hap_id"].map(hid_extreme_vard)
        total_record_df["scatter_hap"] = total_record_df["hap_id"].map(scatter_hid_dict).fillna(False)
        total_record_df["hap_var_count"] = total_record_df["hap_id"].map(hid_var_count)
        total_record_df["ref_genome_similarities"] = total_record_df["hap_id"].map(ref_genome_similarities).fillna(0)
        # total_record_df.loc[:, "coefficient"] = total_record_df["coefficient"] * 100 + total_record_df.loc[:, "var_count"]
        total_record_df.to_csv(compare_haplotype_meta_tab.replace(".tsv", ".raw.tsv"), sep = "\t", index = False)
        logger.info(f"Successfully saved the raw haplotype comparison meta table to {compare_haplotype_meta_tab.replace('.tsv', '.raw.tsv')}. And it looks like \n{total_record_df[:10].to_string(index=False)}\n")
        total_record_df = total_record_df.loc[np.logical_not(total_record_df["extreme_vard"]) & \
                                              np.logical_not(total_record_df["scatter_hap"]) & \
                                              (total_record_df["total_depth"] >= 5) & \
                                              (total_record_df["ref_genome_similarities"] <= 5), :]
        # total_record_df.loc[:, "coefficient"] = np.clip(total_record_df["coefficient"] + total_record_df["ref_genome_similarities"], 10e-3, None)
        if total_record_df.shape[0] == 0:
            failed_lp = True

        if not failed_lp:
            total_record_df.loc[total_record_df["hap_var_count"] == 0, "coefficient"] = -1
            total_record_df = total_record_df.groupby(["chrom", "start", "end"]).filter(lambda x: len(x) > 5)

            if total_record_df.shape[0] == 0:
                failed_lp = True
    else:
        failed_lp = True

    if not failed_lp:
        by_region = total_record_df.groupby(["chrom", "start", "end"], group_keys=False)
        total_record_df = by_region.apply(calulate_coefficient_per_group, logger = logger).reset_index(drop=True)
        total_record_df.to_csv(compare_haplotype_meta_tab, sep = "\t", index = False)
        logger.info(f"Successfully saved the haplotype comparison meta table to {compare_haplotype_meta_tab}. And it looks like \n{total_record_df[:10].to_string(index=False)}\n")
        correct_map_hids, mismap_hids, model_status = lp_solve_remained_haplotypes(total_record_df, logger = logger)
        if model_status == "Infeasible":
            failed_lp = True

    if not failed_lp:
        mismap_hids.update(set([hid for hid in hid_extreme_vard if hid_extreme_vard[hid]]))
        mismap_hids.update(set([hid for hid in ref_genome_similarities if ref_genome_similarities[hid] > 5]))
        mismap_hids.update(set([hid for hid in scatter_hid_dict if scatter_hid_dict[hid] and hid_var_count[hid] > 1]))
        correct_map_hids = correct_map_hids - mismap_hids

        mismap_qnames = set([qname for hid in mismap_hids for qname in hap_qname_info[hid]])
        correct_map_qnames = set([qname for hid in correct_map_hids for qname in hap_qname_info[hid]])
        logger.info(f"Identify {len(mismap_hids)} mismapped haplotypes, {len(mismap_qnames)} mismapped qnames by solving the BINARY INTEGER LINEAR PROGRAMMING.\n")
        logger.info(f"Here is the qnames seemed mismapped: \n{mismap_qnames}\nAnd the mismap haplotype IDs: \n{mismap_hids}\n")
        logger.info(f"Here is the qnames seemed correctly mapped: \n{correct_map_qnames}\nAnd the correct map haplotype IDs: \n{correct_map_hids}\n")
        return correct_map_qnames, mismap_qnames

    if failed_lp:
        mismap_hids = set([hid for hid in hid_extreme_vard if hid_extreme_vard[hid]])
        mismap_hids.update(set([hid for hid in scatter_hid_dict if scatter_hid_dict[hid] and hid_var_count[hid] > 1]))
        mismap_hids.update(set([hid for hid in ref_genome_similarities if ref_genome_similarities[hid] > 5]))
        mismap_qnames = set([qname for hid in mismap_hids for qname in hap_qname_info[hid]])
        total_qnames = set([qname for qnames in hap_qname_info.values() for qname in qnames])
        correct_map_qnames = total_qnames - mismap_qnames
        return correct_map_qnames, mismap_qnames



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
    total_cliques = find_cliques_in_every_component(phased_graph,
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




