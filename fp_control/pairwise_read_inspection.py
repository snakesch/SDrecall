import logging
import numba
import numpy as np
from numba import types

from numba_operators import numba_diff_indices, \
                            numba_sum, \
                            numba_not, \
                            numba_and, \
                            numba_slicing, \
                            numba_indexing_int32, \
                            numba_indexing_int8

logger = logging.getLogger('SDrecall')


@numba.njit(types.bool_[:](types.int32[:]), fastmath=True)
def get_indel_bools(seq_arr):
    return (seq_arr > 1) | (seq_arr == -6)

@numba.njit(types.int32(types.bool_[:]), fastmath=True)
def count_continuous_blocks(arr):
    """
    This function counts the number of blocks of True values in a 1D boolean array.
    """
    if arr.size == 0:
        return 0

    # Pad one False to the beginning and end of arr
    extended_arr = np.empty(arr.size + 2, dtype=np.bool_)
    extended_arr[1:-1] = arr

    shifted_sequence = numba_not(extended_arr[1:])
    block_bools = numba_and(extended_arr[:-1], shifted_sequence)

    return numba_sum(block_bools)


@numba.njit(types.int32(types.int16[:]), fastmath=True)
def count_snv(array):
    snv_bools = array == -4
    return count_continuous_blocks(snv_bools)


@numba.njit(types.int32(types.int16[:]), fastmath=True)
def count_continuous_indel_blocks(array):
    """
    This function counts the number of blocks of -6 and positive integers in a numpy array.
    """
    is_var = (array == -6) | (array > 1)

    return count_continuous_blocks(is_var)


@numba.njit(types.int32(types.int16[:]), fastmath=True)
def count_var(array):
    return count_snv(array) + count_continuous_indel_blocks(array)


@numba.njit(types.bool_(types.int8[:], types.int8[:], types.int8), fastmath=True)
def compare_sequences(read_seq, other_seq, except_char):
    """
    Compares two sequences, return True if they are equivalent. Ambiguous bases ("N") as ignored as except_char.
    """
    if read_seq.shape != other_seq.shape:
        return False

    for i in range(read_seq.size):
        if read_seq[i] != other_seq[i]:
            if read_seq[i] != except_char and other_seq[i] != except_char:
                return False

    return True


def get_read_id(read) -> str:
    """
    Generate a unique identifier for each read in a pair, using also alignment information as BAM flags
    """
    return f"{read.query_name}:{read.flag}"


@numba.njit(types.int32[:](types.int32[:], types.int32), fastmath=True)
def prepare_ref_query_idx_map(qseq_ref_pos_arr, read_ref_start):
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
    numba_arr_2d = np.full((qseq_ref_pos_arr.max() - read_ref_start + 1), -1, dtype=np.int32)
    for qi in range(qseq_ref_pos_arr.size):
        ri = qseq_ref_pos_arr[qi]
        if ri >= 0:
            numba_arr_2d[ri - read_ref_start] = qi
    return numba_arr_2d




@numba.njit(types.Tuple((types.int8[:], types.int8[:], types.int32[:]))(types.int32, types.int32, types.int32, types.int32[:], types.int32[:], types.int8[:], types.int8[:]), fastmath=True)
def get_interval_seq_qual(read_start, 
                          interval_start, 
                          interval_end, 
                          ref_positions,
                          qseq_ref_positions,
                          query_sequence_encoded,
                          query_sequence_qualities):
    interval_start_qidx = ref_positions[interval_start - read_start]
    if interval_start_qidx > -1:
        # logger.debug(f"Found the interval start {interval_start} is not in the reference positions: \n{ref_positions} of read {read.query_name}. Might locate in the middle of a deletion event.")
        while ref_positions[interval_start - read_start] < 0 and interval_start <= interval_end:
            interval_start += 1
        if interval_start > interval_end:
            # logger.debug(f"The whole specified interval (size {size}) is in a deletion event for read {read.query_name}. Now the returned seq is: {return_seq}\n")
            return np.array([np.int8(x) for x in range(0)], dtype=np.int8), np.array([np.int8(x) for x in range(0)], dtype=np.int8), np.array([np.int32(x) for x in range(0)], dtype=np.int32)
        interval_start_qidx = ref_positions[interval_start - read_start]

    interval_end_qidx = ref_positions[interval_end - read_start]
    if not interval_end_qidx:
        # logger.debug(f"Found the interval end {interval_end} is not in the reference positions: \n{ref_positions} of read {read.query_name}. Might locate in the middle of a deletion event.")
        while ref_positions[interval_end - read_start] < 0:
            interval_end -= 1
        interval_end_qidx = ref_positions[interval_end - read_start] + 1
    else:
        interval_end_qidx += 1

    interval_read_seq = query_sequence_encoded[interval_start_qidx:interval_end_qidx]
    interval_qual_seq = query_sequence_qualities[interval_start_qidx:interval_end_qidx]
    qidx_ridx_arr = qseq_ref_positions[interval_start_qidx:interval_end_qidx]

    # interval_read_seq = identify_discontinuity(interval_ref_pos, interval_read_seq)
    # return_seq = interval_read_seq
    # if len(return_seq) < size:
        # logger.debug(f"There are deletion event inside the read. The final return_seq should cover {size}bp, but it turns out to be {len(return_seq)}bp. The interval start is {interval_start}, interval end is {interval_end}. Take a look:\n{return_seq}\nThe preseq is {preseq}\nThe reference_positions are {ref_positions}\n")
    # elif len(return_seq) > size:
        # logger.debug(f"There are insertion event inside the read. The final return_seq should cover {size}bp, but it turns out to be {len(return_seq)}bp. The interval start is {interval_start}, interval end is {interval_end}. Take a look:\n{return_seq}\nThe preseq is {preseq}\nThe reference_positions are {ref_positions}\n")
    return interval_read_seq, interval_qual_seq, qidx_ridx_arr



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



@numba.njit(types.bool_(types.int32, types.int8[:], types.int8[:], types.int32[:], types.int32, types.DictType(types.int8, types.int16)), fastmath=True)
def seq_err_det_stacked_bases(target_read_start,
                              target_read_qseq_encoded,
                              target_read_qseq_qualities,
                              target_read_ref_positions,
                              position0,
                              nested_dict):
    ridx = position0 - target_read_start
    target_read_qidx = numba_indexing_int32(target_read_ref_positions, ridx)
    target_base_encoded = numba_indexing_int8(target_read_qseq_encoded, target_read_qidx)
    base_qual = numba_indexing_int8(target_read_qseq_qualities, target_read_qidx)
    
    if base_qual >= 10:
        return False
    
    ad = nested_dict.get(target_base_encoded, 0)
    dp = nested_dict.get(np.int8(5), 0)

    if int(dp) == 0:
        return False

    af = int(ad) / int(dp)

    if (af <= 0.02 or (int(ad) == 1 and int(dp) >= 10)) and base_qual < 10:
        return True
    else:
        return False



def tolerate_mismatches_two_seq(read1_package, read2_package,
                                abs_diff_indices,
                                nested_ad_dict):
    '''
    This function is used to compare two sequences and identify the mismatch positions
    And to decide whether we are safe to determine the mismatches are originated from sequencing error or not.

    Use seq_err_det_stacked_bases() to determine if the mismatches can be explained by sequencing artifacts, details can be found in the docstring of the function
    '''
    read1_start, read1_qseq_encoded, read1_qseq_qualities, read1_ref_positions = read1_package
    read2_start, read2_qseq_encoded, read2_qseq_qualities, read2_ref_positions = read2_package

    tolerate_mismatches = []
    for diff_ind in abs_diff_indices:
        seq_err1 = seq_err_det_stacked_bases(read1_start,
                                             read1_qseq_encoded,
                                             read1_qseq_qualities,
                                             read1_ref_positions,
                                             diff_ind,
                                             nested_ad_dict[diff_ind])
        # if added_stack_base_dict:
        #     stack_base_dict.update(added_stack_base_dict)

        seq_err2 = seq_err_det_stacked_bases(read2_start,
                                             read2_qseq_encoded,
                                             read2_qseq_qualities,
                                             read2_ref_positions,
                                             diff_ind,
                                             nested_ad_dict[diff_ind])
        tolerate_mismatches.append(seq_err1 or seq_err2)

    if all(tolerate_mismatches):
        return True, len(tolerate_mismatches)
    else:
        return False, 0


def extract_read_qseqs(read, read_ref_pos_dict):
    read_id = get_read_id(read)
    read_start = read.reference_start

    if read_id in read_ref_pos_dict:
        ref_positions, qseq_ref_positions, query_sequence_encoded, query_sequence_qualities = read_ref_pos_dict[read_id]
    else:
        qseq_ref_positions = np.array([ idx if idx is not None else -1 for idx in read.get_reference_positions(full_length=True) ], dtype=np.int32)  # An numpy array mapping query sequence index to reference genome index
        ref_positions = prepare_ref_query_idx_map(qseq_ref_positions, read_start) # An numpy array mapping reference genome index to query sequence index
        base_dict = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
        query_sequence_encoded = np.array([ base_dict[base] for base in read.query_sequence ], dtype=np.int8)
        query_sequence_qualities = np.array(read.query_qualities, dtype=np.int8)
        read_ref_pos_dict[read_id] = (ref_positions, qseq_ref_positions, query_sequence_encoded, query_sequence_qualities)

    return ref_positions, qseq_ref_positions, query_sequence_encoded, query_sequence_qualities, read_ref_pos_dict


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

    ref_positions, qseq_ref_positions, query_sequence_encoded, query_sequence_qualities, read_ref_pos_dict = extract_read_qseqs(read, read_ref_pos_dict)
    other_ref_positions, other_qseq_ref_positions, other_query_sequence_encoded, other_query_sequence_qualities, read_ref_pos_dict = extract_read_qseqs(other_read, read_ref_pos_dict)

    # Below we extract the haplotype vector for both reads, which are an array of int values indicating the alignment status of each base within the reference genome span aligned by the read
    # Detailed explanation can be found in the docstring of the function get_hapvector_from_cigar
    if read_id in read_hap_vectors:
        read_hap_vector = read_hap_vectors[read_id]
    else:
        cigar_arr = np.array(read.cigartuples, dtype=np.int32)
        read_hap_vector = get_hapvector_from_cigar(cigar_arr, query_sequence_encoded)
        read_hap_vectors[read_id] = read_hap_vector

    # We need a unique identifier for every read, which is a string composed by query name and the alignment flag.
    other_read_id = get_read_id(other_read)
    other_start = other_read.reference_start

    if other_read_id in read_hap_vectors:
        other_read_hap_vector = read_hap_vectors[other_read_id]
    else:
        other_cigar_arr = np.array(other_read.cigartuples, dtype=np.int32)
        other_read_hap_vector = get_hapvector_from_cigar(other_cigar_arr, other_query_sequence_encoded)
        read_hap_vectors[other_read_id] = other_read_hap_vector

    # Extract the haplotype vector within the overlapping region
    interval_hap_vector = numba_slicing(read_hap_vector, overlap_start, overlap_end, start)
    interval_other_hap_vector = numba_slicing(other_read_hap_vector, overlap_start, overlap_end, other_start)

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
    read_seq, read_qualities, qr_idx_arr = get_interval_seq_qual(np.int32(read.reference_start), 
                                                                 np.int32(overlap_start), 
                                                                 np.int32(overlap_end - 1), 
                                                                 ref_positions, 
                                                                 qseq_ref_positions, 
                                                                 query_sequence_encoded,
                                                                 query_sequence_qualities)
    other_seq, other_qualities, other_qr_idx_arr = get_interval_seq_qual(np.int32(other_read.reference_start), 
                                                                 np.int32(overlap_start), 
                                                                 np.int32(overlap_end - 1), 
                                                                 other_ref_positions, 
                                                                 other_qseq_ref_positions, 
                                                                 other_query_sequence_encoded,
                                                                 other_query_sequence_qualities)

    # We need to ignore N during the comparison for read_seq and other_seq
    total_match = compare_sequences(read_seq, other_seq, np.int8(4))

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
        # Use numba_diff_indices to find the location indices where two sequences differ
        q_diff_indices = numba_diff_indices(read_seq, other_seq)
        # Map the q_diff_indices to the reference positions
        r_diff_indices = qr_idx_arr[q_diff_indices]
        # If -1 in r_diff_indices, it means the two reads are not aligned properly and we cannot make a conclusion on whether they belong to the same haplotype
        if -1 in r_diff_indices:
            return False, read_ref_pos_dict, read_hap_vectors, None
        read_package = (np.int32(read.reference_start), query_sequence_encoded, query_sequence_qualities, ref_positions)
        other_read_package = (np.int32(other_read.reference_start), other_query_sequence_encoded, other_query_sequence_qualities, other_ref_positions)
        # See if the mismatches can be explained by sequencing artifacts
        tolerate, tolerated_count = tolerate_mismatches_two_seq( read_package,
                                                                 other_read_package,
                                                                 r_diff_indices,
                                                                 nested_ad_dict )

        if tolerate:
            edge_weight = edge_weight - tolerated_count * 20
            if edge_weight >= mean_read_length - 50:
                return True, read_ref_pos_dict, read_hap_vectors, edge_weight
            else:
                return np.nan, read_ref_pos_dict, read_hap_vectors, None
        else:
            return False, read_ref_pos_dict, read_hap_vectors, None