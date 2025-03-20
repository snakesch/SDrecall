import logging
import numba
import numpy as np
from numba import types

from src.log import logger
from src.suppress_warning import *
from fp_control.numba_operators import numba_diff_indices, \
                                        numba_sum, \
                                        numba_not, \
                                        numba_and, \
                                        numba_slicing, \
                                        numba_indexing_int32, \
                                        numba_indexing_int8, \
                                        numba_compare, \
                                        numba_bool_indexing, \
                                        numba_contain


@numba.njit(types.bool_[:](types.int16[:]), fastmath=True)
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
    extended_arr[0] = False  # Explicitly set beginning padding
    extended_arr[-1] = False  # Explicitly set end padding
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
    The output array will be: [0, 1, 4, 5, 6, 7]

    ri stands for reference index, qi stands for query index
    '''
    numba_arr_2d = np.full((qseq_ref_pos_arr.max() - read_ref_start + 1), -1, dtype=np.int32)
    for qi in range(qseq_ref_pos_arr.size):
        ri = qseq_ref_pos_arr[qi]
        if ri >= 0:
            numba_arr_2d[ri - read_ref_start] = qi
    return numba_arr_2d



@numba.njit(types.Tuple((types.int8[:], types.int32[:]))(types.int32, types.int32, types.int32, types.int32[:], types.int32[:], types.int8[:]), fastmath=True)
def get_interval_seq(read_start, 
                     interval_start, 
                     interval_end, 
                     ref_positions,
                     qseq_ref_positions,
                     query_sequence_encoded):
    '''
    This function is used to get the interval sequence and the corresponding reference positions
    The interval_start and interval_end are both 0-indexed, interval_end is exclusive
    '''

    interval_start_qidx = ref_positions[interval_start - read_start]
    while interval_start_qidx < 0 and interval_start <= interval_end:
        interval_start += 1
        interval_start_qidx = ref_positions[interval_start - read_start]
    if interval_start > interval_end:
        # logger.warning(f"The whole specified interval (size {size}) is in a deletion event for read {read_id}. Now the returned seq is: {return_seq}\n")
        return np.array([np.int8(x) for x in range(0)], dtype=np.int8), np.array([np.int32(x) for x in range(0)], dtype=np.int32)

    interval_end_qidx = ref_positions[interval_end - 1 - read_start]
    while interval_end_qidx < 0:
        interval_end -= 1
        interval_end_qidx = ref_positions[interval_end - 1 - read_start]

    interval_read_seq = query_sequence_encoded[interval_start_qidx:interval_end_qidx + 1]
    qidx_ridx_arr = qseq_ref_positions[interval_start_qidx:interval_end_qidx + 1]

    return interval_read_seq, qidx_ridx_arr




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
    assert ref_length > 0, f"The reference length is 0 for the cigar string: \n{cigar_tuples}\n"

    # Create a haplotype vector of zeros
    hapvector = np.empty(ref_length, dtype=np.int16)

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
                if len(base) > 1:
                    # There are more than 1 consecutive bases with mismatches in the cigar string: {cigar_tuples}
                    pass
                elif base == "N":
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



def tolerate_mismatches_two_seq(read1_package, 
                                read2_package,
                                abs_diff_indices,
                                nested_ad_dict,
                                empty_dict):
    '''
    This function is used to compare two sequences and identify the mismatch positions
    And to decide whether we are safe to determine the mismatches are originated from sequencing error or not.

    Use seq_err_det_stacked_bases() to determine if the mismatches can be explained by sequencing artifacts, details can be found in the docstring of the function

    the value at diff_ind in nested_ad_dict should be a dictionary with:
    {0: ad, 1: ad, 2: ad, 3: ad, 4: ad, 5: dp} (ad is allele depth, dp is depth)
    also:
    {A: 0, T: 1, C: 2, G: 3, N: 4}
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
                                             nested_ad_dict.get(diff_ind, empty_dict))

        seq_err2 = seq_err_det_stacked_bases(read2_start,
                                             read2_qseq_encoded,
                                             read2_qseq_qualities,
                                             read2_ref_positions,
                                             diff_ind,
                                             nested_ad_dict.get(diff_ind, empty_dict))
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
        query_sequence_qualities = np.array(read.query_qualities, dtype=np.int8) if read.query_qualities is not None else np.array([np.int8(x) for x in range(0)], dtype=np.int8)
        read_ref_pos_dict[read_id] = (ref_positions, qseq_ref_positions, query_sequence_encoded, query_sequence_qualities)

    return ref_positions, qseq_ref_positions, query_sequence_encoded, query_sequence_qualities, read_ref_pos_dict



def determine_same_haplotype(read, other_read,
                             overlap_start, overlap_end,
                             score_arr,
                             read_hap_vectors = {},
                             nested_ad_dict = {},
                             read_ref_pos_dict = {},
                             mean_read_length = 148,
                             empty_dict = {},
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
    start = np.int32(read.reference_start)

    ref_positions, qseq_ref_positions, query_sequence_encoded, query_sequence_qualities, read_ref_pos_dict = extract_read_qseqs(read, read_ref_pos_dict)
    other_ref_positions, other_qseq_ref_positions, other_query_sequence_encoded, other_query_sequence_qualities, read_ref_pos_dict = extract_read_qseqs(other_read, read_ref_pos_dict)

    # qseq_ref_positions is a numpy array mapping query sequence index to reference genome index, it looks like: [1000, 1001, 1002, 1003, -1, -1, 1004, 1005, 1007, 1008] for a CIGAR str 4=2I2=1D2=
    # qseq_seq_encoded is just a numpy array of integers. A:0, C:1, G:2, T:3, N:4
    # ref_positions is a numpy array mapping reference genome index back to query sequence index, it looks like: [0, 1, 2, 3, 6, 7, -1, 8, 9] for a CIGAR str 4=2I2=1D2=, all integers are the remains of index on ref genome subtracted by start ref pos.

    if read_id in read_hap_vectors:
        read_hap_vector = read_hap_vectors[read_id]
    else:
        cigar_arr = np.array(read.cigartuples, dtype=np.int32)
        read_hap_vector = get_hapvector_from_cigar(cigar_arr, query_sequence_encoded)
        read_hap_vectors[read_id] = read_hap_vector

    other_read_id = get_read_id(other_read)
    other_start = np.int32(other_read.reference_start)

    if other_read_id in read_hap_vectors:
        other_read_hap_vector = read_hap_vectors[other_read_id]
    else:
        other_read_hap_vector = get_hapvector_from_cigar(np.array(other_read.cigartuples, dtype=np.int32), other_query_sequence_encoded)
        read_hap_vectors[other_read_id] = other_read_hap_vector

    interval_hap_vector = numba_slicing(read_hap_vector, overlap_start, overlap_end, start)
    interval_other_hap_vector = numba_slicing(other_read_hap_vector, overlap_start, overlap_end, other_start)

    overlap_span = overlap_end - overlap_start
   
    diff_indices = numba_diff_indices(interval_hap_vector, interval_other_hap_vector)

    # First test whether there are too many mismatches between two reads that we wont tolerate
    if len(diff_indices) >= 3:
        # Cannot tolerate such mismatches
        return False, read_ref_pos_dict, read_hap_vectors, None

    diff_pos_read = interval_hap_vector[diff_indices]
    diff_pos_oread = interval_other_hap_vector[diff_indices]

    read_indel_overlap = get_indel_bools(diff_pos_read)
    oread_indel_overlap = get_indel_bools(diff_pos_oread)

    if numba_sum(read_indel_overlap) > 0:
        return False, read_ref_pos_dict, read_hap_vectors, None
    if numba_sum(oread_indel_overlap) > 0:
        return False, read_ref_pos_dict, read_hap_vectors, None

    # Now use overlap_start and overlap_end to extract the sequence
    read_seq, qr_idx_arr = get_interval_seq(np.int32(read.reference_start), 
                                            np.int32(overlap_start), 
                                            np.int32(overlap_end - 1), 
                                            ref_positions, 
                                            qseq_ref_positions, 
                                            query_sequence_encoded)
    
    other_seq, other_qr_idx_arr = get_interval_seq( np.int32(other_read.reference_start), 
                                                    np.int32(overlap_start), 
                                                    np.int32(overlap_end - 1), 
                                                    other_ref_positions, 
                                                    other_qseq_ref_positions, 
                                                    other_query_sequence_encoded )

    # We need to cut N output comparison for read_seq and other_seq
    
    total_match = compare_sequences(read_seq, other_seq, np.int8(4))

    identical_idx = numba_compare(interval_hap_vector, interval_other_hap_vector)
    identical_part = numba_bool_indexing(interval_hap_vector, identical_idx)
    overlap_span = identical_part.size
    var_count = count_var(identical_part)
    indel_num = count_continuous_indel_blocks(identical_part)

    weight = overlap_span + numba_sum(score_arr[:var_count])
    weight = weight + mean_read_length * 3 * indel_num

    if total_match:
        # logger.debug(f"The two reads {read_id} and {other_read_id} are identical in the overlapping region. The overlap span is {read.reference_name}:{overlap_start}-{overlap_end}. With a weight of {weight}, shared variant count is {var_count}, and shared indel count is {indel_num}. The identical part looks like {identical_part.tolist()}\nThe interval hapvector within overlap region for {read_id} is {interval_hap_vector.tolist()}\nAnd the interval hapvector within overlap region for {other_read_id} is {interval_other_hap_vector.tolist()}")

        # logger.debug(f"The qseq_ref_positions of read {read_id} starting from {read.reference_start} is {qseq_ref_positions}, and the ref_positions is {ref_positions}")
        # logger.debug(f"The qseq_ref_positions of read {other_read_id} starting from {other_read.reference_start} is {other_qseq_ref_positions}, and the ref_positions is {other_ref_positions}")
        return True, read_ref_pos_dict, read_hap_vectors, weight
    else:
        if len(read_seq) != len(other_seq) or interval_hap_vector.size != interval_other_hap_vector.size:
            # logger.debug(f"The two reads {read_id} and {other_read_id} are aligned to the same reference span {read.reference_name}:{overlap_start}-{overlap_end} with different size of sequence.")
            # logger.debug(f"The qseq_ref_positions of read {read_id} starting from {read.reference_start} is {qseq_ref_positions.tolist()}, and the ref_positions is {ref_positions.tolist()}")
            # logger.debug(f"The qseq_ref_positions of read {other_read_id} starting from {other_read.reference_start} is {other_qseq_ref_positions.tolist()}, and the ref_positions is {other_ref_positions.tolist()}")
            return False, read_ref_pos_dict, read_hap_vectors, None

        # Here we know that there are only equal to or smaller than 2 mismatches between the two reads
        # Then we need to know whether the two mismatches are in-trans variants or in-cis variants, if its in-trans we cannot tolerate the two mismatches
        # if read_seq_arr is None:
        #     base_dict = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
        #     read_seq_arr = np.array(list(map(base_dict.get, read_seq)), dtype=np.int8)
        #     other_seq_arr = np.array(list(map(base_dict.get, other_seq)), dtype=np.int8)

        q_diff_indices = numba_diff_indices(read_seq, other_seq)
        # logger.debug(f"The mismatch indices between {read_id} and {other_read_id} within {read.reference_name}:{overlap_start}-{overlap_end} relative to the query sequence are {q_diff_indices.tolist()}")
        r_diff_indices = qr_idx_arr[q_diff_indices]
        # logger.debug(f"The mismatch indices between {read_id} and {other_read_id} within {read.reference_name}:{overlap_start}-{overlap_end} relative to the reference genome are {r_diff_indices.tolist()}")
        if numba_contain(r_diff_indices, np.int32(-1)):
            # logger.debug(f"Within in region {read.reference_name}:{overlap_start}-{overlap_end}, read {read_id} and read {other_read_id} have mismatches on indels, which make this two read pairs not intolerant to the mismatches")
            # logger.debug(f"The query sequence within this region for read {read_id} is {read_seq.tolist()}, and the mismatch indices relative to the query sequence are {q_diff_indices.tolist()}")
            # logger.debug(f"The query sequence within this region for read {other_read_id} is {other_seq.tolist()}, and the mismatch indices relative to the reference genome are {r_diff_indices.tolist()}")
            return False, read_ref_pos_dict, read_hap_vectors, None
        read_package = (np.int32(read.reference_start), query_sequence_encoded, query_sequence_qualities, ref_positions)
        other_read_package = (np.int32(other_read.reference_start), other_query_sequence_encoded, other_query_sequence_qualities, other_ref_positions)
        # logger.info(f"The ref positions that two overlap reads are different are {r_diff_indices.tolist()}, while the read ref start is {read.reference_start} and the other read ref start is {other_read.reference_start}")
        tolerate, tolerated_count = tolerate_mismatches_two_seq(read_package,
                                                                other_read_package,
                                                                r_diff_indices,
                                                                nested_ad_dict,
                                                                empty_dict )

        if tolerate:
            # logger.debug(f"The two reads {read_id} and {other_read_id} within {read.reference_name}:{overlap_start}-{overlap_end} are tolerant to {tolerated_count} mismatches, with a weight of {weight}")
            weight = weight - tolerated_count * 20
            return True, read_ref_pos_dict, read_hap_vectors, weight
        else:
            return False, read_ref_pos_dict, read_hap_vectors, None