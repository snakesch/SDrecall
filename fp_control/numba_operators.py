import numba
# numba.config.THREADING_LAYER = 'omp'
# numba.set_num_threads(4)
from numba import types, prange, get_num_threads


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