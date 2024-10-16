#! /usr/bin/env python3
import numba
import logging
import numpy as np
from scipy import sparse
from numba import prange, types

from numba_operators import numba_and, numba_sum


logger = logging.getLogger("SDrecall")


@numba.njit(types.boolean[:](types.int32), fastmath=True)
def create_default_true_mask(size):
    return np.ones(size, dtype=numba.bool_)


@numba.njit(types.boolean[:](types.int32), fastmath=True, parallel=True)
def para_create_default_true_mask(size):
    return np.ones(size, dtype=numba.bool_)



@numba.njit(types.boolean[:](types.int32, types.int32[:]), fastmath=True)
def reverse_boolean_mask(total_rows, row_indices):
    boolean_mask = np.ones(total_rows, dtype=np.bool_)
    for index in row_indices:
        boolean_mask[index] = False
    return boolean_mask



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



@numba.njit(types.Tuple((types.int32[:], types.boolean[:]))(types.float32[:], types.int32[:], types.int32[:], types.int32, types.boolean[:], types.float32), fastmath=True)
def heuristic_find_largest_edge_weight_clique_sparse(matrix_data,
                                                     matrix_indices,
                                                     matrix_indptr,
                                                     matrix_size,
                                                     initial_index_mask,
                                                     cutoff=0.1):
    # assert sparse.isspmatrix_csr(weight_matrix), "Input matrix must be in CSR format"
    '''
    This function is to find the largest clique in a sparse matrix.
    '''
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



def bk_algorithm(selected_indices,
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