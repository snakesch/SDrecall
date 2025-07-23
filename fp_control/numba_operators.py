import numpy as np
import numba
from numba import types, prange


@numba.njit
def fast_median(data):
    return np.median(data)


@numba.njit(fastmath=True, parallel=True)
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
    for value in array:
        if value == -1:
            return True
    return False


@numba.njit(parallel=False)
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



@numba.njit(types.float32[:,:](types.float32[:,:], types.boolean[:]), fastmath=True)
def numba_ix(matrix, mask):
    return matrix[mask, :][:, mask]



@numba.njit(types.boolean[:](types.int32), fastmath=True)
def create_default_false_mask(size):
    return np.zeros(size, dtype=numba.bool_)



@numba.njit(types.int32[:](types.int32, types.boolean[:]), fastmath=True, parallel=False)
def apply_index_mask(size, initial_index_mask):
    return np.arange(size, dtype=np.int32)[initial_index_mask]


@numba.njit(fastmath=True)
def numba_diff_indices(arr1, arr2):
    """Optimized sequential implementation"""  
    # Fill result array
    result = np.empty(arr1.size, dtype=np.int16)
    idx = 0
    for i in range(arr1.size):
        if arr1[i] != arr2[i]:
            result[idx] = i
            idx += 1
    
    return result[:idx]


@numba.njit(types.int16[:](types.int16[:], types.int32, types.int32, types.int32), fastmath=True)
def numba_slicing(arr, overlap_start, overlap_end, start):
    return arr[overlap_start - start:overlap_end - start]


@numba.njit(types.int32(types.int32[:], types.int32), fastmath=True)
def numba_indexing_int32(arr, idx):
    return arr[idx]


@numba.njit(types.int8(types.int8[:], types.int32), fastmath=True)
def numba_indexing_int8(arr, idx):
    return arr[idx]


@numba.njit(fastmath=True)
def numba_compare(arr1, arr2):
    return arr1 == arr2



@numba.njit(types.int16[:](types.int16[:], types.bool_[:]), fastmath=True)
def numba_bool_indexing(arr, idx_arr):
    return arr[idx_arr]



@numba.njit(types.bool_(types.int32[:], types.int32), fastmath=True)
def numba_contain(arr, scalar):
    return (arr == scalar).any()