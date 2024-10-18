import numpy as np
import numba
from numba import types, prange


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



@numba.njit(types.float32[:,:](types.float32[:,:], types.boolean[:]), fastmath=True)
def numba_ix(matrix, mask):
    return matrix[mask, :][:, mask]



@numba.njit(types.boolean[:](types.int32), fastmath=True)
def create_default_false_mask(size):
    return np.zeros(size, dtype=numba.bool_)



@numba.njit(types.int32[:](types.int32, types.boolean[:]), fastmath=True, parallel=True)
def apply_index_mask(size, initial_index_mask):
    return np.arange(size, dtype=np.int32)[initial_index_mask]




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

