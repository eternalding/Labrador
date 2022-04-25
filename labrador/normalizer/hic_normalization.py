import logging
import scipy.sparse as sparse
from scipy.sparse import coo_array
from scipy.sparse import linalg
import gc
import numpy as np
from tqdm import tqdm


def ICE_normalization(arr: sparse._arrays.coo_array, max_iter: int = 2000, threshold: float = 1e-5):
    """
    ICE normalization.
    Detailed information could be found in https://www.nature.com/articles/nmeth.2148.
    :param arr: Sparse array to be normalized
    :param max_iter: Max iteration for convergence
    :param threshold: Threshold for bias difference between iterations
    :return: normalized sparse array
    """
    logger = logging.getLogger(__name__)

    max_bins = max(arr.row.max(), arr.col.max())

    # Initialize variables
    bias = np.ones(max_bins + 1)
    arr_iter = arr.copy()
    arr_iter.eliminate_zeros()
    axis = 1 if arr.row.max() >= arr.col.max() else 0

    # Start iteration
    for iteration in range(1, max_iter+1):
        # Calculate bin coverage
        bin_coverage = arr_iter.sum(axis=axis)

        # Calculate delta bias
        delta_bias = bin_coverage / bin_coverage[np.nonzero(bin_coverage)].mean()
        delta_bias[delta_bias == 0.0] = 1.0

        # Calculate weights
        row_bias = delta_bias[arr_iter.row]
        col_bias = delta_bias[arr_iter.col]
        weights = row_bias * col_bias

        # Update arr_iter
        arr_iter.data /= weights

        new_bias = bias * delta_bias
        diff = abs(new_bias - bias).sum()
        if diff < threshold:
            logger.info("Reached convergence criteria.")
            break
        bias = new_bias
        if iteration % 50 == 0:
            logger.info(f"Iteration {iteration}: Bias difference = {diff}")

    # Normalize by learnt biases
    normalized_val = [val/(bias[row]*bias[col]) for val, row, col in zip(arr.data, arr.row, arr.col)]
    normalized_arr = coo_array((normalized_val, (arr.row, arr.col)))
    return normalized_arr


def SK_normalization(arr: sparse._arrays.coo_array):
    """
    TODO:
    Sinkhorn & Knopp normalization.
    Balance matrix to a doubly-stochastic matrix.
    Non-negative matrix is required.
    :param arr:
    :return:
    """
    if arr.min() < 0:
        raise RuntimeError("SK normalization requires input matrix >=0!")
    pass


def KR_normalization(arr: sparse._arrays.coo_array):
    """
    TODO:
    :param arr:
    :return:
    """
    pass


def SCN_normalization(arr: sparse._arrays.coo_array, max_iter: int = 5000, threshold: float = 1e-3):
    """
    Sequential Component Normalization.
    Detailed information could be found in https://bmcgenomics.biomedcentral.com/track/pdf/10.1186/1471-2164-13-436.pdf
    :param arr: Sparse array to be normalized
    :param max_iter: Max iteration for convergence
    :param threshold: Minimum requirement of difference between iterations
    :return:
    """
    logger = logging.getLogger(__name__)
    arr_iter = arr.copy()
    arr_iter.eliminate_zeros()

    def l2_norm(arr: sparse._arrays.coo_array, axis: int):
        """scipy norm currently does not support np.ndarray type."""
        return np.sqrt(arr.power(2).sum(axis=axis))

    for iteration in range(max_iter):
        """
        TODO: Fix euc norm problem on scipy
        """
        # Normalize by col
        col_norm = l2_norm(arr_iter, axis=0)
        col_biases = col_norm[arr_iter.col]
        arr_iter.data /= col_biases

        # Normalize by row
        row_norm = l2_norm(arr_iter, axis=1)
        row_biases = row_norm[arr_iter.row]
        arr_iter.data /= row_biases

        row_sum = arr_iter.sum(axis=1)
        col_sum = arr_iter.sum(axis=0)

        if iteration % 50 == 0:
            logger.info(f"Iteration {iteration}: Average sum of rows: {row_sum.mean()}")
    return arr_iter

def HiCNorm_normalization(arr: sparse._arrays.coo_array):
    """
    TODO:
    :param arr:
    :return:
    """
    pass
