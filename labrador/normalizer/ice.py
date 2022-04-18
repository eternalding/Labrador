import logging
import scipy.sparse as sparse
from scipy.sparse import coo_array
import gc
import numpy as np
from tqdm import tqdm


def ICE_normalization(arr: sparse._arrays.coo_array, max_iter: int = 50, threshold: float = 1e-3):
    """
    ICE normalization.
    Detailed information could be found in https://www.nature.com/articles/nmeth.2148.
    :param arr: Sparse array to be normalized
    :param max_iter: Max iteration for convergence
    :param threshold: Threshold for bias difference between iterations
    :return: normalized sparse array
    """
    logger = logging.getLogger(__name__)
    logger.info("Perform ICE normalization.")

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
        weights = [delta_bias[i] * delta_bias[j]
                   for i, j in zip(arr_iter.row, arr_iter.col)]
        # Update arr_iter
        arr_iter.data /= weights

        new_bias = bias * delta_bias
        diff = abs(new_bias - bias).sum()
        if diff < threshold:
            break
        bias = new_bias
        if iteration % 10 == 0:
            logger.info(f"Iteration {iteration}: Bias difference = {diff}")

    # Normalize by learnt biases
    normalized_val = [val/(bias[row]*bias[col]) for val, row, col in zip(arr.data, arr.row, arr.col)]
    normalized_arr = coo_array((normalized_val, (arr.row, arr.col)))
    return normalized_arr



