import logging
import scipy.sparse as sparse
from scipy.sparse import coo_array
import numpy as np
from tqdm import tqdm


def minmax_normalization(arr: sparse._arrays.coo_array):
    logger = logging.getLogger(__name__)
    logger.info("Perform min-max normalization")

    min_val = arr.min()
    max_val = arr.max()

    norm_arr = arr.copy()
    norm_arr.data = (norm_arr.data - min_val) / (max_val - min_val)
    return norm_arr









