import scipy.sparse as sparse
import numpy as np
import multiprocessing
from functools import partial
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon, mannwhitneyu
import logging


def get_binsignal(i: int, arr, w: int):
    row_idx = np.logical_and(i-w <= arr.row, arr.row <= i-1)
    col_idx = np.logical_and(i+1 <= arr.col, arr.col <= i+w)
    valid_idx = row_idx * col_idx
    signal = arr.data[valid_idx].sum() / (w*w)
    return signal


def get_within_area(i: int, arr, w: int):
    row_idx = np.logical_and(i-w <= arr.row, arr.row <= i-1)
    col_idx = np.logical_and( i+1<= arr.col, arr.col <= i+w)
    valid_idx = np.logical_and(row_idx, col_idx)
    signal = arr.data[valid_idx]
    return signal


def get_between_area(i: int, arr, w: int):
    # Upper area
    row_idx = np.logical_and(i-w <= arr.row, arr.row <= i-1)
    col_idx = np.logical_and(i-w <= arr.col, arr.col <= i-1)
    valid_idx = row_idx * col_idx
    upper_area = arr.data[valid_idx]
    # Down area
    row_idx = np.logical_and(i+1 <= arr.row, arr.row <= i+w)
    col_idx = np.logical_and(i+1 <= arr.col, arr.col <= i+w)
    valid_idx = row_idx * col_idx
    down_area = arr.data[valid_idx]
    return np.concatenate([upper_area,down_area])


def TopDom(arr: sparse._arrays.coo_array, window_size: int, multiprocess: bool = True,
           resolution: int = 5000, wilcoxon_test: bool = True):
    """
    TopDom: an efficient and deterministic method for
    identifying topological domains in genomes.
    https://academic.oup.com/nar/article/44/7/e70/2467818
    :return:
    """
    logger = logging.getLogger(__name__)
    logger.info("Perform TopDom TAD caller")
    logger.info("Step 1. Calculate bin signal.")
    # Step 1. Calculate bin signal
    if multiprocess:
        query_func = partial(get_binsignal, arr=arr, w=window_size // resolution)
        cpus = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=cpus)
        bin_signal = pool.map_async(query_func, range(arr.col.max() + 1)).get()
        pool.close()
        pool.join()
    else:
        bin_signal = [get_binsignal(i, arr, window_size // resolution) for i in tqdm(range(arr.col.max()+1))]

    # Step 2. Curve fitting
    logger.info("Step 2. Curve fitting.")
    fitness, prev_fitness, sig_start, sig_end = 0, 0, 0, len(bin_signal) - 1
    sig_cur = 0

    turning_points = []
    while sig_start <= sig_cur <= sig_end:
        line_length = sig_cur - sig_start + 1
        cur_line = np.linspace(bin_signal[sig_start], bin_signal[sig_cur], line_length)
        error = abs(cur_line - bin_signal[sig_start:sig_cur + 1]).sum()
        fitness = line_length - error
        if fitness < prev_fitness:
            turning_points.append(sig_cur - 1)
            sig_start = sig_cur
            prev_fitness = 0
        else:
            prev_fitness = fitness
            sig_cur += 1

    local_minimum = []
    for i in range(1, len(turning_points)-1):
        prev = bin_signal[turning_points[i-1]]
        cur = bin_signal[turning_points[i]]
        next = bin_signal[turning_points[i+1]]
        if prev > cur and next > cur:
            local_minimum.append(turning_points[i])

    if wilcoxon_test:
        logger.info("Step 3. Wilcoxon ranksum test.")
        wilcoxon_significant = []
        w = window_size // resolution
        for point in tqdm(local_minimum):
            between_area = np.pad(get_between_area(point, arr, w), (w * w), constant_values=0)
            within_area = np.pad(get_within_area(point, arr, w), (w * w), constant_values=0)
            pvalue = \
            mannwhitneyu(within_area, between_area, use_continuity=True, alternative='less', nan_policy="omit")[1]
            if pvalue < 0.05:
                wilcoxon_significant.append([point, pvalue])
        wilcoxon_significant_points = [point[0] for point in wilcoxon_significant]
        wilcoxon_significant = np.array(wilcoxon_significant)
        local_minimum = wilcoxon_significant_points

    local_minimum = [val * resolution for val in local_minimum]
    return local_minimum

