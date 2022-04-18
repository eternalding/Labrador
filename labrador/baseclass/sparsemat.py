from scipy.sparse import coo_array
import numpy as np
import logging
from labrador.normalizer.ice import ICE_normalization

SUPPORTED_NORM_METHODS = ["ICE"]

class SparseMat:
    def __init__(self, row_idxs: np.ndarray, col_idxs: np.ndarray, vals: dict,
                 chrom: str, resolution: int):

        self.logger = logging.getLogger(__name__)
        self.chrom = chrom
        self.row_idxs = row_idxs
        self.col_idxs = col_idxs
        self.matrices = {k: coo_array((v, (row_idxs, col_idxs)))
                            for k, v in vals.items()}
        self.resolution = resolution
        self.stats = {}


    def get_stats(self):
        for k, v in self.matrices.items():
            self.stats[k] = {}
            self.stats[k]['min'] = v.min()
            self.stats[k]['mean'] = v.mean()
            self.stats[k]['max'] = v.max()
            stored_vals = v.getnnz()
            self.stats[k]['std'] = np.std(stored_vals)

    def __str__(self):
        if not self.stats:
            self.get_stats()
        msg = ""
        msg += f"Sparse matrix for chromosome {self.chrom}:\n"
        msg += f"Number of records: {len(self.row_idxs)}\n"
        msg += f"Resolution: {self.resolution}\n"
        msg += f"Stats:\n"
        for k, v in self.stats.items():
            msg += f"[{k}]\n"
            for stat, val in self.stats[k].items():
                msg += f"{stat}: {val:.3f}\n"
        return msg

    def get_avail_values(self):
        return list(self.matrices.keys())

    def query(self, start: int, end: int, value: str = "main_value_field", to_symmetric: bool = False):
        if value not in self.matrices:
            self.logger.error(f"Assigned value field {value} doesn't exist. You must manually provide it, or calculate "
                              f"via provided normalization methods. Current value fields: {self.get_avail_values()}")
            raise NotImplementedError
        start_idx = start // self.resolution
        end_idx = end // self.resolution
        arr = self.matrices[value].toarray()[start_idx:end_idx, start_idx:end_idx]
        if to_symmetric:
            arr = arr + arr.T
        return arr

    def normalize(self, method: str, norm_by_value: str = "RawCount", *args, **kargs):
        if method not in SUPPORTED_NORM_METHODS:
            self.logger.error(f"Assigned normalization method {method} is not supported!"
                                "Currently supported methods: {SUPPORTED_NORM_METHODS}")
            raise NotImplementedError

        self.logger.info(f"Normalize {norm_by_value} with {method} normalization.")
        if method == "ICE":
            max_iter = kargs['max_iter'] if 'max_iter' in kargs else 100
            self.matrices[method] = ICE_normalization(self.matrices[norm_by_value], max_iter=max_iter)
        self.logger.info(f"Normalization is done. Value is stored in matrices[{method}]")


