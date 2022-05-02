import pandas as pd
import numpy as np
from typing import List
from pathlib import Path
from labrador.baseclass.specs import DelimiterSpec, CoolerSpec
from labrador.formatter.delimited import read_delimited_file
from labrador.visualizer.plot import plot_interactive_matrix, gen_heatmap, plot_multiple_arrays
from labrador.baseclass.cooler_loader import read_cooler_file
import logging
import gc

SUPPORTED_FORMATS = [".csv", ".tsv", ".cool", ".mcool"]


class Experiment:
    def __init__(self, exp_config: dict = None):
        self.logger = logging.getLogger(__name__)
        self.exp_config = exp_config

        self.spec = {}
        self.sparse_mats = {}

    def create_spec(self, spec_type:str):
        if spec_type == "delimiter":
            self.spec = DelimiterSpec(self.exp_config['Format']['Delimited'])
            for k, v in self.exp_config.items():
                if k != "Format":
                    self.spec.metadata[k] = v

        elif spec_type == "cooler":
            self.spec = CoolerSpec(self.exp_config)

    def read_file(self):
        # Clear sparse matrices if exist and release mem
        self.sparse_mats = None
        gc.collect()

        # Read file
        filename = self.exp_config['Metadata']['data']
        ext = Path(filename).suffix
        if ext in [".gz", ".bz2", ".zip", ".xz", ".csv", ".tsv"]:
            # Create spec
            self.create_spec(spec_type="delimiter")
            self.sparse_mats = read_delimited_file(ext, filename, self.spec)
        elif ext == ".mcool" or ext == ".cool":
            # TODO: Support cooler files
            resolution = self.spec.metadata['resolution']
            fields = list(self.spec.metadata['main_value_field']) + self.spec.metadata['optional_field']
            chroms = self.spec.metadata['chromosomes']
            self.sparse_mats = read_cooler_file(filename, resolution, fields, chroms)
        else:
            raise NotImplementedError(f"Labrador doesn't support file format {ext}! "
                                       "Supported formats: {SUPPORTED_FORMATS}")

    def get_avail_chroms(self):
        return list(self.sparse_mats.keys())

    def query(self, chrom: str, start: int, end: int, value: str = "main_value_field", to_symmetric: bool = False):
        if chrom not in self.sparse_mats:
            self.logger.error(
                f"Assigned chromosome {chrom} doesn't exist. "
                f"Current available chromosome: {self.get_avail_chroms()}")
            raise NotImplementedError
        return self.sparse_mats[chrom].query(start, end, value, to_symmetric)

    def visualize(self, chrom: str, start: int, end: int, value: str,
                         to_symmetric: bool = False, *args, **kwargs):
        resolution = self.spec.metadata['resolution']

        arr = self.query(chrom, start, end, value, to_symmetric=to_symmetric)
        x_labels = np.arange(start, end, resolution)
        y_labels = np.arange(start, end, resolution)
        kwargs['title'] = kwargs['title'] if "title" in kwargs else f"{chrom}: {start}-{end} resolution = {resolution} bp"

        vb = plot_interactive_matrix(arr, x_labels, y_labels, *args, **kwargs)
        return vb

    def normalize(self, chrom: str, method: str, field: str, *args, **kwargs):
        if chrom not in self.sparse_mats:
            raise ValueError(f"Requested chromosome {chrom} does not exist.")
        self.sparse_mats[chrom].normalize(method=method, norm_by_field=field, *args, **kwargs)

    def compare_fields(self, fields: List[str], chrom: str, start: int, end: int):
        num_arrs = len(fields)
        current_fields = self.sparse_mats[chrom].get_available_fields()

        # Check if all requested fields exist
        for field in fields:
            if field not in current_fields:
                raise ValueError(f"Requested field {field} not exist! Current fields are: {current_fields}")

        # Query all fields
        target_arrs = {field: self.query(chrom, start, end, field, to_symmetric=True)
                       for field in fields}

        # Generate heatmaps
        resolution = self.spec.metadata['resolution']
        target_arrs = {field: gen_heatmap(idx, arr, start, end, field, resolution)
                       for idx, (field, arr) in enumerate(target_arrs.items(), 1)}

        fig = plot_multiple_arrays(target_arrs, title=f"{chrom}: {start}-{end}")
        return fig


