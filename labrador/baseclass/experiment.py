import pandas as pd
import numpy as np
from pathlib import Path
from labrador.baseclass.specs import DelimiterSpec, CoolerSpec
from labrador.formatter.delimited import read_delimited_file
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
            for k, v in self.sparse_mats.items():
                print(k," ",v)
        elif ext == ".mcool" or ext == ".cool":
            # TODO: Support cooler files
            pass
        else:
            raise NotImplementedError(f"Labrador doesn't support file format {ext}! "
                                       "Supported formats: {SUPPORTED_FORMATS}")




