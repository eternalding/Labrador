import pandas as pd
import numpy as np
from pathlib import Path
from labrador.baseclass.specs import DelimiterSpec, CoolerSpec
from labrador.formatter.delimited import read_delimited_file
import logging

SUPPORTED_FORMATS = [".csv", ".tsv", ".cool", ".mcool"]


class Experiment:
    def __init__(self, filename: str, delimiter: str = ",", exp_config: dict = None):
        self.logger = logging.getLogger(__name__)
        self.filename = filename
        self.delimiter = delimiter
        self.exp_config = exp_config

        self.spec = {}
        self.sparse_matrices = {}

        # Read file

    def create_spec(self, spec_type:str):
        if spec_type == "delimiter":
            self.spec = DelimiterSpec(self.exp_config)
        elif spec_type == "cooler":
            self.spec = CoolerSpec(self.exp_config)

    def read_file(self):
        ext = Path(self.filename).suffix
        if ext in [".gz", ".bz2", ".zip", ".xz", ".csv", ".tsv"]:
            # Create spec
            self.create_spec(spec_type="delimiter")
            self.sparse_matrices = read_delimited_file(ext, self.filename, self.spec)
        elif ext == ".mcool" or ext == ".cool":
            # TODO: Support cooler files
            pass
        else:
            raise NotImplementedError(f"Labrador doesn't support file format {ext}! "
                                       "Supported formats: {SUPPORTED_FORMATS}")




