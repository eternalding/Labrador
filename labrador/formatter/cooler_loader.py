import logging
from tqdm import tqdm
import cooler
from labrador.baseclass.sparsemat import SparseMat
import h5py
from pathlib import Path
from typing import Optional
from functools import partial


def get_avail_resolutions(hf: h5py._hl.files.File):
    if 'resolutions' not in hf:
        raise RuntimeError("Designated cooler file does not have resolution information.")
    resolutions = sorted(list(hf['resolutions'].keys()), key = lambda x:int(x))
    return resolutions


def read_single_chrom(cooler_file, chrom: str, fields: list, resolution: int):
    logger = logging.getLogger(__name__)
    logger.info(f"Now loading chromosome {chrom}")

    values = {}
    matrix = None
    for field in fields:
        if field == "Count":
            field_name = False
        elif field == "Default":
            field_name = True
        else:
            field_name = field
        matrix = cooler_file.matrix(balance=field_name, sparse=True).fetch(chrom)
        matrix.eliminate_zeros()
        values[field] = matrix.data

    return SparseMat(matrix.row, matrix.col, values, chrom=chrom, resolution=resolution)


def read_cooler_file(data_path: str, resolution: int, fields: list, chroms: Optional[list, str] = None):

    logger = logging.getLogger(__name__)

    hf = h5py.File(data_path, 'r')

    avail_res = get_avail_resolutions(hf)
    if resolution not in avail_res:
        raise RuntimeError(f"Requested resolution {resolution} does not exist. "
                           f"Available resolutions: {avail_res}")

    # Read file
    ext = Path(data_path).suffix
    if ext == ".mcool":
        data_path = f"{data_path}::/resolutions/{resolution}"

    c = cooler.Cooler(data_path)
    # Extract sparse matrices and convert into sparse array
    sparse_mats = {}
    if isinstance(chroms, str):
        sparse_mats[chroms] = read_single_chrom(c, chroms, fields, resolution)
    elif isinstance(chroms, list):
        for chrom in chroms:
            sparse_mats[chrom] = read_single_chrom(c, chrom, fields, resolution)
    else:
        raise NotImplementedError("Chromosome in config file must be list, e.g. ['chr1', 'chr2' ...] "
                                  "or str, e.g. 'chr1'.")