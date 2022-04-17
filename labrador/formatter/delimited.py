import logging
import pandas as pd
from tqdm import tqdm

from labrador.baseclass.sparsemat import SparseMat

CHUNK_SIZE = 16_384


def read_delimited_file(ext: str, file_name: str, spec: dict):
    logger = logging.getLogger(__name__)
    # Get compression format
    if ext in [".csv", ".tsv"]:
        compression = None
    else:
        compression = "gzip" if ext == ".gz" else ext[1:]

    # Get delimiter
    delimiter = spec.required['delimiter']

    # Read delimited data
    logger.info("Now loading Hi-C dataset:")
    df = pd.concat([chunk for chunk in
                    tqdm(pd.read_csv(file_name, compression=compression, sep=delimiter, chunksize=CHUNK_SIZE))])
    logger.info(f"Target dataframe {file_name} loaded. Top 5 records: {df.head()}")

    # Retain only required columns
    if spec.required['chrom1_field'] == spec.required['chrom2_field']:
        logger.info(f"chrom1_field is the same as chrom2_field. Treated as intra-chromosomal data.")
        chrom_field_name = spec.required['chrom1_field']
        df[f"{chrom_field_name}_2"] = df[chrom_field_name]
        spec.required['chrom2_field'] = f"{chrom_field_name}_2"

    df.rename(columns={v: k for k, v in spec.required.items()}, inplace=True)
    df = df[['chrom1_field', 'chrom2_field', 'bin1_field', 'bin2_field', 'main_value_field'] +
             spec.metadata['optional_field']]

    # Drop nan records
    logger.info("Drop na records:")
    df.dropna(inplace=True)

    logger.info(f"Total records: {len(df)}")

    sparse_mats = {}
    # Create sparse matrices
    # TODO: Support inter-chromosomal interactions
    for chrom, group in df.groupby('chrom1_field'):
        row_idxs = (group['bin1_field'] / spec.metadata['resolution']).to_numpy(dtype='int')
        col_idxs = (group['bin2_field'] / spec.metadata['resolution']).to_numpy(dtype='int')
        value_fields = ['main_value_field'] + spec.metadata['optional_field']
        values = {k: group[k].values for k in value_fields}

        sparse_mats[chrom] = SparseMat(row_idxs, col_idxs, values, chrom=chrom,
                                       resolution=spec.metadata['resolution'])

    return sparse_mats


