import logging
import pandas as pd


def read_delimited_file(ext:str, file_name:str, spec:dict):
    logger = logging.getLogger(__name__)

    # Get compression format
    if ext in [".csv", ".tsv"]:
        compression = None
    else:
        compression = "gzip" if ext == ".gz" else ext[1:]

    # Get delimiter
    try:
        delimiter = spec['delimiter']
    except Exception as e:
        logger.error("Delimiter (e.g. ',', '\t', ...) must be specified for compressed/delimited format.")
        raise RuntimeError

    df = pd.read_csv(file_name, compression=compression, sep=delimiter)
    logger.info(f"Target dataframe {file_name} loaded. Top 5 records: {df.head()}")
    return df


