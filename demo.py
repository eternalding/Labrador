import yaml
import argparse
import logging

from labrador.baseclass.experiment import Experiment
from labrador.utils.logger import init_logger

# arg parser
parser = argparse.ArgumentParser(description='Labrador, an integrated tool-kit for Hi-C data analysis.')
parser.add_argument('--configs', type=str, help='Configuration for given Hi-C experiment')
parser.add_argument('--level', type=str, default="INFO", help='Level for logger messages.')
args = parser.parse_args()

# Init logger
init_logger(level=args.level)


def main():
    # Read configuration
    with open(args.configs, 'r') as fp:
        configs = yaml.safe_load(fp)

    # Read Hi-C experiment
    HiC_exp = Experiment(configs)
    HiC_exp.read_file()

if __name__ == "__main__":
    main()