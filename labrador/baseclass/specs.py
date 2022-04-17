from dataclasses import dataclass
import logging


class BaseSpec:
    def __init__(self):
        self.required = {}
        self.metadata = {}


class DelimiterSpec(BaseSpec):
    def __init__(self, exp_config: dict):
        super().__init__()
        self.logger = logging.getLogger(__name__)
        required_fields = ['chrom1_field', 'chrom2_field', 'bin1_field',
                           'bin2_field', 'main_value_field', 'delimiter']

        # Fill-in required fields
        for field in required_fields:
            try:
                self.required[field] = exp_config[field]
            except Exception as E:
                self.logger.error(f"Missing some of the required field. Exception: {E}")
                raise RuntimeError

        # Fill-in optional field
        for k, v in exp_config.items():
            if k not in required_fields:
                self.metadata[k] = v

        # Overwrite resolution
        self.metadata['resolution'] = exp_config['resolution']


@dataclass
class CoolerSpec(BaseSpec):
    def __init__(self, exp_config):
        # TODO
        self.required = {}
        self.optional = {}
        self.required_fields = []
