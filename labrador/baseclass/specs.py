from dataclasses import dataclass
import logging

@dataclass
class BaseSpec:
    required: dict
    optional: dict
    required_fields: list


@dataclass
class DelimiterSpec(BaseSpec):
    def __init__(self, exp_config:dict):
        self.logger = logging.getLogger(__name__)
        self.required_fields = ['chrom1_field', 'chrom2_field', 'bin1_field', 'bin2_field', 'main_value_field']

        # Fill-in required fields
        for field in self.required_fields:
            try:
                self.required[field] = exp_config[field]
            except Exception as E:
                self.logger.error("Missing some of the required field. Exception: {E}")
                raise RuntimeError

        # Fill-in optional field
        for k,v in exp_config:
            if k not in self.required_fields:
                self.optional[k] = v


@dataclass
class CoolerSpec(BaseSpec):
    def __init__(self, exp_config):
        # TODO
        self.required = {}
        self.optional = {}
        self.required_fields = []
