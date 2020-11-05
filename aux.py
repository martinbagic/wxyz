"""."""
import yaml
import numpy as np
import collections
import logging

from xhelp.recorder import Recorder
from config import Config


class Aux:
    def __init__(self, paths_config, cmd_params, recpath):
        def _get_params(path):
            with open(path, "r") as ifile:
                return yaml.safe_load(ifile)

        def _check_params(params):
            assert params["BITS_PER_LOCUS"] % 2 == 0 or params["REPR_MODE"] == "asexual"

        def _add_cmd_params(params, cmd_params):
            """Add GENOME_STRUCT and GENOME_CONST modifications such as GENOME_STRUCT_surv_1"""
            for key, val in cmd_params.items():
                if key.startswith("GENOME_STRUCT_"):
                    trait = key.split("_")[2]
                    index = int(key.split("_")[3])
                    params["GENOME_STRUCT"][trait][index] = val
                elif key.startswith("GENOME_CONST_"):
                    trait = key.split("_")[2]
                    params["GENOME_CONST"][trait] = val
                else:
                    params[key] = val
            return params

        config_params = {}
        for path_config in paths_config:
            config_params = {**_get_params(path_config), **config_params}

        params = _add_cmd_params(config_params, cmd_params)
        _check_params(params)

        logging.info(f"Final parameters: {repr(params)}")

        # simulation constants
        self.CYCLE_NUM = params["CYCLE_NUM"]
        self.LOGGING_RATE = params["LOGGING_RATE"]
        self.REC_EVERY_NTH = params["REC_EVERY_NTH"]
        self.FLUSH_RATE = params["FLUSH_RATE"]
        self.PICKLE_RATE = params["PICKLE_RATE"]

        # simulation variables
        self.stage = 0
        self.max_popid = 0

        # configs
        self.configs = [Config(params)]
        config = self.configs[0]

        # recorder
        self.recorder = Recorder(
            recpath,
            params,
            self.FLUSH_RATE,
            params["MAX_LIFESPAN"],
            config.loci_pos,
            config.BITS_PER_LOCUS,
            config.MATURATION_AGE,
        )

    def get_popid(self):
        self.max_popid += 1
        return self.max_popid

