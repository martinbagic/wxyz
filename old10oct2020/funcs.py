"""."""
import yaml
import pathlib
import logging
import numpy as np
import collections

from xhelp.recorder import Recorder
from xhelp.config import Config

project_path = pathlib.Path(__file__).absolute().parent


class Aux:
    def __init__(self, path_default, params_extra, recpath):
        # def initialize(self, path_default, params_extra, recpath):
        def _get_params():
            with open(path_default, "r") as ifile:
                params_default = yaml.safe_load(ifile)
            return {**params_default, **params_extra}

        def _check_params(params):
            assert params["BITS_PER_LOCUS"] % 2 == 0

        params = _get_params()
        _check_params(params)

        # simulation constants
        self.CYCLE_NUM = params["CYCLE_NUM"]
        self.LOGGING_RATE = params["LOGGING_RATE"]
        self.REC_EVERY_NTH = params["REC_EVERY_NTH"]
        self.ENTRY_LIMIT = params["ENTRY_LIMIT"]
        self.RECORD_PHENOTYPE = params["RECORD_PHENOTYPE"]
        self.PICKLE_RATE = params["PICKLE_RATE"]

        # simulation variables
        self.stage = 0
        self.max_popid = 0

        # recorder
        self.recorder = Recorder(recpath, self.ENTRY_LIMIT)

        # configs
        self.configs = [Config(params)]

    def get_popid(self):
        self.max_popid += 1
        return self.max_popid

