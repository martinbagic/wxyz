"""."""
import yaml
import pathlib
import logging
import numpy as np
import collections

from recorder import Recorder
from overshoot import Overshoot
from interpreter import Interpreter
from phenomap import Phenomap

project_path = pathlib.Path(__file__).absolute().parent


class Config:
    def __init__(self, params):

        # population
        self.max_lifespan = params["max_lifespan"]
        self.maturation_age = params["maturation_age"]
        self.n_neutral_loci = params["n_neutral_loci"]
        self.n_mutation_loci = params["n_mutation_loci"]
        self.bits_per_locus = params["bits_per_locus"]
        self.age_guarantee = params["age_guarantee"]
        self.genome_distribution = params["genome_distribution"]
        self.mutrate_distribution = params["mutrate_distribution"]
        self.surv_bound_lo = params["surv_bound_lo"]
        self.surv_bound_hi = params["surv_bound_hi"]
        self.repr_bound_lo = params["repr_bound_lo"]
        self.repr_bound_hi = params["repr_bound_hi"]
        self.muta_bound_lo = params["muta_bound_lo"]
        self.muta_bound_hi = params["muta_bound_hi"]
        self.interpreter_surv = params["interpreter_surv"]
        self.interpreter_repr = params["interpreter_repr"]
        self.interpreter_neut = params["interpreter_neut"]
        self.interpreter_muta = params["interpreter_muta"]
        self.repr_mode = params["repr_mode"]
        self.recombination_rate = params["recombination_rate"]
        self.mutation_ratio = params["mutation_ratio"]
        self.mutation_probability = params["mutation_probability"]

        # ecology
        self.max_population_size = params["max_population_size"]
        self.overshoot_handling = params["overshoot_handling"]
        self.discrete_generations = params["discrete_generations"]
        self.envmap_rate = params["envmap_rate"]
        self.phenomap_plus = params["phenomap_plus"]
        self.bottleneck_size = params["bottleneck_size"]

        # derived
        self.pos_surv = 0
        self.pos_repr = self.pos_surv + self.max_lifespan
        self.pos_neut = self.pos_repr + self.max_lifespan - self.maturation_age
        self.pos_muta = self.pos_neut + self.n_neutral_loci
        self.pos_end = self.pos_muta + self.n_mutation_loci
        
        self.genome_shape = (
            self.max_population_size,  # ? individuals
            self.max_lifespan * 2
            - self.maturation_age
            + self.n_neutral_loci
            + self.n_mutation_loci,  # ? loci
            self.bits_per_locus,  # ? bits
        )

        # changing
        self.season_countdown = float("inf")
        self.reset_season_countdown()
        self.max_uid = 0

        # interpreter: genotype -> probabilities
        self.interpreter = Interpreter(self.bits_per_locus)

        # phenomap: pleiotropic genotype -> age-linear genotype
        self.phenomap = (
            Phenomap(phenomap_plus=self.phenomap_plus, pos_end=self.pos_end)
            if self.phenomap_plus
            else None
        )

        # overshoot: detects overshoots and returns killing masks
        self.overshoot = Overshoot(
            overshoot_handling=self.overshoot_handling,
            max_population_size=self.max_population_size,
            bottleneck_size=self.bottleneck_size,
        )

        # envmap: genotype -> genotype contextualized in a changing environment
        self.envmap = None
        # self.envmap = (
        #     Envmap() if self.envmap_rate and self.mutation_probability == -1 else None
        # )

    def pheno(self, loci, loci_kind):
        interpreter_kind = {
            "surv": self.interpreter_surv,
            "repr": self.interpreter_repr,
            "neut": self.interpreter_neut,
            "muta": self.interpreter_muta,
        }[loci_kind]
        return self.interpreter(loci, interpreter_kind)

    def reset_season_countdown(self):
        self.season_countdown = (
            self.discrete_generations if self.discrete_generations else float("inf")
        )
        
    def get_uids(self, n):
        """Get an array of unique origin identifiers."""
        uids = np.arange(n) + self.max_uid
        self.max_uid += n
        return uids
        
        


# class EnvSettings:
#     attrs = (
#         "max_population_size",
#     )

#     def __init__(self, params):
#         for attr in self.attrs:
#             setattr(self, attr, params[attr])

#         self.season_countdown = float("inf")

#         self.overshoot = Overshoot(
#             overshoot_handling=self.overshoot_handling,
#             max_population_size=self.max_population_size,
#             bottleneck_size=self.bottleneck_size,
#         )

#         self.reset_season_countdown()

#         # self.envmap = (
#         #     Envmap() if aux.envmap_rate and aux.mutation_probability == -1 else None
#         # )

#     def reset_season_countdown(self):
#         self.season_countdown = (
#             self.discrete_generations if self.discrete_generations else float("inf")
#         )


class Aux:
    def __init__(self, path_default, params_extra, recpath):
    # def initialize(self, path_default, params_extra, recpath):
        def _get_params():
            with open(path_default, "r") as ifile:
                params_default = yaml.safe_load(ifile)
            return {**params_default, **params_extra}

        def _check_params(params):
            assert params["bits_per_locus"] % 2 == 0

        params = _get_params()
        _check_params(params)

        # simulation constants
        # self.project_path = pathlib.Path(__file__).absolute().parent
        self.cycle_num = params["cycle_num"]
        self.logging_rate = params["logging_rate"]
        self.rec_every_nth = params["rec_every_nth"]
        self.entry_limit = params["entry_limit"]
        self.record_phenotype = params["record_phenotype"]
        self.pickle_rate = params["pickle_rate"]

        # simulation variables
        self.stage = 0
        self.max_popid = 0

        # recorder
        self.recorder = Recorder(recpath, self.entry_limit)

        # configs
        self.configs = [Config(params)]
        
    def get_popid(self):
        self.max_popid += 1
        return self.max_popid

# aux = Aux()

# class Aux:
#     attrs = (
#         "cycle_num",
#         "logging_rate",
#         "rec_every_nth",
#         "entry_limit",
#         "record_phenotype",
#         "pickle_rate",
#     )

#     def __init__(self):
#         self.stage = 0
#         self.max_uid = 0

#     def init(self, path_default, params_extra, recpath):

#         # get params

#         # check params

#         # set params for self
#         for attr in self.attrs:
#             setattr(self, attr, params[attr])

#         # make configs
#         configs = [Config(params)]

#         self.recorder = Recorder(recpath, self.entry_limit)


# aux = Aux()


# def init(self, path_default, params_extra, recpath):
#     # default config + job-specific input
#     with open(path_default, "r") as ifile:
#         params_default = yaml.safe_load(ifile)
#     params = {**params_default, **params_extra}

#     # check params
#     assert params["bits_per_locus"] % 2 == 0

#     global POP, ENV, SIM

#     POP = PopSettings(params)
#     ENV = EnvSettings(params)
#     SIM = SimSettings(params, recpath)


# POP, ENV, SIM = None, None, None

# class Aux:
# def __init__(self):
#     self.max_uid = 0
#     self.stage = 0
#     self.season_countdown = float("inf")

# def reset_season_countdown(self):
#     self.season_countdown = (
#         aux.discrete_generations if aux.discrete_generations else float("inf")
#     )

# def init_on_start(self, path_default, params_extra, recpath):
#     ### INPUT PARAMS ###
#     # default config + job-specific input
#     with open(path_default, "r") as ifile:
#         params_default = yaml.safe_load(ifile)
#     params = {**params_default, **params_extra}

#     for key, val in params.items():
#         setattr(self, key, val)

#     assert (
#         params["bits_per_locus"] % 2 == 0
#     )  # check that you have 2 conceptual chromosomes

#     ###
#     ### DERIVED PARAMS ###

#     # positions of first functional locus
#     self.pos_surv = 0
#     self.pos_repr = self.pos_surv + self.max_lifespan
#     self.pos_neut = self.pos_repr + self.max_lifespan - self.maturation_age
#     self.pos_muta = self.pos_neut + self.n_neutral_loci
#     self.pos_end = self.pos_muta + self.n_mutation_loci

#     # genome shape
#     self.genome_shape = (
#         aux.max_population_size,  # ? individuals
#         aux.max_lifespan * 2
#         - aux.maturation_age
#         + aux.n_neutral_loci
#         + aux.n_mutation_loci,  # ? loci
#         aux.bits_per_locus,  # ? bits
#     )

#     # interpreters genotype to phenotype
#     self.phenotype = Phenotype()

#     # overshoot
#     self.overshoot = Overshoot(self.overshoot_handling)

#     # recorder
#     self.recorder = Recorder(recpath, self.entry_limit)

#     if self.discrete_generations:
#         self.reset_season_countdown()

#     self.phenomap = Phenomap() if self.phenomap_plus else None

#     self.envmap = (
#         Envmap() if aux.envmap_rate and aux.mutation_probability == -1 else None
#     )


# aux = Aux()
