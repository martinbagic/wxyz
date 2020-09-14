"""."""
import yaml
import pathlib
import logging
import numpy as np

from recorder import Recorder

# logging settings
logging.basicConfig(
    format="%(asctime)s : %(message)s",
    datefmt="%d/%m/%Y %I:%M:%S",
    level=logging.DEBUG,
)


class Aux:
    # constants
    PATH = pathlib.Path(__file__).absolute().parent

    # variables
    max_uid = None
    stage = None
    season_countdown = None
    resources = None

    # objects
    phenotype = None
    overshoot = None
    recorder = None

    # derived values
    genome_pos = None
    genome_shape = None

    # params_default + params_extra

    def __init__(self):
        self.max_uid = 0
        self.stage = 0
        self.season_countdown = float("inf")

    def reset_season_countdown(self):
        self.season_countdown = (
            aux.discrete_generations if aux.discrete_generations else float("inf")
        )

    def init_on_start(self, path_default, params_extra, recpath):
        ### INPUT PARAMS ###
        # default config + job-specific input
        with open(path_default, "r") as ifile:
            params_default = yaml.safe_load(ifile)
        params = {**params_default, **params_extra}

        for key, val in params.items():
            setattr(self, key, val)

        assert (
            params["bits_per_locus"] % 2 == 0
        )  # check that you have 2 conceptual chromosomes

        ###
        ### DERIVED PARAMS ###

        # positions of first functional locus
        self.pos_surv = 0
        self.pos_repr = self.pos_surv + self.max_lifespan
        self.pos_neut = self.pos_repr + self.max_lifespan - self.maturation_age
        self.pos_muta = self.pos_neut + self.n_neutral_loci
        self.pos_end = self.pos_muta + self.n_mutation_loci

        # genome shape
        self.genome_shape = (
            aux.max_population_size,  # ? individuals
            aux.max_lifespan * 2
            - aux.maturation_age
            + aux.n_neutral_loci
            + aux.n_mutation_loci,  # ? loci
            aux.bits_per_locus,  # ? bits
        )

        self.envmap = (
            np.zeros(self.genome_shape[1:], bool) if self.envmap_rate else None
        )

        # interpreters genotype to phenotype
        self.phenotype = Phenotype()

        # overshoot
        self.overshoot = Overshoot(self.overshoot_handling)

        # recorder
        self.recorder = Recorder(recpath, self.entry_limit)

        self.reset_season_countdown()

    def evolve_envmap(self):
        bit = np.random.choice(range(self.bits_per_locus))
        locus = np.random.choice(range(self.pos_end))
        self.envmap[locus, bit] = False if self.envmap[locus, bit] else True


aux = Aux()


class Phenotype:
    def __init__(self):
        bits_per_chromosome = aux.bits_per_locus // 2
        self.binary_weights = np.repeat(np.arange(1, bits_per_chromosome + 1) ** 2, 2)
        self.binary_max = self.binary_weights.sum()

    def decimal(self, loci):  # applicable to surv, repr, neut and muta
        return loci.sum(1) / aux.bits_per_locus

    def exp(self, loci):  # applicable to muta
        return 1 / 2 ** np.sum(loci, axis=1)

    def binary(self, loci):  # applicable to surv, repr
        return loci.dot(self.binary_weights) / self.binary_max

    def __call__(self, loci, loci_kind):
        interpreter_kind = {
            "surv": aux.interpreter_surv,
            "repr": aux.interpreter_repr,
            "neut": aux.interpreter_neut,
            "muta": aux.interpreter_muta,
        }[loci_kind]

        interpreter = {
            "decimal": self.decimal,
            "exp": self.exp,
            "binary": self.binary,
        }[interpreter_kind]

        phenotype = interpreter(loci)
        return phenotype


class Overshoot:
    def __init__(self, mode):
        self.func = {
            "treadmill_random": self.treadmill_random,
            "bottleneck": self.bottleneck,
            "treadmill_ageist": self.treadmill_ageist,
            "treadmill_boomer": self.treadmill_boomer,
            "limitless": self.limitless,
            "starvation": self.starvation,
        }[mode]

        if mode == "starvation":
            self.starvation_time = 0

    def starvation(self, n):
        """Let starved individuals die."""

        print(n)

        if n > aux.max_population_size:
            self.starvation_time += 1
            surv_probability = 0.95 ** self.starvation_time
            random_probabilities = np.random.random(n)
            mask = random_probabilities > surv_probability
            return mask
        else:
            self.starvation_time = 0
            mask = np.zeros(n, bool)
            return mask

    def limitless(self, n):
        """Do not kill any individuals."""
        mask = np.zeros(n, bool)
        return mask

    def treadmill_random(self, n):
        """Kill some extra individuals randomly."""
        indices = np.random.choice(n, n - aux.max_population_size, replace=False)
        mask = np.zeros(n, bool)
        mask[indices] = True
        return mask

    # TODO: reduce to treadmill_random by passing CONF.bottleneck_size to aux function
    def bottleneck(self, n):
        """Kill all but random few."""
        indices = np.random.choice(n, aux.bottleneck_size, replace=False)
        mask = np.ones(n, bool)
        mask[indices] = False
        return mask

    def treadmill_ageist(self, n):
        """Kill the oldest."""
        # kill the population tail
        mask = np.ones(n, bool)
        mask[: aux.max_population_size] = False
        return mask

    def treadmill_boomer(self, n):
        "Kill the youngest."
        # kill the population head
        mask = np.ones(n, bool)
        mask[-aux.max_population_size :] = False
        return mask

    def __call__(self, n):
        if n <= aux.max_population_size:
            return np.zeros(n, bool)
        else:
            return self.func(n)
