import numpy as np
from xhelp.interpreter import Interpreter
from xhelp.phenomap import Phenomap, PhenomapFake
from xhelp.overshoot import Overshoot
from xhelp.envmap import Envmap, EnvmapFake
from xhelp.season import Season


class Config:
    def __init__(self, params):

        # population
        self.MAX_LIFESPAN = params["MAX_LIFESPAN"]
        self.MATURATION_AGE = params["MATURATION_AGE"]
        self.BITS_PER_LOCUS = params["BITS_PER_LOCUS"]
        self.HEADSUP = params["HEADSUP"]
        self.GENOME_STRUCT = params["GENOME_STRUCT"]
        self.GENOME_CONST = params["GENOME_CONST"]
        self.REPR_MODE = params["REPR_MODE"]
        self.MUTATION_RATIO = params["MUTATION_RATIO"]
        self.RECOMBINATION_RATE = params["RECOMBINATION_RATE"]

        # ecology
        self.MAX_POPULATION_SIZE = params["MAX_POPULATION_SIZE"]
        self.OVERSHOOT_EVENT = params["OVERSHOOT_EVENT"]
        self.DISCRETE_GENERATIONS = params["DISCRETE_GENERATIONS"]
        self.ENVMAP_RATE = params["ENVMAP_RATE"]
        self.PHENOMAP_PLUS = params["PHENOMAP_PLUS"]
        self.CLIFF_SURVIVORSHIP = params["CLIFF_SURVIVORSHIP"]

        self.loci_n = {
            attr: [1, self.MAX_LIFESPAN][vals[0]]
            for attr, vals in self.GENOME_STRUCT.items()
        }
        loci_pos = [0] + np.cumsum(list(self.loci_n.values())).tolist()
        self.loci_pos = {
            attr: (loci_pos[i], loci_pos[i + 1])
            for i, attr in enumerate(self.GENOME_STRUCT)
        }
        self.total_loci = loci_pos[-1]

        # changing
        # self.season_countdown = float("inf")
        # self.reset_season_countdown()
        self.max_uid = 0

        # seasonality
        self.season = Season(DISCRETE_GENERATIONS=self.DISCRETE_GENERATIONS)

        # interpreter: genotype -> probabilities
        self.interpreter = Interpreter(
            BITS_PER_LOCUS=self.BITS_PER_LOCUS,
            REPR_MODE=self.REPR_MODE,
        )

        # phenomap: pleiotropic genotype -> age-linear genotype
        self.phenomap = (
            Phenomap(PHENOMAP_PLUS=self.PHENOMAP_PLUS, pos_end=self.total_loci)
            if self.PHENOMAP_PLUS != []
            else PhenomapFake()
        )

        # overshoot: detects overshoots and returns killing masks
        self.overshoot = Overshoot(
            OVERSHOOT_EVENT=self.OVERSHOOT_EVENT,
            MAX_POPULATION_SIZE=self.MAX_POPULATION_SIZE,
            CLIFF_SURVIVORSHIP=self.CLIFF_SURVIVORSHIP,
        )

        # envmap: genotype -> genotype contextualized in a changing environment
        self.envmap = (
            Envmap(
                BITS_PER_LOCUS=self.BITS_PER_LOCUS,
                total_loci=self.total_loci,
                ENVMAP_RATE=self.ENVMAP_RATE,
            )
            if self.ENVMAP_RATE > 0
            else EnvmapFake()
        )

    # def reset_season_countdown(self):
    #     self.season_countdown = (
    #         self.DISCRETE_GENERATIONS if self.DISCRETE_GENERATIONS else float("inf")
    #     )

    def get_uids(self, n):
        """Get an array of unique origin identifiers."""
        uids = np.arange(n) + self.max_uid
        self.max_uid += n
        return uids
