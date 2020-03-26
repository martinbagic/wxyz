import numpy as np
from record import Record

# starting_genome = np.ones(
#     shape=(CONFIG.population_size, n_total_loci, CONFIG.bits_per_locus), dtype=int
# )
# starting_genome[:, -5:, :] = 0
# starting_genome[:, -1, 0] = 1


class Population:
    def __init__(self, identifier, conf):
        global CONFIG
        CONFIG = conf

        init_popsize = int(CONFIG.population_size_q * CONFIG.max_population_size)

        def set_genomes():
            genomes = (
                np.random.random(
                    size=(init_popsize, self.n_total_loci, CONFIG.bits_per_locus)
                )
                <= CONFIG.genome_distribution
            ).astype(int)

            genomes[:, self.i3 : self.i4] = (
                np.random.random(
                    size=(init_popsize, self.i4 - self.i3, CONFIG.bits_per_locus)
                )
                <= CONFIG.mutrate_distribution
            ).astype(int)

            self.genomes = genomes

        def set_ages():
            self.ages = np.zeros(init_popsize, dtype=int)

        def set_ranges():
            self.i0 = 0
            self.i1 = CONFIG.max_lifespan
            self.i2 = self.i1 + CONFIG.max_lifespan - CONFIG.maturation_age
            self.i3 = self.i2 + CONFIG.n_neutral_loci
            self.i4 = self.i3 + CONFIG.mutation_loci
            self.n_total_loci = self.i4

        self.powers_of_2 = 2 ** np.arange(CONFIG.mutation_loci * CONFIG.bits_per_locus)
        self.sum_powers_of_2 = sum(self.powers_of_2)

        set_ranges()
        set_genomes()
        set_ages()

        self.record = Record(identifier, conf._asdict())

    def is_extinct(self):
        return len(self.genomes) == 0

    def genomelen(self):
        return np.arange(len(self.genomes))

    def survive(self):
        loci = self.genomes[self.genomelen(), self.ages]
        survival_probs = (
            loci.sum(1)
            / CONFIG.bits_per_locus
            * (CONFIG.surv_bound_hi - CONFIG.surv_bound_lo)
            + CONFIG.surv_bound_lo
        )
        random_probs = np.random.random(size=len(self.genomes))
        survived = random_probs < survival_probs

        self.genomes = self.genomes[survived]
        self.ages = self.ages[survived]

    def _get_mutation_probs(self, genomes):
        # extract loci
        loci = genomes[np.arange(len(genomes)), self.i3 : self.i4]

        # flatten
        loci = loci.reshape(len(genomes), (self.i4 - self.i3) * CONFIG.bits_per_locus)

        # interpret
        if CONFIG.mutation_locus_interpreter == "exp":
            mutation_probs = 1 / 2 ** np.sum(loci == 0, axis=1)
        elif CONFIG.mutation_locus_interpreter == "binary":
            mutation_probs = np.dot(loci, self.powers_of_2) / self.sum_powers_of_2

        return mutation_probs

    def reproduce(self):
        def get_success():
            loci = self.genomes[
                self.genomelen(), self.ages - CONFIG.maturation_age + self.i1
            ]
            reproduction_probs = (
                loci.sum(1)
                / CONFIG.bits_per_locus
                * (CONFIG.repr_bound_hi - CONFIG.repr_bound_lo)
                + CONFIG.repr_bound_lo
            )
            random_probs = np.random.random(size=len(self.genomes))
            success = (random_probs < reproduction_probs) & (
                self.ages >= CONFIG.maturation_age
            )
            return success

        def get_new_genomes(parents):

            # get
            mutation_probs = self._get_mutation_probs(parents)

            # broadcast
            mutation_probs = mutation_probs[:, None, None]

            # generate random
            random_probs = np.random.random(
                size=[len(parents), self.n_total_loci, CONFIG.bits_per_locus]
            )

            # condlist
            mutate_0 = (parents == 0) & (
                random_probs < mutation_probs * CONFIG.mutation_ratio
            )
            mutate_1 = (parents == 1) & (random_probs < mutation_probs)
            mutate_no = (~mutate_0) & (~mutate_1)

            # select
            return np.select(
                choicelist=[1, 0, parents], condlist=[mutate_0, mutate_1, mutate_no]
            )

        success = get_success()

        if sum(success) == 0:
            return

        # generate
        new_genomes = get_new_genomes(parents=self.genomes[success])
        new_ages = np.zeros(success.sum(), dtype=int)

        # append
        self.genomes = np.append(self.genomes, new_genomes, axis=0)
        self.ages = np.append(self.ages, new_ages, axis=0)

    def age(self):
        self.ages += 1

        survived = self.ages <= CONFIG.max_lifespan
        self.genomes = self.genomes[survived]
        self.ages = self.ages[survived]

    def handle_overflow(self):
        def bottleneck():
            indices = np.random.choice(
                len(self.genomes), CONFIG.bottleneck_size
            )  # bottleneck
            self.genomes = self.genomes[indices]
            self.ages = self.ages[indices]

        def treadmill_real():
            self.genomes = self.genomes[: -CONFIG.max_population_size]
            self.ages = self.ages[: -CONFIG.max_population_size]

        def treadmill_random():
            pass

        if len(self.genomes) > CONFIG.max_population_size:
            {
                "bottleneck": bottleneck,
                "treadmill_real": treadmill_real,
                "treadmill_random": treadmill_random,
            }[CONFIG.overflow_handling]()

    def cycle(self):
        if len(self.genomes) > 0:
            self.survive()
            self.age()
            self.reproduce()
            self.handle_overflow()
            if len(self.genomes) > 0:
                self.record(self)
