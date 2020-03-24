import numpy as np
from config import CONFIG

powers_of_2 = 2 ** np.arange(CONFIG.mutation_loci * CONFIG.bits_per_locus)
sum_powers_of_2 = sum(powers_of_2)


# starting_genome = np.ones(
#     shape=(CONFIG.population_size, n_total_loci, CONFIG.bits_per_locus), dtype=int
# )
# starting_genome[:, -5:, :] = 0
# starting_genome[:, -1, 0] = 1


class Population:
    i0 = 0
    i1 = CONFIG.max_lifespan
    i2 = i1 + CONFIG.max_lifespan - CONFIG.maturation_age
    i3 = i2 + CONFIG.n_neutral_loci
    i4 = i3 + CONFIG.mutation_loci
    n_total_loci = i4

    def __init__(self):
        def set_genomes():
            self.genomes = (
                np.random.random(
                    size=(CONFIG.population_size, n_total_loci, CONFIG.bits_per_locus)
                )
                <= CONFIG.genome_distribution
            ).astype(int)

        def set_ages():
            # self.ages = np.random.randint(0, 1, size=(CONFIG.population_size))
            self.ages = np.zeros(CONFIG.population_size, dtype=int)

        set_genomes()
        set_ages()

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

    def reproduce(self):
        def get_success():
            loci = self.genomes[
                self.genomelen(), self.ages - CONFIG.maturation_age + Population.i1
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

            # extract loci
            loci = parents[np.arange(len(parents)), i3:i4]

            # flatten
            loci = loci.reshape(len(parents), (i4 - i3) * CONFIG.bits_per_locus)

            # interpret
            if CONFIG.mutation_locus_interpreter == "exp":
                mutation_probs = 1 / 2 ** np.sum(loci, axis=1)
            elif CONFIG.mutation_locus_interpreter == "binary":
                mutation_probs = np.dot(loci, powers_of_2) / sum_powers_of_2

            # broadcast
            mutation_probs = mutation_probs[:, None, None]

            # generate random
            random_probs = np.random.random(
                size=[len(parents), n_total_loci, CONFIG.bits_per_locus]
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
        if len(self.genomes) <= CONFIG.max_population_size:
            return

        def bottleneck():
            indices = np.random.choice(
                len(self.genomes), CONFIG.bottleneck_size
            )  # bottleneck
            self.genomes = self.genomes[indices]
            self.ages = self.ages[indices]

        def treadmill_real():
            pass

        def treadmill_random():
            pass


    # def replenish(self):

    #     # double survivors

    #     self.genomes = np.append(self.genomes, self.genomes, axis=0)
    #     self.ages = np.append(self.ages, self.ages, axis=0)

    #     # -1 for the while break condition
    #     size = int((CONFIG.max_population_size - len(self.genomes)) * 0.2)

    #     # generate
    #     new_genomes = np.random.randint(
    #         0, 2, size=(size, n_total_loci, CONFIG.bits_per_locus)
    #     )
    #     # new_genomes = np.ones(
    #     #     shape=(CONFIG.population_size, n_total_loci, CONFIG.bits_per_locus), dtype=int
    #     # )
    #     new_ages = np.random.randint(0, CONFIG.maturation_age, size=(size))

    #     # append
    #     self.genomes = np.append(self.genomes, new_genomes, axis=0)
    #     self.ages = np.append(self.ages, new_ages, axis=0)

    def cycle(self):
        self.survive()
        self.age()
        self.reproduce()
        self.handle_overflow()
