import numpy as np

import funcs

config = funcs.load_config("config.yml")


# make ranges
i0 = 0
i1 = config.max_lifespan
i2 = i1 + config.max_lifespan - config.maturation_age
i3 = i2 + config.n_neutral_loci
i4 = i3 + config.mutation_loci
n_total_loci = i4


powers_of_2 = 2 ** np.arange(config.mutation_loci * config.bits_per_locus)
sum_powers_of_2 = sum(powers_of_2)

starting_genome = np.ones(
    shape=(config.population_size, n_total_loci, config.bits_per_locus), dtype=int
)
starting_genome[:, -5:, :] = 0
starting_genome[:, -1, 0] = 1


class Population:
    def __init__(self, config):

        self.genomes = np.random.randint(
            0, 2, size=(config.population_size, n_total_loci, config.bits_per_locus)
        )
        self.genomes = np.copy(starting_genome)
        self.ages = np.random.randint(
            0, config.maturation_age, size=(config.population_size)
        )

    def is_extinct(self):
        return len(self.genomes) == 0

    def genomelen(self):
        return np.arange(len(self.genomes))

    def survive(self):
        loci = self.genomes[self.genomelen(), self.ages]
        survival_probs = loci.sum(1) / config.bits_per_locus  # * 0.5 + 0.5
        # print(np.median(survival_probs))
        # print(np.random.choice(survival_probs, 5))
        random_probs = np.random.random(size=len(self.genomes))
        survived = random_probs < survival_probs

        # print("len genomes before survival", len(self.genomes))
        self.genomes = self.genomes[survived]
        # print("len genomes after survival", len(self.genomes))
        self.ages = self.ages[survived]

    def reproduce(self):
        loci = self.genomes[self.genomelen(), self.ages - config.maturation_age + i1]
        reproduction_probs = loci.sum(1) / config.bits_per_locus  # * 0.5 + 0
        random_probs = np.random.random(size=len(self.genomes))
        success = (random_probs < reproduction_probs) & (
            self.ages >= config.maturation_age
        )

        def get_new_genomes(parents):

            # print(parents[-1])
            loci = parents[np.arange(len(parents)), i3:i4]
            # print(loci[-2])
            # print(loci.shape)
            loci = loci.reshape(len(parents), (i4 - i3) * config.bits_per_locus)
            # print(loci[-2])
            # print(loci.shape)

            mutation_probs = np.dot(loci, powers_of_2) / sum_powers_of_2
            # mutation_probs = loci.sum(1) / (config.bits_per_locus * config.mutation_loci)
            # print(mutation_probs[-1], self.genomes[-1, -5:], np.mean(mutation_probs))
            print(np.mean(mutation_probs))
            # print(mutation_probs[-2])
            # print(mutation_probs.shape)

            mutation_probs = mutation_probs[:, None, None]
            # print(mutation_probs.shape)

            random_probs = np.random.random(
                size=[len(parents), n_total_loci, config.bits_per_locus]
            )

            mutate_0 = (parents == 0) & (
                random_probs < mutation_probs * config.mutation_ratio
            )
            mutate_1 = (parents == 1) & (random_probs < mutation_probs)
            mutate_no = (~mutate_0) & (~mutate_1)

            # mutated_q = 1 - sum(mutate_no) / (parents.size)
            # print(mutated_q[:10])

            return np.select(
                choicelist=[1, 0, parents], condlist=[mutate_0, mutate_1, mutate_no]
            )

        if sum(success) > 0:
            # generate
            new_genomes = get_new_genomes(parents=self.genomes[success])
            new_ages = np.zeros(success.sum(), dtype=int)
            # print("len new genomes", len(new_genomes))

            # append
            self.genomes = np.append(self.genomes, new_genomes, axis=0)
            # print("len genomes", len(self.genomes))
            self.ages = np.append(self.ages, new_ages, axis=0)

        # kill
        if len(self.genomes) > config.max_population_size:
            print("bottleneck")
            indices = np.random.choice(
                len(self.genomes), config.bottleneck_size
            )  # bottleneck
            self.genomes = self.genomes[indices]
            self.ages = self.ages[indices]

    def age(self):
        self.ages += 1

        survived = self.ages <= config.max_lifespan
        self.genomes = self.genomes[survived]
        self.ages = self.ages[survived]

    def replenish(self):

        # double survivors

        self.genomes = np.append(self.genomes, self.genomes, axis=0)
        self.ages = np.append(self.ages, self.ages, axis=0)

        # -1 for the while break condition
        size = int((config.max_population_size - len(self.genomes)) * 0.2)

        # generate
        new_genomes = np.random.randint(
            0, 2, size=(size, n_total_loci, config.bits_per_locus)
        )
        # new_genomes = np.ones(
        #     shape=(config.population_size, n_total_loci, config.bits_per_locus), dtype=int
        # )
        new_ages = np.random.randint(0, config.maturation_age, size=(size))

        # append
        self.genomes = np.append(self.genomes, new_genomes, axis=0)
        self.ages = np.append(self.ages, new_ages, axis=0)

    def cycle(self):
        self.survive()
        self.age()
        self.reproduce()
        # if len(self.genomes) < 100:
        #     self.replenish()
