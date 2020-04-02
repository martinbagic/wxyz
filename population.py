import numpy as np
from record import Record


class Population:
    def __init__(self, identifier, conf):
        global CONFIG
        CONFIG = conf

        def set_genomes():
            init_popsize = int(CONFIG.population_size_q * CONFIG.max_population_size)
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

        # genome cutoffs and other interpretation help
        self.i0 = 0
        self.i1 = CONFIG.max_lifespan
        self.i2 = self.i1 + CONFIG.max_lifespan - CONFIG.maturation_age
        self.i3 = self.i2 + CONFIG.n_neutral_loci
        self.i4 = self.i3 + CONFIG.mutation_loci
        self.n_total_loci = self.i4
        self.powers_of_2 = 2 ** np.arange(CONFIG.mutation_loci * CONFIG.bits_per_locus)
        self.sum_powers_of_2 = sum(self.powers_of_2)

        # initialize genomes and metadata
        set_genomes()
        n = len(self.genomes)
        self.ages = np.zeros(n, dtype=int)
        self.origins = np.arange(n)
        self.births = np.zeros(n, dtype=int)
        self.birthdays = np.zeros(n, dtype=int)

        # initialize record
        self.record = Record(identifier, conf._asdict())

        # initialize time
        self.stage = 0

    def is_extinct(self):
        return len(self.genomes) == 0

    def genomelen(self):
        return np.arange(len(self.genomes))

    def survive(self):
        def get_mask():
            loci = self.genomes[self.genomelen(), self.ages]
            survival_probs = (
                loci.sum(1)
                / CONFIG.bits_per_locus
                * (CONFIG.surv_bound_hi - CONFIG.surv_bound_lo)
                + CONFIG.surv_bound_lo
            )
            random_probs = np.random.random(size=len(self.genomes))
            survived = random_probs < survival_probs
            return survived

        mask = get_mask()
        self.kill(np.invert(mask))

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
        def get_mask():
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

        mask = get_mask()

        if sum(mask) == 0:
            return

        # generate
        n = mask.sum()
        new_genomes = get_new_genomes(parents=self.genomes[mask])
        new_ages = np.zeros(n, dtype=int)
        new_origins = self.origins[mask]
        new_births = np.zeros(n, dtype=int)
        new_birthdays = np.zeros(n, dtype=int) + self.stage

        # append
        self.genomes = np.append(self.genomes, new_genomes, axis=0)
        self.ages = np.append(self.ages, new_ages, axis=0)
        self.origins = np.append(self.origins, new_origins, axis=0)
        self.births = np.append(self.births, new_births, axis=0)
        self.birthdays = np.append(self.birthdays, new_birthdays, axis=0)

    def age(self):
        self.ages += 1
        mask = self.ages <= CONFIG.max_lifespan
        self.kill(np.invert(mask))

    def handle_overflow(self):
        def bottleneck():
            boolmask = np.random.choice(len(self.genomes), CONFIG.bottleneck_size)
            return boolmask

        def treadmill_real():
            boolmask = np.zeros(shape=len(self.genomes), dtype=bool)
            boolmask[: -CONFIG.max_population_size] = True
            return boolmask

        if len(self.genomes) > CONFIG.max_population_size:
            funcs = {
                "bottleneck": bottleneck,
                "treadmill_real": treadmill_real,
            }
            func = funcs[CONFIG.overflow_handling]
            boolmask = func()
            self.kill(np.invert(boolmask))

    def kill(self, boolmask, dying=True):
        # print(boolmask.shape)

        attrs = {"genomes", "ages", "births", "birthdays", "origins"}

        # record killed
        for attr in attrs - {"genomes", "ages"}:
            vals = getattr(self, attr)[boolmask]
            self.record.__dict__[attr].extend(vals.tolist())

        if dying:
            vals = self.ages[boolmask]
            self.record.ages.extend(vals)
        else:
            self.record.ages.extend([-1] * sum(boolmask))

        # split and record genomes
        genomes = self.genomes[boolmask]
        d = {
            "mutrates": self._get_mutation_probs(genomes),
            "survloci": genomes[:, : self.i1].sum(2),
            "reprloci": genomes[:, self.i1 : self.i2].sum(2),
            "neutloci": genomes[:, self.i2 : self.i3].sum(2),
        }
        for k, v in d.items():
            self.record.genomes[k].append(v)

        # remove killed
        for attr in attrs:
            new = getattr(self, attr)[np.invert(boolmask)]
            setattr(self, attr, new)

    def killall(self):
        boolmask = np.ones(shape=len(self.genomes), dtype=bool)
        self.kill(boolmask, dying=False)

    def cycle(self):
        self.stage += 1
        if len(self.genomes) > 0:
            self.survive()
            self.age()
            self.reproduce()
            self.handle_overflow()
            # if len(self.genomes) > 0:
            #     self.record(self)
