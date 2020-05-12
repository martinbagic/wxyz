import numpy as np
from record import Record
from collections import Counter


class Nextgen:
    def __init__(self, genomes_shape):
        self.genomes = np.empty((0, genomes_shape[1], genomes_shape[2]), dtype=int)
        self.ages = np.empty(0, dtype=int)
        self.origins = np.empty(0, dtype=int)
        self.uids = np.empty(0, dtype=int)
        self.births = np.empty(0, dtype=int)
        self.birthdays = np.empty(0, dtype=int)


class Population:
    def __init__(self, identifier, conf, opath):
        global CONFIG
        CONFIG = conf

        self.opath = opath

        def set_genomes():
            init_popsize = int(CONFIG.population_size_q * CONFIG.max_population_size)
            self.genome_shape = (init_popsize, self.n_total_loci, CONFIG.bits_per_locus)
            genomes = (
                np.random.random(size=self.genome_shape) <= CONFIG.genome_distribution
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
        self.max_uid = 0

        self.ages = np.zeros(n, dtype=int)
        self.origins = np.zeros(n, dtype=int) - 1
        self.uids = self.get_uids(n)
        self.births = np.zeros(n, dtype=int)
        self.birthdays = np.zeros(n, dtype=int)

        # initialize record
        self.record = Record(identifier, conf._asdict())

        # initialize time
        self.stage = 0

        # initialize nextgen if generations are not overlapping
        if CONFIG.split_generations:
            self.nextgen = Nextgen(self.genome_shape)

    def is_extinct(self):
        return len(self.genomes) == 0

    def genomelen(self):
        return np.arange(len(self.genomes))

    def get_uids(self, n):
        uids = np.arange(n) + self.max_uid
        self.max_uid += n
        return uids

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
        self.kill(np.invert(mask), "genetic")

    def _get_mutation_probs(self, genomes):

        if CONFIG.mutation_probability != -1:
            return np.ones(len(genomes)) * CONFIG.mutation_probability

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
        self.births += mask

        if sum(mask) == 0:
            return

        # generate
        n = mask.sum()
        new_genomes = get_new_genomes(parents=self.genomes[mask])
        new_ages = np.zeros(n, dtype=int)
        new_origins = self.uids[mask].copy()
        new_uids = self.get_uids(n)
        new_births = np.zeros(n, dtype=int)
        new_birthdays = np.zeros(n, dtype=int) + self.stage

        # append

        if CONFIG.split_generations:
            obj = self.nextgen
        else:
            obj = self

        obj.genomes = np.append(obj.genomes, new_genomes, axis=0)
        obj.ages = np.append(obj.ages, new_ages, axis=0)
        obj.origins = np.append(obj.origins, new_origins, axis=0)
        obj.uids = np.append(obj.uids, new_uids, axis=0)
        obj.births = np.append(obj.births, new_births, axis=0)
        obj.birthdays = np.append(obj.birthdays, new_birthdays, axis=0)

    def age(self):
        self.ages += 1
        mask = self.ages <= CONFIG.max_lifespan
        self.kill(np.invert(mask), "max_lifespan")

    def handle_overflow(self):
        def bottleneck():
            # kill all but chosen few
            indices = np.random.choice(
                len(self.genomes), CONFIG.bottleneck_size, replace=False
            )
            boolmask = np.ones(shape=len(self.genomes), dtype=bool)
            boolmask[indices] = False
            return boolmask

        def treadmill_ageist():
            # kill the population tail
            boolmask = np.ones(shape=len(self.genomes), dtype=bool)
            boolmask[: CONFIG.max_population_size] = False
            return boolmask

        def treadmill_boomer():
            # kill the population head
            boolmask = np.ones(shape=len(self.genomes), dtype=bool)
            boolmask[-CONFIG.max_population_size :] = False
            return boolmask

        def treadmill_youngandold():
            # kill the population head and tail
            boolmask = np.zeros(shape=len(self.genomes), dtype=bool)
            killn = (len(self.genomes) - CONFIG.max_population_size) // 2 + 1
            boolmask[-killn:] = True
            boolmask[:killn] = True
            return boolmask

        def treadmill_adults():
            # kill the population middle
            boolmask = np.ones(shape=len(self.genomes), dtype=bool)
            killn = CONFIG.max_population_size // 2
            boolmask[-killn:] = False
            boolmask[:killn] = False
            return boolmask

        def treadmill_random():
            # kill chosen few
            indices = np.random.choice(
                len(self.genomes),
                len(self.genomes) - CONFIG.max_population_size,
                replace=False,
            )
            boolmask = np.zeros(shape=len(self.genomes), dtype=bool)
            boolmask[indices] = True
            return boolmask

        if len(self.genomes) > CONFIG.max_population_size:
            funcs = {
                "bottleneck": bottleneck,
                "treadmill_ageist": treadmill_ageist,
                "treadmill_random": treadmill_random,
                "treadmill_boomer": treadmill_boomer,
                "treadmill_adults": treadmill_adults,
                "treadmill_youngandold": treadmill_youngandold,
            }
            func = funcs[CONFIG.overflow_handling]
            boolmask = func()
            self.kill(boolmask, "overflow")

    def kill(self, boolmask, causeofdeath):
        # print(boolmask.shape)

        attrs = ["genomes", "ages", "births", "birthdays", "origins", "uids"]

        killmask = np.copy(boolmask)

        if not CONFIG.record_immature:
            boolmask = (boolmask) & (self.ages >= CONFIG.maturation_age)

        # record killed
        for attr in attrs[2:]:
            vals = getattr(self, attr)[boolmask]
            self.record.__dict__[attr].extend(vals.tolist())

        if causeofdeath == "sim_end":
            self.record.ages.extend([-1] * sum(boolmask))
        else:
            vals = self.ages[boolmask]
            self.record.ages.extend(vals)

        self.record.causeofdeath.extend([causeofdeath] * sum(boolmask))

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

        self.record.flush(self.opath)

        # remove killed
        for attr in attrs:
            new = getattr(self, attr)[np.invert(killmask)]
            setattr(self, attr, new)

    def killall(self):
        boolmask = np.ones(shape=len(self.genomes), dtype=bool)
        self.kill(boolmask, "sim_end")

    def cycle(self):
        self.stage += 1
        if len(self.genomes) > 0:
            self.survive()
            self.age()
            self.reproduce()
            # print(self.stage, len(self.genomes),end=" ")
            self.handle_overflow()
            # print(len(self.genomes))

            if CONFIG.split_generations and self.stage % 25 == 0 and self.stage > 0:
                self.killall()
                self.bring_nextgen()
        else:
            if CONFIG.split_generations and len(self.nextgen.genomes) > 0:
                self.bring_nextgen()

    def bring_nextgen(self):
        # print(self.stage)
        for attr in ("genomes", "ages", "origins", "uids", "births", "birthdays"):
            val = getattr(self.nextgen, attr)
            setattr(self, attr, val)

        self.nextgen = Nextgen(self.genome_shape)
