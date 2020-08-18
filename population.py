import numpy as np
from record import Record


class Nextgen:
    def __init__(self, genomes_shape):
        self.genomes = np.empty((0, genomes_shape[1], genomes_shape[2]), dtype=int)
        self.ages = np.empty(0, dtype=int)
        self.origins = np.empty(0, dtype=int)
        self.uids = np.empty(0, dtype=int)
        self.births = np.empty(0, dtype=int)
        self.birthdays = np.empty(0, dtype=int)


class Population:
    def __init__(self, conf, opath):
        global CONFIG
        CONFIG = conf

        def set_genomes():
            init_popsize = int(CONFIG.population_size_q * CONFIG.max_population_size)
            self.genome_shape = (init_popsize, self.n_total_loci, CONFIG.bits_per_locus)

            if isinstance(CONFIG.genome_distribution, (int, float)):
                genomes = (
                    np.random.random(size=self.genome_shape)
                    <= CONFIG.genome_distribution
                ).astype(int)

                genomes[:, self.i3 : self.i4] = (
                    np.random.random(
                        size=(init_popsize, self.i4 - self.i3, CONFIG.bits_per_locus)
                    )
                    <= CONFIG.mutrate_distribution
                ).astype(int)
            else:
                genomes = np.array(
                    [CONFIG.genome_distribution for _ in range(init_popsize)]
                )

            assert genomes.shape == self.genome_shape

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
        self.record = Record(opath)

        # initialize time
        self.stage = 0

        # initialize nextgen if generations are not overlapping
        if CONFIG.split_generations:
            self.nextgen = Nextgen(self.genome_shape)
            self.split_in = CONFIG.split_generations
        else:
            self.split_in = float("inf")

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

        def mutate_genomes(genomes):

            # get
            mutation_probs = self._get_mutation_probs(genomes)

            # broadcast
            mutation_probs = mutation_probs[:, None, None]

            # generate random
            random_probs = np.random.random(
                size=[len(genomes), self.n_total_loci, CONFIG.bits_per_locus]
            )

            # condlist
            mutate_0 = (genomes == 0) & (
                random_probs < mutation_probs * CONFIG.mutation_ratio
            )
            mutate_1 = (genomes == 1) & (random_probs < mutation_probs)
            mutate_no = (~mutate_0) & (~mutate_1)

            # select
            return np.select(
                choicelist=[1, 0, genomes], condlist=[mutate_0, mutate_1, mutate_no]
            )

        def recombine_genomes(genomes):
            # make recombined genomes
            rgs = np.zeros(genomes.shape, dtype=int)
            rgs[:, :, ::2] = genomes[:, :, 1::2]
            rgs[:, :, 1::2] = genomes[:, :, ::2]

            # make choice array: when to take recombined and when to take original loci
            reco_fwd = (
                np.random.random(size=genomes.shape[:-1])
                < CONFIG.recombination_rate / 2
            ) * (-2) + 1
            reco_bkd = (
                np.random.random(size=genomes.shape[:-1])
                < CONFIG.recombination_rate / 2
            ) * (-2) + 1
            reco_bkd = reco_bkd[:, ::-1]

            reco_fwd_cum = np.cumprod(reco_fwd, axis=1)
            reco_bkd_cum = np.cumprod(reco_bkd, axis=1)
            reco_bkd_cum = reco_bkd_cum[::-1]

            reco_final = reco_fwd_cum * reco_bkd_cum
            reco_final = (reco_final * (-1) + 1) // 2

            # choose from original and recombined

            recombined = [
                [
                    [genomes[i], rgs[i]][if_reco][locus]
                    for locus, if_reco in enumerate(reco_finali)
                ]
                for i, reco_finali in enumerate(reco_final)
            ]

            return np.array(recombined)

        def assort_genomes(genomes):
            # shuffle parents
            np.random.shuffle(genomes)

            # pair chromosomes
            # also makes sure that no odd parent
            pairs = zip(genomes[::2], genomes[1::2])
            assorted = []

            # a contributes 1st chromosome, b contributes 2nd chromosome
            for a, b in pairs:
                a[:, 1::2] = b[:, 1::2]
                assorted.append(a)

            return np.array(assorted)

        mask = get_mask()
        self.births += mask  # wrong! it can happen that one parent won't find a pair

        if mask.sum() == 0:
            return

        # generate children genomes
        new_genomes = self.genomes[mask]

        # if reproduction mode is sexual
        if CONFIG.repr_mode == "sexual" and mask.sum() >= 2:
            # recombination
            # z = np.copy(new_genomes)
            new_genomes = recombine_genomes(new_genomes)
            # print(len(np.where(z != new_genomes)[0]))

            # assortment
            new_genomes = assort_genomes(new_genomes)

        # mutation
        new_genomes = mutate_genomes(new_genomes)

        # calc number of new individuals
        n = len(new_genomes)

        # generate data
        new_ages = np.zeros(n, dtype=int)
        new_uids = self.get_uids(n)
        new_births = np.zeros(n, dtype=int)
        new_birthdays = np.zeros(n, dtype=int) + self.stage

        # calculate origins
        if CONFIG.repr_mode == "asexual":
            new_origins = self.uids[mask].copy()
        else:
            new_origins = np.zeros(n, dtype=int)

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
            boolmask[:capacity] = False
            return boolmask

        def treadmill_boomer():
            # kill the population head
            boolmask = np.ones(shape=len(self.genomes), dtype=bool)
            boolmask[-capacity:] = False
            return boolmask

        def treadmill_youngandold():
            # kill the population head and tail
            boolmask = np.zeros(shape=len(self.genomes), dtype=bool)
            killn = (len(self.genomes) - capacity) // 2 + 1
            boolmask[-killn:] = True
            boolmask[:killn] = True
            return boolmask

        def treadmill_adults():
            # kill the population middle
            boolmask = np.ones(shape=len(self.genomes), dtype=bool)
            killn = capacity // 2
            boolmask[-killn:] = False
            boolmask[:killn] = False
            return boolmask

        def treadmill_random():
            # kill chosen few
            indices = np.random.choice(
                len(self.genomes), len(self.genomes) - capacity, replace=False,
            )
            boolmask = np.zeros(shape=len(self.genomes), dtype=bool)
            boolmask[indices] = True
            return boolmask

        capacity = CONFIG.max_population_size
        if CONFIG.cyclicity != False and CONFIG.cyclicity != "False":
            cyclic_stages, cyclic_size = CONFIG.cyclicity.split("|")
            capacity += np.sin(2 * np.pi / int(cyclic_stages) * self.stage) * int(
                cyclic_size
            )
            capacity = int(capacity)

        if len(self.genomes) > capacity:
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

    def kill(self, mask_kill, causeofdeath):
        def record_killed():
            # record data of all individuals that are killed

            mask_record = (
                mask_kill
                if CONFIG.record_immature
                else (mask_kill) & (self.ages >= CONFIG.maturation_age)
            )

            ### record every nth individual
            if CONFIG.rec_every_nth > 1:
                indices = mask_record.nonzero()[0][:: CONFIG.rec_every_nth]
                mask_record = np.zeros(mask_record.shape, dtype=bool)
                mask_record[indices] = True

            if not mask_record.sum():
                return

            # add demography data to record instance
            self.record.sid.extend(self.uids[mask_record])  # self id
            self.record.pid.extend(self.origins[mask_record])  # parental id
            self.record.bday.extend(self.birthdays[mask_record])
            self.record.age.extend(self.ages[mask_record])
            self.record.causeofdeath.extend([causeofdeath] * sum(mask_record))

            # add genomes data to record instance
            genomes = self.genomes[mask_record]
            self.record.add_genomes(genomes)

            # write data
            self.record.record()

        def retain_survivors():
            # retain data of all individuals that are not killed

            mask_alive = np.invert(mask_kill)
            for attr in ["genomes", "ages", "births", "birthdays", "origins", "uids"]:
                data_alive = getattr(self, attr)[mask_alive]
                setattr(self, attr, data_alive)

        if mask_kill.sum():
            record_killed()
            retain_survivors()

    def killall(self):
        boolmask = np.ones(shape=len(self.genomes), dtype=bool)
        self.kill(boolmask, "killall")

    def cycle(self):
        def renew_split_in():
            """Renews split generations counter."""
            if CONFIG.split_generations_soft:
                rand = np.random.normal(
                    loc=CONFIG.split_generations, scale=CONFIG.split_generations_soft,
                )
                self.split_in = int(abs(rand) + 1)
            else:
                self.split_in = CONFIG.split_generations

        self.stage += 1
        self.split_in -= 1
        # print(self.stage, self.split_in, len(self.genomes),len(self.nextgen.genomes))
        if len(self.genomes) > 0:
            self.survive()
            self.age()
            self.reproduce()
            self.handle_overflow()

        if CONFIG.split_generations:
            if self.split_in == 0:
                # print("killall", self.split_in, len(self.genomes), self.stage)
                self.killall()
                self.bring_nextgen()
                renew_split_in()

    def bring_nextgen(self):
        for attr in ("genomes", "ages", "origins", "uids", "births", "birthdays"):
            val = getattr(self.nextgen, attr)
            setattr(self, attr, val)

        self.nextgen = Nextgen(self.genome_shape)
