import numpy as np
from funcs import aux


class Population:

    attrs = ("genomes", "ages", "origins", "uids", "births", "birthdays")

    def __init__(
        self,
        genomes=None,
        ages=None,
        origins=None,
        uids=None,
        births=None,
        birthdays=None,
    ):
        def _get_genomes():

            # set values throughout genome
            genomes = (
                np.random.random(size=aux.genome_shape) <= aux.genome_distribution
            ).astype(int)

            # set values for mutation loci separately
            genomes[:, -aux.n_mutation_loci :] = (
                np.random.random(
                    size=(
                        aux.max_population_size,
                        aux.n_mutation_loci,
                        aux.bits_per_locus,
                    )
                )
                <= aux.mutrate_distribution
            ).astype(int)

            return genomes

        # pop attributes
        self.genomes = genomes
        self.ages = ages
        self.origins = origins
        self.uids = uids
        self.births = births
        self.birthdays = birthdays

        # split generations
        self.eggs = None

        if any(getattr(self, attr) is None for attr in self.attrs):
            num = aux.max_population_size
            self.genomes = _get_genomes()
            self.ages = np.zeros(num, int)
            self.origins = np.zeros(num, int) - 1
            self.uids = self._get_uids(num)
            self.births = np.zeros(num, int)
            self.birthdays = np.zeros(num, int)

    ###############
    # CYCLE FUNCS #
    ###############

    def age(self):
        """Increase age of all by one and kill those that surpass max lifespan."""
        self.ages += 1
        mask_kill = self.ages >= aux.max_lifespan
        self._kill(mask_kill=mask_kill, causeofdeath="max_lifespan")

    def eco_survival(self):
        # TODO: implement cyclicity
        # TODO: implement starvation

        # resources updated if overshoot is starvation
        mask_kill = aux.overshoot(n=len(self))
        self._kill(mask_kill=mask_kill, causeofdeath="overshoot")

    def gen_survival(self):
        def _get_mask():
            indices = self.ages
            loci = self.genomes[np.arange(len(self)), indices]
            if not aux.envmap is None:
                envloci = aux.envmap[indices]
                loci = np.logical_xor(envloci, loci)
            survival_probabilities = self._get_survival_probabilities(loci)
            random_probabilities = np.random.random(len(self))
            return random_probabilities < survival_probabilities

        mask_surv = _get_mask()
        self._kill(mask_kill=np.invert(mask_surv), causeofdeath="genetic")

    def reproduction(self):
        def _get_mask():
            indices = self.ages - aux.maturation_age + aux.max_lifespan
            # in some cases the loci will be wrong because the individuals are immature but they will be filtered after
            loci = self.genomes[np.arange(len(self)), indices]
            if not aux.envmap is None:
                envloci = aux.envmap[indices]
                loci = np.logical_xor(envloci, loci)
            reproduction_probabilities = self._get_reproduction_probabilities(loci)
            random_probabilities = np.random.random(len(self))
            return (random_probabilities < reproduction_probabilities) & (
                self.ages >= aux.maturation_age  # filtering; reproduce only if mature
            )

        def _recombine(genomes):
            # make recombined genomes
            rgs = np.zeros(genomes.shape, dtype=int)
            rgs[:, :, ::2] = genomes[:, :, 1::2]
            rgs[:, :, 1::2] = genomes[:, :, ::2]

            # make choice array: when to take recombined and when to take original loci
            reco_fwd = (
                np.random.random(size=genomes.shape[:-1]) < aux.recombination_rate / 2
            ) * (-2) + 1
            reco_bkd = (
                np.random.random(size=genomes.shape[:-1]) < aux.recombination_rate / 2
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

        def _assort(genomes):
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

        def _mutate(genomes):

            if aux.mutation_probability == -1:  # if mutation evolvable
                loci = genomes[np.arange(len(genomes)), aux.pos_muta :]
                loci = (loci == 0).sum(1) / aux.n_mutation_loci
                mutation_probabilities = self._get_mutation_probabilities(loci)
            else:
                mutation_probabilities = (
                    np.ones(len(genomes)) * aux.mutation_probability
                )

            # random
            random_probabilities = np.random.random(genomes.shape)

            # broadcast
            mutation_probabilities = mutation_probabilities[:, None, None]

            # create new
            mutate_0to1 = (genomes == 0) & (
                random_probabilities < (mutation_probabilities * aux.mutation_ratio)
            )
            mutate_1to0 = (genomes == 1) & (
                random_probabilities < mutation_probabilities
            )

            mutate = mutate_0to1 + mutate_1to0

            return np.logical_xor(genomes, mutate)

            # mutate_no = (~mutate_0to1) & (~mutate_1to0)

            # return np.select(
            #     choicelist=[1, 0, genomes],
            #     condlist=[mutate_0to1, mutate_1to0, mutate_no],
            # )

        mask_repr = _get_mask()
        n = mask_repr.sum()

        if not n:  # if nobody reproduces
            return

        # TODO: it can happen that one parent won't find a pair
        self.births += mask_repr

        # get genomes
        genomes = self.genomes[mask_repr]
        if aux.repr_mode == "sexual":
            genomes = _recombine(genomes)
            genomes = _assort(genomes)
        genomes = _mutate(genomes)

        # get origins
        if aux.repr_mode == "asexual":
            origins = self.uids[mask_repr]
        elif aux.repr_mode == "sexual":
            # TODO: set dual origin
            origins = np.zeros(n, int)

        # get eggs
        eggs = Population(
            genomes=genomes,
            ages=np.zeros(n, int),
            origins=origins,
            uids=self._get_uids(n),
            births=np.zeros(n, int),
            birthdays=np.zeros(n, int) + aux.stage,
        )

        # save as eggs if generations are nonoverlapping
        if aux.discrete_generations:
            if self.eggs is None:
                self.eggs = eggs
            else:
                self.eggs += eggs
        else:
            self += eggs

    def __call__(self):

        if len(self):
            # TODO rethink ordering
            self.gen_survival()
            self.eco_survival()
            self.reproduction()
            self.age()

        if aux.season_countdown == 0:
            self._season_shift()
            aux.reset_season_countdown()

    ################
    # HELPER FUNCS #
    ################

    # def _apply_envmap(self, )

    def _get_survival_probabilities(self, loci):
        return (
            aux.phenotype(loci, "surv") * (aux.repr_bound_hi - aux.repr_bound_lo)
            + aux.repr_bound_lo
        )

    def _get_reproduction_probabilities(self, loci):
        return (
            aux.phenotype(loci, "repr") * (aux.repr_bound_hi - aux.repr_bound_lo)
            + aux.repr_bound_lo
        )

    def _get_mutation_probabilities(self, loci):
        return (
            aux.phenotype(loci, "muta") * (aux.muta_bound_hi - aux.muta_bound_lo)
            + aux.muta_bound_lo
        )

    def _kill(self, mask_kill, causeofdeath):

        if not any(mask_kill):
            return

        # record (some or all) of killed individuals
        mask_record = (
            mask_kill.nonzero()[0][:: aux.rec_every_nth]
            if aux.rec_every_nth > 1
            else mask_kill
        )
        # TODO recorder
        aux.recorder.rec(self[mask_record], causeofdeath)

        # retain survivors
        mask_survive = np.invert(mask_kill)
        self *= mask_survive

    def _get_uids(self, n):
        uids = np.arange(n) + aux.max_uid
        aux.max_uid += n
        return uids

    def _season_shift(self):
        # kill all living
        mask_kill = np.ones(len(self), bool)
        self._kill(mask_kill, "season_shift")

        # hatch eggs
        if not self.eggs is None:  # if there are eggs
            self += self.eggs
            self.eggs = None

    ###############
    # MAGIC FUNCS #
    ###############

    def __len__(self):
        n = len(self.genomes)
        assert all(len(getattr(self, attr)) == n for attr in self.attrs)
        return n

    def __getitem__(self, index):
        return Population(
            genomes=self.genomes[index],
            ages=self.ages[index],
            origins=self.origins[index],
            uids=self.uids[index],
            births=self.births[index],
            birthdays=self.birthdays[index],
        )

    def __imul__(self, index):
        for attr in self.attrs:
            setattr(self, attr, getattr(self, attr)[index])
        return self

    def __iadd__(self, pop):
        for attr in self.attrs:
            val = np.concatenate([getattr(self, attr), getattr(pop, attr)])
            setattr(self, attr, val)
        return self
