"""
Main method that is called from the outside is the cycle() method.
"""

import numpy as np
import logging

class Pop:
    def __init__(self, genomes, ages, origins, uids, births, birthdays):
        self.genomes = genomes
        self.ages = ages
        self.origins = origins
        self.uids = uids
        self.births = births
        self.birthdays = birthdays
        
    

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
        #
        popid=None,
        aux=None,
    ):
        
        self.aux = aux
        self.conf = aux.configs[0] 
        
        def _get_genomes():

            # set values throughout genome
            genomes = np.random.random(self.conf.genome_shape) <= self.conf.genome_distribution

            # guarantee survival through prematurity
            # guarantee survival and reproduction of few first mature ages
            if self.conf.age_guarantee > -1:
                genomes[:, : self.conf.maturation_age + self.conf.age_guarantee] = True
                genomes[:, self.conf.pos_repr : self.conf.pos_repr + self.conf.age_guarantee] = True

            # set values for mutation loci separately
            genomes[:, -self.conf.n_mutation_loci :] = (
                np.random.random(
                    size=(
                        self.conf.max_population_size,
                        self.conf.n_mutation_loci,
                        self.conf.bits_per_locus,
                    )
                )
                <= self.conf.mutrate_distribution
            )

            return genomes

        # pop attributes
        self.genomes = genomes
        self.ages = ages
        self.origins = origins
        self.uids = uids
        self.births = births
        self.birthdays = birthdays
        self.popid = popid

        # split generations
        self.eggs = None

        if any(getattr(self, attr) is None for attr in self.attrs):
            logging.info("Initialize population attributes")
            num = self.conf.max_population_size
            self.genomes = _get_genomes()
            self.ages = np.zeros(num, int)
            self.origins = np.zeros(num, int) - 1
            self.uids = self.conf.get_uids(num)
            self.births = np.zeros(num, int)
            self.birthdays = np.zeros(num, int)
            self.popid = self.aux.get_popid()

        self._set_envgenomes()

    ###############
    # CYCLE FUNCS #
    ###############

    def cycle(self):
        """Do one simulation cycle."""

        def _season_shift():
            # kill all living
            mask_kill = np.ones(len(self), bool)
            self._kill(mask_kill, "season_shift")

            # hatch eggs
            if not self.eggs is None:
                self += self.eggs
                self.eggs = None

        if len(self):
            # TODO rethink ordering
            self.gen_survival()
            self.eco_survival()
            self.reproduction()
            self.age()

        # if season over, kill pop and hatch eggs
        if self.conf.season_countdown == 0:
            _season_shift()
            self.conf.reset_season_countdown()
        
        self.conf.season_countdown -= 1
            
        if not self.conf.envmap is None and self.aux.stage % self.conf.envmap_rate == 0:
            self.conf.envmap.evolve()
            self._set_envgenomes()

    def age(self):
        """Increase age of all by one and kill those that surpass max lifespan."""
        self.ages += 1
        mask_kill = self.ages >= self.conf.max_lifespan
        self._kill(mask_kill=mask_kill, causeofdeath="max_lifespan")

    def eco_survival(self):
        """Impose ecological death, i.e. death that arises due to resource scarcity."""
        mask_kill = self.conf.overshoot(n=len(self))
        self._kill(mask_kill=mask_kill, causeofdeath="overshoot")

    def gen_survival(self):
        """Impose genomic death, i.e. death that arises with probability encoded in the genome."""
        mask_surv = self._get_mask(
            indices=self.ages, get_probabilities=self._get_survival_probabilities,
        )
        self._kill(mask_kill=np.invert(mask_surv), causeofdeath="genetic")

    def reproduction(self):
        """Make individuals reproduce."""

        def _recombine(genomes):
            # make recombined genomes
            flat_genomes = genomes.reshape(len(genomes), -1)
            chromosomes1 = flat_genomes[:, ::2]
            chromosomes2 = flat_genomes[:, 1::2]

            # make choice array: when to take recombined and when to take original loci
            # -1 means synapse; +1 means clear
            rr = self.conf.recombination_rate / 2
            reco_fwd = (np.random.random(chromosomes1.shape) < rr) * -2 + 1
            reco_bkd = (np.random.random(chromosomes2.shape) < rr) * -2 + 1

            # propagate synapse
            reco_fwd_cum = np.cumprod(reco_fwd, axis=1)
            reco_bkd_cum = np.cumprod(reco_bkd[:, ::-1], axis=1)[:, ::-1]

            # recombine if both recombining
            reco_final = reco_fwd_cum * reco_bkd_cum
            reco_final = reco_final == -1  # True if both recombining

            # choose bits from first or second chromosome
            recombined = np.empty(flat_genomes.shape, dtype=bool)
            recombined[:, ::2] = np.choose(reco_final, [chromosomes1, chromosomes2])
            recombined[:, 1::2] = np.choose(reco_final, [chromosomes2, chromosomes1])
            recombined = recombined.reshape(genomes.shape)

            # assure on one example that bits are recombining
            if reco_final[0, 0]:
                assert (
                    recombined[0, 0, 0] == chromosomes2[0, 0]
                    and recombined[0, 0, 1] == chromosomes1[0, 0]
                )
            else:
                assert (
                    recombined[0, 0, 0] == chromosomes1[0, 0]
                    and recombined[0, 0, 1] == chromosomes2[0, 0]
                )

            return recombined

        def _assort(genomes):
            # extract parent indices twice, and shuffle
            # TODO: make sure no selfing
            # TODO: redefine interpreters so they consider chromosomes
            order = list(range(len(genomes))) * 2  # every parent has 2 children
            np.random.shuffle(order)

            # locus structure on an example with 10 bits
            # [1 2 1 2 1 2 1 2 1 2]
            # [x   x   x   x   x  ] => 1st chromosome
            # [  x   x   x   x   x] => 2nd chromosome
            assorted = genomes[order]
            assorted[::2, :, 1::2] = assorted[1::2, :, 1::2]
            assorted = assorted[::2]
            return assorted, order

        def _mutate(genomes):
            if self.conf.mutation_probability == -1:  # if mutation evolvable
                loci = genomes[np.arange(len(genomes)), self.conf.pos_muta :]
                loci = (loci == 0).sum(1) / self.conf.n_mutation_loci
                mutation_probabilities = self._get_mutation_probabilities(loci)
            else:
                mutation_probabilities = (
                    np.ones(len(genomes)) * self.conf.mutation_probability
                )

            # random
            random_probabilities = np.random.random(genomes.shape)

            # broadcast
            mutation_probabilities = mutation_probabilities[:, None, None]

            # create new
            mutate_0to1 = (genomes == 0) & (
                random_probabilities < (mutation_probabilities * self.conf.mutation_ratio)
            )
            mutate_1to0 = (genomes == 1) & (
                random_probabilities < mutation_probabilities
            )

            mutate = mutate_0to1 + mutate_1to0

            return np.logical_xor(genomes, mutate)

        mask_repr = self._get_mask(
            indices=self.ages - self.conf.maturation_age + self.conf.max_lifespan,
            get_probabilities=self._get_reproduction_probabilities,
        )
        mask_repr = mask_repr & (
            self.ages >= self.conf.maturation_age
        )  # filtering; ensure only mature individuals can reproduce

        n = mask_repr.sum()
        if not n:  # if nobody reproduces
            return

        # get genomes
        genomes = self.envgenomes[mask_repr]
        if self.conf.repr_mode == "sexual":
            genomes = _recombine(genomes)
            genomes, order = _assort(genomes)
        genomes = _mutate(genomes)

        # get origins
        if self.conf.repr_mode == "asexual":
            origins = self.uids[mask_repr]
        elif self.conf.repr_mode == "sexual":
            origins = np.array(
                [
                    f"{self.uids[order[2*i]]}-{self.uids[order[2*i+1]]}"
                    for i in range(len(order) // 2)
                ]
            )

        # get eggs
        eggs = Population(
            genomes=genomes,
            ages=np.zeros(n, int),
            origins=origins,
            uids=self.conf.get_uids(n),
            births=np.zeros(n, int),
            birthdays=np.zeros(n, int) + self.aux.stage,
            popid="eggs",
            aux=self.aux,
            
        )

        # save as eggs if generations are nonoverlapping
        # otherwise add directly to population
        if self.conf.discrete_generations:
            if self.eggs is None:
                self.eggs = eggs
            else:
                self.eggs += eggs
        else:
            self += eggs

        self.births += len(eggs)

    ################
    # HELPER FUNCS #
    ################

    def _get_survival_probabilities(self, loci):
        return (
            self.conf.pheno(loci, "surv") * (self.conf.surv_bound_hi - self.conf.surv_bound_lo)
            + self.conf.surv_bound_lo
        )

    def _get_reproduction_probabilities(self, loci):
        return (
            self.conf.pheno(loci, "repr") * (self.conf.repr_bound_hi - self.conf.repr_bound_lo)
            + self.conf.repr_bound_lo
        )

    def _get_mutation_probabilities(self, loci):
        return (
            self.conf.pheno(loci, "muta") * (self.conf.muta_bound_hi - self.conf.muta_bound_lo)
            + self.conf.muta_bound_lo
        )

    def _get_mask(self, indices, get_probabilities):
        """
        Extract loci, interpret and return boolean mask.
        The boolean masks says for every individual whether an event (such as mutation, reproduction or genomic death) takes place or not.
        """

        # extract loci and apply phenomap if necessary
        loci = (
            self.envgenomes[np.arange(len(self)), indices]
            if self.conf.phenomap is None
            else self.conf.phenomap(self.envgenomes, indices)
        )

        return np.random.random(len(self)) < get_probabilities(loci)

    def _kill(self, mask_kill, causeofdeath):
        """
        Kill individuals and record their data.
        Killing can occur due to age, genomic death, ecological death, and season shift.
        """

        # skip if no one to kill
        if not any(mask_kill):
            return

        # record (some or all) of killed individuals
        mask_record = (
            mask_kill.nonzero()[0][:: self.aux.rec_every_nth]
            if self.aux.rec_every_nth > 1
            else mask_kill
        )
        # TODO recorder
        self.aux.recorder.rec(self[mask_record], causeofdeath, self.popid)

        # retain survivors
        mask_survive = np.invert(mask_kill)
        self *= mask_survive


    def _set_envgenomes(self):
        """
        Calculate environmental genomes.
        Environmental genomes are genomes that are interpreted given a certain environmental map.
        """
        self.envgenomes = (
            self.genomes if self.conf.envmap is None else self.conf.envmap(self.genomes)
        )

    ###############
    # MAGIC FUNCS #
    ###############

    def __len__(self):
        """Return the number of living individuals."""
        n = len(self.genomes)
        assert all(len(getattr(self, attr)) == n for attr in self.attrs), " ".join(
            (str(len(getattr(self, attr))) for attr in self.attrs)
        )
        return n

    def __getitem__(self, index):
        """Return a subpopulation."""
        return Population(
            genomes=self.genomes[index],
            ages=self.ages[index],
            origins=self.origins[index],
            uids=self.uids[index],
            births=self.births[index],
            birthdays=self.birthdays[index],
            popid=self.popid,
            aux=self.aux,
        )

    def __imul__(self, index):
        """Redefine itself as its own subpopulation."""
        setattr(self, "envgenomes", getattr(self, "envgenomes")[index])
        for attr in self.attrs:
            setattr(self, attr, getattr(self, attr)[index])
        return self

    def __iadd__(self, pop):
        """Add another population to itself."""
        val = np.concatenate([self.envgenomes, pop.envgenomes])
        setattr(self, "envgenomes", val)
        for attr in self.attrs:
            val = np.concatenate([getattr(self, attr), getattr(pop, attr)])
            setattr(self, attr, val)
        return self
