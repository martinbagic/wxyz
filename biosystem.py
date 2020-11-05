"""
Main method that is called from the outside is the cycle() method.
"""

import numpy as np
import logging

from pop import Pop


class Biosystem:
    def __init__(self, pop=None, popid=None, aux=None):

        self.aux = aux
        self.conf = aux.configs[0]
        self.popid = popid if not popid is None else aux.get_popid()
        self.eggs = None  # for split generations

        def _initialize_genomes():
            # make a genome array with random values
            genomes = np.random.random(
                size=(
                    self.conf.MAX_POPULATION_SIZE,
                    self.conf.total_loci,
                    self.conf.BITS_PER_LOCUS,
                )
            )

            # make a bool map using the initial values from GENOME_STRUCT
            for attr, pos in self.conf.loci_pos.items():
                initial_value = self.conf.GENOME_STRUCT[attr][4]
                genomes[:, pos[0] : pos[1]] = (
                    genomes[:, pos[0] : pos[1]] <= initial_value
                )

            genomes = genomes.astype(bool)

            # guarantee survival and reproduction values up to first few mature ages
            if self.conf.HEADSUP > -1:
                headsup = self.conf.MATURATION_AGE + self.conf.HEADSUP
                surv_start = self.conf.loci_pos["surv"][0]
                repr_start = self.conf.loci_pos["repr"][0]
                genomes[:, surv_start : surv_start + headsup,] = True
                genomes[:, repr_start : repr_start + headsup,] = True

            return genomes

        def _get_pop():
            if pop is None:
                # logging.info("Initialize population attributes")
                num = self.conf.MAX_POPULATION_SIZE

                genomes = _initialize_genomes()
                ages = np.zeros(num, int)
                # ages = np.random.randint(0, self.conf.MAX_LIFESPAN, size=num)
                origins = np.zeros(num, int) - 1
                uids = self.conf.get_uids(num)
                births = np.zeros(num, int)
                birthdays = np.zeros(num, int)
                phenomes = self._get_phenome(genomes)
                return Pop(genomes, ages, origins, uids, births, birthdays, phenomes)
            else:
                return pop

        self.pop = _get_pop()

    ###############
    # CYCLE FUNCS #
    ###############

    def cycle(self):
        """Do one simulation cycle."""

        if len(self.pop) + (0 if self.eggs is None else len(self.eggs)) == 0:
            return

        def _hatch_eggs(self):
            if self.eggs is not None:
                self.pop += self.eggs
                self.eggs = None

        if len(self.pop):
            self.eco_survival()
            self.gen_survival()
            self.reproduction()
            self.age()

        # self.conf.season_countdown -= 1
        self.conf.season.countdown -= 1

        # if self.conf.season_countdown == 0:
        if self.conf.season.countdown == 0:
            # kill all living
            mask_kill = np.ones(len(self.pop), bool)
            self._kill(mask_kill, "season_shift")
            # hatch eggs
            _hatch_eggs(self)
            # restart season
            # self.conf.reset_season_countdown()
            self.conf.season.reset()

        # elif self.conf.season_countdown == float("inf"):
        elif self.conf.season.countdown == float("inf"):
            # add newborns to population
            _hatch_eggs(self)

        # evolve envmap if necessary
        self.conf.envmap.evolve(stage=self.aux.stage)

        # pickle self
        if self.aux.PICKLE_RATE and self.aux.stage % self.aux.PICKLE_RATE == 0:
            self.aux.recorder.pickle_pop(self, self.aux.stage)

    def age(self):
        """Increase age of all by one and kill those that surpass max lifespan."""
        self.pop.ages += 1
        mask_kill = self.pop.ages >= self.conf.MAX_LIFESPAN
        self._kill(mask_kill=mask_kill, causeofdeath="max_lifespan")

    def eco_survival(self):
        """Impose ecological death, i.e. death that arises due to resource scarcity."""
        mask_kill = self.conf.overshoot(n=len(self.pop))
        self._kill(mask_kill=mask_kill, causeofdeath="overshoot")

    def gen_survival(self):
        """Impose genomic death, i.e. death that arises with probability encoded in the genome."""
        # mask_surv = self._get_mask(indices=self.pop.ages, attr="surv",)
        probs_surv = self._get_evaluation("surv")
        mask_surv = np.random.random(len(probs_surv)) < probs_surv
        self._kill(mask_kill=~mask_surv, causeofdeath="genetic")

    def reproduction(self):
        """Make individuals reproduce."""

        def _recombine(genomes):
            # make recombined genomes
            flat_genomes = genomes.reshape(len(genomes), -1)
            chromosomes1 = flat_genomes[:, ::2]
            chromosomes2 = flat_genomes[:, 1::2]

            # make choice array: when to take recombined and when to take original loci
            # -1 means synapse; +1 means clear
            rr = self.conf.RECOMBINATION_RATE / 2
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
            # locus structure on an example with 8 bits
            # [1 2 1 2 1 2 1 2]
            # [x   x   x   x  ] => 1st chromosome
            # [  x   x   x   x] => 2nd chromosome

            # extract parent indices twice, and shuffle
            order = np.repeat(
                np.arange(len(genomes)), 2
            )  # every parent sends off two chromosomes into two children
            np.random.shuffle(order)

            # check for selfing
            selfed = np.where(order[::2] == order[1::2])[0] * 2
            if len(selfed) == 1:
                # switch first selfed chromosome with the first chromosome of the previous or next pair
                offset = -2 if selfed[0] > 0 else 2
                order[selfed], order[selfed + offset] = (
                    order[selfed + offset],
                    order[selfed],
                )
            elif len(selfed) > 1:
                # shift first chromosomes of selfed pairs
                order[selfed] = order[np.roll(selfed, 1)]

            assorted = genomes[order]

            # take first chromosomes of first parents and second chromosomes of second parents
            assorted[::2, :, 1::2] = assorted[1::2, :, 1::2]

            # children have chromosomes at positions [::2]
            assorted = assorted[::2]

            return assorted, order

        def _mutate(genomes):
            muta_prob = self._get_evaluation("muta", part=mask_repr)
            muta_prob = muta_prob[mask_repr]

            # random
            random_probabilities = np.random.random(genomes.shape)

            # broadcast to fit [individual, locus, bit] shape
            mutation_probabilities = muta_prob[:, None, None]

            # create new
            mutate_0to1 = (genomes == 0) & (
                random_probabilities
                < (mutation_probabilities * self.conf.MUTATION_RATIO)
            )
            mutate_1to0 = (genomes == 1) & (
                random_probabilities < mutation_probabilities
            )

            mutate = mutate_0to1 + mutate_1to0

            return np.logical_xor(genomes, mutate)

        # check if mature
        mask_mature = self.pop.ages >= self.conf.MATURATION_AGE
        if not any(mask_mature):
            return

        # check if reproducing
        probs_repr = self._get_evaluation("repr", part=mask_mature)
        mask_repr = np.random.random(len(probs_repr)) < probs_repr
        if sum(mask_repr) < 2:  # forgo if not at least two available parents
            return

        # increase births statistics
        self.pop.births += mask_repr

        # copy genomes of parents and modify
        genomes = self.pop.genomes[mask_repr]
        if self.conf.REPR_MODE == "sexual":
            genomes = _recombine(genomes)
            genomes, order = _assort(genomes)
        genomes = _mutate(genomes)
        # assert not genomes is self.pop.genomes[mask_repr]

        # get origins
        if self.conf.REPR_MODE in ("asexual", "diasexual"):
            origins = self.pop.uids[mask_repr]
        elif self.conf.REPR_MODE == "sexual":
            origins = np.array(
                [
                    f"{self.pop.uids[order[2*i]]}.{self.pop.uids[order[2*i+1]]}"
                    for i in range(len(order) // 2)
                ]
            )

        # get eggs
        n = len(genomes)
        eggs = Pop(
            genomes=genomes,
            ages=np.zeros(n, int),
            origins=origins,
            uids=self.conf.get_uids(n),
            births=np.zeros(n, int),
            birthdays=np.zeros(n, int) + self.aux.stage,
            phenomes=self._get_phenome(genomes),
        )

        # save as eggs if generations are nonoverlapping/discrete
        # otherwise add directly to population
        # if self.conf.DISCRETE_GENERATIONS:
        #     if self.eggs is None:
        #         self.eggs = eggs
        #     else:
        #         self.eggs += eggs
        # else:
        #     self.pop += eggs
        if self.eggs is None:
            self.eggs = eggs
        else:
            self.eggs += eggs

    ################
    # HELPER FUNCS #
    ################

    def _get_phenome(self, genomes):
        def _get_interpretome(omes):
            interpretome = np.zeros(shape=omes.shape[:2])

            for attr, pos in self.conf.loci_pos.items():
                # fetch
                loci = omes[:, pos[0] : pos[1]]
                interpreter = self.conf.GENOME_STRUCT[attr][1]

                # interpret
                probs = self.conf.interpreter(loci, interpreter)  # * (hi - lo) + lo

                # add back
                interpretome[:, pos[0] : pos[1]] += probs

            return interpretome

        def _bound(omes):
            """Impose lower and upper bounds for genetically encodable attributes."""
            for attr, pos in self.conf.loci_pos.items():
                lo, hi = self.conf.GENOME_STRUCT[attr][2:4]
                omes[:, pos[0] : pos[1]] = omes[:, pos[0] : pos[1]] * (hi - lo) + lo
            return omes

        envgenomes = self.conf.envmap(genomes)
        interpretome = _get_interpretome(envgenomes)
        phenomes = self.conf.phenomap(interpretome)
        bounded_phenome = _bound(phenomes)

        # if not self.aux.stage % 100:
        #     print(
        #         np.round(interpretome[0][:10], 2), np.round(bounded_phenome[0][:10], 2)
        #     )

        return bounded_phenome

    def _get_evaluation(self, attr, part=None):

        which_individuals = np.arange(len(self.pop))
        if part is not None:
            which_individuals = which_individuals[part]

        # first scenario
        if attr in self.conf.GENOME_CONST:
            probs = self.conf.GENOME_CONST[attr]

        # second and third scenario
        if attr in self.conf.GENOME_STRUCT:
            agespec = self.conf.GENOME_STRUCT[attr][0]

            which_loci = self.conf.loci_pos[attr][0]
            if agespec:
                which_loci += self.pop.ages[which_individuals]

            probs = self.pop.phenomes[which_individuals, which_loci]

            # envgenomes = self.conf.envmap(self.pop.genomes)
            # loci = envgenomes[which_individuals, which_loci]
            # probs = self.conf.interpreter(loci, interpreter) * (hi - lo) + lo

        # expand values back into an array with shape of whole pop
        final_probs = np.zeros(len(self.pop))
        final_probs[which_individuals] += probs

        return final_probs

    # def _get_locus_phenome(self, attr, part=None):
    #     """
    #     Fetch probability values of certain genes at a certain loci.
    #     Three scenarios:
    #         1) value is not individual-specific nor age-specific (constant)
    #         2) value is individual-specific but not age-specific
    #         3) value is individual-specific and age-specific
    #     Use attr to specify the attribute to fetch.
    #     Use part to limit yourself to a part of the population.
    #     """

    #     # probabilities of which individuals do we calculate
    #     which_individuals = np.arange(len(self.pop))
    #     if not part is None:
    #         which_individuals = which_individuals[part]

    #     # first scenario
    #     if attr in self.conf.GENOME_CONST:
    #         probs = self.conf.GENOME_CONST[attr]

    #     # second and third scenario
    #     if attr in self.conf.GENOME_STRUCT:
    #         agespec, interpreter, lo, hi = self.conf.GENOME_STRUCT[attr][0:4]

    #         which_loci = (
    #             self.pop.ages[which_individuals]
    #             if agespec
    #             else self.conf.loci_pos[attr][0]
    #         )

    #         envgenomes = self.conf.envmap(self.pop.genomes)
    #         loci = envgenomes[which_individuals, which_loci]
    #         probs = self.conf.interpreter(loci, interpreter) * (hi - lo) + lo

    #     # expand values back into an array with shape of whole pop
    #     final_probs = np.zeros(len(self.pop))
    #     final_probs[which_individuals] += probs

    #     return final_probs

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
            mask_kill.nonzero()[0][:: self.aux.REC_EVERY_NTH]
            if self.aux.REC_EVERY_NTH > 1
            else mask_kill
        )
        self.aux.recorder.rec(self.pop[mask_record], causeofdeath, self.popid)

        # retain survivors
        self.pop *= ~mask_kill

    def __len__(self):
        return (
            len(self.pop) + len(self.eggs) if self.eggs is not None else len(self.pop)
        )

