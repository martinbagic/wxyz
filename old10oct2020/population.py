"""
Main method that is called from the outside is the cycle() method.
"""

import numpy as np
import logging

from xhelp.pop import Pop


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

            # guarantee survival and reproduction values up to first few mature ages
            if self.conf.AGE_GUARANTEE > -1:
                headsup = self.conf.MATURATION_AGE + self.conf.AGE_GUARANTEE
                surv_start = self.conf.loci_pos["surv"][0]
                repr_start = self.conf.loci_pos["repr"][0]
                genomes[:, surv_start : surv_start + headsup,] = True
                genomes[:, repr_start : repr_start + headsup,] = True

            return genomes

        def _get_pop():
            if pop is None:
                logging.info("Initialize population attributes")
                num = self.conf.MAX_POPULATION_SIZE

                genomes = _initialize_genomes()
                ages = np.zeros(num, int)
                origins = np.zeros(num, int) - 1
                uids = self.conf.get_uids(num)
                births = np.zeros(num, int)
                birthdays = np.zeros(num, int)
                envgenomes = self._get_envgenomes(genomes)
                return Pop(genomes, ages, origins, uids, births, birthdays, envgenomes)
            else:
                return pop

        self.pop = _get_pop()

    ###############
    # CYCLE FUNCS #
    ###############

    def cycle(self):
        """Do one simulation cycle."""

        def _season_shift():
            # kill all living
            mask_kill = np.ones(len(self.pop), bool)
            self._kill(mask_kill, "season_shift")

            # hatch eggs
            if not self.eggs is None:
                self += self.eggs
                self.eggs = None

        if len(self.pop):
            # TODO rethink ordering
            self.gen_survival()
            self.eco_survival()
            self.reproduction()
            self.age()

        self.conf.season_countdown -= 1

        # if season over, kill pop and hatch eggs
        if self.conf.season_countdown == 0:
            _season_shift()
            self.conf.reset_season_countdown()

        if not self.conf.envmap is None and self.aux.stage % self.conf.ENVMAP_RATE == 0:
            self.conf.envmap.evolve()
            self.pop.envgenomes = self._get_envgenomes(self.pop.envgenomes)

    def age(self):
        """Increase age of all by one and kill those that surpass max lifespan."""
        self.pop.ages += 1
        mask_kill = self.pop.ages >= self.conf.MAX_LIFESPAN
        self._kill(mask_kill=mask_kill, causeofdeath="MAX_LIFESPAN")

    def eco_survival(self):
        """Impose ecological death, i.e. death that arises due to resource scarcity."""
        mask_kill = self.conf.overshoot(n=len(self.pop))
        self._kill(mask_kill=mask_kill, causeofdeath="overshoot")

    def gen_survival(self):
        """Impose genomic death, i.e. death that arises with probability encoded in the genome."""
        # mask_surv = self._get_mask(indices=self.pop.ages, attr="surv",)
        probs_surv = self._get_gen_probs("surv")
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
            muta_prob = self._get_gen_probs("muta", part=mask_repr)
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
        probs_repr = self._get_gen_probs("repr", part=mask_mature)
        mask_repr = np.random.random(len(probs_repr)) < probs_repr
        if not any(mask_repr):
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
        if self.conf.REPR_MODE == "asexual":
            origins = self.pop.uids[mask_repr]
        elif self.conf.REPR_MODE == "sexual":
            # TODO: think about the best way to encode double origination
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
            envgenomes=self._get_envgenomes(genomes),
        )

        # save as eggs if generations are nonoverlapping/discrete
        # otherwise add directly to population
        if self.conf.DISCRETE_GENERATIONS:
            if self.eggs is None:
                self.eggs = eggs
            else:
                self.eggs += eggs
        else:
            self.pop += eggs

    ################
    # HELPER FUNCS #
    ################

    def _get_envgenomes(self, genomes):
        """Get genomes that are contextualized in an env."""

        print(self._get_loc_probs(genomes))

        return genomes if self.conf.envmap is None else self.conf.envmap(genomes)
        # if self.conf.envmap is not None:
        #     genomes = self.conf.envmap(genomes)
        # if self.conf.phenomap is not None:
        #     genomes = self.conf.phenomap(genomes)
        # return genomes

    def _get_loc_probs(self, genomes):
        phenome = np.zeros(shape=genomes.shape[:2])

        for attr, pos in self.conf.loci_pos.items():
            loci = genomes[:, pos[0]:pos[1]]
            agespec, interpreter, lo, hi = self.conf.GENOME_STRUCT[attr][0:4]
            probs = self.conf.interpreter(loci, interpreter) * (hi - lo) + lo
            phenome[:, pos[0]:pos[1]] += probs

        t = self.conf.phenomap.transform(phenome)
        print(phenome[0])
        print(t[0])
        exit()

        return phenome

    def _get_gen_probs(self, attr, part=None):
        """
        Fetch probability values of certain genes at a certain loci.
        Three scenarios:
            1) value is not individual-specific nor age-specific (constant)
            2) value is individual-specific but not age-specific
            3) value is individual-specific and age-specific
        Use attr to specify the attribute to fetch.
        Use part to limit yourself to a part of the population.
        """

        # probabilities of which individuals do we calculate
        which_individuals = np.arange(len(self.pop))
        if not part is None:
            which_individuals = which_individuals[part]

        # first scenario
        if attr in self.conf.GENOME_CONST:
            probs = self.conf.GENOME_CONST[attr]

        # second and third scenario
        if attr in self.conf.GENOME_STRUCT:
            agespec, interpreter, lo, hi = self.conf.GENOME_STRUCT[attr][0:4]

            which_loci = (
                self.pop.ages[which_individuals]
                if agespec
                else self.conf.loci_pos[attr][0]
            )

            loci = self.pop.envgenomes[which_individuals, which_loci]
            probs = self.conf.interpreter(loci, interpreter) * (hi - lo) + lo

        # expand values back into an array with shape of whole pop
        final_probs = np.zeros(len(self.pop))
        final_probs[which_individuals] += probs

        return final_probs

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
        # TODO recorder
        self.aux.recorder.rec(self.pop[mask_record], causeofdeath, self.popid)

        # retain survivors
        self.pop *= ~mask_kill

