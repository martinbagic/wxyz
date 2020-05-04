import json
import numpy as np
import time
import logging
import sys
import time
import pandas

import funcs

# add position of parameter set (index in itertools.product) so you can trace it back later more easily


class Record:

    time_format = "%H:%M:%S"

    def __init__(self, identifier, config):
        self.d = {
            "identifier": identifier,
            "config": config,
            "time_start": time.strftime(Record.time_format),
        }

        self.genomes = {"mutrates": [], "neutloci": [], "survloci": [], "reprloci": []}
        self.ages = []
        self.births = []
        self.origins = []
        self.birthdays = []
        self.uids = []

        self.time0 = time.time()

    # def __call__(self, pop):
    #     quantiles = [0.25, 0.5, 0.75]

    #     def get_qs(values):
    #         return [np.quantile(values / 10, q) for q in quantiles]

    #     mutrates = pop._get_mutation_probs(pop.genomes)
    #     ages = pop.ages

    #     tail = {
    #         "genomes": [
    #             list(np.quantile(pop.genomes.sum(2), q, axis=1)) for q in quantiles
    #         ],
    #         "ages": get_qs(ages),
    #         "mutrates": get_qs(mutrates),
    #         "popsize": len(pop.genomes),
    #     }

    #     for k, v in tail.items():
    #         self.d["data"][k].append(v)

    def save(self, opath):
        time_difference = time.time() - self.time0
        self.d["time_run"]: time.strftime("%H:%M:%S", time.gmtime(time_difference))

        self.d["data"] = {
            attr: [x for x in getattr(self, attr)]
            for attr in ("ages", "births", "birthdays", "origins", "uids",)
        }

        self.d["data"]["ages"] = [int(x) for x in self.d["data"]["ages"]]

        # flatten
        for attr in ("neutloci", "survloci", "reprloci"):
            self.genomes[attr] = [x for l in self.genomes[attr] for x in l]
        for q in (0.5, 0.9, 0.1):
            self.d["data"][f"surv{q}"] = [
                funcs.calc_survX(loci, q) for loci in self.genomes["survloci"]
            ]

        for q in (0.5, 0.9, 0.1):
            self.d["data"][f"repr{q}"] = [
                funcs.calc_reprX(loci, age - self.d["config"]["maturation_age"])
                for loci, age in zip(
                    self.genomes["reprloci"], self.d["data"][f"surv{q}"]
                )
            ]

        self.d["data"]["mutrates"] = [
            mutrate for mutrates in self.genomes["mutrates"] for mutrate in mutrates
        ]

        self.d["data"]["neutloci"] = [
            np.mean(neutloci) for neutloci in self.genomes["neutloci"]
        ]

        # with open(opath, "w") as f:
        #     json.dump(self.d, f)

        df = pandas.DataFrame(self.d["data"])
        df.to_csv(opath, index=False)
