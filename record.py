import json
import numpy as np
import time
import logging
import sys
import time

# add position of parameter set (index in itertools.product) so you can trace it back later more easily

class Record:

    time_format = "%H:%M:%S"

    def __init__(self, identifier, config):
        self.d = {
            "identifier": identifier,
            "config": config,
            "time_start": time.strftime(Record.time_format),
            "data": {"genomes": [], "mutrates": [], "popsize": [], "ages": []},
        }

        self.time0 = time.time()

    def __call__(self, pop):
        quantiles = [0.25, 0.5, 0.75]

        def get_qs(values):
            return [np.quantile(values / 10, q) for q in quantiles]

        mutrates = pop._get_mutation_probs(pop.genomes)
        ages = pop.ages

        tail = {
            "genomes": [
                list(np.quantile(pop.genomes.sum(2), q, axis=1)) for q in quantiles
            ],
            "ages": get_qs(ages),
            "mutrates": get_qs(mutrates),
            "popsize": len(pop.genomes),
        }

        for k, v in tail.items():
            self.d["data"][k].append(v)

    def save(self, jobid, opath):
        # print(vars(self))
        time_difference = time.time() - self.time0
        self.d["time_run"]: time.strftime("%H:%M:%S", time.gmtime(time_difference))

        print(sys.getsizeof(self.d["data"]))

        # with open(f"records-cluster/{self.d['identifier']}.json", "w") as f:
        with open(f"records-cluster/{str(time.time())}.json", "w") as f:
            json.dump(self.d, f)
