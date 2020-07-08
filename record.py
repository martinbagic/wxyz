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

    def __init__(self, identifier, config, opath):
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
        self.causeofdeath = []
        self.fullgenomes = []

        self.time0 = time.time()
        self.opath = opath

        with open(self.opath, "w") as f:
            self.quantiles = (0.1, 0.25, 0.5, 0.75, 0.9)
            headers = (
                [
                    "ages",
                    "births",
                    "birthdays",
                    "origins",
                    "uids",
                    "causeofdeath",
                    "fullgenomes",
                ]
                + [f"repr{q}" for q in self.quantiles]
                + [f"surv{q}" for q in self.quantiles]
                + ["mutrates", "neutloci"]
            )
            f.write(",".join(headers) + "\n")

    def write_genomes(self, genomes, uids):
        d = {"genomes": genomes.tolist(), "uids": uids.tolist()}
        # print(d)
        path = self.opath.with_suffix(".genomes")
        with open(path, "w") as f:
            json.dump(d, f)

    def flush(self):
        def record():

            self.d["data"] = {
                attr: [x for x in getattr(self, attr)]
                for attr in (
                    "ages",
                    "births",
                    "birthdays",
                    "origins",
                    "uids",
                    "causeofdeath",
                    "fullgenomes",
                )
            }

            if self.d["data"]["fullgenomes"] == []:
                del self.d["data"]["fullgenomes"]

            self.d["data"]["ages"] = [int(x) for x in self.d["data"]["ages"]]

            # flatten
            for attr in ("neutloci", "survloci", "reprloci"):
                self.genomes[attr] = [x for l in self.genomes[attr] for x in l]
            for q in self.quantiles:
                self.d["data"][f"surv{q}"] = [
                    funcs.calc_survX(loci, q) for loci in self.genomes["survloci"]
                ]

            for q in self.quantiles:
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

#             print(self.d["data"])

            df = pandas.DataFrame(self.d["data"])

            df.to_csv(self.opath, mode="a", index=False, header=False)  # different

        def clean():
            self.genomes = {
                "mutrates": [],
                "neutloci": [],
                "survloci": [],
                "reprloci": [],
            }
            self.ages = []
            self.births = []
            self.origins = []
            self.birthdays = []
            self.uids = []
            self.causeofdeath = []
            self.fullgenomes = []

        record()
        clean()
