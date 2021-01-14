import pandas as pd
import numpy as np
import pickle
import logging
import shutil
import json

# from functools import wraps
# from time import time
# def measure(func):
#     @wraps(func)
#     def _time_it(*args, **kwargs):
#         # start = int(round(time() * 1000))
#         start = time()
#         try:
#             return func(*args, **kwargs)
#         finally:
#             # end_ = int(round(time() * 1000)) - start
#             end_ = time() - start
#             print(f"Total execution time: {end_ if end_ > 0 else 0} ms")
#     return _time_it


class Recorder:
    """Temporarily saves and writes data from killed individuals to files."""

    attrs = {
        "sid",
        "pid",
        "bday",
        "age",
        "causeofdeath",
        "genomes",
        "popid",
        "phenotypes",
        "births",
    }

    def __init__(
        self, opath, params, LOCI_POS, BITS_PER_LOCUS, MATURATION_AGE,
    ):
        self.sid = []
        self.pid = []
        self.bday = []
        self.age = []
        self.causeofdeath = []
        self.genomes = []
        self.phenotypes = []
        self.births = []
        self.popid = []

        self.batch_number = 0
        self.FLUSH_RATE = params["FLUSH_RATE"]
        self.MAX_LIFESPAN = params["MAX_LIFESPAN"]
        self.JSON_RATE = params["JSON_RATE"]
        self.REC_RATE = params["REC_RATE"]

        self.LOCI_POS = LOCI_POS
        self.BITS_PER_LOCUS = BITS_PER_LOCUS
        self.opath = opath

        self.visor_paths = [
            self.opath / f"{self.opath.stem}.json",
            self.opath.parents[2] / "visor" / "csvs" / f"{self.opath.stem}.json",
        ]
        # self.visor_paths.append(self.visor_paths[-1].parent / "visor.json")

        self.rec_json_flag = True
        self.rec_flush_flag = True

        self.visor_data = {
            "bitsperlocus": self.BITS_PER_LOCUS,
            "survloc": self.LOCI_POS["surv"],
            "reprloc": self.LOCI_POS["repr"],
            "lifespan": self.MAX_LIFESPAN,
            "maturationage": MATURATION_AGE,
            "gensurv": [],
            "genrepr": [],
            "phesurv": [],
            "pherepr": [],
            "death_eco": [],
            "death_gen": [],
            "death_end": [],
            # "births": [],
            "params": params,
            "extinct": False,
        }

    # @measure
    def rec(self, pop, causeofdeath, popid):
        """Add population data to self. Flush if too many individuals recorded."""

        if self.rec_flush_flag or self.rec_json_flag:
            # add values to self
            self.sid.extend(pop.uids)
            self.pid.extend(pop.origins)
            self.bday.extend(pop.birthdays)
            self.age.extend(pop.ages)
            self.genomes.extend(
                pop.genomes.reshape(len(pop), -1)
            )  # flatten each genome
            self.causeofdeath.extend([causeofdeath] * len(pop))
            self.popid.extend([popid] * len(pop))
            self.phenotypes.extend(pop.phenotypes)
            self.births.extend(pop.births)

            if len(self) > self.FLUSH_RATE:
                self.save()

    def save(self, force=False):
        def dfize():
            """Rewrite data into three pandas dataframes."""
            df_gen = pd.DataFrame(np.array(self.genomes))
            df_gen.reset_index(drop=True, inplace=True)
            df_gen.columns = [str(c) for c in df_gen.columns]

            df_phe = pd.DataFrame(np.array(self.phenotypes))
            df_phe.reset_index(drop=True, inplace=True)
            df_phe.columns = [str(c) for c in df_phe.columns]

            dem_attrs = self.attrs - {"genomes", "phenotypes"}
            demo = {attr: getattr(self, attr) for attr in dem_attrs}
            df_dem = pd.DataFrame(demo, columns=dem_attrs)
            df_dem.reset_index(drop=True, inplace=True)
            df_dem["pid"] = df_dem.pid.astype(float)
            return df_gen, df_phe, df_dem

        df_gen, df_phe, df_dem = dfize()

        if self.rec_flush_flag or force:
            self.flush(df_gen, df_phe, df_dem)

        if self.rec_json_flag or force:
            self.record_for_visor(df_gen, df_phe, df_dem)

        self.rec_flush_flag = False
        self.rec_json_flag = False

        self._reinit()  # empty attrs
        self.batch_number += 1  # progress batch

    def flush(self, df_gen, df_phe, df_dem):
        """Write data to *.gen and *.dem files and erase all data from self."""
        path = self.opath / str(self.batch_number)
        df_gen.to_feather(path.with_suffix(".gen"))
        df_phe.to_feather(path.with_suffix(".phe"))
        df_dem.to_feather(path.with_suffix(".dem"))

    def record_for_visor(self, gen, phe, dem):
        def get_bits(loci_kind, array, bitsperlocus=self.BITS_PER_LOCUS):
            pos = self.LOCI_POS[loci_kind]
            return array.iloc[:, pos[0] * bitsperlocus : pos[1] * bitsperlocus]

        def get_deaths(death_kind):
            deaths = dem[dem.causeofdeath == death_kind].age.value_counts()
            return [int(deaths.get(age, 0)) for age in range(self.MAX_LIFESPAN + 1)]

        # def get_births():
        #     births = dem.births
        #     return [int(births.get(age, 0)) for age in range(self.MAX_LIFESPAN + 1)]

        data = {
            "gensurv": get_bits("surv", gen).mean(0).astype(float).tolist(),
            "genrepr": get_bits("repr", gen).mean(0).astype(float).tolist(),
            "phesurv": get_bits("surv", phe, 1).median(0).astype(float).tolist(),
            "pherepr": get_bits("repr", phe, 1).median(0).astype(float).tolist(),
            "death_eco": get_deaths("overshoot"),
            "death_gen": get_deaths("genetic"),
            "death_end": get_deaths("max_lifespan"),
            # "births": get_births(),
        }

        for k, v in data.items():
            self.visor_data[k].append(v)

        self.write_to_visor()

    def write_to_visor(self):
        for path in self.visor_paths:
            with open(path, "w") as ofile:
                json.dump(self.visor_data, ofile)

    def pickle_pop(self, obj, stage):
        """Pickle given population."""
        logging.info(f"Pickling the population at stage {stage}.")
        path = self.opath / f"{stage}.pickle"
        with open(path, "wb") as ofile:
            pickle.dump(obj, ofile)

    def __len__(self):
        """Return number of saved individuals."""
        num = len(self.genomes)
        assert all(len(getattr(self, attr)) == num for attr in self.attrs)
        return num

    def _reinit(self):
        """Erase saved data from self."""
        for attr in self.attrs:
            setattr(self, attr, [])
