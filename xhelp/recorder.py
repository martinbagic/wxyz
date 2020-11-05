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
        "phenomes",
        "births",
    }

    def __init__(
        self,
        opath,
        params,
        FLUSH_RATE,
        MAX_LIFESPAN,
        LOCI_POS,
        BITS_PER_LOCUS,
        MATURATION_AGE,
    ):
        self.sid = []
        self.pid = []
        self.bday = []
        self.age = []
        self.causeofdeath = []
        self.genomes = []
        self.phenomes = []
        self.births = []
        self.popid = []

        self.batch_number = 0
        self.FLUSH_RATE = FLUSH_RATE
        self.MAX_LIFESPAN = MAX_LIFESPAN
        self.LOCI_POS = LOCI_POS
        self.BITS_PER_LOCUS = BITS_PER_LOCUS
        self.opath = opath

        self.vizport_paths = [
            self.opath / f"{self.opath.stem}.json",
            self.opath.parents[2] / "vizport" / "csvs" / f"{self.opath.stem}.json",
        ]
        # self.vizport_paths.append(self.vizport_paths[-1].parent / "vizport.json")

        self.vizport_data = {
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
        }

    # @measure
    def rec(self, pop, causeofdeath, popid):
        """Add population data to self. Flush if too many individuals recorded."""
        # add values to self
        self.sid.extend(pop.uids)
        self.pid.extend(pop.origins)
        self.bday.extend(pop.birthdays)
        self.age.extend(pop.ages)
        self.genomes.extend(pop.genomes.reshape(len(pop), -1))  # flatten each genome
        self.causeofdeath.extend([causeofdeath] * len(pop))
        self.popid.extend([popid] * len(pop))
        self.phenomes.extend(pop.phenomes)
        self.births.extend(pop.births)

        # if record limit reached, flush
        if len(self) > self.FLUSH_RATE:
            self.flush()

    # @measure
    def pickle_pop(self, obj, stage):
        """Pickle given population."""
        logging.info(f"Pickling the population at stage {stage}.")
        path = self.opath / f"{stage}.pickle"
        with open(path, "wb") as ofile:
            pickle.dump(obj, ofile)

    def record_for_vizport(self, gen, phe, dem):
        def get_deaths(death_kind):
            deaths = dem[dem.causeofdeath == death_kind].age.value_counts()
            return [int(deaths.get(age, 0)) for age in range(self.MAX_LIFESPAN + 1)]

        # def get_births():
        #     births = dem.births
        #     return [int(births.get(age, 0)) for age in range(self.MAX_LIFESPAN + 1)]

        def get_bits(loci_kind, array, bitsperlocus=self.BITS_PER_LOCUS):
            pos = self.LOCI_POS[loci_kind]
            return array.iloc[:, pos[0] * bitsperlocus : pos[1] * bitsperlocus]

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
            self.vizport_data[k].append(v)

        for path in self.vizport_paths:
            with open(path, "w") as ofile:
                json.dump(self.vizport_data, ofile)

    # @measure
    def flush(self):
        """Write data to *.gen and *.dem files and erase all data from self."""
        # logging.info(f"Flushing {len(self.genomes)} records")

        # write .gen and .dem
        path = self.opath / str(self.batch_number)

        df_gen = pd.DataFrame(np.array(self.genomes))
        df_gen.reset_index(drop=True, inplace=True)
        df_gen.columns = [str(c) for c in df_gen.columns]
        df_gen.to_feather(path.with_suffix(".gen"))

        df_phe = pd.DataFrame(np.array(self.phenomes))
        df_phe.reset_index(drop=True, inplace=True)
        df_phe.columns = [str(c) for c in df_phe.columns]
        df_phe.to_feather(path.with_suffix(".phe"))

        dem_attrs = self.attrs - {"genomes", "phenomes"}
        demo = {attr: getattr(self, attr) for attr in dem_attrs}
        df_dem = pd.DataFrame(demo, columns=dem_attrs)
        df_dem.reset_index(drop=True, inplace=True)
        df_dem["pid"] = df_dem.pid.astype(float)
        df_dem.to_feather(path.with_suffix(".dem"))

        self.record_for_vizport(df_gen, df_phe, df_dem)

        # empty attrs
        self._reinit()

        # progress batch
        self.batch_number += 1

    def __len__(self):
        """Return number of saved individuals."""
        num = len(self.genomes)
        assert all(len(getattr(self, attr)) == num for attr in self.attrs)
        return num

    def _reinit(self):
        """Erase saved data from self."""
        for attr in self.attrs:
            setattr(self, attr, [])
