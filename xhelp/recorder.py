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
        "evaluomes",
    }

    def __init__(self, opath, FLUSH_RATE, MAX_LIFESPAN):
        self.sid = []
        self.pid = []
        self.bday = []
        self.age = []
        self.causeofdeath = []
        self.genomes = []
        self.evaluomes = []
        self.popid = []

        self.batch_number = 0
        self.FLUSH_RATE = FLUSH_RATE
        self.MAX_LIFESPAN = MAX_LIFESPAN

        self.opath = opath

        self.opath.mkdir(exist_ok=True)  # make folder in which data will be saved

        self.path_vizport = self.opath / "vizport.csv"
        self.vizport_header_recorded = False
        # self.vizport_data = {"data": []}

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
        self.evaluomes.extend(pop.evaluomes)

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

    def record_for_vizport(self, gen, eva, dem):

        if not self.vizport_header_recorded:
            header = ",".join(
                [f"gen_{i}" for i in range(gen.shape[1])]
                + [f"eva_{i}" for i in range(eva.shape[1])]
                + [
                    f"age_{i}" for i in range(self.MAX_LIFESPAN + 1)
                ]  # ignoring deaths "at age -1", when simulation end comes
                + ["cod_max_lifespan", "cod_overshoot", "cod_genetic"]
                + ["max_sid"]
            )

            with open(self.path_vizport, "w") as ofile:
                ofile.write(header + "\n")
            self.vizport_header_recorded = True

        agecounts = dem.age.value_counts()
        codcounts = dem.causeofdeath.value_counts()
        lis = (
            gen.mean(0).tolist()
            + eva.median(0).tolist()
            + [agecounts.get(i, 0) for i in range(self.MAX_LIFESPAN + 1)]
            + [codcounts.get(i, 0) for i in ["max_lifespan", "overshoot", "genetic"]]
            + [dem.sid.max()]
        )
        s = ",".join(str(x) for x in lis)

        with open(self.path_vizport, "a") as ofile:
            ofile.write(s + "\n")

        # self.vizport_data["data"].append([float(x) for x in lis])

        # with open(self.path_vizport.with_suffix(".json"), "w") as ofile:
        #     json.dump(self.vizport_data, ofile)

        shutil.copy(
            self.path_vizport,
            str(
                self.opath.parent.parent.parent
                / "wxyz"
                / "vizport"
                / "csvs"
                / self.opath.stem
            )
            + ".csv",
        )

        shutil.copy(
            self.path_vizport,
            self.opath.parent.parent.parent
            / "wxyz"
            / "vizport"
            / "csvs"
            / "vizport.csv",
        )

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

        df_eva = pd.DataFrame(np.array(self.evaluomes))
        df_eva.reset_index(drop=True, inplace=True)
        df_eva.columns = [str(c) for c in df_eva.columns]
        df_eva.to_feather(path.with_suffix(".eva"))

        dem_attrs = self.attrs - {"genomes"}
        demo = {attr: getattr(self, attr) for attr in dem_attrs}
        df_dem = pd.DataFrame(demo, columns=dem_attrs)
        df_dem.reset_index(drop=True, inplace=True)
        df_dem["pid"] = df_dem.pid.astype(
            float
        )  # fix double origination, TODO: can i simplify
        df_dem.to_feather(path.with_suffix(".dem"))

        self.record_for_vizport(df_gen, df_eva, df_dem)

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
