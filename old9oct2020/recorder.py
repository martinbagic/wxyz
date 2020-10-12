import pandas as pd
import numpy as np
import pickle

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

    attrs = {"sid", "pid", "bday", "age", "causeofdeath", "genomes", "popid"}

    def __init__(self, opath, entry_limit):
        self.sid = []
        self.pid = []
        self.bday = []
        self.age = []
        self.causeofdeath = []
        self.genomes = []
        self.popid = []

        self.batch_number = 0
        self.entry_limit = entry_limit

        self.opath = opath

        self.opath.mkdir(exist_ok=True)  # make folder in which data will be saved

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
        self.popid.extend([popid]*len(pop))

        # if record limit reached, flush
        if len(self) > self.entry_limit:
            self.flush()

    # @measure
    def pickle_pop(self, pop, stage):
        """Pickle given population."""
        path = self.opath / f"{stage}.pop"
        with open(path, "wb") as ofile:
            pickle.dump(pop, ofile)

    # @measure
    def flush(self):
        """Write data to *.gen and *.dem files and erase all data from self."""

        # write .gen and .dem
        path = self.opath / str(self.batch_number)

        df_gen = pd.DataFrame(np.array(self.genomes))
        df_gen.reset_index(drop=True, inplace=True)
        df_gen.columns = [str(c) for c in df_gen.columns]
        df_gen.to_feather(path.with_suffix(".gen"))

        dem_attrs = self.attrs - {"genomes"}
        demo = {attr: getattr(self, attr) for attr in dem_attrs}
        df_dem = pd.DataFrame(demo, columns=dem_attrs)
        df_dem.reset_index(drop=True, inplace=True)
        df_dem.to_feather(path.with_suffix(".dem"))

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
