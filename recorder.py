import pandas as pd
import numpy as np
import pickle


class Recorder:
    """Temporarily saves and writes data from killed individuals to files."""

    attrs = {"sid", "pid", "bday", "age", "causeofdeath", "genomes"}

    def __init__(self, opath, entry_limit):
        self._reinit()
        self.batch_number = 0
        self.entry_limit = entry_limit

        self.opath = opath
        self.opath.mkdir(exist_ok=True)  # make folder in which data will be saved

    def rec(self, pop, causeofdeath):
        """Add population data to self. Flush if too many individuals recorded."""
        # add values to self
        self.sid.extend(pop.uids)
        self.pid.extend(pop.origins)
        self.bday.extend(pop.birthdays)
        self.age.extend(pop.ages)
        self.genomes.extend(pop.genomes.reshape(len(pop), -1))  # flatten each genome
        self.causeofdeath.extend([causeofdeath] * len(pop))

        # if overshooting, flush
        if len(self) > self.entry_limit:
            self.flush()

    def pickle_pop(self, pop):
        """Pickle population."""
        path = self.opath / f"{self.stage}.pop"
        with open(path, "wb") as ofile:
            pickle.dump(pop, ofile)

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
        self.sid = []
        self.pid = []
        self.bday = []
        self.age = []
        self.causeofdeath = []
        self.genomes = []
