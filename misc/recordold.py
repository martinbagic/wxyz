import pandas as pd
import numpy as np
import pathlib


class Record:
    def __init__(self, dpath):

        self._reinit_containers()

        dpath.mkdir(exist_ok=True)
        self.opath = dpath

        self.batch_number = 0

        self.entry_limit = 10 ** 4

        self.demography_attrs = ("sid", "pid", "bday", "age", "causeofdeath")

    def _reinit_containers(self):
        self.sid, self.pid, self.bday, self.age, self.causeofdeath = [], [], [], [], []
        self.genomes = []

    def rec(self, pop, causeofdeath):
        self.sid.extend(pop.uids)
        self.pid.extend(pop.origins)
        self.bday.extend(pop.birthdays)
        self.age.extend(pop.ages)
        self.genomes.extend(pop.genomes.reshape(len(pop), -1)) # flatten each genome
        self.causeofdeath.extend([causeofdeath] * len(pop))

        self.flush()


    def record(self, sid, pid, bday, age, causeofdeath, genomes):

        self.sid.extend(sid)
        self.pid.extend(pid)
        self.bday.extend(bday)
        self.age.extend(age)
        self.causeofdeath.extend(causeofdeath)
        self.genomes.extend(genomes)

        self.flush()

    def flush(self, force=False):

        if (len(self.sid) > self.entry_limit) or force:

            path = self.opath / str(self.batch_number)

            print(np.concatenate(self.genomes).shape)

            df_geno = pd.DataFrame(np.array(self.genomes))
            df_geno.reset_index(drop=True, inplace=True)
            df_geno.columns = [str(c) for c in df_geno.columns]
            df_geno.to_feather(path.with_suffix(".gen"))

            data = {attr: getattr(self, attr) for attr in self.demography_attrs}
            df_demo = pd.DataFrame(data, columns=self.demography_attrs)
            df_demo.reset_index(drop=True, inplace=True)
            df_demo.to_feather(path.with_suffix(".dem"))

            self._reinit_containers()

            self.batch_number += 1


# class Contracter:
#     def __init__(self, dpath):

#         file_i = lambda s: int(s.stem)
#         dems = sorted(dpath.glob("*.dem"), key=file_i)
#         gens = sorted(dpath.glob("*.gen"), key=file_i)

#         def agg_feathers(paths):
#             main_df = pd.read_feather(paths[0])
#             for path in paths[1:]:
#                 df = pd.read_feather(path)
#                 main_df = main_df.append(df)
#             main_df.reset_index(drop=True, inplace=True)
#             return main_df

#         dem_df = agg_feathers(dems)
#         gen_df = agg_feathers(gens)

#         dem_df.to_feather(dpath / "dem")
#         gen_df.to_feather(dpath / "gen")

