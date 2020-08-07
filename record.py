import pandas as pd
import numpy as np
import pathlib


class Record:
    def __init__(self, dpath):

        self.init_containers()

        dpath.mkdir(exist_ok=True)
        self.opath = dpath

        # self.opath_demo = opath.with_suffix(".demo.csv")  # demography
        # self.opath_geno = opath.with_suffix(".geno.csv")  # genomes
        # self.opath_demo_compressed = opath.with_suffix(".demo.feather")
        # self.opath_geno_compressed = opath.with_suffix(".geno.feather")

        self.batch_number = 0

        self.entry_limit = 10 ** 4

        self.demography_attrs = ("sid", "pid", "bday", "age", "causeofdeath")

        # self.init_files()

    # def init_files(self):
    #     with open(self.opath_demo, "w") as f:
    #         f.write(",".join(self.demography_attrs) + "\n")
    #     with open(self.opath_geno, "w") as f:
    #         f.write("genome\n")  # overwrite old
    #         # f.write("")

    def init_containers(self):
        self.sid, self.pid, self.bday, self.age, self.causeofdeath = [], [], [], [], []
        self.genomes = pd.DataFrame()

    def add_genomes(self, genomes):
        df = pd.DataFrame(genomes.sum(2))
        self.genomes = self.genomes.append(df)

    def record(self):

        df_geno = self.genomes

        if len(df_geno) > self.entry_limit:

            path = self.opath / str(self.batch_number)

            df_geno.reset_index(drop=True, inplace=True)
            df_geno.columns = [str(c) for c in df_geno.columns]
            df_geno.to_feather(path.with_suffix(".gen"))

            data = {attr: getattr(self, attr) for attr in self.demography_attrs}
            df_demo = pd.DataFrame(data, columns=self.demography_attrs)
            df_demo.reset_index(drop=True, inplace=True)
            df_demo.to_feather(path.with_suffix(".dem"))

            self.init_containers()

            self.batch_number += 1

    # def compress_output(self):
    #     # get a list of unique genomes and a list of genome positions

    #     genomes = list()
    #     for i in range(self.batch_number):
    #         path = self.opath.with_suffix(f".geno.feather.{i}")
    #         df = pd.read_feather(path)
    #         genomes.extend(df.genome.values)
    #     gnmsset = list(set(genomes))
    #     gnmsdict = {g: gid for gid, g in enumerate(gnmsset)}
    #     gid = [gnmsdict[g] for g in genomes]  # genome => genome id

    #     # compress genome data
    #     try:
    #         pd.DataFrame(gnmsset, columns=["genome"]).to_feather(
    #             self.opath_geno_compressed
    #         )
    #         for i in range(self.batch_number):
    #             path = self.opath.with_suffix(f".geno.feather.{i}")
    #             path.unlink()

    #     except:
    #         print(f"Writing to '{self.opath_geno_compressed}' failed.")

    #     # compress demography data
    #     df = pd.read_feather(self.opath.with_suffix(".demo.feather.0"))

    #     for i in range(1, self.batch_number):
    #         path = self.opath.with_suffix(f".demo.feather.{i}")
    #         dfi = pd.read_feather(path)
    #         df = pd.concat([df, dfi])

    #     df["gid"] = gid  # add genome id's to demography data
    #     df = df.reset_index(drop=True)

    #     try:
    #         df.to_feather(self.opath_demo_compressed)
    #         for i in range(self.batch_number):
    #             path = self.opath.with_suffix(f".demo.feather.{i}")
    #             path.unlink()
    #     except:
    #         print(f"Writing to '{self.opath_demo_compressed}' failed.")


class Contracter:
    def __init__(self, dpath):

        file_i = lambda s: int(s.stem)
        dems = sorted(dpath.glob("*.dem"), key=file_i)
        gens = sorted(dpath.glob("*.gen"), key=file_i)

        def agg_feathers(paths):
            main_df = pd.read_feather(paths[0])
            for path in paths[1:]:
                df = pd.read_feather(path)
                main_df = main_df.append(df)
            main_df.reset_index(drop=True, inplace=True)
            return main_df

        dem_df = agg_feathers(dems)
        gen_df = agg_feathers(gens)

        dem_df.to_feather(dpath / "dem")
        gen_df.to_feather(dpath / "gen")        