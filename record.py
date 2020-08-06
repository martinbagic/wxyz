import pandas as pd
import numpy as np
import pathlib


class Record:
    def __init__(self, opath):

        self.init_containers()

        self.opath = opath
        self.opath_demo = opath.with_suffix(".demo.csv")  # demography
        self.opath_geno = opath.with_suffix(".geno.csv")  # genomes
        self.opath_demo_compressed = opath.with_suffix(".demo.feather")
        self.opath_geno_compressed = opath.with_suffix(".geno.feather")

        self.batch_number = 0

        self.entry_limit = 10 ** 6
        self.entry_count = 0

        self.demography_attrs = ("sid", "pid", "bday", "age", "causeofdeath")

        self.init_files()

    def init_files(self):
        with open(self.opath_demo, "w") as f:
            f.write(",".join(self.demography_attrs) + "\n")
        with open(self.opath_geno, "w") as f:
            f.write("genome\n")  # overwrite old

    def init_containers(self):
        self.sid, self.pid, self.bday, self.age, self.causeofdeath = [], [], [], [], []
        self.genomes = []

    def record(self):
        # demo
        data = {attr: getattr(self, attr) for attr in self.demography_attrs}
        df_demo = pd.DataFrame(data)
        df_demo.to_csv(self.opath_demo, mode="a", index=False, header=False)

        # geno
        df_geno = pd.DataFrame(self.genomes, columns=["genome"])
        df_geno.to_csv(self.opath_geno, mode="a", index=False, header=False)

        # check if batch needs to be compressed
        self.entry_count += len(df_geno)
        if self.entry_count > self.entry_limit:
            self.compress_batch()
            self.init_files()
            self.entry_count = 0

        # remove added data
        self.init_containers() 

    def compress_batch(self):
        # geno
        df = pd.read_csv(self.opath_geno)
        df = df.reset_index(drop=True)
        opath = self.opath.with_suffix(f".geno.feather.{self.batch_number}")
        df.to_feather(opath)
        self.opath_geno.unlink()

        # demo
        df = pd.read_csv(self.opath_demo)
        df = df.reset_index(drop=True)
        opath = self.opath.with_suffix(f".demo.feather.{self.batch_number}")
        df.to_feather(opath)
        self.opath_demo.unlink()

        self.batch_number += 1

    def compress_output(self):
        # get a list of unique genomes and a list of genome positions

        genomes = list()
        for i in range(self.batch_number):
            path = self.opath.with_suffix(f".geno.feather.{i}")
            df = pd.read_feather(path)
            genomes.extend(df.genome.values)
        gnmsset = list(set(genomes))
        gnmsdict = {g: gid for gid, g in enumerate(gnmsset)}
        gid = [gnmsdict[g] for g in genomes]  # genome => genome id

        # compress genome data
        try:
            pd.DataFrame(gnmsset, columns=["genome"]).to_feather(
                self.opath_geno_compressed
            )
            for i in range(self.batch_number):
                path = self.opath.with_suffix(f".geno.feather.{i}")
                path.unlink()

        except:
            print(f"Writing to '{self.opath_geno_compressed}' failed.")

        # compress demography data
        df = pd.read_feather(self.opath.with_suffix(".demo.feather.0"))

        for i in range(1, self.batch_number):
            path = self.opath.with_suffix(f".demo.feather.{i}")
            dfi = pd.read_feather(path)
            df = pd.concat([df, dfi])

        df["gid"] = gid  # add genome id's to demography data
        df = df.reset_index(drop=True)

        try:
            df.to_feather(self.opath_demo_compressed)
            for i in range(self.batch_number):
                path = self.opath.with_suffix(f".demo.feather.{i}")
                path.unlink()
        except:
            print(f"Writing to '{self.opath_demo_compressed}' failed.")

