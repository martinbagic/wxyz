import pandas as pd


class Record:
    def __init__(self, opath):

        self.init_containers()

        self.opath = opath
        self.opath_demo = opath.with_suffix(".demo.csv")  # demography
        self.opath_geno = opath.with_suffix(".geno.csv")  # genomes
        self.opath_demo_compressed = opath.with_suffix(".demo.feather")
        self.opath_geno_compressed = opath.with_suffix(".geno.feather")

        self.demography_attrs = ("sid", "pid", "bday", "age", "causeofdeath")

        with open(self.opath_demo, "w") as f:
            f.write(",".join(self.demography_attrs) + "\n")
        with open(self.opath_geno, "w") as f:
            f.write("genome\n")  # overwrite old

    def init_containers(self):
        self.sid, self.pid, self.bday, self.age, self.causeofdeath = [], [], [], [], []
        self.genomes = []

    def record_demography(self):
        data = {attr: getattr(self, attr) for attr in self.demography_attrs}
        df_demo = pd.DataFrame(data)
        df_demo.to_csv(self.opath_demo, mode="a", index=False, header=False)

    def record_genomes(self):
        # data = [
        #     ["".join(locus) for locus in genome.astype(str)]
        #     for genome in self.genomes
        # ]
        data = [
            "".join("".join(locus) for locus in genome.astype(str))
            for genome in self.genomes
        ]
        df_geno = pd.DataFrame(data, dtype=str)
        df_geno.to_csv(self.opath_geno, mode="a", index=False, header=False)

    def compress_output(self):

        # get a list of unique genomes and a list of genome positions
        df = pd.read_csv(self.opath_geno)
        genomes = list(set(df.genome)) # get unique genomes
        gid = [genomes.index(genome) for genome in df.genome]  # genome => genome id

        # compress genome data
        try:
            pd.DataFrame(genomes, columns=["genome"]).to_feather(
                self.opath_geno_compressed
            )
            self.opath_geno.unlink()
        except:
            print(f"Writing to '{self.opath_geno_compressed}' failed.")

        # compress demography data
        df = pd.read_csv(self.opath_demo)
        df["gid"] = gid # add genome id's to demography data
        try:
            df.to_feather(self.opath_demo_compressed)
            self.opath_demo.unlink()
        except:
            print(f"Writing to '{self.opath_demo_compressed}' failed.")

