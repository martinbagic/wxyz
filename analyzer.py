import pathlib
import funcs
import pandas
import pickle
import yaml
import numpy

import functools
import collections

import tqdm


def measure(df, attr, q):
    return df[df.birthdays > 2000][attr].quantile(q)


class Measurer:
    def extinction_time(df):
        if any(df.ages == -1):
            return -1
        else:
            return max(df.ages + df.birthdays)

    def quantile(df, attr, q):
        return df[df.birthdays > 2000][attr].quantile(q)

    def surv05_05(df):
        return df[df.birthdays > 2000]["surv0.5"].quantile(0.5)

    def repr05_05(df):
        return df[df.birthdays > 2000]["repr0.5"].quantile(0.5)

    def mutrates_05(df):
        return df[df.birthdays > 2000]["mutrates"].quantile(0.5)


class Analyzer:
    def __init__(self, exp_name):
        self.dir = funcs.path.parents[0] / exp_name
        # self.mnames = ("extinction_time", "surv05_05", "repr05_05", "mutrates_05")
        self.config = self.load_config()

        self.msr_attrs = [
            "ages",
            "births",
            "surv0.5",
            "surv0.9",
            "surv0.1",
            "repr0.5",
            "repr0.9",
            "repr0.1",
            "mutrates",
        ]

        self.msr_ncycles = [
            100,
            1000,
            10000,
        ]

        self.msr_quantile = [0.1, 0.25, 0.5, 0.75, 0.9]

    def run(self):
        self.analyze()
        self.order()
        self.save_pickle()

    def load_config(self):
        path = self.dir / "cartesian.yml"
        with open(path, "r") as f:
            return yaml.load(f)

    def analyze(self):
        def measure(df, mname):
            mfunc = getattr(Measurer, mname)
            mval = mfunc(df)
            return mval

        def msr(df, attr, quantile, ncycles):
            return df[df.birthdays > ncycles][attr].quantile(quantile)

        # get all paths
        paths = [x for x in self.dir.iterdir() if x.suffix == ".csv"]
        paths.sort(key=lambda x: int(x.name.split(".")[0]))

        # initialize measure dict
        mdict = collections.defaultdict(list)
        # mdict = {mname: [] for mname in self.mnames}

        for path in tqdm.tqdm(paths):

            # read file
            df = pandas.read_csv(path)

            # measure
            for attr in self.msr_attrs:
                for ncycles in self.msr_ncycles:
                    for quantile in self.msr_quantile:
                        mval = msr(df,attr,quantile,ncycles)
                        mdict[f'{attr}|{ncycles}|{quantile}'].append(mval)
            # for mname in mdict:
            #     mval = measure(df, mname)
            #     mdict[mname].append(mval)

        self.mdict = mdict

    def order(self):
        paramlens = [len(l) for l in self.config.values() if len(l) != 1]

        if paramlens:
            marrays = {
                mname: numpy.array(mvals).reshape(paramlens)
                for mname, mvals in self.mdict.items()
            }
        else:
            marrays = {mname: numpy.array(mvals) for mname, mvals in self.mdict.items()}

        self.marrays = marrays

    def save_pickle(self):
        path = self.dir / ("cartesian.analyzerII")
        with open(path, "wb") as f:
            pickle.dump(self, f)

    def load_pickle(self):
        path = self.dir / ("cartesian.analyzerII")
        with open(path, "rb") as f:
            return pickle.load(f)
