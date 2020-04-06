import pathlib
import funcs
import pandas
import pickle
import yaml
import numpy


class Measurer:
    def extinction_time(df):
        if any(df.ages == -1):
            return -1
        else:
            return max(df.ages + df.birthdays)

    


class Analyzer:
    def __init__(self, exp_name):
        self.dir = funcs.path.parents[0] / exp_name
        self.mnames = ("extinction_time",)
        self.config = self.load_config()

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

        # get all paths
        paths = [x for x in self.dir.iterdir() if x.suffix == ".csv"]
        paths.sort(key=lambda x: int(x.name.split(".")[0]))

        print(paths)

        # initialize measure dict
        mdict = {mname: [] for mname in self.mnames}

        for path in paths:

            # read file
            df = pandas.read_csv(path)

            # measure
            for mname in mdict:
                mval = measure(df, mname)
                mdict[mname].append(mval)

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
        path = self.dir / ("cartesian.analyzer")
        with open(path, "wb") as f:
            pickle.dump(self, f)

    def load_pickle(self):
        path = self.dir / ("cartesian.analyzer")
        with open(path, "rb") as f:
            return pickle.load(f)
