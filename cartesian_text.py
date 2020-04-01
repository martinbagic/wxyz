import argparse
import yaml
import json
import itertools

import funcs


class CartesianText:
    def __init__(self):
        self.set_args()
        self.set_paths()
        self.set_yml()

    def set_args(self):
        parser = argparse.ArgumentParser("CartesianText")
        parser.add_argument("dirname", type=str)
        self.args = parser.parse_args()

    def set_paths(self):
        dirpath = funcs.path.parents[0] / self.args.dirname
        self.ymlpath = dirpath / "cartesian.yml"
        self.txtpath = dirpath / "cartesian.txt"

    def set_yml(self):
        with open(self.ymlpath, "r") as f:
            self.yml = yaml.load(f)

    def write_cartesian_txt(self):
        keys = tuple(self.yml.keys())
        valuess = tuple(itertools.product(*self.yml.values()))

        s = ""

        for values in valuess:
            d = dict(zip(keys, values))
            s += json.dumps(d) + "\n"

        with open(self.txtpath, "w") as f:
            f.write(s)


if __name__ == "__main__":
    CartesianText().write_cartesian_txt()
