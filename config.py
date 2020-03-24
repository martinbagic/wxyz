import collections
import yaml


def load_config(path):
    global CONFIG
    with open(path, "r") as f:
        d = yaml.load(f)
        CONFIG = collections.namedtuple("config", d.keys())(*d.values())
