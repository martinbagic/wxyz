import collections
import yaml


def init(config={}):

    with open("config.yml", "r") as f:
        d = yaml.load(f, Loader=yaml.FullLoader)

    for k, v in config.items():
        d[k] = v

    assert d["bits_per_locus"] % 2 == 0 # check that you have 2 conceptual chromosomes

    return collections.namedtuple("Config", d.keys())(*d.values())
