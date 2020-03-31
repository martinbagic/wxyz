import collections
import yaml


def init(config={}):

    with open("config.yml", "r") as f:
        d = yaml.load(f, Loader=yaml.FullLoader)

    for k, v in config.items():
        d[k] = v

    return collections.namedtuple("Config", d.keys())(*d.values())
