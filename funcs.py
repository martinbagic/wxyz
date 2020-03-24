import collections
import yaml

def load_config(path):
    with open(path, "r") as f:
        d = yaml.load(f)
        config = collections.namedtuple("config", d.keys())(*d.values())
    return config
