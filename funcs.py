import yaml
import json
import itertools
import collections

import pathlib
path = pathlib.Path(__file__).parent.absolute()
       
        
def init_config(config={}):

    with open("config.yml", "r") as f:
        d = yaml.load(f, Loader=yaml.FullLoader)

    for k, v in config.items():
        d[k] = v

    return collections.namedtuple("Config", d.keys())(*d.values())

