import yaml
import json
import itertools
import collections

import pathlib

path = pathlib.Path(__file__).parent.absolute()

import logging
logging.basicConfig(
    format="%(asctime)s | %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S",
    level=logging.DEBUG,
)


def init_config(config={}):

    with open("config.yml", "r") as f:
        d = yaml.load(f, Loader=yaml.FullLoader)

    for k, v in config.items():
        d[k] = v

    return collections.namedtuple("Config", d.keys())(*d.values())
