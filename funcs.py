import yaml
import json
import itertools
import collections
import functools

import pathlib

path = pathlib.Path(__file__).parent.absolute()

import logging

logging.basicConfig(
    format="%(asctime)s | %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S",
    level=logging.DEBUG,
)


def init_config(config={}):

    with open(path / "config.yml", "r") as f:
        d = yaml.load(f, Loader=yaml.FullLoader)

    for k, v in config.items():
        d[k] = v

    return collections.namedtuple("Config", d.keys())(*d.values())


def calc_survX(loci, x):
    p = 1

    iloc = 0
    for locus in loci:
        p *= locus / 10
        if p < x:
            return iloc
        iloc += 1

    return iloc


def calc_reprX(loci, x):
    if x < 0:
        x = 0
    return sum(loci[:x]) / 10


def survXfunc(x):
    return functools.partial(calc_survX, x=x)


def reprXfunc(x):
    return functools.partial(calc_reprX, x=x)
