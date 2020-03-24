import time
import pickle
import os
import argparse

from population import Population
from plot import plot


def save_pop(pop):
    filenames = [int(x) for x in os.listdir("pops")] + [0]
    filename = max(filenames) + 1
    print(f"Saving to filename '{filename}'.")

    with open(f"pops/{filename}", "wb") as f:
        pickle.dump(pop, f)
    return filename


def load_pop(filename):
    with open(f"pops/{filename}", "rb") as f:
        pop = pickle.load(f)
    return pop


def get_stable_pop(n_stages=20):
    stable = False

    while not stable:
        pop = Population()
        for i in range(n_stages):
            pop.cycle()
            if pop.is_extinct():
                break
        stable = not pop.is_extinct()

    return pop


def run(cycles=1000, filename=None):

    if filename:
        pop = load_pop(filename)
    else:
        # pop = get_stable_pop(100)
        pop = Population()

    print('population initialized.')
    print('simulation started.')

    for i in range(cycles + 1):
        pop.cycle()
        if pop.is_extinct():
            break
        if not i % 10:
            print(i, len(pop.genomes), sep="\t")

    return pop


def parse_args():
    parser = argparse.ArgumentParser("AEGIS")
    parser.add_argument("cycles", type=int)
    parser.add_argument("-f", "--filename", type=str)
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()

    pop = run(cycles=args.cycles, filename=args.filename)
    if not pop.is_extinct():
        pickle_name = save_pop(pop)
        plot(pop, pickle_name)
