"""."""
import argparse
import logging
import json
import pickle
import time
import numpy as np

import population
from funcs import aux

time_start = time.time()


def get_time_estimations():
    time_diff = time.time() - time_start
    seconds_per_100 = time_diff / aux.stage * 100

    eta = (aux.cycle_num - aux.stage) / 100 * seconds_per_100 / 60

    runtime = time_diff / 3600

    return round(eta, 2), round(seconds_per_100, 2), round(runtime, 2)


def main():

    # set args
    parser = argparse.ArgumentParser("Job")
    parser.add_argument("params_extra", type=str)
    parser.add_argument("jobid", type=str)
    parser.add_argument("dirname", type=str)
    args = parser.parse_args()

    # init params
    aux.init_on_start(
        path_default=aux.PATH / "config_default.yml",
        params_extra=json.loads(args.params_extra),
        recpath=aux.PATH.parent / "experiments" / args.dirname / args.jobid,
    )

    # init population
    if hasattr(aux, "reload_population"):
        with open(aux.reload_population, "rb") as ifile:
            pop = pickle.load(ifile)
    else:
        pop = population.Population()

    # start simulation
    logging.info(args.params_extra)
    logging.info("Simulation started")
    logging.info("stage n | pop cnt | egg cnt | eta/min | s per 100 | runtime")

    for _ in range(aux.cycle_num): # simulation cycle

        pop()

        if not aux.envmap is None and aux.stage % aux.envmap_rate == 0:
            aux.evolve_envmap()

        aux.stage += 1
        aux.season_countdown -= 1

        if aux.stage % aux.pickle_rate == 0:
            aux.pickle_pop(pop)

        pop_n = len(pop)
        eggs_n = len(pop.eggs.genomes) if hasattr(pop.eggs, "genomes") else 0

        if pop_n + eggs_n == 0:  # if extinct
            break

        if aux.stage % aux.logging_rate == 0:
            eta, sper100, runtime = get_time_estimations()
            logging.info(
                f"{aux.stage:>7} | {pop_n:>7} | {eggs_n:>7} | {eta:7} | {sper100:9} | {runtime:7}"
            )

    # kill and record the rest
    mask_kill = np.ones(len(pop), bool)
    pop._kill(mask_kill, "end_of_sim")
    aux.recorder.flush()

    logging.info("Simulation finished")


if __name__ == "__main__":
    main()
