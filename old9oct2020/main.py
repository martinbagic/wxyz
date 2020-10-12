"""."""
import argparse
import logging
import json
import pickle
import time
import numpy as np
import population
import funcs

# logging settings
logging.basicConfig(
    format="%(asctime)s : %(message)s",
    datefmt="%d/%m/%Y %I:%M:%S",
    level=logging.DEBUG,
)


time_start = time.time()


def get_time_estimations(stage, cycle_num):
    time_diff = time.time() - time_start
    seconds_per_100 = time_diff / stage * 100

    eta = (cycle_num - stage) / 100 * seconds_per_100 / 60

    runtime = time.strftime('%H:%M:%S', time.gmtime(time_diff))

    return round(eta, 2), round(seconds_per_100, 2), runtime


def main():

    # set args
    parser = argparse.ArgumentParser("Job")
    parser.add_argument("params_extra", type=str)
    parser.add_argument("jobid", type=str)
    parser.add_argument("dirname", type=str)
    args = parser.parse_args()

    # init params
    aux = funcs.Aux(
    # aux.initialize(
        path_default=funcs.project_path / "config_default.yml",
        params_extra=json.loads(args.params_extra),
        recpath=funcs.project_path.parent / "experiments" / args.dirname / args.jobid,
    )

    # init population
    if hasattr(aux, "reload_population"):
        with open(aux.reload_population, "rb") as ifile:
            pop = pickle.load(ifile)
    else:
        pop = population.Population(aux=aux)

    # run simulation
    logging.info(args.params_extra)
    logging.info("Simulation started")
    logging.info("stage n | pop cnt | egg cnt | eta/min | s per 100 |  runtime")

    for _ in range(aux.cycle_num):

        pop.cycle()


        aux.stage += 1

        if aux.stage % aux.pickle_rate == 0:
            aux.recorder.pickle_pop(pop, aux.stage)

        pop_n = len(pop)
        eggs_n = len(pop.eggs.genomes) if hasattr(pop.eggs, "genomes") else 0

        if pop_n + eggs_n == 0:  # if extinct
            break

        if aux.stage % aux.logging_rate == 0:
            eta, sper100, runtime = get_time_estimations(aux.stage, aux.cycle_num)
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
