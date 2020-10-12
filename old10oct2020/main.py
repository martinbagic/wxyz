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


def get_time_estimations(stage, CYCLE_NUM):
    def get_dhm(timediff):
        d = int(timediff / 86400)
        timediff %= 86400
        h = int(timediff / 3600)
        timediff %= 3600
        m = int(timediff / 60)
        return f"{d}`{h:02}:{m:02}"

    time_diff = time.time() - time_start

    seconds_per_100 = time_diff / stage * 100
    eta = (CYCLE_NUM - stage) / 100 * seconds_per_100

    runtime = get_dhm(time_diff)
    time_per_1M = get_dhm(time_diff / stage * 1000000)
    eta = get_dhm(eta)

    return eta, time_per_1M, runtime


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
        pass
        # with open(aux.reload_population, "rb") as ifile:
        #     biosystem = pickle.load(ifile)
    else:
        biosystem = population.Biosystem(aux=aux)

    # run simulation
    logging.info(f"Extra parameters: {repr(args.params_extra)}")
    logging.info("Simulation started")
    logging.info("|  stage  |  # pop  |  # egg  |   ETA   |   t1M   | runtime |")

    for _ in range(aux.CYCLE_NUM):

        biosystem.cycle()

        aux.stage += 1

        if aux.stage % aux.PICKLE_RATE == 0:
            aux.recorder.pickle_pop(biosystem, aux.stage)

        pop_n = len(biosystem.pop)
        eggs_n = len(biosystem.eggs) if hasattr(biosystem.eggs, "genomes") else 0

        if pop_n + eggs_n == 0:  # if extinct
            break

        if aux.stage % aux.LOGGING_RATE == 0:
            eta, sper1M, runtime = get_time_estimations(aux.stage, aux.CYCLE_NUM)
            logging.info(
                f"| {aux.stage:>7} | {pop_n:>7} | {eggs_n:>7} | {eta:7} | {sper1M:>7} | {runtime:>6} |"
            )

    # kill and record the rest
    mask_kill = np.ones(len(biosystem.pop), bool)
    biosystem._kill(mask_kill, "end_of_sim")
    aux.recorder.flush()
    logging.info("Simulation finished")


if __name__ == "__main__":
    main()
