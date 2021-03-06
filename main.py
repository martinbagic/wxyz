"""."""
import argparse
import logging
import json
import pickle
import time
import numpy as np
import os
import pathlib

from biosystem import Biosystem
from aux import Aux


# logging settings
logging.basicConfig(
    format="%(asctime)s : %(message)s",
    datefmt="%d/%m/%Y %I:%M:%S",
    level=logging.DEBUG,
)


time_start = time.time()

project_path = pathlib.Path(__file__).absolute().parent


def get_memory_usage(stage, path):
    total_memory_usage = (
        sum(f.stat().st_size for f in path.glob("*") if f.is_file())
        / 1024
        / 1024
        / 1024
    )  # GB
    return total_memory_usage


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

    stages_per_min = int(stage / (time_diff / 60))

    runtime = get_dhm(time_diff)
    time_per_1M = get_dhm(time_diff / stage * 1000000)
    eta = get_dhm(eta)

    return eta, time_per_1M, runtime, stages_per_min


def main():

    # set args
    parser = argparse.ArgumentParser("Job")
    parser.add_argument("params_extra", type=str)
    parser.add_argument("jobid", type=str)
    parser.add_argument("dirname", type=str)
    parser.add_argument("-r", "--reload_biosystem", type=str)
    parser.add_argument("-c", "--config_files", nargs="*", type=str, default=[])
    args = parser.parse_args()

    # init params
    recpath = project_path / "experiments" / args.dirname / args.jobid
    recpath.mkdir(parents=True, exist_ok=True)  # make folder for experiment
    aux = Aux(
        paths_config=[
            project_path / "configs" / cfile
            for cfile in args.config_files + ["config_default.yml"]
        ],
        cmd_params=json.loads(args.params_extra),
        recpath=recpath,
    )

    # init population
    if args.reload_biosystem:
        logging.info(f"Reloading biosystem from: {args.reload_biosystem}")
        with open(project_path / "experiments" / args.reload_biosystem, "rb") as ifile:
            pop = pickle.load(ifile).pop
            biosys = Biosystem(pop=pop, aux=aux)
    else:
        logging.info("Initializing biosystem afresh")
        biosys = Biosystem(aux=aux)

    # run simulation
    logging.info("Simulation started")
    logging_headers = [
        "|         " * 8 + "|",
        "|  stage  |   ETA   |   t1M   | runtime | stg/min | til mem | end mem | recpend |",
    ]
    logging.info(logging_headers[1])

    for _ in range(aux.CYCLE_NUM):

        aux.stage += 1
        
        if aux.stage % aux.recorder.JSON_RATE == 0:
            aux.recorder.rec_json_flag = True
        
        if aux.stage % (aux.recorder.JSON_RATE * aux.recorder.REC_RATE) == 0:
            aux.recorder.rec_flush_flag = True
        
        biosys.cycle()

        if len(biosys) == 0:
            logging.info("Biosystem went extinct")
            aux.recorder.visor_data["extinct"] = True  # record going extinct
            break

        if aux.stage % (aux.LOGGING_RATE * 10) == 0:
            logging.info(logging_headers[0])
            logging.info(logging_headers[1])

        if aux.stage % aux.LOGGING_RATE == 0:
            eta, sper1M, runtime, stpermin = get_time_estimations(
                aux.stage, aux.CYCLE_NUM
            )
            til_memuse = get_memory_usage(aux.stage, aux.recorder.opath)
            end_memuse = (aux.CYCLE_NUM - aux.stage) / aux.stage * til_memuse

            logging.info(
                f"| {aux.stage:>7} | {eta:7} | {sper1M:>7} | {runtime:>6} | {stpermin:>7} | {til_memuse:>7.2f} | {end_memuse:>7.2f} | {len(aux.recorder.genomes):>7} |"
            )

    aux.recorder.pickle_pop(biosys, aux.stage)
    mask_kill = np.ones(len(biosys.pop), bool)
    biosys._kill(mask_kill, "end_of_sim")
    aux.recorder.save(force=True)
    aux.recorder.write_to_visor()
    logging.info("Simulation finished")


if __name__ == "__main__":
    main()
