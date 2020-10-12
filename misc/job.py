import argparse
import logging
import json
import pickle
import numpy as np

import population
from funcs import aux

def main():

    # set args
    parser = argparse.ArgumentParser("Job")
    parser.add_argument("params_extra", type=str)
    parser.add_argument("jobid", type=str)
    parser.add_argument("dirname", type=str)
    args = parser.parse_args()

    # init params
    aux.init_on_start(
        path_default=aux.project_path / "config_default.yml",
        params_extra=json.loads(args.params_extra),
        recpath=aux.project_path.parent / "experiments" / args.dirname / args.jobid
    )

    # init population
    if hasattr(aux, "reload_population"):
        with open(aux.reload_population, "rb") as ifile:
            self.pop = pickle.load(ifile)
    else:
        self.pop = population.Pop()


def do(self):

    logging.info(f"{aux.cycle_num} stages...")

    for i in range(aux.cycle_num):

        # TODO split generation

        aux.stage += 1

        self.pop.cycle()

        # if aux.split_in == 0:
        #     self.pop.split()
        #     aux.renew_split_in()

        # self.pop.cycle()

        n0 = len(self.pop)
        n1 = len(self.pop.eggs.genomes) if hasattr(self.pop.eggs, "genomes") else 0

        if i % 100 == 0 and n0 + n1 > 0:
            logging.info(f"...stage {i:<6} | indivs {n0:6} | eggs {n1:6}")

    # kill all
    mask_kill = np.ones(len(self.pop), bool)
    self.pop._kill(mask_kill, "end_of_sim")
    aux.recorder.flush()

    logging.info("...done!")

if __name__=='__main__':
    main()


class Job:
    def __init__(self):

        # set args
        parser = argparse.ArgumentParser("Job")
        parser.add_argument("params_extra", type=str)
        parser.add_argument("jobid", type=str)
        parser.add_argument("dirname", type=str)
        args = parser.parse_args()

        # init params
        aux.init_on_start(
            path_default=aux.project_path / "config_default.yml",
            params_extra=json.loads(args.params_extra),
            recpath=aux.project_path.parent / "experiments" / args.dirname / args.jobid
        )

        # init population
        if hasattr(aux, "reload_population"):
            with open(aux.reload_population, "rb") as ifile:
                self.pop = pickle.load(ifile)
        else:
            self.pop = population.Pop()


    def do(self):

        logging.info(f"{aux.cycle_num} stages...")

        for i in range(aux.cycle_num):

            # TODO split generation

            aux.stage += 1

            self.pop.cycle()

            # if aux.split_in == 0:
            #     self.pop.split()
            #     aux.renew_split_in()

            # self.pop.cycle()

            n0 = len(self.pop)
            n1 = len(self.pop.eggs.genomes) if hasattr(self.pop.eggs, "genomes") else 0

            if i % 100 == 0 and n0 + n1 > 0:
                logging.info(f"...stage {i:<6} | indivs {n0:6} | eggs {n1:6}")

        # kill all
        mask_kill = np.ones(len(self.pop), bool)
        self.pop._kill(mask_kill, "end_of_sim")
        aux.recorder.flush()

        logging.info("...done!")


if __name__ == "__main__":
    Job().do()
