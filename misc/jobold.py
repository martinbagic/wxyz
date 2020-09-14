import argparse
import logging
import json
import population
import funcs


class Job:
    def __init__(self):

        # set args
        parser = argparse.ArgumentParser("Job")
        parser.add_argument("params_extra", type=str)
        parser.add_argument("jobid", type=str)
        parser.add_argument("dirname", type=str)
        args = parser.parse_args()

        # define CONF in funcs
        funcs.Conf(
            path_default=funcs.PATH / "config_default.yml",
            params_extra=json.loads(args.params_extra),
        )

        # set opath
        opath = funcs.PATH.parent / "experiments" / args.dirname / args.jobid

        # init population
        self.pop = population.Population(opath)

    def do(self):

        logging.info(f"{funcs.CONF.cycle_num} stages...")

        for i in range(funcs.CONF.cycle_num):
            self.pop.cycle()

            n0 = len(self.pop.genomes)
            n1 = len(self.pop.nextgen.genomes) if hasattr(self.pop, "nextgen") else 0

            if i % 100 == 0 and n0 + n1 > 0:
                logging.info(f"...stage {i:<6} | indivs {n0:6} | eggs {n1:6}")

        self.pop.killall()

        self.pop.record.flush(force=True)

        logging.info("...done!")


if __name__ == "__main__":
    Job().do()
