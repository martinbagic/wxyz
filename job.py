import json
import argparse
import logging

import population
import funcs


class Job:
    def __init__(self):
        """
        Reads cmd arguments.
        Prepares Config instance for Population.
        Sets the output path for the result of this job.
        """
        self.set_args()
        self.set_conf()
        self.set_opath()

    def set_args(self):
        parser = argparse.ArgumentParser("Job")
        parser.add_argument("params", type=str)
        parser.add_argument("jobid", type=str)
        parser.add_argument("dirname", type=str)
        self.args = parser.parse_args()

    def set_conf(self):
        d = json.loads(self.args.params)
        self.conf = funcs.init_config(d)

    def set_opath(self):
        self.opath = funcs.path.parents[0] / self.args.dirname / (self.args.jobid + '.csv')

    def run(self):
        logging.info("calculating...")
        pop = population.Population(self.args.jobid, self.conf)
        for i in range(self.conf.cycle_num):
            pop.cycle()
            if not i % 100:
                logging.info(f"...stage {i:<5} | n {len(pop.genomes)}")
        pop.killall()
        logging.info("...done!")

        logging.info("writing results...")
        pop.record.save(self.opath)
        logging.info("...done!")


if __name__ == "__main__":
    job = Job()
    job.run()
