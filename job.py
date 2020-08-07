import argparse
import logging
import json
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
        self.opath = (
            funcs.path.parents[0] / "experiments" / self.args.dirname / self.args.jobid
        )

    def run(self):
        logging.info("calculating...")
        logging.info(f"{self.conf.cycle_num} stages...")
        pop = population.Population(self.conf, self.opath)
        for i in range(self.conf.cycle_num):
            pop.cycle()
            if i % 100 == 0 and len(pop.genomes) > 0:
                logging.info(f"...stage {i:<6} | n {len(pop.genomes)}")

        pop.killall()

        # pop.record.compress_batch()
        # pop.record.compress_output()

        logging.info("...done!")

        # logging.info("writing results...")
        # pop.record.save(self.opath)
        # logging.info("...done!")


if __name__ == "__main__":
    job = Job()
    job.run()
