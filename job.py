import funcs
import json
import os


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
        parser = argparse.ArgumentParser("AEGIS")
        parser.add_argument("params", type=str)
        parser.add_argument("jobid", type=int)
        parser.add_argument("dirname", type=str)
        self.args = parser.parse_args()

    def set_conf(self):
        d = json.loads(self.args.params)
        self.conf = funcs.init_config(d)

    def set_opath(self):
        self.opath = funcs.path.parents[0] / self.args.dirname / self.args.jobid

    def run(self):
        pop = population.Population(i, conf)
        pop.cycle(conf.cycle_num)
        pop.record.save(self.opath)


if __name__ == "__main__":
    job = Job()
    job.run()
