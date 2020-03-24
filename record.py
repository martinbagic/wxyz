import yaml
import numpy as np
import time


class Record:

    time_format = "%H:%M:%S"

    def __init__(self, identifier, config):
        self.d = {
            "identifier": identifier,
            "config": config,
            "data": {},
            "time_start": time.strftime(Record.time_format),
        }

        self.time0 = time.time()

    def __call__(self, pop):
        pass

    def save(self):

        time_difference = time.time() - self.time0
        self.d["time_run"]: time.strftime("%H:%M:%S", time.gmtime(time_difference))

        with open(f"records/{self.d['identifier']}", "w") as f:
            yaml.dump(self.d, f)
