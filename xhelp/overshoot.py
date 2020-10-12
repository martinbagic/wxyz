import numpy as np


class Overshoot:
    def __init__(self, OVERSHOOT_HANDLING, MAX_POPULATION_SIZE, BOTTLENECK_SIZE):
        self.MAX_POPULATION_SIZE = MAX_POPULATION_SIZE
        self.BOTTLENECK_SIZE = BOTTLENECK_SIZE
        self.func = {
            "treadmill_random": self.treadmill_random,
            "bottleneck": self.bottleneck,
            "treadmill_ageist": self.treadmill_ageist,
            "treadmill_boomer": self.treadmill_boomer,
            "limitless": self.limitless,
            "starvation": self.starvation,
        }[OVERSHOOT_HANDLING]

        self.consecutive_overshoot_n = 0  # for starvation mode

    def __call__(self, n):
        if n <= self.MAX_POPULATION_SIZE:
            self.consecutive_overshoot_n = 0
            return np.zeros(n, bool)
        else:
            self.consecutive_overshoot_n += 1
            return self.func(n)

    def starvation(self, n):
        """Let starved individuals die."""
        if n > self.MAX_POPULATION_SIZE:
            surv_probability = 0.95 ** self.consecutive_overshoot_n
            random_probabilities = np.random.rand(n)
            mask = random_probabilities > surv_probability
        return mask

    def limitless(self, n):
        """Do not kill any individuals."""
        mask = np.zeros(n, bool)
        return mask

    def treadmill_random(self, n):
        """Kill some extra individuals randomly."""
        indices = np.random.choice(n, n - self.MAX_POPULATION_SIZE, replace=False)
        mask = np.zeros(n, bool)
        mask[indices] = True
        return mask

    def bottleneck(self, n):
        # TODO: reduce to treadmill_random by passing CONF.BOTTLENECK_SIZE to aux function
        """Kill all but random few."""
        indices = np.random.choice(n, self.BOTTLENECK_SIZE, replace=False)
        mask = np.ones(n, bool)
        mask[indices] = False
        return mask

    def treadmill_ageist(self, n):
        """Kill the oldest."""
        mask = np.ones(n, bool)
        mask[: self.MAX_POPULATION_SIZE] = False
        return mask

    def treadmill_boomer(self, n):
        "Kill the youngest."
        mask = np.ones(n, bool)
        mask[-self.MAX_POPULATION_SIZE :] = False
        return mask
