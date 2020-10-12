import numpy as np


class Overshoot:
    def __init__(self, overshoot_handling, max_population_size, bottleneck_size):
        self.max_population_size = max_population_size
        self.bottleneck_size = bottleneck_size
        self.func = {
            "treadmill_random": self.treadmill_random,
            "bottleneck": self.bottleneck,
            "treadmill_ageist": self.treadmill_ageist,
            "treadmill_boomer": self.treadmill_boomer,
            "limitless": self.limitless,
            "starvation": self.starvation,
        }[overshoot_handling]

        self.consecutive_overshoot_n = 0  # for starvation mode

    def __call__(self, n):
        if n <= self.max_population_size:
            self.consecutive_overshoot_n = 0
            return np.zeros(n, bool)
        else:
            self.consecutive_overshoot_n += 1
            return self.func(n)

    def starvation(self, n):
        """Let starved individuals die."""
        if n > self.max_population_size:
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
        indices = np.random.choice(n, n - self.max_population_size, replace=False)
        mask = np.zeros(n, bool)
        mask[indices] = True
        return mask

    def bottleneck(self, n):
        # TODO: reduce to treadmill_random by passing CONF.bottleneck_size to aux function
        """Kill all but random few."""
        indices = np.random.choice(n, self.bottleneck_size, replace=False)
        mask = np.ones(n, bool)
        mask[indices] = False
        return mask

    def treadmill_ageist(self, n):
        """Kill the oldest."""
        mask = np.ones(n, bool)
        mask[: self.max_population_size] = False
        return mask

    def treadmill_boomer(self, n):
        "Kill the youngest."
        mask = np.ones(n, bool)
        mask[-self.max_population_size :] = False
        return mask
