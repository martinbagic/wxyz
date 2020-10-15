import numpy as np


class Overshoot:
    def __init__(self, OVERSHOOT_HANDLING, MAX_POPULATION_SIZE, CLIFF_SURVIVABILITY):
        self.MAX_POPULATION_SIZE = MAX_POPULATION_SIZE
        self.CLIFF_SURVIVABILITY = CLIFF_SURVIVABILITY
        self.func = {
            "treadmill_random": self.treadmill_random,
            "treadmill_boomer": self.treadmill_boomer,
            "treadmill_zoomer": self.treadmill_zoomer,
            "cliff": self.cliff,
            # "limitless": self.limitless,
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

    def cliff(self, n):
        """Kill all but random few."""
        indices = np.random.choice(n, int(self.MAX_POPULATION_SIZE*self.CLIFF_SURVIVABILITY), replace=False)
        mask = np.ones(n, bool)
        mask[indices] = False
        return mask

    def treadmill_boomer(self, n):
        """Kill the oldest. Let youngest live."""
        mask = np.ones(n, bool)
        mask[-self.MAX_POPULATION_SIZE :] = False
        return mask

    def treadmill_zoomer(self, n):
        """Kill the youngest. Let oldest live."""
        mask = np.ones(n, bool)
        mask[: self.MAX_POPULATION_SIZE] = False
        return mask
