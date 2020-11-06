import numpy as np


class Interpreter:
    def __init__(self, BITS_PER_LOCUS, REPR_MODE):
        self.ploidy = {"sexual": 2, "asexual": 1, "diasexual": 2}[REPR_MODE]

        self.binary_weights = 2 ** np.arange(BITS_PER_LOCUS // self.ploidy)[::-1]
        self.binary_max = self.binary_weights.sum()

        self.binary_switch_weights = self.binary_weights.copy()
        self.binary_switch_weights[-1:] = 0
        self.binary_switch_max = self.binary_switch_weights.sum()

        self.interpreter_map = {
            "uniform": self._uniform,
            "exp": self._exp,
            "binary": self._binary,
            "binary_exp": self._binary_exp,
            "binary_switch": self._binary_switch,
            "switch": self._switch,
            # "simple_1": self._simple_1,
        }

    # def _simple_1(self, loci):
    #     """Locus is evaluated as 1 if its first bit is 1. Otherwise it is 0."""
    #     print(loci[:,0])
    #     print(loci[:, 0] == True)
    #     return loci[:, 0] == 1

    def _switch(self, loci):
        """Locus is evaluated as 0, 1 or randomly evaluated as 0 or 1."""
        sums = loci.mean(2)
        rand_values = np.random.random(loci.shape[:-1]) < 0.5
        return np.select(
            [sums == 0, (sums > 0) & (sums < 1), sums == 1], [0, rand_values, 1]
        )

    def _binary_switch(self, loci):
        """Fundamentally a binary interpreter with switch loci that turn on or completely turn off the gene."""
        where_on = loci[:, :, -1:] == 1
        where_on = where_on.reshape(loci.shape[:-1])
        # assert all(loci[0, 0, -2:]) == where_on[0, 0]
        values = np.zeros(loci.shape[:-1], float)
        values[where_on] = (
            loci[where_on].dot(self.binary_switch_weights) / self.binary_switch_max
        )
        return values

    def _uniform(self, loci):  # applicable to surv, repr, neut and muta
        return loci.sum(-1) / loci.shape[-1]

    def _exp(self, loci):  # applicable to muta
        return 1 / 2 ** np.sum(~loci, axis=1)

    def _binary(self, loci):  # applicable to surv, repr
        return loci.dot(self.binary_weights) / self.binary_max

    def _binary_exp(self, loci):  # applicable to muta
        binary = self._binary(loci)
        return 0.98 ** binary

    def __call__(self, loci, interpreter_kind):
        interpreter = self.interpreter_map[interpreter_kind]
        if self.ploidy == 2:
            loci = np.logical_and(loci[:, :, ::2], loci[:, :, 1::2])
        interpretome = interpreter(loci)
        return interpretome

