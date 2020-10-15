import numpy as np


class Interpreter:
    def __init__(self, BITS_PER_LOCUS):
        self.BITS_PER_LOCUS = BITS_PER_LOCUS
        bits_per_chromosome = self.BITS_PER_LOCUS // 2
        self.binary_weights = np.repeat(
            2 ** np.arange(bits_per_chromosome), 2
        )  # [ 1  1  2  2  4  4  8  8 ...]
        self.binary_max = self.binary_weights.sum()

        self.binary_switch_weights = self.binary_weights.copy()
        self.binary_switch_weights[-2:] = 0
        self.binary_switch_max = self.binary_switch_weights.sum()

        self.interpreter_map = {
            "decimal": self._decimal,
            "exp": self._exp,
            "binary": self._binary,
            "binary_exp": self._binary_exp,
            "binary_switch": self._binary_switch,
            "switch": self._switch,
        }

    def _switch(self, loci):
        """Gene is evaluated at 0, 1 or randomly evaluated at 0 or 1."""
        sums = loci.mean(2)
        rand_values = np.random.random(loci.shape[:-1]) < 0.5
        return np.select(
            [sums == 0, (sums > 0) & (sums < 1), sums == 1], [0, rand_values, 1]
        )

    def _binary_switch(self, loci):
        """Fundamentally a binary interpreter with switch loci that turn on or completely turn off the gene."""
        where_on = loci[:, :, -2:].sum(2) == 2
        # assert all(loci[0, 0, -2:]) == where_on[0, 0]
        values = np.zeros(loci.shape[:-1], float)
        values[where_on] = (
            loci[where_on].dot(self.binary_switch_weights) / self.binary_switch_max
        )
        return values

    def _decimal(self, loci):  # applicable to surv, repr, neut and muta
        return loci.sum(-1) / self.BITS_PER_LOCUS

    def _exp(self, loci):  # applicable to muta
        return 1 / 2 ** np.sum(loci, axis=1)

    def _binary(self, loci):  # applicable to surv, repr
        return loci.dot(self.binary_weights) / self.binary_max

    def _binary_exp(self, loci):  # applicable to muta
        binary = self._binary(loci)
        return 0.98 ** binary

    def __call__(self, loci, interpreter_kind):
        interpreter = self.interpreter_map[interpreter_kind]
        return interpreter(loci)

