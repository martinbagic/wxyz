import numpy as np


class Interpreter:
    def __init__(self, BITS_PER_LOCUS):
        self.BITS_PER_LOCUS = BITS_PER_LOCUS
        bits_per_chromosome = self.BITS_PER_LOCUS // 2
        self.binary_weights = np.repeat(
            2 ** np.arange(bits_per_chromosome), 2
        )  # [ 1  1  2  2  4  4  8  8 ...]
        self.binary_max = self.binary_weights.sum()

        self.interpreter_map = {
            "decimal": self._decimal,
            "exp": self._exp,
            "binary": self._binary,
            "binary_exp": self._binary_exp,
        }

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

