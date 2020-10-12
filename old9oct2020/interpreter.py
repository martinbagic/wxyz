import numpy as np


class Interpreter:
    def __init__(self, bits_per_locus):
        self.bits_per_locus = bits_per_locus
        bits_per_chromosome = self.bits_per_locus // 2
        self.binary_weights = np.repeat(
            2 ** np.arange(bits_per_chromosome), 2
        )  # [ 1  1  2  2  4  4  8  8 ...]
        self.binary_max = self.binary_weights.sum()

    def _decimal(self, loci):  # applicable to surv, repr, neut and muta
        return loci.sum(1) / self.bits_per_locus

    def _exp(self, loci):  # applicable to muta
        return 1 / 2 ** np.sum(loci, axis=1)

    def _binary(self, loci):  # applicable to surv, repr
        return loci.dot(self.binary_weights) / self.binary_max

    def __call__(self, loci, interpreter_kind):
        # interpreter_kind = {
        #     "surv": aux.interpreter_surv,
        #     "repr": aux.interpreter_repr,
        #     "neut": aux.interpreter_neut,
        #     "muta": aux.interpreter_muta,
        # }[loci_kind]

        interpreter = {
            "decimal": self._decimal,
            "exp": self._exp,
            "binary": self._binary,
        }[interpreter_kind]

        return interpreter(loci)
