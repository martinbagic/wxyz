import numpy as np


class Envmap:
    def __init__(self, BITS_PER_LOCUS, total_loci, ENVMAP_RATE):
        self.map_ = np.zeros((total_loci, BITS_PER_LOCUS), bool)
        self.envmap_rate = ENVMAP_RATE

    def __call__(self, genomes):
        return np.logical_xor(self.map_, genomes)

    def evolve(self, stage):
        if stage % self.envmap_rate == 0:
            locus = np.random.choice(np.arange(self.map_.shape[0]))
            bit = np.random.choice(np.arange(self.map_.shape[1]))
            self.map_[locus, bit] = ~self.map_[locus, bit]


class EnvmapFake:
    def __call__(self, genomes):
        return genomes

    def evolve(self, stage):
        pass
