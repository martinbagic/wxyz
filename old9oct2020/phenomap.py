import numpy as np


class Phenomap:
    def __init__(self, phenomap_plus, pos_end):
        self.m = np.diag([1.0] * pos_end)

        for phenoage, genoage, val in phenomap_plus:
            self.m[phenoage, genoage] = val

    def __call__(self, genomes, indices):
        weighted_genes = genomes * self.m[indices, :, None]
        summed_genes = weighted_genes.sum(1)
        clipped_genes = np.clip(summed_genes, 0, 1)
        return clipped_genes
