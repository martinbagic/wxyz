import numpy as np


class Phenomap:
    def __init__(self, PHENOMAP_PLUS, pos_end):
        self.map_ = np.diag([1.0] * pos_end)
        for geno_i, pheno_i, weight in PHENOMAP_PLUS:
            self.map_[geno_i, pheno_i] = weight

    def __call__(self, probs):
        return np.clip(probs.dot(self.map_), 0, 1)

    # def __call__(self, genomes):
    #     # values = genomes.dot(self.m)
    #     # print(genomes[0])
    #     # print(values[0])
    #     return genomes

    # def __call__(self, genomes, indices):
    #     weighted_genes = genomes * self.m[indices, :, None]
    #     summed_genes = weighted_genes.sum(1)
    #     clipped_genes = np.clip(summed_genes, 0, 1)
    #     return clipped_genes


class PhenomapFake:
    def __call__(self, probs):
        return probs
