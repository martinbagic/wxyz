import numpy as np


class Phenomap:
    def __init__(self, PHENOMAP_PLUS, pos_end):
        self.m = np.diag([1.0] * pos_end)
        for geno_i, pheno_i, weight in PHENOMAP_PLUS:
            self.m[geno_i, pheno_i] = weight

    def transform(self, probs):
        # NOTE: I am not clipping values because probabilities
        #       greater than 1 are guaranteed (as is 1) and
        #       probabilities lesser than 1 are impossible (as is 0)
        return probs.dot(self.m)

    def __call__(self, genomes):
        # values = genomes.dot(self.m)
        # print(genomes[0])
        # print(values[0])
        return genomes

    # def __call__(self, genomes, indices):
    #     weighted_genes = genomes * self.m[indices, :, None]
    #     summed_genes = weighted_genes.sum(1)
    #     clipped_genes = np.clip(summed_genes, 0, 1)
    #     return clipped_genes
