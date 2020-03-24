import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pickle
import time
import numpy as np

from population import i1, i2, i3, i4
from population import Population



def plot(self, filename):
    plt.figure(figsize=(8, 4))

    size = len(self.genomes)
    size = size if size < 100 else 100
    indices = np.random.choice(np.arange(len(self.genomes)), size=size, replace=False)
    genomes = self.genomes[indices]
    ages = self.ages[indices]

    array = genomes.sum(2)
    sort_i = np.argsort(ages)
    array = array[sort_i].transpose()
    plt.imshow(array, cmap="bwr", interpolation="nearest")

    for i, c in zip((i4, i3, i2, i1), "rgbm"):
        plt.axvline(-3, 1 - i / i4, 1, lw=5, c=c, solid_capstyle="butt")
    plt.colorbar()
    plt.yticks([])
    plt.xticks([])
    plt.xlim(-5.7, size - 0.5)
    plt.ylim(i4 - 0.5, -0.5)
    plt.xlabel("genome")
    plt.ylabel("locus")
    plt.savefig(f"plots/{filename}.genomes.png", dpi=300, bbox_inches="tight")

    plt.clf()
    plt.figure(figsize=(5, 5))
    plt.hist(ages, bins=len(ages), cumulative=True)
    plt.savefig(f"plots/{filename}.ages.png", dpi=300, bbox_inches="tight")

    plt.clf()
    plt.figure(figsize=(5, 5))
    plt.plot(np.mean(array, axis=1))
    plt.savefig(f"plots/{filename}.quality.png", dpi=300, bbox_inches="tight")

    plt.close()
