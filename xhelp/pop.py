import numpy as np


class Pop:
    attrs = (
        "genomes",
        "ages",
        "origins",
        "uids",
        "births",
        "birthdays",
        "evalomes",
    )

    def __init__(self, genomes, ages, origins, uids, births, birthdays, evalomes):
        self.genomes = genomes
        self.ages = ages
        self.origins = origins
        self.uids = uids
        self.births = births
        self.birthdays = birthdays
        self.evalomes = evalomes

    def __len__(self):
        """Return the number of living individuals."""
        n = len(self.genomes)
        assert all(len(getattr(self, attr)) == n for attr in self.attrs), " ".join(
            (str(len(getattr(self, attr))) for attr in self.attrs)
        )
        return n

    def __getitem__(self, index):
        """Return a subpopulation."""
        return Pop(
            genomes=self.genomes[index],
            ages=self.ages[index],
            origins=self.origins[index],
            uids=self.uids[index],
            births=self.births[index],
            birthdays=self.birthdays[index],
            evalomes=self.evalomes[index],
        )

    def __imul__(self, index):
        """Redefine itself as its own subpopulation."""
        for attr in self.attrs:
            setattr(self, attr, getattr(self, attr)[index])
        return self

    def __iadd__(self, pop):
        """Add another population to itself."""
        for attr in self.attrs:
            val = np.concatenate([getattr(self, attr), getattr(pop, attr)])
            setattr(self, attr, val)
        return self
    
        
