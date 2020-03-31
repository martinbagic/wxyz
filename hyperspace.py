import config
import population
import itertools
import os

Cartesian = {
    "max_population_size": [1000, 2000, 5000, 10000],
    "surv_bound_lo": [0, 0.5, 0.9],
    "repr_bound_hi": [0.5, 0.9],
    "mutation_ratio": [1.5, 1, 0.1, 0.01],
    "bottleneck_size": [100, 900],
}

Cartesian = {
    "max_population_size": [2000],
    "surv_bound_lo": [0.5],
    "repr_bound_hi": [0.5],
    "mutation_ratio": [1],
}

cycle_num = 3000

os.mkdir('records-cluster')

carts = tuple(itertools.product(*Cartesian.values()))

for i, cart in enumerate(carts):
    d = dict(zip(Cartesian.keys(), cart))
    conf = config.init(d)
    print(conf)
    pop = population.Population(i, conf)

    for cycle in range(cycle_num):
        if not cycle % 200:
            print(i, len(carts), cycle)
        pop.cycle()
        if pop.is_extinct():
            break

    pop.record.save()
