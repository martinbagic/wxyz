# Instructions
1. Prepare experiment
   1. make sure wxyz is updated
   2. mkdir experiment folder parallel to wxyz
   3. prepare cartesian.yml and put into folder
2. Run experiment
   1. run `orderall.sh`
   2. zip data
   3. run analyzer `Analyzer(dirname).run()`
   4. download cartesian.analyzer
3. Explore


# To-dos
- [ ] do few experiments with long runtimes to check when values stabilize
- [ ] do few repeated experiments to check the stochasticity



# how to run locally
```
python3 job.py '{"max_population_size": 2000, "surv_bound_lo": 0, "repr_bound_hi": 1, "overflow_handling": "treadmill_random", "cycle_num": 2000, "mutrate_distribution": 0.95, "genome_distribution": 1}' 1 experiment2
```