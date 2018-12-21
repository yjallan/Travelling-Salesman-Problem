All codes are implemented in Python 3.

The folders and files should be located as shown below:

OUTPUT
	-{to store all output files: .sol and .trace}

DATA
	-{all data files in .tsp format}
	
Run instructions:

General Format:

python -inst <city_name> -alg [BnB | Approx | LS1 | LS2] -time <cutoff_in_seconds>[-seed <random_seed for LS>]

Example: 

1.Branch and Bound (BnB) or approx run for 600 seconds time cut-off for Cincinnati:

1.a. BnB: python -inst Cincinnati -alg BnB -time 600

1.b. Approx: python -inst Cincinnati -alg Approx -time 600


2.Local Search algorithms: Simulated Annealing (LS1) and Neighborhood 2-opt Exchange (LS2) for random seed 10

2.a. LS1: python -inst Cincinnati -alg LS1 -time 600 -seed 10

2.b. LS2: python -inst Cincinnati -alg LS2 -time 600 -seed 10