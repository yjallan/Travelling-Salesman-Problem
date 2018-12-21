import argparse, sys

parser=argparse.ArgumentParser()

parser.add_argument('-inst', '--fname', help='File instance to run')
parser.add_argument('-alg', '--algo', help='Algorithm to be used')
parser.add_argument('-time', '--cutoff', type=int, help='Cutoff time in seconds')
parser.add_argument('-seed', '--seed', type=int, help='Random seed for local search algorithm')

args=parser.parse_args()

Algos=['BnB', 'Approx', 'LS1', 'LS2']

# Raise all possible argument errors
if (args.fname is None) or (args.algo is None) or (args.cutoff is None):
    print('Desired format: python -inst <city_name> -alg [BnB | Approx | LS1 | LS2] -time <cutoff_in_seconds>[-seed <random_seed for LS>]')
    raise IOError('Not enough argument!')
if (args.algo=='LS1' or args.algo=='LS2') and args.seed is None:
    raise IOError('Random seed missing for Local Search algorithm')
if not args.algo in Algos:
    raise TypeError('Algorithm should be chosen from: BnB | Approx | LS1 | LS2')
if not isinstance(args.cutoff, int):
    raise TypeError('Time Cutoff should be an integer')

if args.algo=='Approx':
    # from Approx import *
    from Approx import *
    Approx(args.fname, args.cutoff)
elif args.algo=="BnB":
    from BnB import *
    BranchAndBound(args.fname, args.cutoff)
elif args.algo=='LS1':
    from LS1 import *
    LS1(args.fname, args.cutoff, args.seed)
elif args.algo=='LS2':
    from LS2 import *
    LS2(args.fname, args.cutoff, args.seed)
