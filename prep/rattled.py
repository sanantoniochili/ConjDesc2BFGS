from ase import *
import sys
import argparse
import numpy as np
from ase.io import read, write
from ase.visualize import view
from ase.calculators.gulp import GULP

DATAPATH = "../../../Data/"

''' Get stdev from user '''
parser = argparse.ArgumentParser(
    description='Define input.')
parser.add_argument(
    'stdev', metavar='--stdev', type=float,
    help='Standard deviation')
parser.add_argument(
    'draw', metavar='--draw',
    help='Number of same derivation trial')
parser.add_argument(
    'out', metavar='--outfolder',
    help='Name of ouput folder')
args = parser.parse_args()

# print("stdev: "+str(args.stdev))
print("draw: "+str(args.draw))

# np.random.randint(low=1, high=12, size=(2, 4))
atoms = Atoms("SrTiO3",

              cell=[[4.00, 0.00, 0.00],
                    [0.00, 4.00, 0.00],
                    [0.00, 0.00, 4.00]],

              positions=[[0, 0, 0],
                         [2, 2, 2],
                         [0, 2, 2],
                         [2, 0, 2],
                         [2, 2, 0]],
              pbc=True)
atoms = atoms.repeat((3, 1, 1))

''' Perturb atoms '''
atoms.rattle(stdev=args.stdev, rng=np.random.RandomState(np.random.seed()))

import os
if not os.path.isdir(DATAPATH+args.out):
  os.mkdir(DATAPATH+args.out)
  print("Created "+DATAPATH+args.out+" directory.")
if not os.path.isdir(DATAPATH+args.out+"/dev"+str(args.stdev)):
  os.mkdir(DATAPATH+args.out+"/dev"+str(args.stdev))
  print("Created "+DATAPATH+args.out+"/dev"+str(args.stdev)+" directory.")

''' Write to .cif files '''
filename = DATAPATH+args.out+"/dev"+str(args.stdev)+"/dev"+str(args.stdev)+"draw"+args.draw+".cif"
if not os.path.isfile(filename):
  write(filename, atoms, format='cif')
else:
  print("File exists.")
