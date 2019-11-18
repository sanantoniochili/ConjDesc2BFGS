import sys
import argparse
import fileinput
from ase import *
from ase.io import read
from ase.visualize import view
from ase.calculators.gulp import GULP
from ase.optimize.sciopt import SciPyFminBFGS, SciPyFminCG


class Method:
	'''
		Initialisation
		name:       Name of method
		keywords:   Keywords for GULP calculator
		library:    Library of potential
	'''

	def __init__(self, name, keywords=[], options=[], library=None):
		self.name = name
		self.keywords = keywords
		self.options = options
		self.library = library

	''' Crystal structure positions and forces '''

	def set_atoms(self, atoms):
		self.atoms = atoms

	''' Set GULP calculator '''

	def set_calc(self):
		keywords_ = ['opti', 'c6']
		keywords_.append(self.name)
		if len(self.keywords):
			self.keywords = keywords_ + self.keywords
		self.keywords = ' '.join(self.keywords)
		calc = GULP(keywords=self.keywords,
					options=self.options,
					library=self.library)
		self.atoms.set_calculator(calc)

	''' Calculate energy '''

	def calc(self):
		return self.atoms.get_potential_energy()


if __name__ == "__main__":
	''' Get input file and method to use from user '''
	parser = argparse.ArgumentParser(
		description='Define input')
	parser.add_argument(
		'method', metavar='--method', type=str,
		help='Algorithm to use')
	parser.add_argument(
		'ifilename', metavar='--input', type=str,
		help='Name of input file')
	parser.add_argument(
		'-s', '--switch',
		action='store_true',
		help='Switch methods')
	parser.add_argument(
		'-t', type=int, default=-1,
		help='Set time_out')
	args = parser.parse_args()

	''' File with structure '''
	print("-----About to use "+args.method +
		  " with input file:" + str(args.ifilename))

	''' Set options '''
	options = []
	if args.switch:
		print("-----Switching methods.")
		# inuser = input("-----Options: ")  # get options from user
		options = ['switch_minimiser bfgs gnorm 0.1']
		# if inuser:
		# options = [inuser]
		print("-----Using options: "+options[0])
	if args.t != -1:
		print("-----Using time_out.")
		options = ['time ' + str(args.t)]
		print("-----Using options: "+options[0])

	with open('temp.txt', 'w') as f:
		for item in options:
			f.write("%s\n" % item)

	''' Set GULP parameters and calculate energy '''
	m = Method(args.method, ['conp', 'full', 'nosymm'], options, 'buck.lib')
	m.set_atoms(read(args.ifilename))
	m.set_calc()  # set GULP
	e = m.calc()  # structure energy
