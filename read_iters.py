import os
import sys
import math
import mmap
import argparse
import collections
import pandas as pd
from ase.io import read as aread
from ase.units import eV, Ang
import numpy as np
from ase import Atoms
sys.path.append('/home/sanantoniochili/Desktop/PhD/Scripts/Switch_Implementation/gradients_implementation/cysrc')

LIBFILE = "../../Data/Libraries/buck.lib"
RADFILE = "../../Data/Libraries/radii.lib"

class Info:
	def __init__(self, file, catg):  # create dict
		self.file = file
		self.catg = catg

import re
def tryint(s):
	try:
		return int(s)
	except:
		return s

def alphanum_key(s):
	""" Turn a string into a list of string and number chunks.
		"z23a" -> ["z", 23, "a"]
	"""
	return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def to_cartesian(positions, a, b, c, alpha, beta, gamma,
							   angle_in_degrees=True):
	"""
	Convert scaled to Cartesian coordinates.

	"""
	atoms = Atoms(scaled_positions=positions, cell=[a,b,c,alpha,beta,gamma])
	return atoms.get_positions()

def read_iters(args):
	"""
	Record information about the structural relaxation
	from GULP output files into a csv file.

	"""
	flist = []
	dirs = [d for d in os.listdir(args.test_dir) # find all directory objects
			if os.path.isdir(os.path.join(args.test_dir,d))] # check if is directory
	for d in dirs: # for every method folder
		dpath = os.path.join(args.test_dir,d)
		flist += [str(path) for path in Path(dpath).rglob('*.got')] # list of every got file
	flist = sorted(flist, key=alphanum_key)

	count=0
	ofilename = args.test_dir+'/'+args.ofilename # output file
	structs_dict = {}

	# Get energies and gnorms for all files 
	#	in method directory in one dataframe
	for filename in flist:
		with open(filename, 'r') as file:
			try:
				filename_str = str(filename)
				info = Info(file, {})
				info.catg['structure'] = [
					filename_str.split('/')[-1].split('.')[0]]
				info.catg['method'] = filename_str.split('/')[2]
				info.catg['folder'] = filename_str.split('/')[4]
				info.catg['opt_succ'] = False
				df = pd.DataFrame.from_dict(info.catg, orient='columns')
			except:
				continue

			# lists for each csv file
			Es = []
			Gns = []
			Steps = []
			Inters = []
			Times = []

			symbols = []
			iflag = True # keep only first interatomic potential
			pos_flaf = False
			struct_dict = {}
			positions = []

			# cell parameters
			a, b, c = .0,.0,.0
			alpha, beta, gamma = .0,.0,.0

			lines = file.readlines()
			for i, line in enumerate(lines):
				if args.step: # keeping records of step sizes 
							  # and respective energy levels
					step = None
					energy = None
					if "new    ..." in line:
						line = line.split()[2:]
						# take care of overlapping values
						if len(line) == 1: 
							line = line[0]
							if ("*" in line):
								step = line.split('*')[0]
								energy = -math.inf
							elif ("-" in line): 
								both = line.lstrip('-')
								step = both.split('-')[0]
								# check if stepsize has negative value
								if line[0] == "-": 
									step = "-"+step
								energy = "-"+both.split('-')[-1]
						else:
							step = line[0]
							energy= line[1]

						Es.append(energy)
						Steps.append(step)
						Gns.append(gnorm)
						Times.append(time)

				if "Cycle" in line:
					energy = line.split(':')[2].split()[0]
					gnorm = line.split(':')[3].split()
					time = line.split(':')[4].split()
					if "**" in energy:
						energy = -math.inf
					if "**" in gnorm:
						gnorm = -math.inf

					# to avoid duplicates and initialise lists
					if not args.step or (len(Es)==0): 
						if args.energy:
							Es.append(energy)
						if args.gnorm:
							Gns.append(gnorm)
						if args.time:
							Times.append(time)
					else:
						if args.gnorm:
							Gns[-1] = gnorm
						if args.time:
							Times[-1] = time

				if "Interatomic potentials     =" in line:
					if args.interatomic & iflag:
						Inters.append(line.split()[-2])
						iflag = False

				# read cell parameters
				if ("a =" in line) & ("alpha" in line):
					import re
					a_line = line.replace(" ", "").rstrip('\n').split("alpha")
					a = float(a_line[0].split("=")[-1])
					alpha = float(a_line[1].split("=")[-1])

				elif ("b =" in line) & ("beta" in line):
					import re
					b_line = line.replace(" ", "").rstrip('\n').split("beta")
					b = float(b_line[0].split("=")[-1])
					beta = float(b_line[1].split("=")[-1])

				elif ("c =" in line) & ("gamma" in line):
					import re
					c_line = line.replace(" ", "").rstrip('\n').split("gamma")
					c = float(c_line[0].split("=")[-1])
					gamma = float(c_line[1].split("=")[-1])

				if "Optimisation achieved" in line:
					# Check if optimisation succeeded
					info.catg['opt_succ'] = True

				if ("Final" in line) & ("coordinates" in line):
					iline = i + 5
					scaled_positions = []
					while(True):
						# read positions table
						iline = iline + 1
						if lines[iline].find("------------") != -1:
							break
						xyz = lines[iline].split()[3:6]
						symbol = lines[iline].split()[1]+lines[iline].split()[0]
						XYZ = [float(x) * Ang for x in xyz]
						scaled_positions.append(XYZ)
						symbols.append(symbol)
						# checks for cell vectors
						assert(a*b*c)
						assert(alpha*beta*gamma)
						positions = to_cartesian(np.array(scaled_positions),
							a=a,b=b,c=c,
							alpha=alpha,beta=beta,gamma=gamma)
						# fill dataframe
						i = 0
						for symbol in symbols:
							struct_dict[symbol] = tuple(positions[i,])
							i += 1

			if args.positions:
				key = tuple([info.catg['structure'][0],
					info.catg['method'],
					info.catg['folder']])
				struct_dict['opt_succ'] = info.catg['opt_succ']
				structs_dict[key] = struct_dict

			if args.energy:
				dfe = pd.DataFrame(Es).T
				dfe_ = df.join(dfe)
				dfe_ = dfe_.set_index(['structure', 'method'])

			if args.gnorm:
				dfg = pd.DataFrame(Gns).T
				dfg_ = df.join(dfg)
				dfg_ = dfg_.set_index(['structure', 'method'])

			if args.step:
				dfs = pd.DataFrame(Steps).T
				dfs_ = df.join(dfs)
				dfs_ = dfs_.set_index(['structure', 'method'])

			if args.time:
				dfst = pd.DataFrame(Times).T
				dfst_ = df.join(dfst)
				dfst_ = dfst_.set_index(['structure', 'method'])

			if args.interatomic:
				dfsei = pd.DataFrame(Inters).T
				dfsei_ = df.join(dfsei)
				dfsei_ = dfsei_.set_index(['structure', 'method'])

			# Check if Buckingham catastrophe happened
			if symbols==[]:
				continue
			import potential
			Bpot = potential.Buckingham()
			chemical_symbols = np.array([re.split('([0-9]+)', s)[0] for s in symbols])
			Bpot.set_parameters(
				filename=LIBFILE,
				chemical_symbols=chemical_symbols,
				radius_lib=RADFILE)
			thres_value = 0.5
			check = Bpot.catastrophe_check(positions, thres_value)
			# if catastrophe happened keep the threshold that was used to check
			# the ions' distance
			struct_dict['catastrophe'] = check*thres_value 

		''' Merge dataframes '''
		if count:
			if args.energy:
				dfes = pd.concat([dfes,dfe_], axis=0, sort=False)
			if args.gnorm:
				dfgs = pd.concat([dfgs,dfg_], axis=0, sort=False)
			if args.step:
				dfss = pd.concat([dfss,dfs_], axis=0, sort=False)
			if args.time:
				dfsts = pd.concat([dfsts,dfst_], axis=0, sort=False)
			if args.interatomic:
				dfseis = pd.concat([dfseis,dfsei_], axis=0, sort=False)

		else: # initialise
			if args.energy:
				dfes = dfe_
			if args.gnorm:
				dfgs = dfg_
			if args.step:
				dfss = dfs_
			if args.time:
				dfsts = dfst_
			if args.interatomic:
				dfseis = dfsei_

		count += 1

	# for d in dirs:
	if args.energy:
		with open(args.test_dir+'/'+args.ofilename+'_energy.csv', 'w') as f:
			dfes.to_csv(f, header=True)

	if args.gnorm:
		with open(args.test_dir+'/'+args.ofilename+'_gnorm.csv', 'w') as f:
			dfgs.to_csv(f, header=True)

	if args.step:
		with open(args.test_dir+'/'+args.ofilename+'_step.csv', 'w') as f:
			dfss.to_csv(f, header=True)

	if args.time:
		with open(args.test_dir+'/'+args.ofilename+'_time.csv', 'w') as f:
			dfsts.to_csv(f, header=True)

	if args.interatomic:
		with open(args.test_dir+'/'+args.ofilename+'_interatomic.csv', 'w') as f:
			dfseis.to_csv(f, header=True)

	if args.positions:
		# print(structs_dict)
		import pickle
		pickle.dump( structs_dict, 
			open( args.test_dir+'/'+args.ofilename+"_positions.p", "wb" ) )
		
def read_initial_positions(args):
	# Initial cif files
	print("Give data directory path:")
	DATAPATH = input()
	print("Give directories with datafiles:")
	data_flist = [folder for folder in input().split(',')]

	init_dict = {}
	print("Creating dictionary to keep initial positions...")
	for folder in data_flist:
		folder_path = DATAPATH+'/'+folder
		for path in Path(folder_path).rglob('*.cif'): # every structure
			struct_dict = {}
			atoms = aread(path)
			positions = atoms.get_positions()
			count_ions = 0
			for ion in atoms.get_chemical_symbols():
				struct_dict[ion+str(count_ions)] = tuple(list(positions[count_ions,:]))
				count_ions += 1
			init_dict[(path.name,folder)] = struct_dict
			# print("File "+path.name+" is done.")

	print("Dumped "+str(len(init_dict.keys()))+" structures\' initial \
	positions into "+args.test_dir+"/temp.p")
	import pickle
	pickle.dump( init_dict, open( args.test_dir+"/temp.p", "wb" ) )

import csv
from pathlib import Path
def read_trajectory_positions(args):
	import sys
	import numpy as np
	sys.path.append('/home/sanantoniochili/Desktop/PhD/Scripts/GULP_Python')
	from gulp import read_gulp
	"""This function gathers the positions of the ions from every trajectory file 
	produced along with the initial and final positions from the respective
	.cif files found. 

	The positions are saved into a list preceded by the name 
	of the structure, the kind of the structure ('random' or 'rattled') and the name
	of the relaxation method ('conj','bfgs','switch' or 'switch_o'). For every ion
	in the list there is its chemical symbol with the number of the trajectory file, 
	or 'init' or 'final' for initial and final structures, followed by a 3-tuple 
	marking its 3D position in Euclidean space. 

	It is taken for granted that the order of the ions in the files
	remains tha same and the positions are listed in the same order every time.

	The trajectory files are read in ascending order w.r.t. their numbering.

	"""
	ofilename = args.test_dir+'/'+args.ofilename
	flist = []

	import pickle
	init_dict = pickle.load( open( args.test_dir+"/temp.p", "rb" ) )
	print("Searching for trajectory and final files..")
	dirs = [d for d in os.listdir(args.test_dir) # find all directory objects
			if os.path.isdir(os.path.join(args.test_dir,d))] # check if is directory

	for method in dirs: # for every method folder
		print("For method "+method+":")
		MAP = args.test_dir+'/'+method+'/map_files.txt'
		print("Opened map file to match initial files to output files.")
		# Find trajectory and final files
		dout = method+'/output'
		if not os.path.exists(os.path.join(args.test_dir,dout)):
			continue
		path = args.test_dir+'/'+dout
		sampledirs = [d_ for d_ in os.listdir(path) # find all directory objects
			if os.path.isdir(os.path.join(path,d_))] # check if is directory
		structs_dict = {}

		for r in sampledirs: # rattled or random
			rpath = os.path.join(path,r)
			sdirs = [d for d in os.listdir(rpath) # find all directory objects
			if os.path.isdir(rpath)] # check if is directory
			for d in sdirs: # every structure
				with open(MAP,'r') as mapfile:
					# print("Looking through trajectories of "+d+".",end=" ")

					# Find correct pairs of cif and grs/gin/got files
					# d contains the structure name and r is the kind (random or rattled)
					for line in mapfile:
						if (r+'/'+d+'\n' in line):
							map_from = line.split(":")[0].strip(' ').split('/')[-1]

					# Define a dict for each structure
					fpath = rpath+'/'+d
					name,method,folder = fpath.split('/')[-1].split('.')[0],fpath.split('/')[-4],fpath.split('/')[-2]
					struct_dict = {}

					# Prepare dictionary
					traj_list = sorted([f for f in os.listdir(fpath) if "grs" in f], key=alphanum_key)

					# Add initial positions
					init_struct = init_dict[(map_from,folder)]
					struct_dict['initial'] = {ion : init_struct[ion] for ion in init_struct}
					print(struct_dict)

					# Add intermediate positions
					count = 0
					for file in traj_list:
						filename = fpath+'/'+file
						atoms = read_gulp(filename)
						positions = atoms.get_positions()
						struct_dict[file] = \
						{atoms.get_chemical_symbols()[i]+str(i) : tuple(list(positions[i,:])) \
																			for i in range(len(atoms.positions))}
					# Add structure to all as a dictionary of files
					structs_dict[(name, method, folder)] = struct_dict

		# Print to pickle for every method
		import pickle
		outfile = args.test_dir+'/'+method+'/'+args.ofilename+'_positions.p'
		pickle.dump( structs_dict, open( outfile, "wb" ) )
	print("Done.")					


if __name__ == "__main__":
	''' Get input file and method to use from user '''
	parser = argparse.ArgumentParser(
		description='Define input')
	parser.add_argument(
		'ofilename', metavar='--output', type=str,
		help='.csv file to produce')
	parser.add_argument(
		'test_dir', metavar='--test_folder', type=str,
		help='Define test environment folder')
	parser.add_argument(
		'-t', '--trajectory', action='store_true',
		help='Define if there whole trajectory is to be found.')
	parser.add_argument(
		'-e', '--energy', action='store_true',
		help='Records of final energy')
	parser.add_argument(
		'-g', '--gnorm', action='store_true',
		help='Records of final gradient norm')
	parser.add_argument(
		'-s', '--step', action='store_true',
		help='Records of step sizes')
	parser.add_argument(
		'-ti', '--time', action='store_true',
		help='Records of cpu time')
	parser.add_argument(
		'-i', '--interatomic', action='store_true',
		help='Records of interatomic energy values')
	parser.add_argument(
		'-p', '--positions', action='store_true',
		help='Records of final positions.')
	parser.add_argument(
		'--all', action='store_true',
		help='Records of all values')
	args = parser.parse_args()

	if args.all:
		args.energy = True
		args.gnorm  = True
		args.step   = True
		args.time = True
		args.interatomic  = True
		args.positions = True

	if args.trajectory:
		# read_initial_positions(args)
		read_trajectory_positions(args)
	else:
		read_iters(args)