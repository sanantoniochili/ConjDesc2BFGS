import os
import sys
import math
import mmap
import argparse
import collections
import pandas as pd
from ase.io import read as aread

class Info:
	def __init__(self, file, catg):  # create dict
		self.file = file
		self.catg = catg

def read_iters(args):
	flist = []
	dirs = [d for d in os.listdir(args.test_dir) # find all directory objects
			if os.path.isdir(os.path.join(args.test_dir,d))] # check if is directory
	for d in dirs: # for every method folder
		d += '/output'
		if not os.path.exists(os.path.join(args.test_dir,d)):
			continue
		path = args.test_dir+'/'+d
		sampledirs = [d_ for d_ in os.listdir(path) # find all directory objects
			if os.path.isdir(os.path.join(path,d_))] # check if is directory
		for r in sampledirs: # rattled or random
			rpath = os.path.join(path,r)
			flist += [os.path.join(rpath,file) # list all files
						for file in os.listdir(rpath) if file.endswith(".got")]
	count=0
	ofilename = args.test_dir+'/'+args.ofilename
	# Get energies and gnorms for all files 
	#	in method directory in one dataframe
	for filename in flist:
		count_stepc=0 # count step size iters
		with open(filename, 'r') as file:
			info = Info(file, {})
			info.catg['structure'] = [
				filename.split('/')[-1].split('.')[0]]
			info.catg['method'] = filename.split('/')[-4]
			info.catg['folder'] = filename.split('/')[-2]
			df = pd.DataFrame.from_dict(info.catg, orient='columns')

			Es = [] # lists for each file
			Gns = []
			Steps = []
			Inters = []
			iflag = True # keep only first interatomic potential

			for line in file:
				if args.step: # keeping records of step sizes and respective energy levels
					if "new    ..." in line:
						line_ = line.split('...')[1].rstrip().lstrip(' ').split(' ') # remove first part of line
						nline_ = list(filter(None, line_)) # remove blanks
						step = nline_[0]
						energy = nline_[-1]
						if ("*" in energy[1:]):
							step = energy.split('*')[0]
							energy = -math.inf
						elif ("-" in energy[1:]): # take care of overlapping values
							both = energy.lstrip('-')
							step = both.split('-')[0]
							if energy[0] == "-":
								step = "-"+step
							energy = "-"+both.split('-')[-1]

						Es.append(energy)
						Steps.append(step)
						count_stepc += 1

				if "Cycle" in line:
					if args.energy:
						energy = line.split(':')[2].split(' ')[-3]
					if args.gnorm:
						gnorm = line.split(':')[3].split(' ')[-3]
					# if args.energy:
						# if "**" not in energy:
							# Es.append(energy)
					if args.gnorm:
						if "**" not in gnorm:
							Gns.append(gnorm)
					count_stepc=0 # step is stabilized

				if "Interatomic potentials     =" in line:
					if args.interatomic & iflag:
						Inters.append(line.split()[-2])
						iflag = False

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

			if args.interatomic:
				dfsei = pd.DataFrame(Inters).T
				dfsei_ = df.join(dfsei)
				dfsei_ = dfsei_.set_index(['structure', 'method'])

			# dfsi = pd.DataFrame(Stepi).T
			# dfsi_ = df.join(dfsi)
			# dfsi_ = dfsi_.set_index(['structure', 'method'])

		''' Merge dataframes '''
		if count:
			if args.energy:
				dfes = pd.concat([dfes,dfe_], axis=0, sort=False)
			if args.gnorm:
				dfgs = pd.concat([dfgs,dfg_], axis=0, sort=False)
			if args.step:
				dfss = pd.concat([dfss,dfs_], axis=0, sort=False)
			if args.interatomic:
				dfseis = pd.concat([dfseis,dfsei_], axis=0, sort=False)
			# dfsis = pd.concat([dfsis,dfsi_], axis=0, sort=False)												
		else: # initialise
			if args.energy:
				dfes = dfe_
			if args.gnorm:
				dfgs = dfg_
			if args.step:
				dfss = dfs_
			if args.interatomic:
				dfseis = dfsei_
			# dfsis = dfsi_
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

	if args.interatomic:
		with open(args.test_dir+'/'+args.ofilename+'_interatomic.csv', 'w') as f:
			dfseis.to_csv(f, header=True)

		# with open(args.test_dir+'/'+args.ofilename+'_stepi.csv', 'w') as f:
		# 	dfsis.to_csv(f, header=True)

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

import csv
from pathlib import Path
def read_positions(args):
	import sys
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

	# # Initial cif files
	# print("Give data directory path:")
	# DATAPATH = input()
	# print("Give directories with datafiles:")
	# data_flist = [folder for folder in input().split(',')]

	# df_init = pd.DataFrame()
	# print("Creating dataframes to keep initial positions...")
	# for folder in data_flist:
	# 	folder_path = DATAPATH+'/'+folder
	# 	for path in Path(folder_path).rglob('*.cif'): # every structure
	# 		struct_dict = {'structure':path.name, 
	# 			'folder':folder}
	# 		atoms = aread(path)
	# 		positions = atoms.get_positions()
	# 		count_ions = 0
	# 		for ion in atoms.get_chemical_symbols():
	# 			struct_dict[ion+str(count_ions)] = tuple(list(positions[count_ions,:]))
	# 			count_ions += 1
	# 		df_init = df_init.append([struct_dict], ignore_index=True, sort=False)
	# 		# print("File "+path.name+" is done.")
	# df_init = df_init.set_index('structure')
	# df_init.to_csv(args.test_dir+"/temp.csv")

	df_init = pd.read_csv(args.test_dir+"/temp.csv", index_col=['structure','folder'])
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

		# Empty Dataframe
		df = pd.DataFrame()
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

					# Define a Dataframe for each structure
					fpath = rpath+'/'+d
					name,method,folder = fpath.split('/')[-1].split('.')[0],fpath.split('/')[-4],fpath.split('/')[-2]
					# index = [[name], # structure name
					# 		[folder], # random or rattled
					# 		[method]] # method e.g. bfgs

					# Prepare Dataframe
					columns = df_init.columns
					traj_list = sorted([f for f in os.listdir(fpath) if "grs" in f], key=alphanum_key)
					ncolumns = pd.MultiIndex.from_product([['initial']+traj_list,columns]).append(pd.Index(['structure','folder','method']))
					df_struct = pd.DataFrame(columns=ncolumns)
					df_struct = df_struct.set_index(['structure','folder','method'])

					# with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
					    # print(df_struct)

					# Add initial positions
					for col in df_init:
						df_struct.at[(name,folder,method),('initial',col)] = df_init.loc[map_from,folder][col]

					# Add intermediate positions
					count = 0
					for file in traj_list:
						filename = fpath+'/'+file
						atoms = read_gulp(filename)
						positions = atoms.get_positions()
						count_ions = 0
						for ion in atoms.get_chemical_symbols():
							col = ion+str(count_ions) 
							df_struct.loc[(name,folder,method),(file,col)] = \
											tuple(list(positions[count_ions,:]))
							count_ions+=1

					# with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
					#     print(df_struct)
					
					# Add to method dataframe
					if df.empty:
						df = df_struct
					else:
						df = df.append(df_struct)
						# df = df.reset_index(drop=True)
		df.to_csv(args.test_dir+'/'+args.ofilename+'_'+method+'.csv')

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
		'-e', '--energy', action='store_true',
		help='Records of final energy')
	parser.add_argument(
		'-g', '--gnorm', action='store_true',
		help='Records of final gradient norm')
	parser.add_argument(
		'-s', '--step', action='store_true',
		help='Records of step sizes')
	parser.add_argument(
		'-i', '--interatomic', action='store_true',
		help='Records of interatomic energy values')
	parser.add_argument(
		'--all', action='store_true',
		help='Records of all values')
	args = parser.parse_args()

	if args.all:
		args.energy = True
		args.gnorm  = True
		args.step   = True
		args.interatomic  = True

	read_positions(args)