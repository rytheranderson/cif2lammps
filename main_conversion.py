from __future__ import print_function
import multiprocessing 
from multiprocessing import Pool
import numpy as np
import glob
import os
import re
import sys
import time
from write_lammps_data import lammps_inputs

from force_field_construction import UFF
# add more force field classes here as they are made

force_fields = ['UFF']

def serial_conversion(directory, force_field=UFF, outdir='unopt_lammps_data', charges=False, parallel=False, replication='1x1x1'):

	try:
		os.mkdir(outdir)
	except OSError:
		pass

	print('conversion running serial on a single core')

	cifs = sorted(glob.glob(directory + os.sep + '*.cif'))
	for cif in cifs:
		print('converting ', cif, '...')
		#try:
		lammps_inputs([cif, force_field, outdir, charges, replication])
		#except:
		#	print('ERROR during converion of', cif)
		#	continue

	print('--- cifs in', directory, 'converted and placed in', outdir, '---')

def parallel_conversion(directory, force_field=UFF, outdir='unopt_lammps_data', charges=False, parallel=True, replication='1x1x1'):

	try:
		os.mkdir(outdir)
	except OSError:
		pass

	print('conversion running on ' + str(multiprocessing.cpu_count()) + ' cores')

	cifs = sorted(glob.glob(directory + os.sep + '*.cif'))
	args = [[cif, force_field, outdir, charges, replication] for cif in cifs]
	pool = Pool(multiprocessing.cpu_count())
	results_par = pool.map_async(lammps_inputs, args) 
	pool.close()
	pool.join()

	print('--- cifs in', directory, 'converted and placed in', outdir, '---')

def run_conversion():

	arguments = sys.argv[1:]
	directory = arguments[0]

	optional_arguments = {}
	for arg in arguments[1:]:
		if '--' in arg and 'parallel':
			parse_arg = re.sub('[--]', '', arg).split('=')
			if parse_arg[1] in force_fields:
				if parse_arg[1] == 'UFF':
					value = UFF
				# other options go here as more forcefields are made
			else:
				if parse_arg[1] == 'True':
					value = True
				elif parse_arg[1] == 'False':
					value = False
				else:
					value = parse_arg[1]

			optional_arguments[parse_arg[0]] = value

	try:
		if optional_arguments['parallel']:
			parallel_conversion(directory, **optional_arguments)
		else:
			serial_conversion(directory, **optional_arguments)
			
	except KeyError:
		serial_conversion(directory, **optional_arguments)

start_time = time.time()
if __name__ == '__main__': 
	run_conversion()
print("conversion took %s seconds " % np.round((time.time() - start_time), 3))
