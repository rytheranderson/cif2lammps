from __future__ import print_function
from cif2system import initialize_system, replication_determination, write_cif_from_system
from small_molecule_construction import add_small_molecules
import atomic_data
import os
import numpy as np
import datetime
import math
import warnings
from ase import Atoms, Atom
from ase.io import write

import UFF4MOF_constants
import UFF_constants
import Dreiding_constants

# add more force field classes here as they are made

UFF4MOF_atom_parameters = UFF4MOF_constants.UFF4MOF_atom_parameters
UFF4MOF_bond_orders_0 = UFF4MOF_constants.UFF4MOF_bond_orders_0

UFF_atom_parameters = UFF_constants.UFF_atom_parameters
UFF_bond_orders_0 = UFF_constants.UFF_bond_orders_0

Dreiding_atom_parameters = Dreiding_constants.Dreiding_atom_parameters
Dreiding_bond_orders_0 = Dreiding_constants.Dreiding_bond_orders_0

mass_key = atomic_data.mass_key

def isfloat(value):
	"""
		determines if a value is a float
	"""
	try:
		float(value)
		return True
	except ValueError:
		return False

def lammps_inputs(args):

	cifname, force_field, ff_string, sm_ff_string, outdir, charges, replication, read_pymatgen = args

	# add more forcefields here as they are created
	if ff_string == 'UFF4MOF':
		FF_args = {'FF_parameters':UFF4MOF_atom_parameters, 'bond_orders':UFF4MOF_bond_orders_0}
		cutoff = 12.50
		mixing_rules='shift yes mix geometric'
	elif ff_string == 'UFF':
		FF_args = {'FF_parameters':UFF_atom_parameters, 'bond_orders':UFF_bond_orders_0}
		cutoff = 12.50
		mixing_rules='shift yes mix geometric'
	elif ff_string == 'Dreiding':
		FF_args = {'FF_parameters':Dreiding_atom_parameters, 'bond_orders':Dreiding_bond_orders_0}
		cutoff = 12.50
		mixing_rules='shift yes mix arithmetic'
	elif ff_string == 'MZHB':
		FF_args = {'FF_parameters':None, 'bond_orders':None}
		cutoff = 12.50
		mixing_rules='shift yes mix arithmetic'
	elif ff_string == 'ZIFFF':
		FF_args = {'FF_parameters':UFF4MOF_atom_parameters, 'bond_orders':UFF4MOF_atom_parameters}
		cutoff = 12.50
		mixing_rules='shift yes mix arithmetic'

	system = initialize_system(cifname, charges=charges, read_pymatgen=read_pymatgen)
	system, replication = replication_determination(system, replication, cutoff)
	FF = force_field(system, cutoff, FF_args)
	FF.compile_force_field(charges=charges)

	if sm_ff_string != None:
		add_small_molecules(FF, sm_ff_string)
	else:
		if len(FF.system['SM_graph'].nodes()) != 0:
			warnings.warn('extra-framework molecules detected, but no small molecule force field is specified!')

	first_line = "Created by Ryther's extremely high quality code on " + str(datetime.datetime.now())

	SG = FF.system['graph']
	N_atoms, ty_atoms = (len(SG.nodes()), len(FF.atom_types))
	N_bonds, ty_bonds = FF.bond_data['count']
	N_angles, ty_angles = FF.angle_data['count']

	try:
		N_dihedrals, ty_dihedrals = FF.dihedral_data['count']
	except AttributeError:
		N_dihedrals = 0
		ty_dihedrals = None
	
	try:
		N_impropers, ty_impropers = FF.improper_data['count']
	except AttributeError:
		N_impropers = 0
		ty_impropers = None

	a,b,c,alpha,beta,gamma = system['box']
	lx = np.round(a, 8)
	xy = np.round(b * np.cos(math.radians(gamma)), 8)
	xz = np.round(c * np.cos(math.radians(beta)), 8)
	ly = np.round(np.sqrt(b**2 - xy**2), 8)
	yz = np.round((b*c*np.cos(math.radians(alpha)) - xy*xz)/ly, 8)
	lz = np.round(np.sqrt(c**2 - xz**2 - yz**2), 8)

	suffix = cifname.split('/')[-1].split('.')[0] + '_' + replication
	data_name = 'data.' + suffix

	#with open('check.xyz', 'w') as out:
	#	out.write(str(len(SG.nodes())))
	#	out.write('\n')
	#	out.write('\n')
	#	for atom, atom_data in SG.nodes(data=True):
	#		pos = [np.round(v,8) for v in atom_data['cartesian_position']]
	#		line = [atom_data['element_symbol']] + pos
	#		out.write('{} {} {} {}'.format(*line))
	#		out.write('\n')

	with open(outdir + os.sep + data_name, 'w') as data:
		data.write(first_line + '\n')
		data.write('\n')
		data.write('    ' + str(N_atoms) + ' atoms\n')
		data.write('    ' + str(N_bonds) + ' bonds\n')
		data.write('    ' + str(N_angles) + ' angles\n')
		data.write('    ' + str(N_dihedrals) + ' dihedrals\n')
		data.write('    ' + str(N_impropers) + ' impropers\n')
		data.write('\n')
		data.write('    ' + str(ty_atoms) + ' atom types\n')
		data.write('    ' + str(ty_bonds) + ' bond types\n')
		data.write('    ' + str(ty_angles) + ' angle types\n')
		data.write('    ' + str(ty_dihedrals) + ' dihedral types\n')
		data.write('    ' + str(ty_impropers) + ' improper types\n')
		data.write('\n')
		data.write('0.00000000 ' + str(lx) + ' xlo xhi\n')
		data.write('0.00000000 ' + str(ly) + ' ylo yhi\n')
		data.write('0.00000000 ' + str(lz) + ' zlo zhi\n')
		data.write(str(xy) + ' ' + str(xz) + ' ' + str(yz) + ' xy xz yz \n')
		data.write('\n')
		data.write('Masses \n')
		data.write('\n')
		for fft in FF.atom_masses:
			mass = FF.atom_masses[fft]
			aty = FF.atom_types[fft]
			data.write('{:>5} {:>10}'.format(aty, mass))
			data.write('\n')
		data.write('\n')
		data.write('Pair Coeffs\n')
		data.write('\n')

		for aty in FF.pair_data['params']:
			params = FF.pair_data['params'][aty]
			params = [np.round(x,6) if isfloat(x) else x for x in params]
			comment = FF.pair_data['comments'][aty]
			style = FF.pair_data['style']

			if 'hybrid' in style:
				# type
				data.write('    {:<3}'.format(aty))
				# style needs to be written for hybrid
				data.write('{:<20}'.format(params[0]))
				format_string = ' '.join(['{:{w}.{p}f}' if not np.issubdtype(x, np.integer) else '{:{w}}' for x in params[1:]])
				data.write(format_string.format(*params[1:], w=12, p=5))
				data.write(' '.join(['#'] + comment))
				data.write('\n')
			else:
				# type
				data.write('    {:<3}'.format(aty))
				format_string = ' '.join(['{:{w}.{p}f}' if not np.issubdtype(x, np.integer) else '{:{w}}' for x in params[1:]])
				data.write(format_string.format(*params[1:], w=12, p=5))
				data.write(' '.join([' #'] + comment))
				data.write('\n')

		data.write('\n')
		data.write('Bond Coeffs\n')
		data.write('\n')

		for bty in FF.bond_data['params']:
			params = FF.bond_data['params'][bty]
			params = [np.round(x,6) if isfloat(x) else x for x in params]
			comment = FF.bond_data['comments'][bty]
			style = FF.bond_data['style']

			if 'hybrid' in style:
				# type
				data.write('    {:<3}'.format(bty))
				# style needs to be written for hybrid
				data.write('{:<20}'.format(params[0]))

				try:
					format_string = ' '.join(['{:{w}.{p}f}' if not np.issubdtype(x, np.integer) else '{:{w}}' for x in params[1:]])
				except TypeError:
					format_string = ' '.join(['{:{w}}' for x in params[1:]])
				
				data.write(format_string.format(*params[1:], w=12, p=5))
				data.write(' '.join([' #'] + comment))
				data.write('\n')
			else:
				# type
				data.write('    {:<3}'.format(bty))

				try:
					format_string = ' '.join(['{:{w}.{p}f}' if not np.issubdtype(x, np.integer) else '{:{w}}' for x in params[1:]])
				except TypeError:
					format_string = ' '.join(['{:{w}}' for x in params[1:]])

				data.write(format_string.format(*params[1:], w=12, p=5))
				data.write(' '.join([' #'] + comment))
				data.write('\n')

		data.write('\n')
		data.write('Angle Coeffs\n')
		data.write('\n')

		for aty in FF.angle_data['params']:
			params = FF.angle_data['params'][aty]
			params = [np.round(x,6) if isfloat(x) else x for x in params]
			comment = FF.angle_data['comments'][aty]
			style = FF.angle_data['style']

			if 'hybrid' in style:
				# type
				data.write('    {:<3}'.format(aty))
				# style needs to be written for hybrid
				data.write('{:<20}'.format(params[0]))
				
				try:
					format_string = ' '.join(['{:{w}.{p}f}' if not np.issubdtype(x, np.integer) else '{:{w}}' for x in params[1:]])
				except TypeError:
					format_string = ' '.join(['{:{w}}' for x in params[1:]])

				data.write(format_string.format(*params[1:], w=12, p=5))
				data.write(' '.join([' #'] + comment))
				data.write('\n')
			else:
				# type
				data.write('    {:<3}'.format(aty))

				try:
					format_string = ' '.join(['{:{w}.{p}f}' if not np.issubdtype(x, np.integer) else '{:{w}}' for x in params[1:]])
				except TypeError:
					format_string = ' '.join(['{:{w}}' for x in params[1:]])

				data.write(format_string.format(*params[1:], w=12, p=5))
				data.write(' '.join([' #'] + comment))
				data.write('\n')

		if N_dihedrals != 0:

			data.write('\n')
			data.write('Dihedral Coeffs\n')
			data.write('\n')

			for dty in FF.dihedral_data['params']:
				params = FF.dihedral_data['params'][dty]
				params = [np.round(x,6) if isfloat(x) else x for x in params]
				comment = FF.dihedral_data['comments'][dty]
				style = FF.dihedral_data['style']

				if 'hybrid' in style:
					# type
					data.write('    {:<3}'.format(dty))
					# style needs to be written for hybrid
					data.write('{:<20}'.format(params[0]))
					format_string = ' '.join(['{:{w}.{p}f}' if not np.issubdtype(x, np.integer) else '{:{w}}' for x in params[1:]])
					data.write(format_string.format(*params[1:], w=12, p=5))
					data.write(' '.join([' #'] + comment))
					data.write('\n')
				else:
					# type
					data.write('    {:<3}'.format(dty))
					format_string = ' '.join(['{:{w}.{p}f}' if not np.issubdtype(x, np.integer) else '{:{w}}' for x in params[1:]])
					data.write(format_string.format(*params[1:], w=12, p=5))
					data.write(' '.join([' #'] + comment))
					data.write('\n')

		if N_impropers != 0:

			data.write('\n')
			data.write('Improper Coeffs\n')
			data.write('\n')
	
			for ity in FF.improper_data['params']:
				params = FF.improper_data['params'][ity]
				params = [np.round(x,6) if isfloat(x) else x for x in params]
				comment = FF.improper_data['comments'][ity]
				style = FF.improper_data['style']
	
				if 'hybrid' in style:
					# type
					data.write('    {:<3}'.format(ity))
					# style needs to be written for hybrid
					data.write('{:<20}'.format(params[0]))
					format_string = ' '.join(['{:{w}.{p}f}' if not np.issubdtype(x, np.integer) else '{:{w}}' for x in params[1:]])
					data.write(format_string.format(*params[1:], w=12, p=5))
					data.write(' '.join([' #'] + comment))
					data.write('\n')
				else:
					# type
					data.write('    {:<3}'.format(ity))
					format_string = ' '.join(['{:{w}.{p}f}' if not np.issubdtype(x, np.integer) else '{:{w}}' for x in params[1:]])
					data.write(format_string.format(*params[1:], w=12, p=5))
					data.write(' '.join([' #'] + comment))
					data.write('\n')

		data.write('\n')
		data.write('Atoms\n')
		data.write('\n')
		total_charge = 0.0

		for a in SG.nodes(data=True):
			atom_data = a[1]
			index = a[0]
			force_field_type = atom_data['force_field_type']
			lammps_type = FF.atom_types[force_field_type]
			charge = atom_data['charge']
			total_charge += charge
			pos = [np.round(v,8) for v in atom_data['cartesian_position']]

			data.write('{:>5} {:<5} {:<5} {:8.5f} {:12.5f} {:12.5f} {:12.5f}'.format(index, atom_data['mol_flag'], lammps_type, charge, pos[0], pos[1], pos[2]))
			data.write('\n')

		if charges and abs(total_charge) > 0.001:
			warnings.warn('There is a potentially significant net charge of', total_charge)
			
		data.write('\n')
		data.write('Bonds\n')
		data.write('\n')

		bond_index = 0
		for bond_type in FF.bond_data['all_bonds']:
			for bond in FF.bond_data['all_bonds'][bond_type]:
				bond_index += 1
				data.write('{:>5} {:<5} {:<6} {:<6}'.format(bond_index, bond_type, bond[0], bond[1]))
				data.write('\n')

		data.write('\n')
		data.write('Angles\n')
		data.write('\n')

		angle_index = 0
		for angle_type in FF.angle_data['all_angles']:
			for angle in FF.angle_data['all_angles'][angle_type]:
				angle_index += 1
				data.write('{:>5} {:<5} {:<6} {:<6} {:<6}'.format(angle_index, angle_type, angle[0], angle[1], angle[2]))
				data.write('\n')

		if N_dihedrals != 0:

			data.write('\n')
			data.write('Dihedrals\n')
			data.write('\n')
	
			dihedral_index = 0
			for dihedral_type in FF.dihedral_data['all_dihedrals']:
				for dihedral in FF.dihedral_data['all_dihedrals'][dihedral_type]:
					dihedral_index += 1
					data.write('{:>5} {:<5} {:<6} {:<6} {:<6} {:<6}'.format(dihedral_index, dihedral_type, dihedral[0], dihedral[1], dihedral[2], dihedral[3]))
					data.write('\n')

		if N_impropers != 0:
			
			data.write('\n')
			data.write('Impropers\n')
			data.write('\n')
	
			improper_index = 0
			for improper_type in FF.improper_data['all_impropers']:
				for improper in FF.improper_data['all_impropers'][improper_type]:
					improper_index += 1
					data.write('{:>5} {:<5} {:<6} {:<6} {:<6} {:<6}'.format(improper_index, improper_type, improper[0], improper[1], improper[2], improper[3]))
					data.write('\n')

	with open(outdir + os.sep + 'in.' + suffix, 'w') as infile:
		infile.write('units           real\n')
		infile.write('atom_style      full\n')
		infile.write('boundary        p p p\n')
		infile.write('\n')
		infile.write('pair_style      ' + FF.pair_data['style'] + ' ' + str(FF.cutoff) + '\n')
		infile.write('bond_style      ' + FF.bond_data['style'] + '\n')
		infile.write('angle_style     ' + FF.angle_data['style'] + '\n')
		try:
			infile.write('dihedral_style  ' + FF.dihedral_data['style'] + '\n')
		except AttributeError:
			pass
		try:
			infile.write('improper_style  ' + FF.improper_data['style'] + '\n')
		except AttributeError:
			pass
		if charges:
			infile.write('kspace_style    ewald 0.000001\n')
		infile.write('\n')
		sb = FF.pair_data['special_bonds']
		infile.write('dielectric      1.0\n')
		infile.write('pair_modify     ' + mixing_rules + '\n')
		infile.write('special_bonds   ' + sb + '\n')
		infile.write('box             tilt large\n')
		infile.write('read_data       ' + data_name + '\n')
		infile.write('\n')

