from __future__ import print_function
from cif2system import initialize_system, duplicate_system
import atomic_data
import sys
import os
import numpy as np
import math
import datetime
import math
import itertools

from force_field_construction import UFF
# add more force field classes here as they are made

def isfloat(value):
	"""
		determines if a value is a float
	"""
	try:
		float(value)
		return True
	except ValueError:
		return False

mass_key = atomic_data.mass_key

def lammps_inputs(args):

	cifname, force_field, outdir, charges, replication = args

	system = initialize_system(cifname, charges=charges)
	
	if 'min_atoms' in replication:
		
		min_atoms = int(replication.split(':')[-1])
		box = system['box']

		pi = np.pi
		a,b,c,alpha,beta,gamma = box
		ax = a
		ay = 0.0
		az = 0.0
		bx = b * np.cos(gamma * pi / 180.0)
		by = b * np.sin(gamma * pi / 180.0)
		bz = 0.0
		cx = c * np.cos(beta * pi / 180.0)
		cy = (c * b * np.cos(alpha * pi /180.0) - bx * cx) / by
		cz = (c ** 2.0 - cx ** 2.0 - cy ** 2.0) ** 0.5
		unit_cell = np.asarray([[ax,ay,az],[bx,by,bz],[cx,cy,cz]]).T
		inv_uc = np.linalg.inv(unit_cell)

		G = system['graph']
		Natoms = float(len(G.nodes()))
		
		duplications = int(math.ceil(min_atoms/Natoms))
		rvals = range(duplications + 1)[1:]
		shapes = itertools.product(rvals, rvals, rvals)
		shapes = [s for s in shapes if reduce((lambda x, y: x * y), s) == duplications]
		shape_deviations = [(i, np.std([shapes[i][0]*a, shapes[i][1]*b, shapes[i][2]*c])) for i in range(len(shapes))]
		shape_deviations.sort(key = lambda x:x[1])
		selected_shape = shapes[shape_deviations[0][0]]
		
		replication = 'x'.join(map(str, selected_shape))
		system = duplicate_system(system, replication)
		replication='ma' + str(min_atoms)

	elif 'min_length' in replication:
		
		min_length = float(replication.split(':')[-1])
		replication=''

	elif 'x' in replication and replication != '1x1x1':
		
		system = duplicate_system(system, replication)

	FF = force_field(system)
	FF.compile_force_field(charges=charges)

	first_line = "Created by Ryther's extremely high quality code on " + str(datetime.datetime.now())

	SG = FF.system['graph']
	N_atoms, ty_atoms = (len(SG.nodes()), len(FF.atom_types))
	N_bonds, ty_bonds = FF.bond_data['count']
	N_angles, ty_angles = FF.angle_data['count']
	N_dihedrals, ty_dihedrals = FF.dihedral_data['count']
	N_impropers, ty_impropers = FF.improper_data['count']

	a,b,c,alpha,beta,gamma = system['box']
	lx = np.round(a, 8)
	xy = np.round(b * np.cos(math.radians(gamma)), 8)
	xz = np.round(c * np.cos(math.radians(beta)), 8)
	ly = np.round(np.sqrt(b**2 - xy**2), 8)
	yz = np.round((b * c*np.cos(math.radians(alpha)) - xy*xz)/ly, 8)
	lz = np.round(np.sqrt(c**2 - xz**2 - yz**2), 8)

	suffix = cifname.split('/')[-1].split('.')[0] + '_' + replication
	data_name = 'data.' + suffix

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
				format_string = ' '.join(['{:{w}.{p}f}' if not np.issubdtype(x, np.integer) else '{:{w}}' for x in params[1:]])
				data.write(format_string.format(*params[1:], w=12, p=5))
				data.write(' '.join([' #'] + comment))
				data.write('\n')
			else:
				# type
				data.write('    {:<3}'.format(bty))
				format_string = ' '.join(['{:{w}.{p}f}' if not np.issubdtype(x, np.integer) else '{:{w}}' for x in params[1:]])
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
				format_string = ' '.join(['{:{w}.{p}f}' if not np.issubdtype(x, np.integer) else '{:{w}}' for x in params[1:]])
				data.write(format_string.format(*params[1:], w=12, p=5))
				data.write(' '.join([' #'] + comment))
				data.write('\n')
			else:
				# type
				data.write('    {:<3}'.format(aty))
				format_string = ' '.join(['{:{w}.{p}f}' if not np.issubdtype(x, np.integer) else '{:{w}}' for x in params[1:]])
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

		for a in SG.nodes(data=True):
			atom_data = a[1]
			index = atom_data['index']
			force_field_type = atom_data['force_field_type']
			lammps_type = FF.atom_types[force_field_type]
			if charges:
				charge = atom_data['charge']
			else:
				charge = 0.0
			pos = [np.round(v,8) for v in atom_data['cartesian_position']]

			data.write('{:>5} {:<5} {:<5} {:8.5f} {:12.5f} {:12.5f} {:12.5f}'.format(index, '444', lammps_type, charge, pos[0], pos[1], pos[2]))
			data.write('\n')

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
		infile.write('pair_style      ' + FF.pair_data['style'] + ' 12.5\n')
		infile.write('bond_style      ' + FF.bond_data['style'] + '\n')
		infile.write('angle_style     ' + FF.angle_data['style'] + '\n')
		infile.write('dihedral_style  ' + FF.dihedral_data['style'] + '\n')
		infile.write('improper_style  ' + FF.improper_data['style'] + '\n')
		if charges:
			infile.write('kspace_style    ewald 0.000001\n')
		infile.write('\n')
		sb = FF.pair_data['special_bonds']
		infile.write('dielectric      1.0\n')
		infile.write('pair_modify     shift yes mix geometric\n')
		infile.write('special_bonds   ' + sb + '\n')
		infile.write('box             tilt large\n')
		infile.write('read_data       ' + data_name + '\n')
		infile.write('\n')

