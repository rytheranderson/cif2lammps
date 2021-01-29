import atomic_data
from textwrap import dedent

def water(last_atom_ID, last_bond_ID, last_angle_ID, model='TIP4P_cutoff'):

	ID_O = last_atom_ID + 1
	ID_H = ID_O + 1

	BT = last_bond_ID + 1
	AT = last_angle_ID + 1

	charge_dict = {
	'TIP4P_cutoff': (-1.0400, 0.5200),
	'TIP4P_2005':   (-1.1128, 0.5564),
	'TIP4P_long':   (-1.0484, 0.5242),
	'TIP3P_long':   (-0.8300, 0.4150)
	}

	LJ_dict = {
	# LAMMPS has a special TIP4P pair_style that automatically adds the M site
	'TIP4P_cutoff': {ID_O: ('lj/cut/tip4p/cut' , 0.15500, 3.15360), ID_H: ('lj/cut/tip4p/cut' , 0.0, 0.0), 'style': 'lj/cut/tip4p/cut'},
	'TIP4P_2005':   {ID_O: ('lj/cut/tip4p/long', 0.15500, 3.15360), ID_H: ('lj/cut/tip4p/long', 0.0, 0.0), 'style': 'lj/cut/tip4p/long'},
	'TIP4P_long':   {ID_O: ('lj/cut/tip4p/long', 0.15500, 3.15360), ID_H: ('lj/cut/tip4p/long', 0.0, 0.0), 'style': 'lj/cut/tip4p/long'},
	'TIP3P_long':   {ID_O: ('lj/cut/coul/long' , 0.15500, 3.15360), ID_H: ('lj/cut/coul/long' , 0.0, 0.0), 'style': 'lj/cut/coul/long'}
	}

	bond_dict = {
	# TIP4P is a rigid model (use fix shake), force constants should just be reasonable values
	# TIP3P has force constants if a flexible model is desired
	'TIP4P_cutoff': {BT: {'style':'harmonic', 'params':(100.0, 0.9572), 'comments':'# O_water H_water'}},
	'TIP4P_2005':   {BT: {'style':'harmonic', 'params':(100.0, 0.9572), 'comments':'# O_water H_water'}},
	'TIP4P_long':   {BT: {'style':'harmonic', 'params':(100.0, 0.9572), 'comments':'# O_water H_water'}},
	'TIP3P_long':   {BT: {'style':'harmonic', 'params':(450.0, 0.9572), 'comments':'# O_water H_water'}} 
	}

	angle_dict = {
	# TIP4P is a rigid model (use fix shake), force constants should just be reasonable values
	# TIP3P has force constants if a flexible model is desired
	'TIP4P_cutoff': {AT: {'style':'harmonic', 'params':(50.0, 104.52), 'comments':'# O_water H_water'}},
	'TIP4P_2005':   {AT: {'style':'harmonic', 'params':(50.0, 104.52), 'comments':'# O_water H_water'}},
	'TIP4P_long':   {AT: {'style':'harmonic', 'params':(50.0, 104.52), 'comments':'# O_water H_water'}},
	'TIP3P_long':   {AT: {'style':'harmonic', 'params':(55.0, 104.52), 'comments':'# O_water H_water'}}
	}

	qO,qH = charge_dict[model]
	LJ_params = LJ_dict[model]
	bond_params = bond_dict[model]
	angle_params = angle_dict[model]

	molfile = dedent("""# Water molecule. useable for TIP3P or TIP4P in LAMMPS.

	                3 atoms
	                2 bonds
	                1 angle

	                Coords

	                1    1.12456   0.09298   1.27452
	                2    1.53683   0.75606   1.89928
	                3    0.49482   0.56390   0.65678

	                Types

	                1    {ID_O}
	                2    {ID_H}
	                3    {ID_H}

	                Charges

	                1    {qO}
	                2    {qH}
	                3    {qH}

	                Bonds

	                1    {BT} 1 2
	                2    {BT} 1 3

	                Angles

	                1    {AT} 2 1 3""".format(**locals())).strip()

	return molfile, LJ_params, bond_params, angle_params
