from __future__ import print_function
import sys
import re
import numpy as np
import networkx as nx
import glob

PT = ['H' , 'He', 'Li', 'Be', 'B' , 'C' , 'N' , 'O' , 'F' , 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar',
	  'K' , 'Ca', 'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
	  'Rb', 'Sr', 'Y' , 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I' , 'Xe',
	  'Cs', 'Ba', 'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 
	  'Ra', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th', 
	  'Pa', 'U' , 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'FG', 'X' ]

def nn(string):
	return re.sub('[^a-zA-Z]','', string)

def nl(string):
	return re.sub('[^0-9]','', string)

def isfloat(value):
	"""
		determines if a value is a float
	"""
	try:
		float(value)
		return True
	except ValueError:
		return False

def iscoord(line):
	"""
		identifies coordinates in CIFs
	"""
	if nn(line[0]) in PT and line[1] in PT and False not in map(isfloat,line[2:5]):
		return True
	else:
		return False
	
def isbond(line):
	"""
		identifies bonding in cifs
	"""
	if nn(line[0]) in PT and nn(line[1]) in PT and isfloat(line[2]) and True not in map(isfloat,line[3:]):
		return True
	else:
		return False

def PBC3DF_sym(vec1, vec2):
	"""
		applies periodic boundary to distance between vec1 and vec2 (fractional coordinates)
	"""
	dist = vec1 - vec2
	sym_dist = [(1.0, dim - 1.0) if dim > 0.5 else (-1.0, dim + 1.0) if dim < -0.5 else (0, dim) for dim in dist]
	sym = np.array([s[0] for s in sym_dist])
	ndist = np.array([s[1] for s in sym_dist])

	return ndist, sym

def cif_read(filename, charges=False):

	with open(filename,'r') as f:
		f = f.read()
		f = filter(None, f.split('\n'))

	names = []
	elems = []
	fcoords = []
	charge_list = []
	bonds = []

	for line in f:
		s = line.split()
		if '_cell_length_a' in line:
			a = s[1]
		if '_cell_length_b' in line:
			b = s[1]
		if '_cell_length_c' in line:
			c = s[1]
		if '_cell_angle_alpha' in line:
			alpha = s[1]
		if '_cell_angle_beta' in line:
			beta = s[1]
		if '_cell_angle_gamma' in line:
			gamma = s[1]
		if iscoord(s):
			names.append(s[0])
			elems.append(s[1])
			fvec = np.array([np.round(float(v),8) for v in s[2:5]])
			fcoords.append(fvec)
			if charges:
				charge_list.append(float(s[-1]))
			else:
				charge_list.append(0.0)
		if isbond(s):
			bonds.append((s[0],s[1],s[3],s[4]))

	pi = np.pi
	a,b,c,alpha,beta,gamma = map(float,(a,b,c,alpha,beta,gamma))
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
	
	norm_vec = fcoords[1]
	ccoords = []

	for l in fcoords:
		vec = l
		dist,sym = PBC3DF_sym(norm_vec,vec)
		vec = vec - sym
		vec = np.dot(unit_cell,vec)
		ccoords.append(vec)

	fcoords = np.asarray(fcoords)
	ccoords = np.asarray(ccoords)
	charges = np.asarray(charges)

	return elems, names, ccoords, fcoords, charge_list, bonds, (a,b,c,alpha,beta,gamma)

def initialize_system(filename, charges=False):

	elems, names, ccoords, fcoords, charge_list, bonds, uc_params = cif_read(filename, charges=charges)
	A,B,C,alpha,beta,gamma = uc_params

	G = nx.Graph()
	index = 0
	index_key = {}
	for e, n, cc, fc, charge in zip(elems, names, ccoords, fcoords, charge_list):
		index += 1
		G.add_node(index, element_symbol=e, index=index, force_field_type='', cartesian_position=cc, fractional_position=fc, charge=charge)
		index_key[n] = index

	for b in bonds:
		G.add_edge(index_key[b[0]], index_key[b[1]], sym_code=b[2], bond_type=b[3])

	return {'box':(A,B,C,alpha,beta,gamma), 'graph':G}

