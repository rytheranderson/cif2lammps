import networkx as nx
import numpy as np
import atomic_data
from itertools import groupby, combinations
import small_molecule_constants
from cif2system import PBC3DF_sym

mass_key = atomic_data.mass_key

def add_small_molecules(FF, ff_string):
	
	if ff_string == 'TraPPE':
		SM_constants = small_molecule_constants.TraPPE
	if ff_string == 'TIP4P':
		SM_constants = small_molecule_constants.TIP4P
	# insert more force fields here if needed
	else:
		raise ValueError('the small molecule force field', ff_string, 'is not defined')

	SG = FF.system['graph']
	SMG = FF.system['SM_graph']
	mol_flag = 1
	max_ind = FF.system['max_ind']
	index = max_ind

	box = FF.system['box']
	a,b,c,alpha,beta,gamma = box
	pi = np.pi
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
	inv_unit_cell = np.linalg.inv(unit_cell)

	add_nodes = []
	add_edges = []
	comps = []

	for comp in nx.connected_components(SMG):

		mol_flag += 1
		comp = sorted(list(comp))
		ID_string = sorted([SMG.node[n]['element_symbol'] for n in comp])
		ID_string = [(key, len(list(group))) for key, group in groupby(ID_string)]
		ID_string = ''.join([str(e) for c in ID_string for e in c])
		comps.append(ID_string)

		for n in comp:

			data = SMG.node[n]

			SMG.node[n]['mol_flag'] = str(mol_flag)

			if ID_string == 'H2O1':
				SMG.node[n]['force_field_type'] = SMG.node[n]['element_symbol'] + '_w' 
			else:
				SMG.node[n]['force_field_type'] = SMG.node[n]['element_symbol'] + '_' + ID_string

		# add COM sites where relevant, extend this list as new types are added
		if ID_string in ('O2', 'N2'):

			coords = []
			anchor = SMG.node[comp[0]]['fractional_position']

			for n in comp:

				data = SMG.node[n]
				data['mol_flag'] = str(mol_flag)
				fcoord = data['fractional_position']
				mic = PBC3DF_sym(fcoord, anchor)
				fcoord += mic[1]
				ccoord = np.dot(unit_cell, fcoord)
				coords.append(ccoord)

			ccom = np.average(coords, axis=0)
			fcom = np.dot(inv_unit_cell, ccom)
			index += 1

			if ID_string == 'O2':
				fft = 'O_com'
			elif ID_string == 'N2':
				fft = 'N_com'

			ndata =  {'element_symbol':'NA', 'mol_flag':mol_flag, 'index':index, 'force_field_type':fft, 'cartesian_position':ccom, 'fractional_position':fcom, 'charge':0.0, 'replication':np.array([0.0,0.0,0.0]), 'duplicated_version_of':None}
			edata =  {'sym_code':None, 'length':None, 'bond_type':None}

			add_nodes.append([index, ndata])
			add_edges.extend([(index, comp[0], edata), (index, comp[1], edata)])

	for n, data in add_nodes:
		SMG.add_node(n, **data)
	for e0, e1, data in add_edges:
		SMG.add_edge(e0, e1, **data)

	ntypes = max([FF.atom_types[ty] for ty in FF.atom_types])
	nbonds = max([i for i in FF.bond_data['params']])
	nangles = max([i for i in FF.angle_data['params']])
	
	try:
		ndihedrals = max([i for i in FF.dihedral_data['params']])
	except ValueError:
		ndihedrals = 0
	try:
		nimpropers = max([i for i in FF.improper_data['params']])
	except ValueError:
		nimpropers = 0

	new_bond_types = {}
	new_angle_types = {}
	new_dihedral_types = {}
	new_improper_types = {}

	for subG, ID_string in zip([SMG.subgraph(c).copy() for c in nx.connected_components(SMG)], comps):

		constants = SM_constants[ID_string]

		# add new atom types
		for name,data in sorted(subG.nodes(data=True), key=lambda x:x[0]):

			fft = data['force_field_type']
			chg = constants['pair']['charges'][fft]
			data['charge'] = chg
			SG.add_node(name, **data)

			try: 

				FF.atom_types[fft] += 0

			except KeyError:

				ntypes += 1
				FF.atom_types[fft] = ntypes
				style = constants['pair']['style']
				vdW = constants['pair']['vdW'][fft]
				FF.pair_data['params'][FF.atom_types[fft]] = (style, vdW[0], vdW[1])
				FF.pair_data['comments'][FF.atom_types[fft]] = [fft, fft]
				FF.atom_masses[fft] = mass_key[data['element_symbol']]

				if 'hybrid' not in FF.pair_data['style'] and style != FF.pair_data['style']:
					FF.pair_data['style'] = ' '.join(['hybrid', FF.pair_data['style'], style])
				elif 'hybrid' in FF.pair_data['style'] and style in FF.pair_data['style']:
					pass
				elif 'hybrid' in FF.pair_data['style'] and style not in FF.pair_data['style']:
					FF.pair_data['style'] += ' ' + style

		# add new bonds
		used_bonds = []
		ty = nbonds
		for e0,e1,data in subG.edges(data=True):

			bonds = constants['bonds']
			fft_i = SG.node[e0]['force_field_type']
			fft_j = SG.node[e1]['force_field_type']
			# make sure the order corresponds to that in the molecule dictionary
			bond = tuple(sorted([fft_i, fft_j]))

			try:

				style = bonds[bond][0]

				if bond not in used_bonds:

					ty = ty + 1
					new_bond_types[bond] = ty
					FF.bond_data['params'][ty] = list(bonds[bond])
					FF.bond_data['comments'][ty] = list(bond)

					used_bonds.append(bond)

				if 'hybrid' not in FF.bond_data['style'] and style != FF.bond_data['style']:
					FF.bond_data['style'] = ' '.join(['hybrid', FF.bond_data['style'], style])
				elif 'hybrid' in FF.bond_data['style'] and style in FF.bond_data['style']:
					pass
				elif 'hybrid' in FF.bond_data['style'] and style not in FF.bond_data['style']:
					FF.bond_data['style'] += ' ' + style

				if ty in FF.bond_data['all_bonds']:
					FF.bond_data['count'] = (FF.bond_data['count'][0] + 1, FF.bond_data['count'][1] + 1)
					FF.bond_data['all_bonds'][ty].append((e0,e1))
				else:
					FF.bond_data['count'] = (FF.bond_data['count'][0] + 1, FF.bond_data['count'][1] + 1)
					FF.bond_data['all_bonds'][ty] = [(e0,e1)]

			except KeyError:
				pass
				
		# add new angles
		used_angles = []
		ty = nangles
		for name,data in subG.nodes(data=True):

			angles = constants['angles']
			nbors = list(subG.neighbors(name))

			for comb in combinations(nbors, 2):

				j = name
				i, k = comb
				fft_i = subG.node[i]['force_field_type']
				fft_j = subG.node[j]['force_field_type']
				fft_k = subG.node[k]['force_field_type']

				angle = sorted((fft_i, fft_k))
				angle = (angle[0], fft_j, angle[1])

				try:
	
					style = angles[angle][0]
					FF.angle_data['count'] = (FF.angle_data['count'][0] + 1, FF.angle_data['count'][1])
	
					if angle not in used_angles:
						
						ty = ty + 1
						new_angle_types[angle] = ty
						FF.angle_data['count'] = (FF.angle_data['count'][0], FF.angle_data['count'][1] + 1)
						FF.angle_data['params'][ty] = list(angles[angle])
						FF.angle_data['comments'][ty] = list(angle)

						used_angles.append(angle)
	
					if 'hybrid' not in FF.angle_data['style'] and style != FF.angle_data['style']:
						FF.angle_data['style'] = ' '.join(['hybrid', FF.angle_data['style'], style])
					elif 'hybrid' in FF.angle_data['style'] and style in FF.angle_data['style']:
						pass
					elif 'hybrid' in FF.angle_data['style'] and style not in FF.angle_data['style']:
						FF.angle_data['style'] += ' ' + style
					
					if ty in FF.angle_data['all_angles']:
						FF.angle_data['all_angles'][ty].append((i,j,k))
					else:
						FF.angle_data['all_angles'][ty] = [(i,j,k)]

				except KeyError:
					pass

		# add new dihedrals

	FF.bond_data['count'] = (FF.bond_data['count'][0], len(FF.bond_data['params']))
	FF.angle_data['count'] = (FF.angle_data['count'][0], len(FF.angle_data['params']))




