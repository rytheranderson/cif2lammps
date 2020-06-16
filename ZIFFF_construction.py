from __future__ import print_function
import itertools
import atomic_data
from force_field_construction import force_field

metals = atomic_data.metals
mass_key = atomic_data.mass_key

class ZIFFF(force_field):

	def __init__(self, system, cutoff, args):

		self.system = system
		self.cutoff = cutoff
		self.args = args

	def type_atoms(self):

		SG = self.system['graph']
		types = []
		imidazole_ring_atoms = []

		for atom,data in SG.nodes(data=True):

			element_symbol = data['element_symbol']

			if element_symbol in ('Zn','Cd'):
				nborhood = nx.ego_graph(G,atom,radius=2)
				print(nborhood)

		# assign Zn, N, and C types
		for atom in SG.nodes(data=True):
			
			name, inf = atom
			element_symbol = inf['element_symbol']
			nbors = [a for a in SG.neighbors(name)]
			nbor_symbols = [SG.nodes[n]['element_symbol'] for n in nbors]
			mass = mass_key[element_symbol]

			# one type of Zn
			if element_symbol == 'Zn':
				ty = 'Zn'
			# one type of N
			if element_symbol == 'N':
				ty = 'N'
			# three types of C
			if element_symbol == 'C':
				if len(nbors) == 3 and 'H' not in nbor_symbols:
					ty = 'C1'
				elif len(nbors) == 3 and 'H' in nbor_symbols:
					ty = 'C2'
				elif len(nbors) == 4 and 'N' not in nbor_symbols:
					ty = 'C3'

			hyb = None

			types.append((ty, element_symbol, mass))
			SG.node[name]['force_field_type'] = ty
			SG.node[name]['hybridization'] = hyb

		# assign H types once all carbons have been typed
		for atom in SG.nodes(data=True):
			
			name, inf = atom
			element_symbol = inf['element_symbol']
			nbors = [a for a in SG.neighbors(name)]
			nbor_types = [SG.nodes[n]['force_field_type'] for n in nbors]
			mass = mass_key[element_symbol]

			# two types of H
			if element_symbol == 'H':
				if 'C2' in nbor_types:
					ty = 'H2'
				elif 'C3' in nbor_types:
					ty = 'H3'

			hyb = None

			types.append((ty, element_symbol, mass))
			SG.node[name]['force_field_type'] = ty
			SG.node[name]['hybridization'] = hyb

		types = set(types)
		Ntypes = len(types)
		atom_types = dict((ty[0],i+1) for i,ty in zip(range(Ntypes), types))
		atom_element_symbols = dict((ty[0], ty[1]) for ty in types)
		atom_masses = dict((ty[0],ty[2]) for ty in types)

		self.system['graph'] = SG
		self.atom_types = atom_types
		self.atom_element_symbols = atom_element_symbols
		self.atom_masses = atom_masses

	def bond_parameters(self, bond, bond_order):
		
		return ('style', *constants)

	def angle_parameters(self, angle, r_ij, r_jk):

		constants = []
		return ('style', *constants)

	def dihedral_parameters(self, bond, hybridization, element_symbols, nodes):

		constants = []
		return ('style', *constants)

	def improper_parameters(self, fft_i, O_2_flag):
		
		constants = []
		return ('style', *constants)

	def pair_parameters(self, charges=False):

		self.pair_data = {'params':params, 'style':style, 'special_bonds':sb, 'comments':comments}

	def enumerate_bonds(self):

		self.bond_data = {'all_bonds':all_bonds, 'params':bond_params, 'style':'harmonic', 'count':(count, len(all_bonds)), 'comments':bond_comments}

	def enumerate_angles(self):

		self.angle_data = {'all_angles':all_angles, 'params':angle_params, 'style':style, 'count':(count, len(all_angles)), 'comments':angle_comments}

	def enumerate_dihedrals(self):
		
		self.dihedral_data = {'all_dihedrals':all_dihedrals, 'params':indexed_dihedral_params, 'style':'harmonic', 'count':(count, len(all_dihedrals)), 'comments':dihedral_comments}

	def enumerate_impropers(self):
		
		self.improper_data = {'all_impropers':all_impropers, 'params':improper_params, 'style':'fourier', 'count':(count, len(all_impropers)), 'comments':improper_comments}

	def compile_force_field(self, charges=False):

		self.type_atoms()
		self.pair_parameters(charges)
		self.enumerate_bonds()
		self.enumerate_angles()
		self.enumerate_dihedrals()
		self.enumerate_impropers()
		

