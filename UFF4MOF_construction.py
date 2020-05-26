from __future__ import print_function
import numpy as np
import math
import itertools
import atomic_data
from force_field_construction import force_field
from cif2system import PBC3DF_sym

metals = atomic_data.metals
mass_key = atomic_data.mass_key

class UFF4MOF(force_field):

	def __init__(self, system, cutoff, args):

		self.system = system
		self.cutoff = cutoff
		self.args = args

		pi = np.pi
		a,b,c,alpha,beta,gamma = system['box']
		ax = a
		ay = 0.0
		az = 0.0
		bx = b * np.cos(gamma * pi / 180.0)
		by = b * np.sin(gamma * pi / 180.0)
		bz = 0.0
		cx = c * np.cos(beta * pi / 180.0)
		cy = (c * b * np.cos(alpha * pi /180.0) - bx * cx) / by
		cz = (c ** 2.0 - cx ** 2.0 - cy ** 2.0) ** 0.5
		self.unit_cell = np.asarray([[ax,ay,az],[bx,by,bz],[cx,cy,cz]]).T

	def type_atoms(self):

		UFF4MOF_atom_parameters = self.args['FF_parameters']
		SG = self.system['graph']
		types = []

		for atom in SG.nodes(data=True):
			
			name, inf = atom
			element_symbol = inf['element_symbol']
			nbors = [a for a in SG.neighbors(name)]
			nbor_symbols = [SG.nodes[n]['element_symbol'] for n in nbors]
			bond_types = [SG.get_edge_data(name, n)['bond_type'] for n in nbors]
			mass = mass_key[element_symbol]

			if len(nbors) > 1:

				dist_j, sym_j = PBC3DF_sym(SG.node[nbors[0]]['fractional_position'], inf['fractional_position'])
				dist_k, sym_k = PBC3DF_sym(SG.node[nbors[1]]['fractional_position'], inf['fractional_position'])
	
				dist_j = np.dot(self.unit_cell, dist_j)
				dist_k = np.dot(self.unit_cell, dist_k)
	
				cosine_angle = np.dot(dist_j, dist_k) / (np.linalg.norm(dist_j) * np.linalg.norm(dist_k))
				angle = (180.0/np.pi) * np.arccos(cosine_angle)

			# Atom typing for UFF4MOF, this can be made much more robust with pattern matching,
			# but this works for most ToBaCCo MOFs, use at your own risk.
			ty = None
			if 'A' in bond_types and element_symbol != 'O':
				ty = element_symbol + '_' + 'R'
				hyb = 'resonant'
			else:
				# Group 1
				if element_symbol == 'H':
					ty = element_symbol + '_'
					hyb = 'sp1'
				# Group 6
				elif element_symbol in ('C', 'Si'):
					if len(element_symbol) == 1:
						ty = element_symbol + '_' + str(len(nbors) - 1)
						hyb = 'sp' + str(len(nbors) - 1)
					else:
						ty = element_symbol + str(len(nbors) - 1)
						hyb = 'sp' + str(len(nbors) - 1)
				# Group 7
				elif element_symbol in ('N'):
					ty = element_symbol + '_' + str(len(nbors))
					hyb = 'sp' + str(len(nbors))
				# Group 8
				elif element_symbol in ('O', 'S'):
					
					### oxygens ###
					
					if element_symbol == 'O':
						# =O for example
						if len(nbors) == 1:
							ty = 'O_1'
							hyb = 'sp1'
						# -OH, for example
						elif len(nbors) == 2 and 'A' not in bond_types and 'D' not in bond_types and not any(i in metals for i in nbor_symbols):
							ty = 'O_3'
							hyb = 'sp3'
						# coordinated solvent, same parameters as O_3, but different name to modulate bond orders
						elif len(nbors) == 2 and len([i for i in nbor_symbols if i in metals]) == 1 and 'H' in nbor_symbols:
							ty = 'O_3_M'
							hyb = 'sp3'
						# furan oxygen, for example
						elif len(nbors) == 2 and 'A' in bond_types and not any(i in metals for i in nbor_symbols):
							ty = 'O_R'
							hyb = 'sp2'
						# carboxyllic oxygen
						elif len(nbors) == 2 and 'D' in bond_types and not any(i in metals for i in nbor_symbols):
							ty = 'O_2'
							hyb = 'sp2'
						# carboxylate oxygen bound to metal node, same parameters as O_2, but different name to modulate bond orders
						elif len(nbors) == 2 and any(i in metals for i in nbor_symbols) and 'C' in nbor_symbols:
							ty = 'O_2_M'
							hyb = 'sp2'
						# 3 connected oxygens
						elif len(nbors) == 3 and any(i in metals for i in nbor_symbols):

							dist_triangle = abs(angle-120.0)
							dist_tetrahedral = abs(angle-109.47)

							if dist_triangle < dist_tetrahedral:
								ty = 'O_2_z'
								hyb = 'sp2'

							else:
								ty = 'O_3_f'
								hyb = 'sp3'

						# 4 connected oxygens
						elif len(nbors) == 4:

							dist_square_linear = min((abs(angle-90.0), abs(angle-180.0)))
							dist_tetrahedral = abs(angle-109.47)

							# update this based on tetrahedral vs. square planar geometry
							if dist_square_linear < dist_tetrahedral:
								ty = 'O_4_f'
								hyb = 'sp3'

							else:
								ty = 'O_3_f'
								hyb = 'sp3'

						# error if no type is identified
						else:
							raise ValueError('Oxygen with neighbors ' + ' '.join(nbor_symbols) + ' is not parametrized')
					# sulfur case is simple
					elif element_symbol == 'S':
						ty = 'S_' + str(len(nbors) + 1)
						hyb = 'sp' + str(len(nbors) + 1)
				# Group 9
				elif element_symbol in ('F', 'Cl', 'Br') and len(nbor_symbols) == 1:
					ty = element_symbol + '_'
					hyb = 'sp1'
				elif element_symbol == 'Cl' and len(nbor_symbols) == 4:
					ty = 'Cl_f'
					hyb = 'sp1'
				
				### metals ###
				
				elif element_symbol in metals:

					hyb = 'NA'

					dist_square_linear = min((abs(angle-90.0), abs(angle-180.0)))
					dist_tetrahedral = abs(angle-109.47)

					if len(element_symbol) == 1:
						add_symbol = element_symbol + '_'
					else:
						add_symbol = element_symbol

					# 2 connected, linear
					if len(nbors) == 2 and abs(angle - 180.0) < 10.0:
						add_symbol + '1f1'

					# 4 connected, square planar or tetrahedral
					elif len(nbors) == 4:

						if dist_square_linear < dist_tetrahedral:
							try:
								UFF4MOF_atom_parameters[add_symbol + '4f2']
								ty = add_symbol + '4f2'
							except KeyError:
								ty = add_symbol + '4+2'
						else:
							try:
								UFF4MOF_atom_parameters[add_symbol + '3f2']
								ty = add_symbol + '3f2'
							except KeyError:
								ty = add_symbol + '3+2'

					# paddlewheels
					elif len(nbors) == 5 and any(i in metals for i in nbor_symbols) and (dist_square_linear < dist_tetrahedral):
						try:
							UFF4MOF_atom_parameters[add_symbol + '4f2']
							ty = add_symbol + '4f2'
						except KeyError:
							ty = add_symbol + '4+2'

					# M3O(CO2H)6 metals, e.g. MIL-100
					elif len(nbors) in (5,6) and not any(i in metals for i in nbor_symbols) and (dist_square_linear < dist_tetrahedral):
						try:
							UFF4MOF_atom_parameters[add_symbol + '6f3']
							ty = add_symbol + '6f3'
						except KeyError:
							try:
								UFF4MOF_atom_parameters[add_symbol + '6+3']
								ty = add_symbol + '6+3'
							except KeyError:
								UFF4MOF_atom_parameters[add_symbol + '6+2']
								ty = add_symbol + '6+2'

					# 7 and 8c metals
					elif len(nbors) in (7,8) and element_symbol in ('Cd', 'Eu', 'Tb', 'Zr'):
						ty = element_symbol + '8f4'

					else:
						raise ValueError('No UFF4MOF type identified for ' + element_symbol + ' with neighbors ' + ' '.join(nbor_symbols))

				# if no type can be identified
				else:
					raise ValueError('No UFF4MOF type identified for ' + element_symbol + ' with neighbors ' + ' '.join(nbor_symbols))
			
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
		
		UFF4MOF_atom_parameters = self.args['FF_parameters']

		i,j = bond
		params_i = UFF4MOF_atom_parameters[i]
		params_j = UFF4MOF_atom_parameters[j]

		r0_i, theta0_i, x1_i, D1_i, zeta_i, Z1_i, V_i, X_i = params_i
		r0_j, theta0_j, x1_j, D1_j, zeta_j, Z1_j, V_j, X_j = params_j

		# bond-order correction
		rbo = -0.1332 * (r0_i+r0_j) * np.log(bond_order)
		# electronegativity correction
		ren = r0_i*r0_j * (((np.sqrt(X_i) - np.sqrt(X_j))**2)) / (X_i*r0_i + X_j*r0_j)
		# equilibrium distance
		r_ij = r0_i + r0_j + rbo - ren
		r_ij3 = r_ij * r_ij * r_ij
		# force constant (1/2 factor should be included here for LAMMPS)
		k_ij = 0.5 * 664.12 * ((Z1_i*Z1_j)/r_ij3)

		return ('harmonic', k_ij, r_ij)

	def angle_parameters(self, angle, r_ij, r_jk):
		
		UFF4MOF_atom_parameters = self.args['FF_parameters']

		i,j,k = angle
		angle_style = 'cosine/periodic'

		params_i = UFF4MOF_atom_parameters[i]
		params_j = UFF4MOF_atom_parameters[j]
		params_k = UFF4MOF_atom_parameters[k]

		r0_i, theta0_i, x1_i, D1_i, zeta_i, Z1_i, V_i, X_i = params_i
		r0_j, theta0_j, x1_j, D1_j, zeta_j, Z1_j, V_j, X_j = params_j
		r0_k, theta0_k, x1_k, D1_k, zeta_k, Z1_k, V_k, X_k = params_k

		# linear
		if theta0_j == 180.0:
			n = 1
			b = 1
		# trigonal planar
		elif theta0_j == 120.0:
			n = 3
			b = -1
		# square planar or octahedral
		elif theta0_j == 90.0:
			n = 4
			b = 1
		# general non-linear
		else:
			n = 'NA'
			b = 'NA'

		cosT0 = np.cos(math.radians(theta0_j))
		sinT0 = np.sin(math.radians(theta0_j))

		r_ik = np.sqrt(r_ij**2.0 + r_jk**2.0 - 2.0*r_ij*r_jk*cosT0)
		# force constant
		K = ((664.12*Z1_i*Z1_k)/(r_ik**5.0)) * (3.0*r_ij*r_jk*(1.0-cosT0**2.0)-r_ik**2.0*cosT0)

		# general non-linear
		if theta0_j not in (90.0, 120.0, 180.0):

			angle_style = 'fourier'
			C2 = 1.0/(4*sinT0**2) 
			C1 = -4*C2*cosT0
			C0 = C2*(2*cosT0**2+1)
			
			return (angle_style, K, C0, C1, C2)

		# this is needed to correct the LAMMPS angle energy calculation
		K *= 0.5

		return (angle_style, K, b, n)

	def dihedral_parameters(self, bond, hybridization, element_symbols, nodes):

		fft_j, fft_k, bond_order = bond
		hyb_j, hyb_k = hybridization
		els_j, els_k = element_symbols
		node_j, node_k = nodes 

		SG = self.system['graph']
		UFF4MOF_atom_parameters = self.args['FF_parameters']

		con_j = SG.degree(node_j) - 1
		con_k = SG.degree(node_k) - 1

		mult = con_j * con_k
		if mult == 0.0:
			return 'NA'

		# cases taken from the DREIDING paper (same cases, different force constants for UFF)
		# they are not done in order to save some lines, I don't know of a better way for doing
		# this besides a bunch of conditionals.
		if hyb_j == 'sp3' and hyb_k == 'sp3':
			# case (a)
			phi0 = 60.0
			n = 3.0
			V_j = UFF4MOF_atom_parameters[fft_j][6]
			V_k = UFF4MOF_atom_parameters[fft_k][6]
			V = np.sqrt(V_j*V_k)
			# case (h)
			if els_j == 'O' and els_k == 'O':
				phi0 = 90.0
				n = 2.0
				V = 2.0
			elif els_j == 'S' and els_k == 'S':
				phi0 = 90.0
				n = 2.0
				V = 6.8

		elif (hyb_j in ('sp2', 'resonant') and hyb_k == 'sp3') or (hyb_k in ('sp2', 'resonant') and hyb_j == 'sp3'):
			# case (b)
			phi0 = 180.0
			n = 6.0
			V = 2.0
			# case (i) 
			if hyb_j == 'sp3' and els_j in ('O','S'):
				phi0 = 180.0
				n = 2.0
				U_j = UFF4MOF_atom_parameters[fft_j][6]
				U_k = UFF4MOF_atom_parameters[fft_k][6]
				V = 5 * np.sqrt(U_j*U_k) * (1.0 + 4.18 * np.log(bond_order))
			elif hyb_k == 'sp3' and els_k in ('O','S'):
				phi0 = 180.0
				n = 2.0
				U_j = UFF4MOF_atom_parameters[fft_j][6]
				U_k = UFF4MOF_atom_parameters[fft_k][6]
				V = 5 * np.sqrt(U_j*U_k) * (1.0 + 4.18 * np.log(bond_order))
			# case (j) not needed for the current ToBaCCo MOFs

		# case (c, d, e, f)
		elif hyb_j in ('sp2', 'resonant') and hyb_k in ('sp2', 'resonant'):
			phi0 = 180.0
			n = 2.0
			U_j = UFF4MOF_atom_parameters[fft_j][6]
			U_k = UFF4MOF_atom_parameters[fft_k][6]
			V = 5 * np.sqrt(U_j*U_k) * (1.0 + 4.18 * np.log(bond_order))

		# case (g)
		elif hyb_j == 'sp1' or hyb_k == 'sp1':
			return 'NA'

		elif hyb_j == 'NA' or hyb_k == 'NA':
			return 'NA'
		
		# divide by multiplicity and halve to match UFF paper
		V /= mult
		V *= 0.5
		d = -1.0 * np.cos(math.radians(n*phi0))

		return ('harmonic', V, int(d), int(n))

	def improper_parameters(self, fft_i, O_2_flag):
		
		if fft_i in ('N_R', 'C_R', 'C_2'):

			# constants for C_R and N_R
			C0 = 1.0
			C1 = -1.0
			C2 = 0.0
			K = 6.0/3.0
			al = 1

			# constants for bound O_2
			if O_2_flag:
				K = 50.0/3.0

		else:
			return None

		return ('fourier', K, C0, C1, C2, al)

	def pair_parameters(self, charges=False):
		
		UFF4MOF_atom_parameters = self.args['FF_parameters']
		atom_types = self.atom_types
		params = {}
		comments = {}

		# determine style and special bonds
		if charges:
			style = 'lj/cut/coul/long'
			sb = 'lj/coul 0.0 0.0 1.0'
		else:
			style = 'lj/cut'
			sb = 'lj 0.0 0.0 1.0'

		for a in atom_types:
			ID = atom_types[a]
			data = UFF4MOF_atom_parameters[a]
			x_i = data[2] * (2**(-1.0/6.0))
			D_i = data[3]
			params[ID] = (style, D_i, x_i)
			comments[ID] = [a,a]

		self.pair_data = {'params':params, 'style':style, 'special_bonds':sb, 'comments':comments}

	def enumerate_bonds(self):

		SG = self.system['graph']
		bond_order_dict = self.args['bond_orders']

		bonds = {}
		for e in SG.edges(data=True):

			i,j,data = e
			fft_i = SG.node[i]['force_field_type']
			fft_j = SG.node[j]['force_field_type']
			bond_type = data['bond_type']

			# look for the bond order, otherwise use the convention based on the bond type
			try:
				bond_order = bond_order_dict[(fft_i,fft_j)]
			except KeyError:
				try:
					bond_order = bond_order_dict[(fft_j,fft_i)]
				except KeyError:
					bond_order = bond_order_dict[bond_type]

			bond = tuple(sorted([fft_i, fft_j]) + [bond_order])

			# add to list if bond type already exists, else add a new type
			try:
				bonds[bond].append((i,j))
			except KeyError:
				bonds[bond] = [(i,j)]

			data['bond_order'] = bond_order

		bond_params = {}
		bond_comments = {}
		all_bonds = {}
		ID = 0
		count = 0
		# index bonds by ID
		for b in bonds:

			ID += 1
			bond_order = b[2]
			bond = (b[0], b[1])
			params = self.bond_parameters(bond, float(bond_order))
			bond_params[ID] = list(params)
			bond_comments[ID] = list(bond) + ['bond order=' + str(bond_order)]
			all_bonds[ID] = bonds[b]
			count += len(bonds[b])

		self.bond_data = {'all_bonds':all_bonds, 'params':bond_params, 'style':'harmonic', 'count':(count, len(all_bonds)), 'comments':bond_comments}

	def enumerate_angles(self):
		
		SG = self.system['graph']
		bonds = self.bond_data['all_bonds']
		bond_params = self.bond_data['params']
		inv_bonds = dict((b,bt) for bt in bonds for b in bonds[bt])
		angles = {}

		for n in SG.nodes(data=True):

			name, data = n
			nbors = list(SG.neighbors(name))

			for comb in itertools.combinations(nbors, 2):

				j = name
				i, k = comb

				fft_i = SG.node[i]['force_field_type']
				fft_j = SG.node[j]['force_field_type']
				fft_k = SG.node[k]['force_field_type']

				sort_ik = sorted([(fft_i,i),(fft_k,k)], key=lambda x:x[0])
				fft_i, i = sort_ik[0]
				fft_k, k = sort_ik[1]

				# look up bond constants (don't need to calculate again, yay!)
				try:
					bond_type_ij = inv_bonds[(i,j)]
				except KeyError:
					bond_type_ij = inv_bonds[(j,i)]
				try:
					bond_type_jk = inv_bonds[(j,k)]
				except KeyError:
					bond_type_jk = inv_bonds[(k,j)]

				r_ij = bond_params[bond_type_ij][2]
				r_jk = bond_params[bond_type_jk][2]

				angle = sorted((fft_i, fft_k))
				angle = (angle[0], fft_j, angle[1], r_ij, r_jk)

				# add to list if angle type already exists, else add a new type
				try:
					angles[angle].append((i,j,k))
				except KeyError:
					angles[angle] = [(i,j,k)]

		angle_params = {}
		angle_comments = {}
		all_angles = {}
		ID = 0
		count = 0
		styles = []

		# index angles by ID
		for a in angles:

			ID += 1
			fft_i, fft_j, fft_k, r_ij, r_jk = a
			angle = (fft_i, fft_j, fft_k)
			params = self.angle_parameters(angle, r_ij, r_jk)
			styles.append(params[0])
			angle_params[ID] = list(params)
			angle_comments[ID] = list(angle)
			all_angles[ID] = angles[a]
			count += len(angles[a])
		
		styles = set(styles)
		if len(styles) == 1:
			style = list(styles)[0]
		else:
			style = 'hybrid ' + ' '.join(styles)

		self.angle_data = {'all_angles':all_angles, 'params':angle_params, 'style':style, 'count':(count, len(all_angles)), 'comments':angle_comments}

	def enumerate_dihedrals(self):
		
		SG = self.system['graph']
		dihedrals = {}
		dihedral_params = {}

		for e in SG.edges(data=True):

			j,k = e[0:2]
			fft_j = SG.node[j]['force_field_type']
			fft_k = SG.node[k]['force_field_type']
			hyb_j = SG.node[j]['hybridization']
			hyb_k = SG.node[k]['hybridization']
			els_j = SG.node[j]['element_symbol']
			els_k = SG.node[k]['element_symbol']
			bond_order = e[2]['bond_order']
			nodes = (j,k)

			nbors_j = [n for n in SG.neighbors(j) if n != k]
			nbors_k = [n for n in SG.neighbors(k) if n != j]

			il_pairs = list(itertools.product(nbors_j, nbors_k))
			dihedral_list = [(p[0],j,k,p[1]) for p in il_pairs]

			bond = sorted([fft_j, fft_k])
			bond = (bond[0], bond[1], bond_order)
			hybridization = (hyb_j, hyb_k)
			element_symbols = (els_j, els_k)

			# here I calculate  parameters for each dihedral (I know) but I prefer identifying
			# those dihedrals before passing to the final dihedral data construction.
			params = self.dihedral_parameters(bond, hybridization, element_symbols, nodes)
			
			if params != 'NA':
				try:
					dihedrals[bond].extend(dihedral_list)
				except KeyError:
					dihedrals[bond] = dihedral_list
					dihedral_params[bond] = params

		all_dihedrals = {}
		dihedral_comments = {}
		indexed_dihedral_params = {}
		ID = 0
		count = 0
		for d in dihedrals:

			ID += 1
			dihedral = ('X', d[0], d[1], 'X')
			params = dihedral_params[d]
			all_dihedrals[ID] = dihedrals[d]
			indexed_dihedral_params[ID] = list(dihedral_params[d])
			dihedral_comments[ID] = list(dihedral) + ['bond order=' + str(d[2])]
			count += len(dihedrals[d])

		self.dihedral_data = {'all_dihedrals':all_dihedrals, 'params':indexed_dihedral_params, 'style':'harmonic', 'count':(count, len(all_dihedrals)), 'comments':dihedral_comments}

	def enumerate_impropers(self):
		
		SG = self.system['graph']
		impropers = {}

		for n in SG.nodes(data=True):
			
			i, data = n
			nbors = list(SG.neighbors(i))

			if len(nbors) == 3:
				
				fft_i = data['force_field_type']
				fft_nbors = tuple(sorted([SG.node[m]['force_field_type'] for m in nbors]))
				O_2_flag = False
				# force constant is much larger if j,k, or l is O_2
				if 'O_2' in fft_nbors or 'O_2_M' in fft_nbors:
					O_2_flag = True
				j,k,l = nbors

				# only need to consider one combination
				imps = [[i, j, k, l]]

				try:
					impropers[(fft_i, O_2_flag)].extend(imps)
				except KeyError:
					impropers[(fft_i, O_2_flag)] = imps

		all_impropers = {}
		improper_params = {}
		improper_comments = {}
		ID = 0
		count = 0
		for i in impropers:
		
			fft_i, O_2_flag = i	

			params = self.improper_parameters(fft_i, O_2_flag)

			if params != None:
				ID += 1
				improper_params[ID] = list(params)
				improper_comments[ID] = [i[0], 'X', 'X', 'X', 'O_2 present=' + str(O_2_flag)]
				all_impropers[ID] = impropers[i]
				count += len(impropers[i])
				
		self.improper_data = {'all_impropers':all_impropers, 'params':improper_params, 'style':'fourier', 'count':(count, len(all_impropers)), 'comments':improper_comments}

	def compile_force_field(self, charges=False):

		self.type_atoms()
		self.pair_parameters(charges)
		self.enumerate_bonds()
		self.enumerate_angles()
		self.enumerate_dihedrals()
		self.enumerate_impropers()
		
