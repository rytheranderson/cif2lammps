TraPPE =  {
	'O2': {
		'pair': {'style': 'lj/cut/coul/long', 'vdW': {'O_O2': (0.0974,3.02), 'O_com': (0.0,0.0)}, 'charges': {'O_O2': -0.113, 'O_com': 0.226}},
		'bonds': {('O_O2','O_com'): ('harmonic',100.0,0.604)}, # molecule should be kept rigid, force constants don't matter
		'angles': {('O_O2','O_com','O_O2'): ('harmonic',100.0,180.0)}, # molecule should be kept rigid, force constants don't matter
		'dihedrals': None,
		'impropers': None
	},
	'N2': {
		'pair': {},
		'bonds': {},
		'angles': {},
		'dihedrals': None,
		'impropers': None
	},
	'H2O1': {
		'pair': {},
		'bonds': {},
		'angles': {},
		'dihedrals': None,
		'impropers': None
	}
}

TIP4P =  {
	'H2O1': {
		'pair': {'style': 'lj/cut/tip4p/long', 'vdW': {'H_w': (0.0,0.0), 'O_w': (0.16275,3.16435)}, 'charges': {'H_w': 0.5242, 'O_w': 0.0}},
		'bonds': {('H_w','O_w'): ('harmonic', 100.0, 0.9572)},
		'angles': {('H_w','O_w','H_w'): ('harmonic', 100.0, 104.52)},
		'dihedrals': None,
		'impropers': None
	},
	'Cl1': {
		'pair': {'style': 'lj/cut/coul/long', 'vdW': {'Cl_Cl1': (0.22700, 3.51638)}, 'charges': {'Cl_Cl1': -1.0}},
		'bonds': None,
		'angles': None,
		'dihedrals': None,
		'impropers': None
	}
}