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
	'H2O': {
		'pair': {},
		'bonds': {},
		'angles': {},
		'dihedrals': None,
		'impropers': None
	}
}

TIP4P =  {
	'H2O': {
		'pair': {'style': 'lj/cut/coul/long', 'vdW': {'H': (0.0,0.0), 'O': (0.16275,3.16435), 'X':(0.0,0.0)}, 'charges': {'H': 0.5242, 'O': 0.0, 'X': -1.0484}},
		'bonds': {('O','H'): ('zero',None,None), ('O','X'): ('zero',None,None)},
		'angles': {('H','O','H'): ('zero',None,None), ('H','O','X'): ('zero',None,None)},
		'dihedrals': None,
		'impropers': None
	}
}