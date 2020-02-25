
TraPPE =  {
	'O2': {
		'pair': {'style':'lj/cut/coul/long', 'vdW': {'O_O2': (0.0974,3.02), 'O_com': (0.0,0.0)}, 'charges': {'O_O2':-0.113, 'O_com':0.226}},
		'bonds': {('O_O2','O_com'): ('harmonic',100.0,0.604)}, # use fix rigid or fix shake
		'angles': {('O_O2','O_com','O_O2'): ('harmonic',100.0,180.0)}, # use fix rigid or fix shake
		'dihedrals': None,
		'impropers': None
	},
	'N2': {
		'pair': {},
		'bonds': {},
		'angles': {},
		'dihedrals': {},
		'impropers': {}
	},
	'H2O': {
		'pair': {},
		'bonds': {},
		'angles': {},
		'dihedrals': {},
		'impropers': {}
	}
}