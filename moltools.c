/***************************************************************************

    moltools-python

    Python module for manipulation of atomic coordinates
    Copyright (C) 2012, Borys Szefczyk

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 ***************************************************************************/


#include "moltools.h"

//#include "periodic_table.c"




/*static PyObject *exposed_energy(PyObject *self, PyObject *args, PyObject *kwds) {
	FFType ff = FF_NONE;
	char *ff_type;
	float box[3] = { 0.0, 0.0, 0.0 };
	unsigned short int pbc = 0;
	double energy;
	int i, type;
	const char *ff_type_map[] = { "12-6", NULL };

	static char *kwlist[] = {
		"coordinates", "types", "ff_type", "ff", "box", NULL };

	PyObject *py_types, *py_ff;
	PyArrayObject *py_coords, *py_box = NULL;
	PyObject *py_result = NULL;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|sO!O!", kwlist,
			&PyArray_Type, &py_coords,
			&PyList_Type, &py_types,
			&ff_type,
			&PyList_Type, &py_ff,
			&PyArray_Type, &py_box))
		return NULL;

	if ( py_box != NULL ) {
		type = PyArray_TYPE(py_box);
		switch(type) {
			case NPY_FLOAT:
				box[0] = *( (float*) PyArray_GETPTR1(py_box, 0) );
				box[1] = *( (float*) PyArray_GETPTR1(py_box, 1) );
				box[2] = *( (float*) PyArray_GETPTR1(py_box, 2) );
				break;
			case NPY_DOUBLE:
				box[0] = *( (double*) PyArray_GETPTR1(py_box, 0) );
				box[1] = *( (double*) PyArray_GETPTR1(py_box, 1) );
				box[2] = *( (double*) PyArray_GETPTR1(py_box, 2) );
				break;
			default:
				PyErr_SetString(PyExc_ValueError, "Incorrect type in box vector");
				return NULL;
		}
		pbc = 1;
	}

	i = 0;
	while ( ff_type_map[i] != NULL ) {
		if (!strcmp(ff_type_map[i], ff_type)) ff = i;
		i += 1;
	}
	if (ff == FF_NONE) {
		PyErr_SetString(PyExc_ValueError, "Force field not implemented");
		return NULL; }

	energy = evaluate_energy(py_coords, py_types, ff, py_ff, box);
	py_result = PyFloat_FromDouble(energy);
	return py_result;
}*/


static PyObject *mass_list(PyObject *self, PyObject *args) {
	int nat, i, j;
	extern Element element_table[];
	char *symbol;
	double *masses;
	npy_intp dims[1];
	PyObject *py_symbols;
	PyObject *py_result = NULL;

	if(!PyArg_ParseTuple(args, "O!", &PyList_Type, &py_symbols))
		return NULL;

	nat = PyList_Size(py_symbols);
	masses = (double*) malloc( nat * sizeof(double) );

	for ( i = 0; i < nat; i++ ) {
		symbol = PyString_AsString(PyList_GetItem(py_symbols, i));
		j = 0;
		while ( strcmp(element_table[j].symbol, symbol) && element_table[j].number != -1 ) j++;
		if(element_table[j].number == -1) {
			PyErr_SetString(PyExc_RuntimeError, "Symbol unrecognized.");
			return NULL;
		}
		masses[i] = element_table[j].mass;
	}

	dims[0] = nat; 
	py_result = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, (double*) masses);

	return py_result;

}


static PyObject *rot_inertia(PyObject *self, PyObject *args, PyObject *kwds) {
	return NULL;
}


static PyMethodDef moltoolsMethods[] = {
    {"read", (PyCFunction)exposed_read, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"dict = read(filename [, unit ] )\n"
		"\n"
		"Reads a coordinate file and returns data as a dictionary.\n"
		"Currently, supports XYZ, GRO, XTC and Molden formats. In some file\n"
		"formats the unit can be chosen between unit = angs, bohr, nm.\n"
		"\n"
		"XYZ:    supports extended files (charges in the fifth column)\n"
		"        as well as standard files; coordinates are assumed to be\n"
		"        in Angstroms. Also supports multiframe files. In that\n"
		"        case, returns a list of dictionaries.\n"
		"\n"
		"Molden: supports groups N_GEO, GEOCONV (energies only),\n"
		"        GEOMETRIES (XYZ only), ATOMS (Angstroms only).\n"
		"\n"
		"Output structure, dictionary including various keys, depending\n"
		"on what is found in the file:\n"
		"\n"
		"nofatoms        - number of atoms (int)\n"
		"comment         - comment line extracted from second line (str)\n"
		"symbols         - atom symbols (list)\n"
		"coordinates     - array of [nofatoms,3] coordinates\n"
		"                  (numpy array)\n"
		"charges         - array of [nofatoms] point charges\n"
		"                  (numpy array)\n"
		"energies        - energies of the subsequent configurations\n"
		"                  (numpy array)\n"
		"geometries      - list of dictionaries, one for each configuration\n"
		"atomic_numbers  - atomic numbers (nuclear charge, numpy array, int)\n"
		"\n" },
    {"write", (PyCFunction)exposed_write, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"write(FILENAME, SYMBOLS, COORDS, [ comment=COMMENT,\n"
		"      residues=RESNAM, residue_numbers=RESNUM, box=PBCBOX,\n"
		"      format=FORMAT, mode=MODE ])\n"
		"\n"
		"Write structure in XYZ or GRO format. Following parameters are\n"
		"obligatory:\n"
		"FILENAME - file name (string)\n"
		"SYMBOLS  - atomic symbols (list)\n"
		"COORDS   - coordinates (numpy array)\n"
		"Following parameters are optional:\n"
		"FORMAT   - file format, 'XYZ' or 'GRO' (default 'XYZ')\n"
		"MODE     - writing mode, 'w' or 'a' (default 'w')\n"
		"COMMENT  - comment text (string), default is empty\n"
		"RESNAM   - residue names (list of strings), default are empty\n"
		"           strings\n"
		"RESNUM   - residue numbers (list if int), default is 1\n"
		"PBCBOX   - Dimentions of the box, used in PBC (numpy array);\n"
		"           required when GRO format is used\n"
		"\n" },
	{"findBonds", (PyCFunction)find_bonds, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"findBonds(symbols, coordinates, factor=1.3, format=list)\n"
		"\n"
		"Finds connections between atoms. 'symbols' is a list of atomic\n"
		"types, 'coordinates' is a numpy array of coordinates and factor\n"
		"is a factor used to multiply the covalent radii. Depending on the\n"
		"optional parameter 'format', returns either a dictionary of\n"
		"indices { <center> : (<center1>, <center2>, ... ) ... } for every\n"
		"atom or a list of unique bonds as tuples [ (center1, center2),\n"
		"(center1, center3), ... ]\n"
		"\n" },
	{"findMolecules", (PyCFunction)find_molecules, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"findMolecules(numberOfAtoms, bonds)\n"
		"\n"
		"Given a list of bonds as tuples, find sets of topologically\n"
		"connected atoms (molecules).\n"
		"\n" },
    /*{"energy", (PyCFunction)exposed_energy, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"energy(coordinates, types, ff_type, ff, box=None)\n"
		"\n"
		"Evaluate the energy of the configuration <coordinates>. <types> are assigned\n"
		"from <ff> and the energy is calculated according to <ff_type>, using PBC only\n"
		"if <box> dimensions are specified.\n"
		"\n" },*/
	{"mass_list", (PyCFunction)mass_list, METH_VARARGS,
		"\n"
		"mass_list(symbols)\n"
		"\n"
		"For a given list of atomic symbols, returns numpy array of masses.\n"
		"\n" },
	{"COM", (PyCFunction)centerofmass, METH_VARARGS,
		"\n"
		"COM(coordinates, masses)\n"
		"\n"
		"For given coordinates (numpy array) evaluate the center of mass\n"
		"(numpy array).\n"
		"\n" },
	{"MEP_distance", (PyCFunction)mep_distance, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"MEP_distance(coordinates1, coordinates2, masses=None)\n"
		"\n"
		"Computes mass-weighted distance of two structures, as used in\n"
		"IRC/MEP calculations. If masses are omitted, no weighting is done.\n"
		"\n" },
	{"distanceMatrix", (PyCFunction)distanceMatrix, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"distanceMatrix(coordinates, box=None, squared=False)\n"
		"\n"
		"Calculates distance matrix, taking PBC into account when specified.\n"
		"box - dimentions of the PBC box (numpy array).\n"
		"squared - use square of distance to avoid calculating roots.\n"
		"\n" },
	{"measureAngleCosine", (PyCFunction)measureAngleCosine, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"measureAngleCosine(atom1, atom2, atom3, box=None)\n"
		"\n"
		"Returns the cosine of the angle between 1-2-3.\n"
		"box - dimentions of the PBC box (numpy array).\n"
		"\n" },
	{"findHBonds", (PyCFunction)findHBonds, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"findHBonds(symbols, coordinates, acceptors,\n"
		"           box=None, cutoff=3.5, carbon_cutoff=cutoff, angle_cutoff=30)\n"
		"\n"
		"Returns a list of hydrogen bonds according to definition similar to Chandler'.\n"
		"donor - acceptor distance must be less than cutoff, unless if donor/acceptor\n"
		"is carbon; in such a case the carbon_cutoff is used. Additionally, the angle\n"
		"H - donor - acceptor is required to be less than angle_cutoff degrees.\n"
		"\n" },
	{"inertia", (PyCFunction)inertia, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"inertia_tensor = inertia(coordinates, masses)\n"
		"\n"
		"Computes the moments of inertia tensor of a molecule.\n"
		"\n" },
	{"rot_inertia", (PyCFunction)rot_inertia, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"coordinates_rot, principal_moments = rot_inertia(coordinates, inertia_tensor)\n"
		"\n"
		"Rotates the coordinates so that the principal moments of inertia\n"
		"are aligned with axes of the coordinate system.\n"
		"\n" },

    {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC initmoltools(void)
{
    PyObject *md;
	PyObject *exposed_atom_symbols, *exposed_atom_names;
	PyObject *exposed_atom_masses, *exposed_symbol2number;
	PyObject *exposed_covalentradii;
	extern PyTypeObject EAMffType;
	extern PyTypeObject MoleculeType;

	/* Use system-wide locale, but make sure that decimal point is a point! */
	setlocale(LC_ALL, "");
	setlocale(LC_NUMERIC, "C");

	if (PyType_Ready(&EAMffType) < 0)
		return;
	if (PyType_Ready(&MoleculeType) < 0)
		return;

    md = Py_InitModule3("moltools", moltoolsMethods,
	     "The moltools module provides some classes and functions related to molecular "
		 "modelling. The idea is to facilitate writing scripts for processing molecular "
		 "data, using standard types. For atomic coordinates, numpy arrays are used, "
		 "since they are fast and implement linear algebra.");
    if(md == NULL) return;

	import_array();

	if( build_tables(&exposed_atom_symbols, &exposed_atom_names,
		&exposed_atom_masses, &exposed_symbol2number,
		&exposed_covalentradii) ) return;

	/* Add the lists to the module */
	PyModule_AddObject(md, "AtomicSymbols", exposed_atom_symbols);
	PyModule_AddObject(md, "AtomicNames", exposed_atom_names);
	PyModule_AddObject(md, "AtomicMasses" , exposed_atom_masses);
	PyModule_AddObject(md, "AtomicNumbers" , exposed_symbol2number);
	PyModule_AddObject(md, "CovalentRadii" , exposed_covalentradii);

	Py_INCREF(&EAMffType);
	PyModule_AddObject(md, "EAMff", (PyObject *)&EAMffType);

	Py_INCREF(&MoleculeType);
	PyModule_AddObject(md, "Molecule", (PyObject *)&MoleculeType);
}

