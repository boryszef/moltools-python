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

#include "periodic_table.c"



static PyObject *exposed_energy(PyObject *self, PyObject *args, PyObject *kwds) {
	char *ff_type;
	float box[3] = { 0.0, 0.0, 0.0 };
	long int pbc = 0;

	static char *kwlist[] = {
		"coordinates", "types", "ff_type", "ff", "box",
		"pbc", NULL };

	PyObject *py_types, *py_coords, *py_ff, *py_box, *py_pbc;
	PyObject *py_result;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|sO!O!O!", kwlist,
			&PyArray_Type, &py_coords,
			&PyList_Type, &py_types,
			&ff_type,
			&PyDict_Type, &py_ff,
			&PyArray_Type, &py_box,
			&PyBool_Type, &py_pbc))
		return NULL;

	pbc = PyInt_AsLong(py_pbc);
	if ( pbc ) {
		box[0] = *( (double*) PyArray_GETPTR1(py_box, 0) );
		box[1] = *( (double*) PyArray_GETPTR1(py_box, 1) );
		box[2] = *( (double*) PyArray_GETPTR1(py_box, 2) );
	}

	py_result = PyFloat_FromDouble(0.0);
	return py_result;
}



static PyObject *exposed_read(PyObject *self, PyObject *args, PyObject *kwds) {

	const char *filename;
	char *line;
	char *unit = NULL;
	enum { GUESS, XYZ, MOLDEN, FRAC, GRO } type = GUESS;
	char *str_type = NULL;
	float factor;
	FILE *fd, *test;
	long int fpos;
	struct stat fst;
	char ext[5];

	static char *kwlist[] = {"file", "format", "unit", NULL};

	PyObject *py_result, *py_list;

	py_result = NULL;

	/* Process the arguments */
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|ss",
		kwlist, &filename, &str_type, &unit)) return NULL;

	/* Get the unit measure for coordinates and set the factor accordingly */
	if (unit == NULL) factor = 1.0;
	else {
		make_lowercase(unit);
		if ( !strcmp(unit, "angs") ) factor = 1.0;
		else if ( !strcmp(unit, "bohr") ) factor = BOHR;
		else if ( !strcmp(unit, "nm") ) factor = 10.0;
		else {
			PyErr_SetString(PyExc_ValueError, "Unrecognized measure unit name");
			return NULL; }
	}

	/* Set the enum symbol of the file format */
	if ( str_type != NULL ) {
		if      ( !strcmp(str_type,    "XYZ") ) type = XYZ;
		else if ( !strcmp(str_type, "MOLDEN") ) type = MOLDEN;
		else if ( !strcmp(str_type,   "FRAC") ) type = FRAC;
		else if ( !strcmp(str_type,    "GRO") ) type = GRO;
	}

	/* Open the coordinate file */
	if ( (fd = fopen(filename, "r")) == NULL ) {
		PyErr_SetFromErrno(PyExc_IOError);
		return NULL; }

	/* Guess the file format, if not given explicitly */
	if ( type == GUESS ) {
		strcpy(ext, filename + strlen(filename) - 4);
		if      ( !strcmp(ext, ".xyz") ) type = XYZ;
		else if ( !strcmp(ext, ".gro") ) type = GRO;
		else {
			/* Extract the first line */
			if ( (test = fopen(filename, "r")) == NULL ) {
				PyErr_SetFromErrno(PyExc_IOError);
				return NULL; }
			if ( (line = readline(test)) == NULL ) {
				PyErr_SetFromErrno(PyExc_IOError);
				return NULL; }
			make_lowercase(line);
			stripline(line);
			fclose(test);

			/* Perhaps it's Molden format? */
			if ( !strcmp(line, "[molden format]") ) type = MOLDEN;

			free(line);
		}
	}

	/* Router */
	switch(type) {
		case XYZ:
			py_result = read_xyz(fd, factor);
			/* Support for multiframe XYZ: check the current position in
 			 * the file; if we are near the end, finish; else, continue
 			 * reading. The criteria is that there should be no more than
 			 * 3 bytes left. */
			fpos = ftell(fd);
			if ( fpos == -1 ) {
				PyErr_SetFromErrno(PyExc_IOError);
				return NULL; }
			if ( fstat(fileno(fd), &fst) == -1 ) {
				PyErr_SetFromErrno(PyExc_IOError);
				return NULL; }
			if ( fst.st_size - fpos > 3 ) {
				py_list = PyList_New(1);
				PyList_SetItem(py_list, 0, py_result);
			} else {
				break;
			}
			while ( fst.st_size - fpos > 3 ) {
				py_result = read_xyz(fd, factor);
				PyList_Append(py_list, py_result);
				fpos = ftell(fd);
				if ( fpos == -1 ) {
					PyErr_SetFromErrno(PyExc_IOError);
					return NULL; }
			}
			py_result = py_list;
			break;
		case MOLDEN:
			py_result = read_molden(fd);
			break;
		case FRAC:
			py_result = read_fractional(fd);
			break;
		case GRO:
			py_result = read_gro(fd);
			break;
		/* If the file format is GUESS or different,
		   it means we've failed to guess :-(        */
		case GUESS:
		default:
			PyErr_SetString(PyExc_ValueError, "Unsupported file format");
			return NULL;
	}

	fclose(fd);
	return py_result;
}


static PyObject *exposed_write(PyObject *self, PyObject *args, PyObject *kwds) {

	char *filename;
	char *comment = NULL;
	char *str_format = NULL;
	char *mode = NULL;
	enum { XYZ, GRO } format = XYZ;
	FILE *fd;
	int nat, i;

	static char *kwlist[] = {
		"file", "symbols", "coordinates", "comment", "residues",
		"residue_numbers", "box", "format", "mode", NULL };

	PyObject *py_symbols, *py_coords, *val;
	PyObject *py_resnam = NULL, *py_resid = NULL, *py_box = NULL;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sO!O!|sO!O!O!ss", kwlist,
			&filename,
			&PyList_Type, &py_symbols,
			&PyArray_Type, &py_coords,
			&comment,
			&PyList_Type, &py_resnam,
			&PyList_Type, &py_resid,
			&PyArray_Type, &py_box,
			&str_format,
			&mode))
		return NULL;

	if ( mode == NULL || !strcmp(mode, "w") ) {
		if( (fd = fopen(filename, "w")) == NULL ) {
			PyErr_SetFromErrno(PyExc_IOError);
			return NULL; }
	} else if ( !strcmp(mode, "a") ) {
		if( (fd = fopen(filename, "a")) == NULL ) {
			PyErr_SetFromErrno(PyExc_IOError);
			return NULL; }
	} else {
		PyErr_SetString(PyExc_ValueError, "Unsupported file mode");
		return NULL; }

	if ( str_format != NULL ) {
		if      ( !strcmp(str_format, "XYZ") ) format = XYZ;
		else if ( !strcmp(str_format, "GRO") ) format = GRO;
	}

	switch(format) {
		case XYZ:
			if ( write_xyz(fd, py_symbols, py_coords, comment) == -1 )
				return NULL;
			break;
		case GRO:
			if ( py_box == NULL ) {
				PyErr_SetString(PyExc_ValueError, "Box vector must be given, when GRO format is used");
				return NULL; }
			nat = PyList_Size(py_symbols);
			if ( py_resnam == NULL ) {
				py_resnam = PyList_New(nat);
				val = Py_BuildValue("s", "");
				for ( i = 0; i < nat; i++ ) {
					PyList_SetItem(py_resnam, i, val);
					Py_INCREF(val);
				}
				Py_DECREF(val);
			}
			if ( py_resid == NULL ) {
				py_resid = PyList_New(nat);
				val = Py_BuildValue("i", 1);
				for ( i = 0; i < nat; i++ ) {
					PyList_SetItem(py_resid, i, val);
					Py_INCREF(val);
				}
				Py_DECREF(val);
			}
			if ( write_gro(fd, py_symbols, py_coords, comment,
			               py_resnam, py_resid, py_box) == -1 )
				return NULL;
			break;
		default:
			PyErr_SetString(PyExc_RuntimeError, "This was unexpected...");
			return NULL;
	}

	fclose(fd);

	return Py_None;
}


static PyObject *find_bonds(PyObject *self, PyObject *args) {
	int i, j, nat;
	double ax, ay, az, bx, by, bz, ar, br, dist;
	npy_intp *numpyint;

	PyObject *val1, *val2, *tmp_list;
	PyObject *py_symbols, *py_coords, *py_types;
	PyObject *py_result = NULL;

	if(!PyArg_ParseTuple(args, "O!O!O!",
			&PyList_Type, &py_symbols,
			&PyArray_Type, &py_coords,
			&PyDict_Type, &py_types))
		return NULL;

	nat = PyList_Size(py_symbols);
	numpyint = PyArray_DIMS(py_coords);
	if (numpyint[0] != nat) {
		PyErr_SetString(PyExc_RuntimeError, "Size of symbol and coordinate lists does not match.");
		return NULL;
	}

	py_result = PyDict_New();
	for (i = 0; i < nat; i++) {
		ax = *( (double*) PyArray_GETPTR2(py_coords, i, 0) );
		ay = *( (double*) PyArray_GETPTR2(py_coords, i, 1) );
		az = *( (double*) PyArray_GETPTR2(py_coords, i, 2) );
		val1 = PyList_GetItem(py_symbols, i); // borrowed
		val2 = PyDict_GetItem(py_types, val1); // borrowed
		ar = PyFloat_AsDouble(val2);
		tmp_list = PyList_New(0); // new
		for (j = 0; j < nat; j++) {
			if (i == j) continue;
			bx = *( (double*) PyArray_GETPTR2(py_coords, j, 0) );
			by = *( (double*) PyArray_GETPTR2(py_coords, j, 1) );
			bz = *( (double*) PyArray_GETPTR2(py_coords, j, 2) );
			val1 = PyList_GetItem(py_symbols, j); // borrowed
			val2 = PyDict_GetItem(py_types, val1); // borrowed
			br = PyFloat_AsDouble(val2);
			dist = sq(bx-ax) + sq(by-ay) + sq(bz-az);
			if (dist < sq(ar+br)) {
				val2 = PyInt_FromLong(j); // new
				PyList_Append(tmp_list, val2);
				Py_DECREF(val2);
			}
		}
		val1 = PyInt_FromLong(i); // new
		if (PyList_Size(tmp_list))
			PyDict_SetItem(py_result, val1, tmp_list);
		Py_DECREF(val1);
		Py_DECREF(tmp_list);
	}

	return py_result;
}


static PyMethodDef moltoolsMethods[] = {
    {"read", (PyCFunction)exposed_read, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"dict = read(filename [, unit ] )\n"
		"\n"
		"Reads a coordinate file and returns data as a dictionary.\n"
		"Currently, supports only XYZ and Molden formats. In some file formats\n"
		"the unit can be chosen between unit = angs, bohr, nm.\n"
		"\n"
		"XYZ:    supports extended files (charges in the fifth column) as well\n"
		"        as standard files; coordinates are assumed to be in Angstroms.\n"
		"        Also supports multiframe files. In that case, returns a list\n"
		"        of dictionaries.\n"
		"\n"
		"Molden: supports groups N_GEO, GEOCONV (energies only), GEOMETRIES\n"
		"        (XYZ only), ATOMS (Angstroms only).\n"
		"\n"
		"Output structure, dictionary including various keys, depending on what\n"
		"is found in the file:\n"
		"\n"
		"number_of_atoms - number of atoms (int)\n"
		"comment         - comment line extracted from second line (str)\n"
		"symbols         - atom symbols (list)\n"
		"coordinates     - array of [number_of_atoms,3] coordinates (numpy array)\n"
		"charges         - array of [number_of_atoms] point charges (numpy array)\n"
		"energies        - energies of the subsequent configurations (numpy array)\n"
		"geometries      - list of dictionaries, one for each configuration\n"
		"atomic_numbers  - atomic numbers (nuclear charge, numpy array, int)\n"
		"\n" },
    {"write", (PyCFunction)exposed_write, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"FIXME"
		"\n" },
	{"find_bonds", (PyCFunction)find_bonds, METH_VARARGS,
		"\n"
		"find_bonds(symbols, coordinates, vdWradii) -> topology\n"
		"\n"
		"Finds connections between atoms. 'symbols' is a list of atomic types,\n"
		"'coordinates' is a numpy array of coordinates and 'vdWradii' is a dictionary\n"
		"of van der Waals radii for atomic types found in 'symbols'. Returns 'topology,'\n"
		"which is a dictionary of indices <center> : (<center1>, <center2>, ... )\n"
		"\n" },
    {"energy", (PyCFunction)exposed_energy, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"energy(coordinates, types, ff_type, ff, box=array(0.0, 0.0, 0.0), pbc=False)\n"
		"\n"
		"Evaluate the energy of the configuration <coordinates>. <types> are assigned\n"
		"from <ff> and the energy is calculated according to <ff_type>, using <box>,\n"
		"only if pbc=True.\n"
		"\n" },
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


int build_tables(PyObject **list_symbols, PyObject **list_names,
                 PyObject **list_masses, PyObject **symbol2number) {

	extern Element element_table[];
	PyObject *key, *val;
	int elements_defined, idx;
	double mass;

	/* Check how many elements we know */
	elements_defined = 0;
	while ( element_table[elements_defined].number >= 0 )
		elements_defined += 1;

	/* Prepare lists of symbols, names and masses */
	*list_symbols  = PyList_New(elements_defined + 1);
	*list_names    = PyList_New(elements_defined + 1);
	*list_masses   = PyList_New(elements_defined + 1);
	*symbol2number = PyDict_New();

	/* Set the first element to None, so that indexing starts at 1 */
	val = Py_BuildValue("s", "");
	PyList_SetItem(*list_symbols, 0, val);
	PyList_SetItem(*list_names, 0, val);
	PyList_SetItem(*list_masses, 0, Py_None);

	/* Fill with data */
	for(idx = 0; idx < elements_defined; idx++) {

		/* Atom symbol */
		val = Py_BuildValue("s", element_table[idx].symbol);
		PyList_SetItem(*list_symbols, element_table[idx].number, val);

		/* Atom name */
		val = Py_BuildValue("s", element_table[idx].name);
		PyList_SetItem(*list_names, element_table[idx].number, val);

		/* Atom mass */
		mass = element_table[idx].mass;
		val = mass >= 0 ? Py_BuildValue("d", mass) : Py_None;
		PyList_SetItem(*list_masses, element_table[idx].number, val);

		/* Hash table symbol->number */
		key = Py_BuildValue("s", element_table[idx].symbol);
		val = Py_BuildValue("i", element_table[idx].number);
		PyDict_SetItem(*symbol2number, key, val);
		Py_DECREF(key);
		Py_DECREF(val);

	}

	return 0;
}


PyMODINIT_FUNC initmoltools(void)
{
    PyObject *md;
	PyObject *exposed_atom_symbols, *exposed_atom_names;
	PyObject *exposed_atom_masses, *exposed_symbol2number;

	/* Use system-wide locale, but make sure that decimal point is a point! */
	setlocale(LC_ALL, "");
	setlocale(LC_NUMERIC, "C");

    md = Py_InitModule("moltools", moltoolsMethods);
    if(md == NULL) return;
	import_array();

	if( build_tables(&exposed_atom_symbols, &exposed_atom_names,
		&exposed_atom_masses, &exposed_symbol2number) ) return;

	/* Add the lists to the module */
	PyModule_AddObject(md, "AtomicSymbols", exposed_atom_symbols);
	PyModule_AddObject(md, "AtomicNames", exposed_atom_names);
	PyModule_AddObject(md, "AtomicMasses" , exposed_atom_masses);
	PyModule_AddObject(md, "AtomicNumbers" , exposed_symbol2number);
}

