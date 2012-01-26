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



#define BUFFER_LENGTH 256
static PyObject *read_xyz(FILE *fd, float factor) {

	int nofatoms, pos;
	char *buffer, *buffpos, *token;
	float *xyz, *charges;
	short int charges_present;

	npy_intp dims[2];

	PyObject *key, *val, *py_result, *py_coord, *py_syms, *py_charges;


	/* Create the dictionary that will be returned */
	py_result = PyDict_New();


	/* Read number of atoms */
	if( (buffer = readline(fd)) == NULL) {
		PyErr_SetFromErrno(PyExc_IOError);
		return NULL; }
	if( sscanf(buffer, "%d", &nofatoms) != 1 ) {
		PyErr_SetString(PyExc_IOError, "Incorrect atom number");
		return NULL; }

    val = Py_BuildValue("i", nofatoms);
	key = PyString_FromString("number_of_atoms");
	PyDict_SetItem(py_result, key, val);
	Py_DECREF(key);
	Py_DECREF(val);


	/* Read the comment line */
	if((buffer = readline(fd)) == NULL) {
		PyErr_SetFromErrno(PyExc_IOError);
		return NULL; }
	buffer[strlen(buffer)-1] = '\0';

	val = Py_BuildValue("s", buffer);
	free(buffer);
	key = PyString_FromString("comment");
	PyDict_SetItem(py_result, key, val);
	Py_DECREF(key);
	Py_DECREF(val);


	/* Set-up the raw arrays for coordinates and charges */
	xyz = (float*) malloc(3 * nofatoms * sizeof(float));
	if(xyz == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return NULL; }
	charges = (float*) malloc(nofatoms * sizeof(float));
	if(charges == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return NULL; }
	charges_present = 0;


	py_syms = PyList_New(nofatoms);

	/* Atom loop */
	for(pos = 0; pos < nofatoms; pos++) {

		/* Get the whole line */
		if((buffer = readline(fd)) == NULL) {
			PyErr_SetFromErrno(PyExc_IOError);
			return NULL; }
		buffer[strlen(buffer)-1] = '\0';
		buffpos = buffer;

		/* Read symbol */
		token = strtok(buffpos, " ");
		val = Py_BuildValue("s", token);
		PyList_SetItem(py_syms, pos, val);

		/* Read coordinates */
		if ( (token = strtok(NULL, " ")) == NULL) {
			PyErr_SetString(PyExc_IOError, "Missing coordinate");
			return NULL; }
		xyz[3*pos + 0] = atof(token) * factor;
		if ( (token = strtok(NULL, " ")) == NULL) {
			PyErr_SetString(PyExc_IOError, "Missing coordinate");
			return NULL; }
		xyz[3*pos + 1] = atof(token) * factor;
		if ( (token = strtok(NULL, " ")) == NULL) {
			PyErr_SetString(PyExc_IOError, "Missing coordinate");
			return NULL; }
		xyz[3*pos + 2] = atof(token) * factor;

		/* Read charge, if present */
		token = strtok(NULL, " ");
		if ( token != NULL ) {

			/* This is bad: until now, there were no charges */
			if ( pos > 0 && !charges_present ) {
				PyErr_SetString(PyExc_IOError, "Unexpected charges found");
				return NULL;
			}

			charges_present = 1;
			charges[pos] = atof(token);

		} else {

			/* This is bad: we were expecting charges here and found nothing */
			if ( pos > 0 && charges_present ) {
				PyErr_SetString(PyExc_IOError, "Missing charges");
				return NULL;
			}
		}

		/* Free the line buffer */
		free(buffer);
	}


	/* Add symbols to the dictionary */
	key = PyString_FromString("symbols");
	PyDict_SetItem(py_result, key, py_syms);
	Py_DECREF(key);
	Py_DECREF(py_syms);


	/* Add coordinates to the dictionary */
	dims[0] = nofatoms;
	dims[1] = 3;
	py_coord = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, (float*) xyz);
	/***************************************************************
	 * Do not free the raw array! It will be still used by Python! *

	free(*xyz);
	free(xyz);

     * when the ver. 1.7 arrives, use PyArray_SetBaseObject        *
     * to prevent memory leaks.                                    *
     ***************************************************************/

	key = PyString_FromString("coordinates");
	PyDict_SetItem(py_result, key, py_coord);
	Py_DECREF(key);
	Py_DECREF(py_coord);


	/* Add charges, if present */
	if ( charges_present ) {
		py_charges = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT, charges);
		key = PyString_FromString("charges");
		PyDict_SetItem(py_result, key, py_charges);
		Py_DECREF(key);
		Py_DECREF(py_charges);
	} else
		/* Free the charges ONLY if the Python object was not created! */
		free(charges);


	return py_result;
}


/* Read molden file */
static PyObject *read_molden(FILE *fd) {

	char *line;
	char *buffpos, *token;
	char **line_store;
	int n_geo, i, nat;
	float *xyz;
	int *anum;

	double *geoconv;
	npy_intp dims[2];
	PyObject *py_anum, *py_syms, *py_geom, *py_geoconv, *py_result, *key, *val;


	/* Prepare dictionary */

	py_result = PyDict_New();

	/* Loop over the lines */

	line = readline(fd);

	while ( strlen(line) ) {

		stripline(line);
		make_lowercase(line);

		/* This is a start of a new block */
		if ( line[0] == '[' ) {

			/* Number of geometries present in the file */
			if ( !strcmp(line, "[n_geo]") ) {
				free(line);
				line = readline(fd);
				stripline(line);
				n_geo = atoi(line);
				key = PyString_FromString("number_of_geometries");
				val = Py_BuildValue("i", n_geo);
				PyDict_SetItem(py_result, key, val);
				Py_DECREF(key);
				Py_DECREF(val);

				/* This will be used by arrays depending on number of geoms */
				dims[0] = n_geo;
				dims[1] = 3;
			}

			/* Energy of the subsequent geometries */
			else if ( !strcmp(line, "[geoconv]" ) ) {
				free(line);
				line = readline(fd);
				geoconv = (double*) malloc(n_geo * sizeof(double));
				if(geoconv == NULL) {
					PyErr_SetFromErrno(PyExc_MemoryError);
					return NULL; }
				for( i = 0; i < n_geo; i++) {
					free(line);
					line = readline(fd);
					stripline(line);
					geoconv[i] = atof(line);
				}
				py_geoconv = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, geoconv);
				key = PyString_FromString("energies");
				PyDict_SetItem(py_result, key, py_geoconv);
				Py_DECREF(key);
				Py_DECREF(py_geoconv);
			}

			/* Energy of the subsequent geometries */
			else if ( !strcmp(line, "[geometries] (xyz)" ) ) {
				py_geom = PyList_New(n_geo);
				for ( i = 0; i < n_geo; i++ )
					PyList_SetItem(py_geom, i, read_xyz(fd, 1.0));
				key = PyString_FromString("geometries");
				PyDict_SetItem(py_result, key, py_geom);
				Py_DECREF(key);
				Py_DECREF(py_geom);
			}

			/* Section 'atoms' present - this is one-geometry file */
			else if ( !strcmp(line, "[atoms] angs" ) ) {

				n_geo = 1;
				key = PyString_FromString("number_of_geometries");
				val = Py_BuildValue("i", n_geo);
				PyDict_SetItem(py_result, key, val);
				Py_DECREF(key);
				Py_DECREF(val);

				free(line);
				line = readline(fd);
				stripline(line);

				/* We don't know how many atoms are there, so we have to
                   store the lines. */
				nat = 0;
				line_store = (char**) malloc(10 * sizeof(char*));
				if ( line_store == NULL ) {
					PyErr_SetFromErrno(PyExc_MemoryError);
					return NULL; }
				while ( line[0] != '[' ) {
					line_store[nat++] = line;
					if ( nat % 10 == 0 ) {
						line_store = realloc(line_store, (nat + 10) * sizeof(char*));
						if( line_store == NULL ) {
							PyErr_SetFromErrno(PyExc_MemoryError);
						return NULL; }
					}
					line = readline(fd);
					stripline(line);
				}

				xyz = (float*) malloc(3 * nat * sizeof(float));
				if(xyz == NULL) {
					PyErr_SetFromErrno(PyExc_MemoryError);
					return NULL; }
				anum = (int*) malloc(nat * sizeof(int));
				if(anum == NULL) {
					PyErr_SetFromErrno(PyExc_MemoryError);
					return NULL; }

				py_syms = PyList_New(nat);

				/* Loop over atoms */
				for ( i = 0; i < nat; i++ ) {

					buffpos = line_store[i];
					token = strtok(buffpos, " ");
					val = Py_BuildValue("s", token);
					PyList_SetItem(py_syms, i, val);

					token = strtok(NULL, " ");
					/* not used */

					token = strtok(NULL, " ");
					anum[i] = atoi(token);

					token = strtok(NULL, " ");
					xyz[3*i+0] = atof(token);

					token = strtok(NULL, " ");
					xyz[3*i+1] = atof(token);

					token = strtok(NULL, " ");
					xyz[3*i+2] = atof(token);

					/* Get rid of the line. */
					free(line_store[i]);

				}

				/* Free the stored line pointers. */
				free(line_store);

				/* Add symbols to the dictionary */
				key = PyString_FromString("symbols");
				PyDict_SetItem(py_result, key, py_syms);
				Py_DECREF(key);
				Py_DECREF(py_syms);

				/* Add coordinates to the dictionary */
				dims[0] = nat;
				dims[1] = 3;
				py_geom = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, (float*) xyz);
				key = PyString_FromString("coordinates");
				PyDict_SetItem(py_result, key, py_geom);
				Py_DECREF(key);
				Py_DECREF(py_geom);

				/* Add atomic numbers to the dictionary */
				py_anum = PyArray_SimpleNewFromData(1, dims, NPY_INT, (int*) anum);
				key = PyString_FromString("atomic_numbers");
				PyDict_SetItem(py_result, key, py_anum);
				Py_DECREF(key);
				Py_DECREF(py_anum);

				/* This is to avoid reading the next line! */
				continue;
			}

		}

		line = readline(fd);

	}

	return py_result;
}



static PyObject *exposed_read(PyObject *self, PyObject *args, PyObject *kwds) {

	const char *filename;
	char *line;
	char *unit = NULL;
	float factor;
	FILE *fd, *test;

	static char *kwlist[] = {"file", "unit", NULL};

	PyObject *py_result;

	py_result = NULL;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|s",
		kwlist, &filename, &unit)) return NULL;

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

	if ( (fd = fopen(filename, "r")) == NULL ) {
		PyErr_SetFromErrno(PyExc_IOError);
		return NULL; }

	/* Guess the file format */

	/* If the filename ends with .xyz, it is xyz file */
	if ( !strcmp(filename + strlen(filename) - 4, ".xyz") ) {

		py_result = read_xyz(fd, factor);

	} else {

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
		if ( !strcmp(line, "[molden format]") ) {
			py_result = read_molden(fd);
		}

		free(line);

	}

	fclose(fd);
	return py_result;
}


static PyObject *exposed_write(PyObject *self, PyObject *args, PyObject *kwds) {

	char *filename;
	char *comment = NULL;
	char *str_format = NULL;
	enum { XYZ, GRO } format = XYZ;
	FILE *fd;
	int nat, i;

	static char *kwlist[] = {"file", "symbols", "coordinates", "comment", "residues", "residue_numbers", "box", "format", NULL};

	PyObject *py_symbols, *py_coords, *val;
	PyObject *py_resnam = NULL, *py_resid = NULL, *py_box = NULL;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sO!O!|sO!O!O!s", kwlist,
			&filename,
			&PyList_Type, &py_symbols,
			&PyArray_Type, &py_coords,
			&comment,
			&PyList_Type, &py_resnam,
			&PyList_Type, &py_resid,
			&PyArray_Type, &py_box,
			&str_format))
		return NULL;

	if( (fd = fopen(filename, "w")) == NULL ) {
		PyErr_SetFromErrno(PyExc_IOError);
		return NULL;
	}

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
		"\nTODO: multi-frame xyz\n"
		"\n" },
    {"write", (PyCFunction)exposed_write, METH_VARARGS | METH_KEYWORDS, "" },
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

