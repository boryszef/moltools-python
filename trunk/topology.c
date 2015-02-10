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


typedef struct {
	int len;
	int *idx;
} Group;


static int groupOverlap(Group a, Group b) {
	int i, j;

	for (i = 0; i < a.len; i++) {
		for (j = 0; j < b.len; j++)
			if (a.idx[i] == b.idx[j]) return 1;
	}
	return 0;
}

static void groupMerge(Group *a, Group b) {
	int i, j, k, newLen;
	int duplicate;

	newLen = a->len + b.len;
	a->idx = (int*)realloc(a->idx, newLen * sizeof(int));
	for (i = a->len, j = 0; j < b.len; j++) {
		duplicate = 0;
		for (k = 0; k < a->len; k++)
			if ( (a->idx)[k] == b.idx[j] ) {
				duplicate = 1;
				break;
			}
		if (!duplicate)
			(a->idx)[i++] = b.idx[j];
	}
	newLen = i;
	a->idx = (int*)realloc(a->idx, newLen * sizeof(int));
	a->len = newLen;
}

static void groupDelete(Group *a) {
	free(a->idx);
	a->idx = NULL;
	a->len = 0;
}

static void groupPurge(int *n, Group *g) {
	int i, j;

	for (i = 0; i < *n;) {
		if(g[i].len == 0) {
			for (j = i; j < *n-1; j++)
				g[j] = g[j+1];
			*n -= 1;
		} else i++;
	}
}

/*static void groupPrint(int n, Group *g) {
	int i, j;

	printf("%d groups:\n", n);
	for (i = 0; i < n; i++) {
		printf("%3d: [", i);
		for (j = 0; j < g[i].len; j++)
			printf(" %3d,", g[i].idx[j]);
		printf("]\n");
	}
}*/

PyObject *find_bonds(PyObject *self, PyObject *args, PyObject *kwds) {
	extern Element element_table[];
	int i, j, nat, type, idx, start;
	float ax, ay, az, bx, by, bz;
	double ar, br, dist;
	npy_intp *numpyint;
	float factor = 1.3;
	char *format = NULL;
	enum Formats { FMT_LIST, FMT_DICT } fmt = FMT_LIST;

	static char *kwlist[] = {
		"coordinates", "types", "factor", "format", NULL };

	PyObject *val1, *val2, *tmp_list, *tmptup;
	PyObject *py_symbols;
	PyArrayObject *py_coords;
	PyObject *py_result = NULL;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|fs", kwlist,
			&PyList_Type, &py_symbols,
			&PyArray_Type, &py_coords,
			&factor,
			&format))
		return NULL;

	if(format != NULL && !strcmp(format, "list"))
		fmt = FMT_LIST;
	else if (format != NULL && !strcmp(format, "dict"))
		fmt = FMT_DICT;
	else if (format != NULL) {
		PyErr_SetString(PyExc_RuntimeError, "Unrecognized format.");
		return NULL;
	}

	nat = PyList_Size(py_symbols);
	numpyint = PyArray_DIMS(py_coords);
	if (numpyint[0] != nat) {
		PyErr_SetString(PyExc_RuntimeError, "Size of symbol and coordinate lists does not match.");
		return NULL;
	}

	type = PyArray_TYPE(py_coords);
	if (type != NPY_FLOAT && type != NPY_DOUBLE) {
		PyErr_SetString(PyExc_ValueError, "Coordinates must be of FLOAT or DOUBLE type.");
		return NULL;
	}

	if (fmt == FMT_DICT) {
		py_result = PyDict_New();
	} else {
		py_result = PyList_New(0);
	}

	for (i = 0; i < nat; i++) {
		if (type == NPY_FLOAT) {
			ax = *( (float*) PyArray_GETPTR2(py_coords, i, 0) );
			ay = *( (float*) PyArray_GETPTR2(py_coords, i, 1) );
			az = *( (float*) PyArray_GETPTR2(py_coords, i, 2) );
		} else {
			ax = *( (double*) PyArray_GETPTR2(py_coords, i, 0) );
			ay = *( (double*) PyArray_GETPTR2(py_coords, i, 1) );
			az = *( (double*) PyArray_GETPTR2(py_coords, i, 2) );
		}
		val1 = PyList_GetItem(py_symbols, i); // borrowed
		//val2 = PyDict_GetItem(py_types, val1); // borrowed
		//ar = PyFloat_AsDouble(val2);
		idx = getElementIndexBySymbol(PyString_AsString(val1));
		if(element_table[idx].number == -1) {
			PyErr_SetString(PyExc_RuntimeError, "Symbol unrecognized.");
			return NULL; }
		ar = element_table[idx].covalent_radius;
		if(ar < 0) {
			PyErr_SetString(PyExc_RuntimeError, "Covalent radius undefined.");
			return NULL; }
		tmp_list = PyList_New(0); // new
		val1 = PyInt_FromLong(i); // new
		if (fmt == FMT_DICT) start = 0;
		else start = i+1;
		for (j = start; j < nat; j++) {
			if (i == j) continue;
			if (type == NPY_FLOAT) {
				bx = *( (float*) PyArray_GETPTR2(py_coords, j, 0) );
				by = *( (float*) PyArray_GETPTR2(py_coords, j, 1) );
				bz = *( (float*) PyArray_GETPTR2(py_coords, j, 2) );
			} else {
				bx = *( (double*) PyArray_GETPTR2(py_coords, j, 0) );
				by = *( (double*) PyArray_GETPTR2(py_coords, j, 1) );
				bz = *( (double*) PyArray_GETPTR2(py_coords, j, 2) );
			}
			val2 = PyList_GetItem(py_symbols, j); // borrowed
			idx = getElementIndexBySymbol(PyString_AsString(val2));
			if(element_table[idx].number == -1) {
				PyErr_SetString(PyExc_RuntimeError, "Symbol unrecognized.");
				return NULL; }
			br = element_table[idx].covalent_radius;
			if(br < 0) {
				PyErr_SetString(PyExc_RuntimeError, "Covalent radius undefined.");
				return NULL; }
			dist = sq(bx-ax) + sq(by-ay) + sq(bz-az);
			//if (dist < sq((ar+br) * factor)) {
			if (sqrt(dist) < (ar+br) * factor) {
				val2 = PyInt_FromLong(j); // new
				if (fmt == FMT_DICT)
					PyList_Append(tmp_list, val2);
				else {
					tmptup = PyTuple_New(2);
					PyTuple_SetItem(tmptup, 0, val1);
					PyTuple_SetItem(tmptup, 1, val2);
					PyList_Append(py_result, tmptup);
					Py_DECREF(tmptup);
				}
				Py_DECREF(val2);
			}
		}
		if (fmt == FMT_DICT)
			if (PyList_Size(tmp_list))
				PyDict_SetItem(py_result, val1, tmp_list);
		Py_DECREF(val1);
		Py_DECREF(tmp_list);
	}

	return py_result;
}


PyObject *find_molecules(PyObject *self, PyObject *args, PyObject *kwds) {
	int merged, nbonds, i, j, nmols, indexFound, natoms;
	unsigned int pyNAtoms;
	long int idx1, idx2, idx3;
	Group *groups;

	static char *kwlist[] = {
		"natoms", "bonds", NULL };

	PyObject *bond, *mol, *atomIdx;
	PyObject *py_bonds;
	PyObject *py_result = NULL;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "IO!", kwlist,
			&pyNAtoms, &PyList_Type, &py_bonds))
		return NULL;

	nbonds = PyList_Size(py_bonds);
	nmols = nbonds;

	groups = (Group*) malloc(nmols * sizeof(Group));
	for (i = 0; i < nmols; i++) {
		bond = PyList_GetItem(py_bonds, i); // borrowed
		atomIdx = PyTuple_GetItem(bond, 0);
		idx1 = PyInt_AsLong(atomIdx);
		atomIdx = PyTuple_GetItem(bond, 1);
		idx2 = PyInt_AsLong(atomIdx);
		if (!PyTuple_Check(bond)) {
			PyErr_SetString(PyExc_RuntimeError, "List of bonds should contain tuples.");
			return NULL; }
		if (PyTuple_Size(bond) != 2) {
			PyErr_SetString(PyExc_RuntimeError, "Bond tuple must have exactly two indices.");
			return NULL; }
		groups[i].len = 2;
		groups[i].idx = (int*) malloc(2 * sizeof(int));
		groups[i].idx[0] = idx1;
		groups[i].idx[1] = idx2;
	}
	//groupPrint(nmols, groups);
		
	py_result = PyList_New(0);

	do {
		merged = 0;
		for (i = 0; i < nmols; i++) {
			for (j = i+1; j < nmols; j++) {
				if (groupOverlap(groups[i], groups[j])) {
					groupMerge(groups+i, groups[j]);
					groupDelete(groups+j);
					merged = 1;
				}
			}
		}
		groupPurge(&nmols, groups);
		//groupPrint(nmols, groups);
	} while(merged);
		
	for (i = 0; i < nmols; i++) {
		mol = PyList_New(groups[i].len);
		for(j = 0; j < groups[i].len; j++)
			PyList_SetItem(mol, j, PyInt_FromLong(groups[i].idx[j]));
		PyList_Append(py_result, mol);
		Py_DECREF(mol);
	}

	return py_result;
}

