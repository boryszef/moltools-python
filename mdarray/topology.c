/***************************************************************************

    mdarray

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


#include "topology.h"
#include "periodic_table.h"
#include "utils.h"


// Check if two groups have at least one atom in common
//
static int groupOverlap(Group a, Group b) {
	int i, j;

	for (i = 0; i < a.len; i++) {
		for (j = 0; j < b.len; j++)
			if (a.idx[i] == b.idx[j]) return 1;
	}
	return 0;
}

// Append atoms from group b to the end of a
//
static void groupMerge(Group *a, Group b) {
	int i, j, k, newLen;
	int duplicate;

	newLen = a->len + b.len;
	a->idx = (int*)realloc(a->idx, newLen * sizeof(int));
	for (i = a->len, j = 0; j < b.len; j++) {
		duplicate = 0;
		// Serch a for duplicates of b[j]
		for (k = 0; k < a->len; k++)
			if ( (a->idx)[k] == b.idx[j] ) {
				duplicate = 1;
				break;
			}
		if (!duplicate)
			(a->idx)[i++] = b.idx[j];
	}
	// Truncate the pre-allocated array (if there were duplicates)
	newLen = i;
	a->idx = (int*)realloc(a->idx, newLen * sizeof(int));
	a->len = newLen;
}

// Clean-up group
//
static void groupDelete(Group *a) {
	free(a->idx);
	a->idx = NULL;
	a->len = 0;
}

// Remove empty groups and recalculate their count
//
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


PyObject *find_bonds(PyObject *self, PyObject *args, PyObject *kwds) {
	extern Element element_table[];
	int i, j, nat, type, idx, start;
	double ax, ay, az, bx, by, bz;
	double ar, br, dist, covsum;
	npy_intp *numpyint;
	float factor = 1.3;
	char *format = NULL;
	enum Formats { FMT_LIST, FMT_DICT } fmt = FMT_LIST;

	static char *kwlist[] = {
		"coordinates", "types", "factor", "format", NULL };

	PyObject *val1, *val2, *tmp_list, *tmptup;
	PyObject *py_symbols;
	PyObject *py_coords;
	PyObject *py_result = NULL;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|fs", kwlist,
			&PyList_Type, &py_symbols,
			&PyArray_Type, &py_coords,
			&factor,
			&format))
		return NULL;

	if(format != NULL && !strcmp(format, "bonds"))
		fmt = FMT_LIST;
	else if (format != NULL && !strcmp(format, "atoms"))
		fmt = FMT_DICT;
	else if (format != NULL) {
		PyErr_SetString(PyExc_RuntimeError, "Unrecognized format.");
		return NULL;
	}

	nat = PyList_Size(py_symbols);
	numpyint = PyArray_DIMS((PyArrayObject*)py_coords);
	if (numpyint[0] != nat) {
		PyErr_SetString(PyExc_RuntimeError, "Size of symbol and coordinate lists does not match.");
		return NULL;
	}

	type = PyArray_TYPE((PyArrayObject*)py_coords);

	py_result = PyList_New(0);

	for (i = 0; i < nat; i++) {

		ax = (double)getFromArray2D(py_coords, type, i, 0);
		ay = (double)getFromArray2D(py_coords, type, i, 1);
		az = (double)getFromArray2D(py_coords, type, i, 2);

		val1 = PyList_GetItem(py_symbols, i); // borrowed
		idx = getElementIndexBySymbol(PyUnicode_AsUTF8(val1));
		if(idx == -1) {
			PyErr_SetString(PyExc_RuntimeError, "Symbol unrecognized.");
			return NULL; }

		ar = element_table[idx].covalent_radius;
		if(ar < 0) {
			PyErr_SetString(PyExc_RuntimeError, "Covalent radius undefined.");
			return NULL; }

		if (fmt == FMT_DICT)
			//tmp_list = PyList_New(0); // new
			tmp_list = PySet_New(NULL); // new
		val1 = PyLong_FromLong(i); // new

		if (fmt == FMT_DICT) start = 0;
		else start = i+1;

		for (j = start; j < nat; j++) {

			if (i == j) continue;

			bx = (double)getFromArray2D(py_coords, type, j, 0);
			by = (double)getFromArray2D(py_coords, type, j, 1);
			bz = (double)getFromArray2D(py_coords, type, j, 2);

			val2 = PyList_GetItem(py_symbols, j); // borrowed
			idx = getElementIndexBySymbol(PyUnicode_AsUTF8(val2));
			if(idx == -1) {
				PyErr_SetString(PyExc_RuntimeError, "Symbol unrecognized.");
				return NULL; }

			br = element_table[idx].covalent_radius;
			if(br < 0) {
				PyErr_SetString(PyExc_RuntimeError, "Covalent radius undefined.");
				return NULL; }

			dist = sq(bx-ax) + sq(by-ay) + sq(bz-az);
			covsum = (ar+br) * factor;

			if (dist < sq(covsum)) {

				val2 = PyLong_FromLong(j); // new

				if (fmt == FMT_DICT)
					//PyList_Append(tmp_list, val2);
					PySet_Add(tmp_list, val2);
				else {
					tmptup = PyTuple_New(2);
					PyTuple_SetItem(tmptup, 0, val1); // steals reference
					Py_INCREF(val1); // decreased after the inner loop
					PyTuple_SetItem(tmptup, 1, val2);
					Py_INCREF(val2); // decreased just below
					PyList_Append(py_result, tmptup);
					Py_DECREF(tmptup);
				}
				Py_DECREF(val2);
			}
		}
		if (fmt == FMT_DICT) {
			PyList_Append(py_result, tmp_list);
			Py_DECREF(tmp_list);
		}
	}

	return py_result;
}


// This function uses data from findBonds to assemble atoms into molecules
//
PyObject *find_molecules(PyObject *self, PyObject *args, PyObject *kwds) {
	int merged, nbonds, i, j, nmols;
	unsigned int natoms;
	long int idx1, idx2;
	int *checkTable;
	Group *groups;

	static char *kwlist[] = {
		"natoms", "bonds", NULL };

	PyObject *bond, *mol, *atomIdx, *num;
	PyObject *py_bonds;
	PyObject *py_result = NULL;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "IO!", kwlist,
			&natoms, &PyList_Type, &py_bonds))
		return NULL;

	nbonds = PyList_Size(py_bonds);
	nmols = nbonds;

	// Makes preliminary single group (molecule) for each bond
	groups = (Group*) malloc(nmols * sizeof(Group));
	for (i = 0; i < nmols; i++) {
		bond = PyList_GetItem(py_bonds, i); // borrowed
		atomIdx = PyTuple_GetItem(bond, 0);
		idx1 = PyLong_AsLong(atomIdx);
		atomIdx = PyTuple_GetItem(bond, 1);
		idx2 = PyLong_AsLong(atomIdx);
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
		
	py_result = PyList_New(0);

	// Loop through the groups until there is nothing to merge
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
	} while(merged);

	checkTable = (int*) malloc(natoms * sizeof(int));
	for (i = 0; i < natoms; i++) checkTable[i] = 0;

	// Make list for each molecule and mark atoms which are part of
	// something bigger
	for (i = 0; i < nmols; i++) {
		//mol = PyList_New(groups[i].len);
		mol = PySet_New(NULL);
		for(j = 0; j < groups[i].len; j++) {
			//PyList_SetItem(mol, j, PyLong_FromLong(groups[i].idx[j]));
			num = PyLong_FromLong(groups[i].idx[j]);
			PySet_Add(mol, num);
			Py_DECREF(num);
			checkTable[groups[i].idx[j]] = 1;
		}
		PyList_Append(py_result, mol);
		Py_DECREF(mol);
	}

	// For each unassigned atom, create a list with this single member
	for (i = 0; i < natoms; i++) {
		if (!checkTable[i]) {
			//mol = PyList_New(1);
			mol = PySet_New(NULL);
			//PyList_SetItem(mol, 0, PyLong_FromLong(i));
			num = PyLong_FromLong(i);
			PySet_Add(mol, num);
			Py_DECREF(num);
			PyList_Append(py_result, mol);
			Py_DECREF(mol);
		}
	}

	return py_result;
}

