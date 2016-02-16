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


#define PY_ARRAY_UNIQUE_SYMBOL MOLTOOLS
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "moltools.h"
#include "periodic_table.h"
//#include "writers.h"
//#include "readers.h"
#include "measure.h"
#include "constants.h"
#include "topology.h"



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


/*static PyObject *rot_inertia(PyObject *self, PyObject *args, PyObject *kwds) {
	return NULL;
}*/


static PyMethodDef moltoolsMethods[] = {

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
		"       __________________________\n"
		"MEP = √ Σ_i (m_i * (δr_i)^2) / M \n"
		"\n"
		"When masses are skipped, the function simply assumes m = 1, M = N.\n"
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
	/*{"rot_inertia", (PyCFunction)rot_inertia, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"coordinates_rot, principal_moments = rot_inertia(coordinates, inertia_tensor)\n"
		"\n"
		"Rotates the coordinates so that the principal moments of inertia\n"
		"are aligned with axes of the coordinate system.\n"
		"\n" },*/
	{"quatfit", (PyCFunction)quatfit, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"fitted = quatfit(reference, fit, masses=None)\n"
		"\n"
		"Performs a quaternion fit of 'fit' ndarray (Nx3) to 'reference'\n"
		"ndarray (Nx3). If masses (list) are given, the fit is mass weighted.\n"
		"Returns an (Nx3) ndarray with fitted coordinates. Remember to translate\n"
		"both molecules to the origin of the coordinate system.\n"
		"\n" },

    {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC initmoltools(void)
{
    PyObject *md;
	PyObject *exposed_atom_symbols, *exposed_atom_names;
	PyObject *exposed_atom_masses, *exposed_symbol2number;
	PyObject *exposed_covalentradii;
	//extern PyTypeObject EAMffType;
	extern PyTypeObject TrajectoryType;
	//extern PyTypeObject FrameType;

	/* Use system-wide locale, but make sure that decimal point is a point! */
	setlocale(LC_ALL, "");
	setlocale(LC_NUMERIC, "C");

	//if (PyType_Ready(&EAMffType) < 0)
	//	return;
	if (PyType_Ready(&TrajectoryType) < 0)
		return;
	//if (PyType_Ready(&FrameType) < 0)
	//	return;

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

	//Py_INCREF(&EAMffType);
	//PyModule_AddObject(md, "EAMff", (PyObject *)&EAMffType);

	Py_INCREF(&TrajectoryType);
	PyModule_AddObject(md, "Trajectory", (PyObject *)&TrajectoryType);

	//Py_INCREF(&FrameType);
	//PyModule_AddObject(md, "Frame", (PyObject *)&FrameType);
}

