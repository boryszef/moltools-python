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


/* This header contains basic declarations, should come first */
#include "mdarray.h"




static PyMethodDef mdarrayMethods[] = {
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


static struct PyModuleDef mdarrayModule = {
    PyModuleDef_HEAD_INIT,
    "mdarray",   /* name of module */
    /* module documentation, may be NULL */
    "The mdarray module provides the Trajectory class and functions related "
	 "to molecular modelling. The idea is to facilitate writing scripts for "
	 "processing molecular data, using standard types. For atomic coordinates, "
	 "numpy arrays are used, since they are fast and implement linear algebra.",
    -1,           /* size of per-interpreter state of the module,
                     or -1 if the module keeps state in global variables. */
    mdarrayMethods,
	 NULL, NULL, NULL, NULL
};


PyMODINIT_FUNC PyInit_mdarray(void)
{
	PyObject *md;
	extern PyTypeObject TrajectoryType;

	/* Use system-wide locale, but make sure that decimal point is a point! */
	setlocale(LC_ALL, "");
	setlocale(LC_NUMERIC, "C");

	if (PyType_Ready(&TrajectoryType) < 0)
		return NULL;

	md = PyModule_Create(&mdarrayModule);
	if (md == NULL) return NULL;

	Py_INCREF(&TrajectoryType);
	PyModule_AddObject(md, "Trajectory", (PyObject *)&TrajectoryType);

	import_array();

	return md;
}

