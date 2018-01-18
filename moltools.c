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




static PyMethodDef moltoolsMethods[] = {

    {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC initmoltools(void)
{
    extern PyTypeObject TrajectoryType;


    /* Use system-wide locale, but make sure that decimal point is a point! */
    setlocale(LC_ALL, "");
    setlocale(LC_NUMERIC, "C");

    if (PyType_Ready(&TrajectoryType) < 0)
        return;

    md = Py_InitModule3("moltools", moltoolsMethods,
         "The moltools module provides some classes and functions related to molecular "
         "modelling. The idea is to facilitate writing scripts for processing molecular "
         "data, using standard types. For atomic coordinates, numpy arrays are used, "
         "since they are fast and implement linear algebra.");
    if(md == NULL) return;

    import_array();

    Py_INCREF(&TrajectoryType);
    PyModule_AddObject(md, "Trajectory", (PyObject *)&TrajectoryType);
}

