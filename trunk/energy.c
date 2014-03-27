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

double evaluate_energy(PyArrayObject *py_coords, PyObject *py_types, FFType ff_type, PyObject *py_ff, float *box) {
	double energy = 0.0;
	double energy_nb = 0.0;
	double d2, d6;
	float ax, ay, az, bx, by, bz;
	double eps, sigma;
	int i, j, natoms;
	npy_intp *numpyint;
	PyObject *tmp;

	numpyint = PyArray_DIMS(py_coords);
	natoms = (int) numpyint[0];
	switch(ff_type) {
		case FF_TWELVESIX:
			for (i = 0; i < natoms; i++) {
				ax = *( (float*) PyArray_GETPTR2(py_coords, i, 0) );
				ay = *( (float*) PyArray_GETPTR2(py_coords, i, 1) );
				az = *( (float*) PyArray_GETPTR2(py_coords, i, 2) );
				for (j = i+1; j < natoms; j++) {
					bx = *( (float*) PyArray_GETPTR2(py_coords, j, 0) );
					by = *( (float*) PyArray_GETPTR2(py_coords, j, 1) );
					bz = *( (float*) PyArray_GETPTR2(py_coords, j, 2) );
					tmp = PyList_GetItem(py_ff, i*natoms+j);
					sigma = PyFloat_AsDouble(PyList_GetItem(tmp, 0));
					eps = PyFloat_AsDouble(PyList_GetItem(tmp, 1));
					d2 = sq(ax-bx) + sq(ay-by) + sq(az-bz);
					d6 = pow(d2,3);
					energy_nb += eps * pow(sigma,6) * (pow(sigma/d2,6) - 1.0/d6);
				}
			}
			break;
		default:
			return 0.0;
	}

	energy = energy_nb;
	return energy;
}
