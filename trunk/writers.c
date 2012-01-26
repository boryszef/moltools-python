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


int write_xyz(FILE *fd, PyObject *py_symbols, PyObject *py_coords, char *comment) {

	int nat, i;
	float x, y, z;
	char *s;

	nat = PyList_Size(py_symbols);
    fprintf(fd, "%d\n", nat);
	if( comment != NULL )
        fprintf(fd, "%s\n", comment);
    else
        fprintf(fd, "\n");

    for ( i = 0; i < nat; i++ ) {
        x = *( (float*) PyArray_GETPTR2(py_coords, i, 0) );
        y = *( (float*) PyArray_GETPTR2(py_coords, i, 1) );
        z = *( (float*) PyArray_GETPTR2(py_coords, i, 2) );
        s = PyString_AsString(PyList_GetItem(py_symbols, i));
        fprintf(fd, "%-3s  %12.8f  %12.8f  %12.8f\n", s, x, y, z);
    }

    return nat;
}

int write_gro(FILE *fd, PyObject *py_symbols, PyObject *py_coords, char *comment,
              PyObject *py_resnam, PyObject *py_resid, PyObject *py_box) {

	int nat, i, type;
	long int resid;
	float x, y, z;
	char *s, *resnam;
	PyObject *val;

	if( comment != NULL )
        fprintf(fd, "%s\n", comment);
    else
        fprintf(fd, "\n");
	nat = PyList_Size(py_symbols);
    fprintf(fd, "%5d\n", nat);

    for ( i = 0; i < nat; i++ ) {
		resid = PyInt_AsLong(PyList_GetItem(py_resid, i));
		resnam = PyString_AsString(PyList_GetItem(py_resnam, i));
        x = *( (float*) PyArray_GETPTR2(py_coords, i, 0) ) / 10.0;
        y = *( (float*) PyArray_GETPTR2(py_coords, i, 1) ) / 10.0;
        z = *( (float*) PyArray_GETPTR2(py_coords, i, 2) ) / 10.0;
        s = PyString_AsString(PyList_GetItem(py_symbols, i));
        fprintf(fd, "%5d%-5s%5s%5ld%8.3f%8.3f%8.3f\n", i+1, resnam, s, resid, x, y, z);
    }

	/* Do some testing on the array */
	if ( PyArray_NDIM(py_box) != 1 ) {
		PyErr_SetString(PyExc_ValueError, "Unsupported box shape");
		return -1; }
	else {
		type = PyArray_TYPE(py_box);
		switch(type) {
			case NPY_FLOAT:
				x = *( (float*) PyArray_GETPTR1(py_box, 0) );
				y = *( (float*) PyArray_GETPTR1(py_box, 1) );
				z = *( (float*) PyArray_GETPTR1(py_box, 2) );
				break;
			case NPY_DOUBLE:
				x = *( (double*) PyArray_GETPTR1(py_box, 0) );
				y = *( (double*) PyArray_GETPTR1(py_box, 1) );
				z = *( (double*) PyArray_GETPTR1(py_box, 2) );
				break;
			default:
				PyErr_SetString(PyExc_ValueError, "Incorrect type in box vector");
				return -1;
		}
		x /= 10.0;
		y /= 10.0;
		z /= 10.0;
		fprintf(fd, "%10.5f%10.5f%10.5f\n", x, y, z);
	}

	return nat;
}
