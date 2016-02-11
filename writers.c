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
#define NO_IMPORT_ARRAY
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "writers.h"
#include "moltools.h"


PyObject *exposed_write(PyObject *self, PyObject *args, PyObject *kwds) {

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

	PyObject *py_symbols, *val;
	PyArrayObject *py_coords, *py_box = NULL;
	PyObject *py_resnam = NULL, *py_resid = NULL;

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

	Py_RETURN_NONE;
}



int write_xyz(FILE *fd, PyObject *py_symbols, PyArrayObject *py_coords, char *comment) {

	int nat, i, type;
	double x, y, z;
	char *s;

	nat = PyList_Size(py_symbols);
    fprintf(fd, "%d\n", nat);
	if( comment != NULL )
        fprintf(fd, "%s\n", comment);
    else
        fprintf(fd, "\n");

	type = PyArray_TYPE(py_coords);
    for ( i = 0; i < nat; i++ ) {
		switch(type) {
			case NPY_FLOAT:
		        x = *( (float*) PyArray_GETPTR2(py_coords, i, 0) );
    		    y = *( (float*) PyArray_GETPTR2(py_coords, i, 1) );
        		z = *( (float*) PyArray_GETPTR2(py_coords, i, 2) );
				break;
			case NPY_DOUBLE:
		        x = *( (double*) PyArray_GETPTR2(py_coords, i, 0) );
    		    y = *( (double*) PyArray_GETPTR2(py_coords, i, 1) );
        		z = *( (double*) PyArray_GETPTR2(py_coords, i, 2) );
				break;
			default:
				PyErr_SetString(PyExc_ValueError, "Incorrect type in box vector");
				return -1;
		}
        s = PyString_AsString(PyList_GetItem(py_symbols, i));
        fprintf(fd, "%-3s  %12.8lf  %12.8lf  %12.8lf\n", s, x, y, z);
    }

    return nat;
}

int write_gro(FILE *fd, PyObject *py_symbols, PyArrayObject *py_coords, char *comment,
              PyObject *py_resnam, PyObject *py_resid, PyArrayObject *py_box) {

	int nat, i, type;
	long int resid;
	float x, y, z;
	char *s, *resnam;

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
