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

#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <locale.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define BOHR 0.529177209

#define sq(a) ((a) * (a))

typedef const struct element__ {
		int number;
		double mass;
		const char *symbol;
		const char *name;
	} Element;

char *readline(FILE *);
int make_lowercase(char *);
int stripline(char *);
PyObject *read_xyz(FILE *fd, float factor);
PyObject *read_molden(FILE *fd);
PyObject *read_fractional(FILE *fd);
int write_xyz(FILE *, PyObject *, PyObject *, char *);
int write_gro(FILE *, PyObject *, PyObject *, char *, PyObject *, PyObject *, PyObject *);

