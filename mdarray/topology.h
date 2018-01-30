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


#ifndef __TOPOLOGY_H__
#define __TOPOLOGY_H__

/* This should come before other Numpy-related declarations in every *
 * file that does not define the module's init function              */
#define NO_IMPORT_ARRAY

/* Make sure the general declarations are made first */
#include "mdarray.h"

typedef struct {
	int len;
	int *idx;
} Group;

PyObject *find_molecules(PyObject *self, PyObject *args, PyObject *kwds);
PyObject *find_bonds(PyObject *self, PyObject *args, PyObject *kwds);

#endif /* __TOPOLOGY_H__ */
