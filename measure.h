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


#ifndef __MEASURE_H__
#define __MEASURE_H__


PyObject *findHBonds(PyObject *self, PyObject *args, PyObject *kwds);
PyObject *measureAngleCosine(PyObject *self, PyObject *args, PyObject *kwds);
PyObject *distanceMatrix(PyObject *self, PyObject *args, PyObject *kwds);
PyObject *centerofmass(PyObject *, PyObject *);
PyObject *inertia(PyObject *, PyObject *, PyObject *);
PyObject *mep_distance(PyObject *, PyObject *, PyObject *);

#endif /* __MEASURE_H__ */
