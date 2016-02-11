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


#ifndef __WRITERS_H__
#define __WRITERS_H__


PyObject *exposed_write(PyObject *, PyObject *, PyObject *);
int write_xyz(FILE *, PyObject *, PyArrayObject *, char *);
int write_gro(FILE *, PyObject *, PyArrayObject *, char *, PyObject *, PyObject *, PyArrayObject *);

#endif /* __WRITERS_H__ */
