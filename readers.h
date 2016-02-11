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


#ifndef __READERS_H__
#define __READERS_H__

PyObject *read_xyz(FILE *fd, float factor);
PyObject *read_molden(FILE *fd);
PyObject *read_fractional(FILE *fd);
PyObject *read_gro(FILE *fd);
#ifdef HAVE_GROMACS
PyObject *read_xtc(const char *filename);
#endif

PyObject *exposed_read(PyObject *, PyObject *, PyObject *);

#endif /* __READERS_H__ */
