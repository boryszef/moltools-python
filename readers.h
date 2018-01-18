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

#ifdef HAVE_GROMACS
	#include <gromacs/utility/smalloc.h>
	#include <gromacs/fileio/xtcio.h>
#endif

static int read_topo_from_xyz(Trajectory *self);
PyObject* read_molden_sections(FILE *fd);
static int read_topo_from_molden(Trajectory *self);
static int read_topo_from_gro(Trajectory *self);
static PyObject *read_frame_from_xyz(Trajectory *self);
static PyObject *read_frame_from_molden_atoms(Trajectory *self);
static PyObject *read_frame_from_molden_geometries(Trajectory *self);
static PyObject *read_frame_from_gro(Trajectory *self);
static PyObject *read_frame_from_xtc(Trajectory *self);


#endif /* __READERS_H__ */
