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

#ifndef __TRAJECTORY_H__
#define __TRAJECTORY_H__

typedef struct {

	PyObject_HEAD

	enum { GUESS, XYZ, MOLDEN, GRO, XTC } type;
	enum { ANGS, BOHR, NM } units;
	char mode;
	char *fileName; /* Used while opening the file and for __repr__ */
	FILE *fd;
	/* Used for keeping track of the position in the file while reading     *
	 * frames. Two variables are needed, because some formats, like Molden, *
	 * store geometries and energies in different parts of the file.        */
#ifdef HAVE_GROMACS
	t_fileio *xd;
	rvec *xtcCoord;
#endif
	long filePosition1;
	long filePosition2;
	enum { MLGEOMETRY, MLATOMS, MLUNKNOWN } moldenStyle;
	int nOfAtoms;
	int nOfFrames;
	int lastFrame;
	PyObject *symbols; /* list of symbols */
	PyObject *atomicNumbers; /* atomic numbers */
	PyObject *resIDs; /* residue numbers */
	PyObject *resNames; /* residue names */

	PyObject *moldenSections; /* Sections in Molden file and offsets */

} Trajectory;

static int read_topo_from_xyz(Trajectory *self);
static int read_topo_from_molden(Trajectory *self);
static int read_topo_from_gro(Trajectory *self);
static PyObject *read_frame_from_xyz(Trajectory *self);
static PyObject *read_frame_from_molden_atoms(Trajectory *self);
static PyObject *read_frame_from_molden_geometries(Trajectory *self);
static PyObject *read_frame_from_gro(Trajectory *self);
#ifdef HAVE_GROMACS
static PyObject *read_frame_from_xtc(Trajectory *self);
#endif

#endif /* __TRAJECTORY_H__ */
