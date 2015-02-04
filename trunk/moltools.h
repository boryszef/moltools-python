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
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include "structmember.h"
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
		float covalent_radius;
	} Element;

typedef enum ff_type__ {
		FF_TWELVESIX = 0,
		FF_NONE,
	} FFType;

char *readline(FILE *);
int make_lowercase(char *);
int stripline(char *);
PyObject *read_xyz(FILE *fd, float factor);
PyObject *read_molden(FILE *fd);
PyObject *read_fractional(FILE *fd);
PyObject *read_gro(FILE *fd);
int write_xyz(FILE *, PyObject *, PyArrayObject *, char *);
int write_gro(FILE *, PyObject *, PyArrayObject *, char *, PyObject *, PyObject *, PyArrayObject *);
double evaluate_energy(PyArrayObject *,PyObject *, FFType, PyObject *, float *);

typedef struct {
	PyObject_HEAD
	PyObject *rscale; /* values of r used by RHO and Z */
	PyObject *rhoscale; /* values of rho used by F */
	PyObject *rho; /* values of rho */
	PyObject *rho2; /* values of the second derivative of rho */
	PyObject *Z; /* values of Z */
	PyObject *Z2; /* values of the second derivative of Z */
	PyObject *F; /* values of F */
	PyObject *F2; /* values of the second derivative of F */
} EAMff;

typedef struct {

	PyObject_HEAD
	
	PyObject *frames; /* 3D array of frames */
	float *frames_raw;
	npy_intp frames_dim[3];
	
	PyObject *symbols; /* list of symbols */
	
	PyObject *comment; /* comment */
	
	PyObject *charges; /* comment */
	float *charges_raw;
	npy_intp charges_dim[2];
	
	PyObject *energies; /* energies */
	float *energies_raw;
	npy_intp energies_dim[1];
	
	PyObject *atomicnumbers; /* atomic numbers */
	
	int natoms; /* number of atoms */
	
	int nframes; /* number of frames */

} Molecule;

void cspline_calculate_drv2(double y2[], int n, double x[], double y[]);
//double cspline_interpolate_y(double v, int n, double x[], double y[], double y2[]);
double cspline_interpolate_y(double v, PyObject *, PyObject *, PyObject *);
double cspline_interpolate_drv(double v, PyObject *, PyObject *, PyObject *);

PyObject *exposed_read(PyObject *, PyObject *, PyObject *);
PyObject *exposed_write(PyObject *, PyObject *, PyObject *);
PyObject *centerofmass(PyObject *, PyObject *);
PyObject *inertia(PyObject *, PyObject *, PyObject *);
PyObject *mep_distance(PyObject *, PyObject *, PyObject *);

int read_topo_from_xyz(FILE *fd, Molecule *self);
int read_topo_from_molden(FILE *fd, Molecule *self);
int read_topo_from_molden(FILE *fd, Molecule *self);
int read_frame_from_xyz(FILE *fd, float factor, Molecule *);
