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
#ifdef HAVE_GROMACS
	#include <gromacs/utility/smalloc.h>
	#include <gromacs/fileio/xtcio.h>
#endif

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
#ifdef HAVE_GROMACS
PyObject *read_xtc(const char *filename);
#endif
int write_xyz(FILE *, PyObject *, PyArrayObject *, char *);
int write_gro(FILE *, PyObject *, PyArrayObject *, char *, PyObject *, PyObject *, PyArrayObject *);
//double evaluate_energy(PyArrayObject *,PyObject *, FFType, PyObject *, float *);

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

	enum { GUESS, XYZ, MOLDEN, FRAC, GRO, XTC } type;

	char mode;

	char *filename; /* Used while opening the file and for __repr__ */

	FILE *fd;
#ifdef HAVE_GROMACS
	t_fileio *xd;
#endif

	int nofatoms;
	
	PyObject *symbols; /* list of symbols */
	
	PyObject *atomicnumbers; /* atomic numbers */

} Trajectory;


	
void cspline_calculate_drv2(double y2[], int n, double x[], double y[]);
//double cspline_interpolate_y(double v, int n, double x[], double y[], double y2[]);
double cspline_interpolate_y(double v, PyObject *, PyObject *, PyObject *);
double cspline_interpolate_drv(double v, PyObject *, PyObject *, PyObject *);

PyObject *exposed_read(PyObject *, PyObject *, PyObject *);
PyObject *exposed_write(PyObject *, PyObject *, PyObject *);
PyObject *centerofmass(PyObject *, PyObject *);
PyObject *inertia(PyObject *, PyObject *, PyObject *);
PyObject *mep_distance(PyObject *, PyObject *, PyObject *);
PyObject *find_molecules(PyObject *self, PyObject *args, PyObject *kwds);
PyObject *find_bonds(PyObject *self, PyObject *args, PyObject *kwds);
PyObject *distanceMatrix(PyObject *self, PyObject *args, PyObject *kwds);
PyObject *measureAngleCosine(PyObject *self, PyObject *args, PyObject *kwds);
PyObject *findHBonds(PyObject *self, PyObject *args, PyObject *kwds);

int read_topo_from_xyz(Trajectory *self);
int read_topo_from_molden(Trajectory *self);
int read_topo_from_gro(Trajectory *self);

double *boxArray2double(double box[], PyArrayObject *arr);
void wrapCartesian(double point[3], double box[3]);
void nearestImage(double center[3], double other[3], double half[3]);
double threePointAngleCosine(double A[3], double B[3], double C[3]);
double distanceSquare(double p[3], double q[3]);
double copyPoint(double dst[3], double src[3]);
