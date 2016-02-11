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

#include "eam.h"


/* Methods of the EAMff object */

static void EAMff_dealloc(EAMff* self)
{
	Py_XDECREF(self->rscale);
	Py_XDECREF(self->rhoscale);
	Py_XDECREF(self->rho);
	Py_XDECREF(self->rho2);
	Py_XDECREF(self->Z);
	Py_XDECREF(self->Z2);
	Py_XDECREF(self->F);
	Py_XDECREF(self->F2);
	self->ob_type->tp_free((PyObject*)self);
}

static PyObject *EAMff_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	EAMff *self;
	npy_intp dims[1];

	dims[0] = 0;
	self = (EAMff *)type->tp_alloc(type, 0);
	if (self != NULL) {
	    self->rscale = PyArray_EMPTY(1, dims, NPY_DOUBLE, 0);
	    if (self->rscale == NULL) {
	        Py_DECREF(self);
	        return NULL; }
	    
	    self->rhoscale = PyArray_EMPTY(1, dims, NPY_DOUBLE, 0);
	    if (self->rhoscale == NULL) {
	        Py_DECREF(self);
	        return NULL; }
	    
	    self->rho = PyArray_EMPTY(1, dims, NPY_DOUBLE, 0);
	    if (self->rho == NULL) {
	        Py_DECREF(self);
	        return NULL; }
	    
	    self->rho2 = PyArray_EMPTY(1, dims, NPY_DOUBLE, 0);
	    if (self->rho2 == NULL) {
	        Py_DECREF(self);
	        return NULL; }
	    
	    self->Z = PyArray_EMPTY(1, dims, NPY_DOUBLE, 0);
	    if (self->Z == NULL) {
	        Py_DECREF(self);
	        return NULL; }
	    
	    self->Z2 = PyArray_EMPTY(1, dims, NPY_DOUBLE, 0);
	    if (self->Z2 == NULL) {
	        Py_DECREF(self);
	        return NULL; }
	    
	    self->F = PyArray_EMPTY(1, dims, NPY_DOUBLE, 0);
	    if (self->F == NULL) {
	        Py_DECREF(self);
	        return NULL; }
	    
	    self->F2 = PyArray_EMPTY(1, dims, NPY_DOUBLE, 0);
	    if (self->F2 == NULL) {
	        Py_DECREF(self);
	        return NULL; }
	    
    }

    return (PyObject *)self;
}

static int EAMff_init(EAMff *self, PyObject *args, PyObject *kwds) {

	FILE *eamfile;
	char *text;
	int nr, nd, i;
	double dr, dd;
	npy_intp dims[1];
	double *rscale, *rhoscale, *rho, *rho2, *Z, *Z2, *F, *F2;

	PyObject *py_filename;

	static char *kwlist[] = {
		"filename", NULL };

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist, &PyString_Type, &py_filename))
		return -1;

	eamfile = fopen(PyString_AsString(py_filename), "r");
	text = readline(eamfile);
	free(text);
	text = readline(eamfile);
	free(text);
	fscanf(eamfile, "%d %lf %d %lf %*f", &nd, &dd, &nr, &dr);

	rscale = malloc(nr * sizeof(double));
	rho = malloc(nr * sizeof(double));
	rho2 = malloc(nr * sizeof(double));
	Z = malloc(nr * sizeof(double));
	Z2 = malloc(nr * sizeof(double));
	if(rscale == NULL || rho == NULL || rho2 == NULL || Z == NULL || Z2 == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return -1; }

	rhoscale = malloc(nd * sizeof(double));
	F = malloc(nd * sizeof(double));
	F2 = malloc(nd * sizeof(double));
	if(rhoscale == NULL || F == NULL || F2 == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return -1; }

	for (i = 0; i < nr; i++)
		rscale[i] = i*dr;
	for (i = 0; i < nd; i++)
		rhoscale[i] = i*dd;
	for (i = 0; i < nd; i++)
		fscanf(eamfile, "%lf", F+i);
	for (i = 0; i < nr; i++)
		fscanf(eamfile, "%lf", Z+i);
	for (i = 0; i < nr; i++)
		fscanf(eamfile, "%lf", rho+i);

	fclose(eamfile);

	dims[0] = nr;
	self->rscale = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, rscale);
	self->rho = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, rho);
	self->Z = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, Z);

	cspline_calculate_drv2(Z2, -nr, rscale, Z);
	cspline_calculate_drv2(rho2, -nr, rscale, rho);
	self->Z2 = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, Z2);
	self->rho2 = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, rho2);

	dims[0] = nd;
	self->rhoscale = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, rhoscale);
	self->F = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, F);

	cspline_calculate_drv2(F2, nr, rhoscale, F);
	self->F2 = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, F2);

	return 0;
}

static PyMemberDef EAMff_members[] = {
    {"rscale", T_OBJECT_EX, offsetof(EAMff, rscale), 0,
     "Values of r used by rho, rho2, Z and Z2."},
    {"rhoscale", T_OBJECT_EX, offsetof(EAMff, rhoscale), 0,
     "Values of rho used by F and F2."},
    {"rho", T_OBJECT_EX, offsetof(EAMff, rho), 0,
     "Values of electron density."},
    {"rho2", T_OBJECT_EX, offsetof(EAMff, rho2), 0,
     "Second derivative of electron density."},
    {"Z", T_OBJECT_EX, offsetof(EAMff, Z), 0,
     "Values of Z."},
    {"Z2", T_OBJECT_EX, offsetof(EAMff, Z2), 0,
     "Second derivative of Z."},
    {"F", T_OBJECT_EX, offsetof(EAMff, F), 0,
     "Values of Z."},
    {"F2", T_OBJECT_EX, offsetof(EAMff, F2), 0,
     "Second derivative of Z."},
    {NULL}  /* Sentinel */
};

static PyObject *EAMff_energy(EAMff *self, PyObject *args) {
	int i, j, p, q, nat, ctype;
	double d, z, f, phi, E, rho;
	double ax, ay, az, bx, by, bz, dist;
	npy_intp *dims;

	PyArrayObject *py_coords;
	PyObject *py_symbols;
	PyObject *py_result = NULL;

	if(!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &py_coords, &PyList_Type, &py_symbols))
		return NULL;

	dims = PyArray_DIMS(py_coords);
	nat = PyList_Size(py_symbols);
	if (nat != dims[0]) {
		PyErr_SetString(PyExc_RuntimeError, "Number of coordinates and symbols does not match");
		return NULL; }
	ctype = PyArray_TYPE(py_coords);
	E = 0.0;
	for (i = 0; i < nat; i++) {
		rho = 0.0;
		if (ctype == NPY_FLOAT) {
			ax = *( (float*) PyArray_GETPTR2(py_coords, i, 0) );
			ay = *( (float*) PyArray_GETPTR2(py_coords, i, 1) );
			az = *( (float*) PyArray_GETPTR2(py_coords, i, 2) );
		} else {
			ax = *( (double*) PyArray_GETPTR2(py_coords, i, 0) );
			ay = *( (double*) PyArray_GETPTR2(py_coords, i, 1) );
			az = *( (double*) PyArray_GETPTR2(py_coords, i, 2) );
		}
		for (j = 0; j < nat; j++) {
			if (i == j) continue;
			if (ctype == NPY_FLOAT) {
				bx = *( (float*) PyArray_GETPTR2(py_coords, j, 0) );
				by = *( (float*) PyArray_GETPTR2(py_coords, j, 1) );
				bz = *( (float*) PyArray_GETPTR2(py_coords, j, 2) );
			} else {
				bx = *( (double*) PyArray_GETPTR2(py_coords, j, 0) );
				by = *( (double*) PyArray_GETPTR2(py_coords, j, 1) );
				bz = *( (double*) PyArray_GETPTR2(py_coords, j, 2) );
			}
			d = sqrt(sq(bx-ax) + sq(by-ay) + sq(bz-az));
			rho += cspline_interpolate_y(d, self->rscale, self->rho, self->rho2);
			if ( j < i ) continue;
			z = cspline_interpolate_y(d, self->rscale, self->Z, self->Z2);
			phi = 27.2 * 0.529 * z*z / d;
			E += phi;
		}
		f = cspline_interpolate_y(rho, self->rhoscale, self->F, self->F2);
		E += f;
	}

	return PyFloat_FromDouble(E);
}


static PyObject *EAMff_gradient(EAMff *self, PyObject *args) {
	int i, j, k, p, q, nat, ctype;
	double ix, iy, iz, jx, jy, jz, kx, ky, kz, dist;
	double dij2, dik2, dkj2, dij, dik, dkj;
	double rhoi, rhok, rho1j[3], rho1k, tmp;
	double F1i, F1k, z, z1, phi;
	float *gradient;
	npy_intp *dims;

	PyArrayObject *py_coords;
	PyObject *py_symbols;
	PyObject *py_result = NULL;
	PyObject *py_grad;

	if(!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &py_coords, &PyList_Type, &py_symbols))
		return NULL;

	dims = PyArray_DIMS(py_coords);
	nat = PyList_Size(py_symbols);
	if (nat != dims[0]) {
		PyErr_SetString(PyExc_RuntimeError, "Number of coordinates and symbols does not match");
		return NULL; }
	ctype = PyArray_TYPE(py_coords);

	/* Set-up the raw arrays for gradients */
	gradient = (float*) malloc(3 * nat * sizeof(float));
	if(gradient == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return NULL; }

    for (k = 0; k < nat; k++) {
        gradient[3*k  ] = 0.0;
        gradient[3*k+1] = 0.0;
        gradient[3*k+2] = 0.0;
		if (ctype == NPY_FLOAT) {
			kx = *( (float*) PyArray_GETPTR2(py_coords, k, 0) );
			ky = *( (float*) PyArray_GETPTR2(py_coords, k, 1) );
			kz = *( (float*) PyArray_GETPTR2(py_coords, k, 2) );
		} else {
			kx = *( (double*) PyArray_GETPTR2(py_coords, k, 0) );
			ky = *( (double*) PyArray_GETPTR2(py_coords, k, 1) );
			kz = *( (double*) PyArray_GETPTR2(py_coords, k, 2) );
		}
        for (i = 0; i < nat; i++) {
            if ( i == k ) continue;
			if (ctype == NPY_FLOAT) {
				ix = *( (float*) PyArray_GETPTR2(py_coords, i, 0) );
				iy = *( (float*) PyArray_GETPTR2(py_coords, i, 1) );
				iz = *( (float*) PyArray_GETPTR2(py_coords, i, 2) );
			} else {
				ix = *( (double*) PyArray_GETPTR2(py_coords, i, 0) );
				iy = *( (double*) PyArray_GETPTR2(py_coords, i, 1) );
				iz = *( (double*) PyArray_GETPTR2(py_coords, i, 2) );
			}
			dik2 = sq(ix-kx) + sq(iy-ky) + sq(iz-kz);
            dik = sqrt(dik2);
            rhoi = 0.0;
            for (j = 0; j < nat; j++) {
                if (i == j) continue;
				if (ctype == NPY_FLOAT) {
					jx = *( (float*) PyArray_GETPTR2(py_coords, j, 0) );
					jy = *( (float*) PyArray_GETPTR2(py_coords, j, 1) );
					jz = *( (float*) PyArray_GETPTR2(py_coords, j, 2) );
				} else {
					jx = *( (double*) PyArray_GETPTR2(py_coords, j, 0) );
					jy = *( (double*) PyArray_GETPTR2(py_coords, j, 1) );
					jz = *( (double*) PyArray_GETPTR2(py_coords, j, 2) );
				}
				dij2 = sq(ix-jx) + sq(iy-jy) + sq(iz-jz);
                dij = sqrt(dij2);
                rhoi += cspline_interpolate_y(dij, self->rscale, self->rho, self->rho2);
            }
            F1i = cspline_interpolate_drv(rhoi, self->rhoscale, self->F, self->F2);
            rho1k = cspline_interpolate_drv(dik, self->rscale, self->rho, self->rho2);
            F1i *= rho1k / dik;
            gradient[3*k  ] -= F1i * (ix - kx);
            gradient[3*k+1] -= F1i * (iy - ky);
            gradient[3*k+2] -= F1i * (iz - kz);
            z = cspline_interpolate_y(dik, self->rscale, self->Z, self->Z2);
            z1 = cspline_interpolate_drv(dik, self->rscale, self->Z, self->Z2);
            phi = z*z1;
            phi -= 0.5 * z*z / dik;
            phi *= 27.2 * 0.529 / dik2;
            gradient[3*k  ] -= phi * (ix - kx);
            gradient[3*k+1] -= phi * (iy - ky);
            gradient[3*k+2] -= phi * (iz - kz);
        }
        rhok = 0.0;
        rho1j[0] = 0.0;
        rho1j[1] = 0.0;
        rho1j[2] = 0.0;
        for (j = 0; j < nat; j++) {
            if (j == k) continue;
			if (ctype == NPY_FLOAT) {
				jx = *( (float*) PyArray_GETPTR2(py_coords, j, 0) );
				jy = *( (float*) PyArray_GETPTR2(py_coords, j, 1) );
				jz = *( (float*) PyArray_GETPTR2(py_coords, j, 2) );
			} else {
				jx = *( (double*) PyArray_GETPTR2(py_coords, j, 0) );
				jy = *( (double*) PyArray_GETPTR2(py_coords, j, 1) );
				jz = *( (double*) PyArray_GETPTR2(py_coords, j, 2) );
			}
			dij2 = sq(kx-jx) + sq(ky-jy) + sq(kz-jz);
            dkj = sqrt(dkj2);
            rhok += cspline_interpolate_y(dkj, self->rscale, self->rho, self->rho2);
            tmp = cspline_interpolate_drv(dkj, self->rscale, self->rho, self->rho2);
            tmp /= dkj;
            rho1j[0] += tmp * (kx - jx);
            rho1j[1] += tmp * (ky - jy);
            rho1j[2] += tmp * (kz - jz);
            z = cspline_interpolate_y(dkj, self->rscale, self->Z, self->Z2);
            z1 = cspline_interpolate_drv(dkj, self->rscale, self->Z, self->Z2);
            phi = z*z1;
            phi -= 0.5 * z*z / dkj;
            phi *= 27.2 * 0.529 / dkj2;
            gradient[3*k  ] += phi * (kx - jx);
            gradient[3*k+1] += phi * (ky - jy);
            gradient[3*k+2] += phi * (kz - jz);
        }
        F1k = cspline_interpolate_drv(rhok, self->rhoscale, self->F, self->F2);
        gradient[3*k  ] += F1k * rho1j[0];
        gradient[3*k+1] += F1k * rho1j[1];
        gradient[3*k+2] += F1k * rho1j[2];
    }

	dims[0] = nat;
	dims[1] = 3;
	py_grad = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, (float*) gradient);

	return py_grad;
}


static PyMethodDef EAMff_methods[] = {
    {"energy", (PyCFunction)EAMff_energy, METH_VARARGS,
		"\n"
		"EAM_energy(coordinates, types)\n"
		"\n"
		"Evaluate the EAM energy of the configuration <coordinates>.\n"
		"\n" },
    {"gradient", (PyCFunction)EAMff_gradient, METH_VARARGS,
		"\n"
		"EAM_gradient(coordinates, types)\n"
		"\n"
		"Evaluate the EAM gradients of the configuration <coordinates>.\n"
		"\n" },
	{NULL}  /* Sentinel */
};

PyTypeObject EAMffType = {
	PyObject_HEAD_INIT(NULL)
	0,                         /*ob_size*/
	"moltools.EAMff",          /*tp_name*/
	sizeof(EAMff),             /*tp_basicsize*/
	0,                         /*tp_itemsize*/
	(destructor)EAMff_dealloc, /*tp_dealloc*/
	0,                         /*tp_print*/
	0,                         /*tp_getattr*/
	0,                         /*tp_setattr*/
	0,                         /*tp_compare*/
	0,                         /*tp_repr*/
	0,                         /*tp_as_number*/
	0,                         /*tp_as_sequence*/
	0,                         /*tp_as_mapping*/
	0,                         /*tp_hash */
	0,                         /*tp_call*/
	0,                         /*tp_str*/
	0,                         /*tp_getattro*/
	0,                         /*tp_setattro*/
	0,                         /*tp_as_buffer*/
	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
	"EAMff objects",           /* tp_doc */
	0,		               /* tp_traverse */
	0,		               /* tp_clear */
	0,		               /* tp_richcompare */
	0,		               /* tp_weaklistoffset */
	0,		               /* tp_iter */
	0,		               /* tp_iternext */
	EAMff_methods,             /* tp_methods */
	EAMff_members,             /* tp_members */
	0,                         /* tp_getset */
	0,                         /* tp_base */
	0,                         /* tp_dict */
	0,                         /* tp_descr_get */
	0,                         /* tp_descr_set */
	0,                         /* tp_dictoffset */
	(initproc)EAMff_init,      /* tp_init */
	0,                         /* tp_alloc */
	EAMff_new,                 /* tp_new */
};

/* End of methods of the EAMff object */



/*
 * For a given set of points (x,y) calculate second derivatives,
 * stored in y2. y2[0] is set to 0 in order to have linear function
 * for x < 0. If n > 0, set y2[n] to 0 as well. If n < 0, set
 * first derivative y1[n] = 0.
 * Adapted from Numerical Recipes in C.
 */
void cspline_calculate_drv2(double y2[], int n, double x[], double y[]) {
	int N, i;
	double *v;
	double sig, p, vn, qn;

	N = n < 0 ? -n : n;
	v = malloc((N-1) * sizeof(double));
	y2[0] = 0.0;
	v[0] = 0.0;

	for (i = 1; i < N-1; i++) {
		sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
		p = sig * y2[i-1] + 2.0;
		y2[i] = (sig - 1.0) / p;
		v[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]) - (y[i] - y[i-1]) / (x[i] - x[i-1]);
		v[i] = (6.0 * v[i] / (x[i+1] - x[i-1]) - sig * v[i-1]) / p;
	}

	if ( n < 0 ) {
		qn = 0.5;
		vn = -3.0 / sq(x[N-1]-x[N-2]) * (y[N-1]-y[N-2]);
	} else {
		qn = 0.0;
		vn = 0.0;
	}
	y2[N-1] = (vn - qn * v[N-2]) / (qn * y2[N-2] + 1.0);
	for (i = N-2; i >= 0; i--)
		y2[i] = y2[i]*y2[i+1] + v[i];

	free(v);
}


/*
 * Interpolate y(v) value using tabulated c-spline function.
 * Uses y2 calculated with spline_secondderv.
 * Adapted from Numerical Recipes in C.
 */
double cspline_interpolate_y(double v, PyObject *x, PyObject *y, PyObject *yp) {
	int k, khi, klo, n;
	double yint, h, a, b, slope, x0, y0, xlast, ylast, xcurr, x1, x2, y1, y2, y21, y22;
	npy_intp *dims;

	dims = PyArray_DIMS((PyArrayObject*)x);
	n = dims[0];
	/* 
	 * If v is out of range, get the derivative of the first/last point
	 * and extrapolate using linear function.
	 */
	x0 = *( (double*) PyArray_GETPTR1((PyArrayObject*)x, 0) );
	y0 = *( (double*) PyArray_GETPTR1((PyArrayObject*)y, 0) );
	xlast = *( (double*) PyArray_GETPTR1((PyArrayObject*)x, n-1) );
	ylast = *( (double*) PyArray_GETPTR1((PyArrayObject*)y, n-1) );
	if ( v < x0 ) {
		slope = cspline_interpolate_drv(x0, x, y, yp);
		yint = y0 - slope * (x0 - v);
		return yint;
	}

	if ( v > xlast ) {
		slope = cspline_interpolate_drv(xlast, x, y, yp);
		yint = ylast + slope * (v - xlast);
		return yint;
	}

	klo = 0;
	khi = n-1;
	while (khi - klo > 1) {
		k = (khi + klo) >> 1;
		xcurr = *( (double*) PyArray_GETPTR1((PyArrayObject*)x, k) );
		if (xcurr > v) khi = k;
		else klo = k;
	}

	x1 = *( (double*) PyArray_GETPTR1((PyArrayObject*)x, klo) );
	x2 = *( (double*) PyArray_GETPTR1((PyArrayObject*)x, khi) );
	h = x2 - x1;
	a = (x2 - v) / h;
	b = (v - x1) / h;
	y1 = *( (double*) PyArray_GETPTR1((PyArrayObject*)y, klo) );
	y2 = *( (double*) PyArray_GETPTR1((PyArrayObject*)y, khi) );
	y21 = *( (double*) PyArray_GETPTR1((PyArrayObject*)yp, klo) );
	y22 = *( (double*) PyArray_GETPTR1((PyArrayObject*)yp, khi) );
	yint = a*y1 + b*y2 + ((a*a*a - a) * y21 + (b*b*b - b) * y22) * (h*h) / 6.0;
	return yint;
}



/*
 * Calculate interpolated/extrapolated y'(v) value using tabulated c-spline function.
 * Uses y2 calculated with spline_secondderv.
 */
double cspline_interpolate_drv(double v, PyObject *x, PyObject *y, PyObject *yp) {
	int k, khi, klo, n;
	double ydrv, h, a, b, c, d, x0, xlast, xcurr, x1, x2, y1, y2, y21, y22;
	npy_intp *dims;

	dims = PyArray_DIMS((PyArrayObject*)x);
	n = dims[0];

	x0 = *( (double*) PyArray_GETPTR1((PyArrayObject*)x, 0) );
	xlast = *( (double*) PyArray_GETPTR1((PyArrayObject*)x, n-1) );
	/*
	 * If v fails outside the range (can be negative!), set v=0, in order to
	 * return derivative that corresponds to the first point
	 */
	if ( v < x0 ) {
		v = 0.0; }

	/*
	 * If v fails outside the range (can be negative!), set v=x[n], in order to
	 * return derivative that corresponds to the last point.
	 */
	if ( v > xlast ) {
		v = xlast; }

	klo = 0;
	khi = n-1;
	while (khi - klo > 1) {
		k = (khi + klo) >> 1;
		xcurr = *( (double*) PyArray_GETPTR1((PyArrayObject*)x, k) );
		if (xcurr > v) khi = k;
		else klo = k;
	}

	x1 = *( (double*) PyArray_GETPTR1((PyArrayObject*)x, klo) );
	x2 = *( (double*) PyArray_GETPTR1((PyArrayObject*)x, khi) );
	y1 = *( (double*) PyArray_GETPTR1((PyArrayObject*)y, klo) );
	y2 = *( (double*) PyArray_GETPTR1((PyArrayObject*)y, khi) );
	y21 = *( (double*) PyArray_GETPTR1((PyArrayObject*)yp, klo) );
	y22 = *( (double*) PyArray_GETPTR1((PyArrayObject*)yp, khi) );
	h = x2 - x1;
	a = (x2 - v) / h;
	b = (v - x1) / h;
	ydrv = (y2 - y1) / (x2 - x1) - (3*a*a - 1) / 6.0 * h * y21
	   + (3*b*b - 1) / 6.0 * h * y22;

	return ydrv;
}

