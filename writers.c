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

#include "trajectory.h"


int traj_write_gro(Trajectory *self, PyObject *py_coords, PyObject *py_vel,
                PyObject *py_box, char *comment) {

    int i, resid, type;
    npy_half hx, hy, hz;
    float x, y, z;
    double dx, dy, dz;
    long double lx, ly, lz;
    char *s, *resnam;
    char empty[1] = "";
    int box_order[9][2] = { {0,0}, {1,1}, {2,2}, {0,1}, {0,2}, {1,0}, {1,2}, {2,0}, {2,1} };
    double box[9];
    npy_intp *dims;

    if( comment != NULL )
        fprintf(self->fd, "%s\n", comment);
    else
        fprintf(self->fd, "\n");
    fprintf(self->fd, "%5d\n", self->nofatoms);

    type = PyArray_TYPE((PyArrayObject*)py_coords);
    for (i = 0; i < self->nofatoms; i++) {
        if (self->resids != Py_None)
            resid = *((int*) PyArray_GETPTR1((PyArrayObject*)self->resids, i));
        else
            resid = 1;
        if (self->resnames != Py_None)
            resnam = PyString_AsString(PyList_GetItem(self->resnames, i));
        else
            resnam = empty;
        s = PyString_AsString(PyList_GetItem(self->symbols, i));
        /* Depending on the dtype used in python, the pointer has to be    *
         * casted to the right C type and the correct C type variable must *
         * be used.                                                        */
        switch(type) {
            case NPY_HALF:
                hx = *( (npy_half*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 0) ) / 10.0;
                hy = *( (npy_half*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 1) ) / 10.0;
                hz = *( (npy_half*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 2) ) / 10.0;
                x = npy_half_to_float(hx);
                y = npy_half_to_float(hy);
                z = npy_half_to_float(hz);
                fprintf(self->fd, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", resid, resnam, s, i+1, x, y, z);
                break;
            case NPY_FLOAT:
                x = *( (float*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 0) ) / 10.0;
                y = *( (float*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 1) ) / 10.0;
                z = *( (float*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 2) ) / 10.0;
                fprintf(self->fd, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", resid, resnam, s, i+1, x, y, z);
                break;
            case NPY_DOUBLE:
                dx = *( (double*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 0) ) / 10.0;
                dy = *( (double*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 1) ) / 10.0;
                dz = *( (double*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 2) ) / 10.0;
                fprintf(self->fd, "%5d%-5s%5s%5d%8.3lf%8.3lf%8.3lf\n", resid, resnam, s, i+1, dx, dy, dz);
                break;
            case NPY_LONGDOUBLE:
                lx = *( (long double*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 0) ) / 10.0;
                ly = *( (long double*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 1) ) / 10.0;
                lz = *( (long double*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 2) ) / 10.0;
                fprintf(self->fd, "%5d%-5s%5s%5d%8.3Lf%8.3Lf%8.3Lf\n", resid, resnam, s, i+1, lx, ly, lz);
                break;
            default:
                PyErr_SetString(PyExc_ValueError, "Incorrect type of coordinate array");
                return -1;
        }
    }

    /* Do some testing on the array */
    if (py_box != NULL) {
        if ( PyArray_NDIM((PyArrayObject*)py_box) != 2 ) {
            PyErr_SetString(PyExc_ValueError, "Incorrect box shape");
            return -1; }
        else {
            dims = PyArray_DIMS((PyArrayObject*)py_box);
            if (dims[0] != 3 || dims[1] != 3) {
                PyErr_SetString(PyExc_ValueError, "Incorrect box shape");
                return -1; }
            type = PyArray_TYPE((PyArrayObject*)py_box);
            for (i = 0; i < 9; i++) {
                switch(type) {
                    case NPY_HALF:
                        hx = *( (npy_half*) PyArray_GETPTR2((PyArrayObject*)py_box, box_order[i][0], box_order[i][1]) );
                        x = npy_half_to_float(hx)/10.0;
                        box[i] = (float)x;
                        break;
                    case NPY_FLOAT:
                        x = *( (float*) PyArray_GETPTR2((PyArrayObject*)py_box, box_order[i][0], box_order[i][1]) );
                        x /= 10.0;
                        box[i] = (float)x;
                        break;
                    case NPY_DOUBLE:
                        dx = *( (double*) PyArray_GETPTR2((PyArrayObject*)py_box, box_order[i][0], box_order[i][1]) );
                        dx /= 10.0;
                        box[i] = (float)dx;
                        break;
                    case NPY_LONGDOUBLE:
                        lx = *( (long double*) PyArray_GETPTR2((PyArrayObject*)py_box, box_order[i][0], box_order[i][1]) );
                        lx /= 10.0;
                        box[i] = (float)lx;
                        break;
                    default:
                        PyErr_Format(PyExc_ValueError, "Incorrect type in box vector (%u)", type);
                        return -1;
                }
            }
            for (i = 0; i < 3; i++)
                fprintf(self->fd, "%10.5f", box[i]);
            if (fabs(box[3]) > 1e-6 || fabs(box[4]) > 1e-6 || fabs(box[5]) > 1e-6 ||
                fabs(box[6]) > 1e-6 || fabs(box[7]) > 1e-6 || fabs(box[8]) > 1e-6)
                for (i = 3; i < 9; i++)
                    fprintf(self->fd, "%10.5f", box[i]);
            fprintf(self->fd, "\n");
        }
    } else
        fprintf(self->fd, "%10.5f%10.5f%10.5f\n", 0.0, 0.0, 0.0);

    return 0;
}


