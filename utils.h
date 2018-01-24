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


#ifndef __UTILS_H__
#define __UTILS_H__

#include <Python.h>
#include <numpy/halffloat.h>

#define sq(a) ((a) * (a))

char *readline(FILE *);
int make_lowercase(char *);
int stripline(char *);
float strPartFloat(char *buf, int pos, int len);
int getElementIndexBySymbol(const char *symbol);
double *vectorToDouble(double dvec[], PyArrayObject *arr);
void wrapCartesian(double point[3], double box[3]);
void nearestImage(double center[3], double other[3], double half[3]);
double threePointAngleCosine(double A[3], double B[3], double C[3]);
double distanceSquare(double p[3], double q[3]);
void copyPoint(double dst[3], double src[3]);
float getFromArray2D(PyObject *arr, int i, int j);

#endif /* __MOLTOOLS_H__ */
