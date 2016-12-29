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

#include "periodic_table.h"
#include "utils.h"


/*****************************************************************************
 *  This is a helper function that reads a line from the file desctiptor fd  *
 *  and returns an allocated string, including the new-line character.       *
 *  If the read operation return no characters an empty string is returned.  *
 *  In any case, the string is terminated with the '\0' character.           *
 *  Make sure, you free the allocated memory after use!                      *
 *****************************************************************************/

#define CHUNK 10
char *readline(FILE *fd) {
	int length, last;
	long int offset;
	char buffer[CHUNK];
	char *line;

	/* Initial offset in the file */
	if ( (offset = ftell(fd)) == -1) {
		PyErr_SetFromErrno(PyExc_IOError);
		return NULL; }

	/* Get the first chunk of data; if nothing read, make it empty */
	if ( fgets(buffer, CHUNK, fd) == NULL ) buffer[0] = '\0';
	last = strlen(buffer) - 1;
	length = last + 1;

	/* Continue reading because this is not the end of the line yet */
	while( buffer[last] != '\n' && last != -1) {

		/* Get the next chunk; if empty - finish here */
		if ( fgets(buffer, CHUNK, fd) == NULL ) buffer[0] = '\0';
		last = strlen(buffer) - 1;
		length += last + 1;
	}

	/* Rewind the file */
	if ( fseek(fd, offset, SEEK_SET) == -1 ) {
		PyErr_SetFromErrno(PyExc_IOError);
		return NULL; }

	/* Allocate the memory */
	line = (char*) malloc( length + 1 );

	/* Really read the data */
	if ( fgets(line, length + 1, fd) == NULL ) return NULL;
	//printf("Line:%s:\n", line);

	return line;
}




int make_lowercase(char *line) {
	int length, i, c, count = 0;

	length = strlen(line);
	for(i = 0; i < length; i++ ) {
		c = tolower(line[i]);
		if ( c != line[i] ) count++;
		line[i] = c;
	}

	return count;
}



int stripline(char *line) {
	int length, i, start, end;
	char c;

	length = strlen(line);
	start = 0;
	while( start < length ) {
		c = line[start++];
		if ( c != ' ' && c != '\t' && c != '\n' && c != '\r' ) break;
	}
	start -= 1;
	end = length-1;
	while( end > 0 ) {
		c = line[end--];
		if ( c != ' ' && c != '\t' && c != '\n' && c != '\r' ) break;
	}
	end += 1;
	length = end-start+1;
	for(i = 0; i < length; i++) {
		line[i] = line[i+start];
	}
	line[i] = '\0';

	return length;
}


/* Extract number from buf (at pos, with length len) and return as float *
 * For example: "abc  45.6 xyz"                                          *
 *               0123456789012                                           *
 *               strPartFloat(buf, 4, 5) -> " 45.6" -> 45.6              */

float strPartFloat(char *buf, int pos, int len) {
	char *number;
	float fret;

	number = (char*) malloc((len+1) * sizeof(char));
	strncpy(number, buf+pos, len);
	number[len] = '\0';
	stripline(number);
	fret = (float)atof(number);
	free(number);
	return fret;
}



int getElementIndexBySymbol(const char *symbol) {
	extern Element element_table[];
	int idx = 0;

	while (element_table[idx].number != -1) {
		if (!strcmp(symbol, element_table[idx].symbol)) return idx;
		idx += 1;
	}

	return -1;
}


double *vectorToDouble(double dvec[], PyArrayObject *arr) {
	int type;
	npy_intp *numpyint;

    numpyint = PyArray_DIMS(arr);
    if (numpyint[0] != 3) {
        PyErr_SetString(PyExc_ValueError, "Cartesian vector should contain exactly 3 numbers");
        return NULL; }

	type = PyArray_TYPE(arr);
    if( type != NPY_FLOAT && type != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "Incorrect type of the array");
        return NULL; }
    switch(type) {
        case NPY_FLOAT:
            dvec[0] = *( (float*) PyArray_GETPTR1(arr, 0) );
            dvec[1] = *( (float*) PyArray_GETPTR1(arr, 1) );
            dvec[2] = *( (float*) PyArray_GETPTR1(arr, 2) );
            break;
        case NPY_DOUBLE:
            dvec[0] = *( (double*) PyArray_GETPTR1(arr, 0) );
            dvec[1] = *( (double*) PyArray_GETPTR1(arr, 1) );
            dvec[2] = *( (double*) PyArray_GETPTR1(arr, 2) );
            break;
    }
    return dvec;
}


void wrapCartesian(double point[3], double box[3]) {
	point[0] = fmod(point[0], box[0]);
	point[1] = fmod(point[1], box[1]);
	point[2] = fmod(point[2], box[2]);
	if (point[0] < 0.0) point[0] += box[0];
	if (point[1] < 0.0) point[1] += box[1];
	if (point[2] < 0.0) point[2] += box[2];
}



void nearestImage(double center[3], double other[3], double half[3]) {
	if (center[0] - other[0] > half[0]) other[0] += half[0]+half[0];
	else if (center[0] - other[0] < -half[0]) other[0] -= half[0]+half[0];
	if (center[1] - other[1] > half[1]) other[1] += half[1]+half[1];
	else if (center[1] - other[1] < -half[1]) other[1] -= half[1]+half[1];
	if (center[2] - other[2] > half[2]) other[2] += half[2]+half[2];
	else if (center[2] - other[2] < -half[2]) other[2] -= half[2]+half[2];
}


double threePointAngleCosine(double A[3], double B[3], double C[3]) {
	double p[3], q[3];
	double lp, lq, cos;

	p[0] = A[0] - B[0];
	p[1] = A[1] - B[1];
	p[2] = A[2] - B[2];

	q[0] = C[0] - B[0];
	q[1] = C[1] - B[1];
	q[2] = C[2] - B[2];

	lp = sq(p[0]) + sq(p[1]) + sq(p[2]);
	lq = sq(q[0]) + sq(q[1]) + sq(q[2]);

	cos = (p[0]*q[0] + p[1]*q[1] + p[2]*q[2])/sqrt(lp*lq);

	return cos;
}


double distanceSquare(double p[3], double q[3]) {
	return sq(p[0]-q[0]) + sq(p[1]-q[1]) + sq(p[2]-q[2]);
}


void copyPoint(double dst[3], double src[3]) {
	dst[0] = src[0];
	dst[1] = src[1];
	dst[2] = src[2];
}

float getFromArray2D(PyObject *arr, int i, int j) {
	int type;
	npy_half hx;
	float fx;
	double dx;
	long double lx;

	type = PyArray_TYPE((PyArrayObject*)arr);
	switch(type) {
		case NPY_HALF:
	        hx = *( (npy_half*) PyArray_GETPTR2((PyArrayObject*)arr, i, j));
			fx = npy_half_to_float(hx);
        	return fx;
			break;
		case NPY_FLOAT:
	        fx = *( (float*) PyArray_GETPTR2((PyArrayObject*)arr, i, j));
        	return fx;
			break;
		case NPY_DOUBLE:
	        dx = *( (double*) PyArray_GETPTR2((PyArrayObject*)arr, i, j));
        	return (float)dx;
			break;
		case NPY_LONGDOUBLE:
	        lx = *( (long double*) PyArray_GETPTR2((PyArrayObject*)arr, i, j));
        	return (float)lx;
			break;
		default:
			PyErr_SetString(PyExc_ValueError, "Incorrect type of coordinate array");
			return NAN;
			break;
    }
}
