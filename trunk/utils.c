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




#include "moltools.h"


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
	if ( (offset = ftell(fd)) == -1) return NULL;

	/* Get the first chunk of data; if nothing read, make it empty */
	if ( fgets(buffer, CHUNK, fd) == NULL ) buffer[0] = '\0';
	last = strlen(buffer) - 1;
	length = last + 1;
	//printf("First chunk:%s:\n", buffer);

	/* Continue reading because this is not the end of the line yet */
	while( buffer[last] != '\n' && last != -1) {

		/* Get the next chunk; if empty - finish here */
		if ( fgets(buffer, CHUNK, fd) == NULL ) buffer[0] = '\0';
		last = strlen(buffer) - 1;
		length += last + 1;
		//printf("Next  chunk:%s:\n", buffer);
	}
	//printf("Last = %d, length = %d\n", (int)last, length);

	/* Rewind the file */
	if ( fseek(fd, offset, SEEK_SET) != 0 ) return NULL;

	/* Allocate the memory */
	line = (char*) malloc( length + 1 );

	/* Realliy read the data */
	if ( fgets(line, length + 1, fd) == NULL ) buffer[0] = '\0';
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
	//printf("%d %d %d\n", length, start, end);
	length = end-start+1;
	for(i = 0; i < length; i++) {
		line[i] = line[i+start];
	}
	line[i] = '\0';

	return length;
}

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
double cspline_interpolate_y(double v, int n, double x[], double y[], double y2[]) {
	int k, khi, klo;
	double yint, h, a, b, slope;

	/* 
	 * If v is out of range, get the derivative of the first/last point
	 * and extrapolate using linear function.
	 */
	if ( v < x[0] ) {
		/* Extrapolate from the slope at x[0] */
		//slope = (y[1] - y[0]) / (x[1] - x[0]) - (x[1] - x[0]) / 3.0 * y2[0] - (x[1] - x[0]) / 6.0 * y2[1];
		slope = cspline_interpolate_drv(x[0], n, x, y, y2);
		yint = y[0] - slope * (x[0] - v);
		return yint;
	}

	if ( v > x[n-1] ) {
		/* Extrapolate from the slope at x[n-1] */
		//slope = (y[n-1] - y[n-2]) / (x[n-1] - x[n-2]) - (x[n-1] - x[n-2]) / 3.0 * y2[n-2] - (x[n-1] - x[n-2]) / 6.0 * y2[n-1];
		slope = cspline_interpolate_drv(x[n-1], n, x, y, y2);
		yint = y[n-1] + slope * (v - x[n-1]);
		return yint;
	}

	klo = 0;
	khi = n-1;
	while (khi - klo > 1) {
		k = (khi + klo) >> 1;
		if (x[k] > v) khi = k;
		else klo = k;
	}

	h = x[khi] - x[klo];
	a = (x[khi] - v) / h;
	b = (v - x[klo]) / h;
	yint = a*y[klo] + b*y[khi] + ((a*a*a - a) * y2[klo] + (b*b*b - b) * y2[khi]) * (h*h) / 6.0;
	return yint;
}



/*
 * Calculate interpolated/extrapolated y'(v) value using tabulated c-spline function.
 * Uses y2 calculated with spline_secondderv.
 */
double cspline_interpolate_drv(double v, int n, double x[], double y[], double y2[]) {
	int k, khi, klo;
	double yp, h, a, b, c, d;

	/*
	 * If v fails outside the range (can be negative!), set v=0, in order to
	 * return derivative that corresponds to the first point
	 */
	if ( v < x[0] ) {
		/* Extrapolate from the slope at x[0] */
		//yp = (y[1] - y[0]) / (x[1] - x[0]) - (x[1] - x[0]) / 3.0 * y2[0] - (x[1] - x[0]) / 6.0 * y2[1];
		//return yp;
		v = 0.0;
	}

	/*
	 * If v fails outside the range (can be negative!), set v=x[n], in order to
	 * return derivative that corresponds to the last point.
	 */
	if ( v > x[n-1] ) {
		/* Extrapolate from the slope at x[n-1] */
		//yp = (y[n-1] - y[n-2]) / (x[n-1] - x[n-2]) - (x[n-1] - x[n-2]) / 3.0 * y2[n-2] - (x[n-1] - x[n-2]) / 6.0 * y2[n-1];
		//return yp;
		v = x[n-1];
	}

	klo = 0;
	khi = n-1;
	while (khi - klo > 1) {
		k = (khi + klo) >> 1;
		if (x[k] > v) khi = k;
		else klo = k;
	}

	h = x[khi] - x[klo];
	a = (x[khi] - v) / h;
	b = (v - x[klo]) / h;
	yp = (y[khi] - y[klo]) / (x[khi] - x[klo]) - (3*a*a - 1) / 6.0 * h * y2[klo]
	   + (3*b*b - 1) / 6.0 * h * y2[khi];

	/*h = x[klo] - x[khi];
	d = (y2[klo] - y2[khi]) / h;
	c = y2[klo] - d*x[klo];
	b = (y[klo] - y[khi] - c * (sq(x[klo]) - sq(x[khi])) + d * (sq(x[klo])*x[klo] - sq(x[khi])*x[khi])) / h;
	yp = b + c*v + d*v*v;*/
	return yp;
}


int getElementIndexBySymbol(const char *symbol) {
	extern Element element_table[];
	int idx = 0;

	while (element_table[idx].number != -1) {
		if (!strcmp(symbol, element_table[idx++].symbol)) break;
	}

	return idx;
}
