/***************************************************************************

    mdarray

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


#ifndef __MDARRAY_H__
#define __MDARRAY_H__

/* These two declarations should be package-wide, so put them here! */
#define PY_ARRAY_UNIQUE_SYMBOL MDARRAY
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/* To make sure that Python and Numpy stuff is declared everywhere, *
 * put this here, in the header                                     */
#include <locale.h>
#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>
#include <numpy/halffloat.h>

#ifdef HAVE_GROMACS
	#include <gromacs/utility/smalloc.h>
	#include <gromacs/fileio/xtcio.h>
#endif

#define BOHRTOANGS 0.529177209

#define ARRAY_REAL double
#define NPY_ARRAY_REAL NPY_DOUBLE
#define TYPE_NAME(t) __TYPE_NAME(t)
#define __TYPE_NAME(t) #t


#endif /* __MDARRAY_H__ */
