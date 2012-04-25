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



PyObject *read_xyz(FILE *fd, float factor) {

	int nofatoms, pos;
	char *buffer, *buffpos, *token;
	float *xyz, *charges;
	short int charges_present;

	npy_intp dims[2];

	PyObject *key, *val, *py_result, *py_coord, *py_syms, *py_charges;


	/* Create the dictionary that will be returned */
	py_result = PyDict_New();


	/* Read number of atoms */
	if( (buffer = readline(fd)) == NULL) {
		PyErr_SetFromErrno(PyExc_IOError);
		return NULL; }
	if( sscanf(buffer, "%d", &nofatoms) != 1 ) {
		PyErr_SetString(PyExc_IOError, "Incorrect atom number");
		return NULL; }

    val = Py_BuildValue("i", nofatoms);
	key = PyString_FromString("number_of_atoms");
	PyDict_SetItem(py_result, key, val);
	Py_DECREF(key);
	Py_DECREF(val);


	/* Read the comment line */
	if((buffer = readline(fd)) == NULL) {
		PyErr_SetFromErrno(PyExc_IOError);
		return NULL; }
	buffer[strlen(buffer)-1] = '\0';

	val = Py_BuildValue("s", buffer);
	free(buffer);
	key = PyString_FromString("comment");
	PyDict_SetItem(py_result, key, val);
	Py_DECREF(key);
	Py_DECREF(val);


	/* Set-up the raw arrays for coordinates and charges */
	xyz = (float*) malloc(3 * nofatoms * sizeof(float));
	if(xyz == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return NULL; }
	charges = (float*) malloc(nofatoms * sizeof(float));
	if(charges == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return NULL; }
	charges_present = 0;


	py_syms = PyList_New(nofatoms);

	/* Atom loop */
	for(pos = 0; pos < nofatoms; pos++) {

		/* Get the whole line */
		if((buffer = readline(fd)) == NULL) {
			PyErr_SetFromErrno(PyExc_IOError);
			return NULL; }
		buffer[strlen(buffer)-1] = '\0';
		buffpos = buffer;

		/* Read symbol */
		token = strtok(buffpos, " ");
		val = Py_BuildValue("s", token);
		PyList_SetItem(py_syms, pos, val);

		/* Read coordinates */
		if ( (token = strtok(NULL, " ")) == NULL) {
			PyErr_SetString(PyExc_IOError, "Missing coordinate");
			return NULL; }
		xyz[3*pos + 0] = atof(token) * factor;
		if ( (token = strtok(NULL, " ")) == NULL) {
			PyErr_SetString(PyExc_IOError, "Missing coordinate");
			return NULL; }
		xyz[3*pos + 1] = atof(token) * factor;
		if ( (token = strtok(NULL, " ")) == NULL) {
			PyErr_SetString(PyExc_IOError, "Missing coordinate");
			return NULL; }
		xyz[3*pos + 2] = atof(token) * factor;

		/* Read charge, if present */
		token = strtok(NULL, " ");
		if ( token != NULL ) {

			/* This is bad: until now, there were no charges */
			if ( pos > 0 && !charges_present ) {
				PyErr_SetString(PyExc_IOError, "Unexpected charges found");
				return NULL;
			}

			charges_present = 1;
			charges[pos] = atof(token);

		} else {

			/* This is bad: we were expecting charges here and found nothing */
			if ( pos > 0 && charges_present ) {
				PyErr_SetString(PyExc_IOError, "Missing charges");
				return NULL;
			}
		}

		/* Free the line buffer */
		free(buffer);
	}


	/* Add symbols to the dictionary */
	key = PyString_FromString("symbols");
	PyDict_SetItem(py_result, key, py_syms);
	Py_DECREF(key);
	Py_DECREF(py_syms);


	/* Add coordinates to the dictionary */
	dims[0] = nofatoms;
	dims[1] = 3;
	py_coord = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, (float*) xyz);
	/***************************************************************
	 * Do not free the raw array! It will be still used by Python! *

	free(*xyz);
	free(xyz);

     * when the ver. 1.7 arrives, use PyArray_SetBaseObject        *
     * to prevent memory leaks.                                    *
     ***************************************************************/

	key = PyString_FromString("coordinates");
	PyDict_SetItem(py_result, key, py_coord);
	Py_DECREF(key);
	Py_DECREF(py_coord);


	/* Add charges, if present */
	if ( charges_present ) {
		py_charges = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT, charges);
		key = PyString_FromString("charges");
		PyDict_SetItem(py_result, key, py_charges);
		Py_DECREF(key);
		Py_DECREF(py_charges);
	} else
		/* Free the charges ONLY if the Python object was not created! */
		free(charges);


	return py_result;
}


/* Read molden file */
PyObject *read_molden(FILE *fd) {

	char *line;
	char *buffpos, *token;
	char **line_store;
	int n_geo, i, nat;
	float *xyz;
	int *anum;

	double *geoconv;
	npy_intp dims[2];
	PyObject *py_anum, *py_syms, *py_geom, *py_geoconv, *py_result, *key, *val;


	/* Prepare dictionary */

	py_result = PyDict_New();

	/* Loop over the lines */

	line = readline(fd);

	while ( strlen(line) ) {

		stripline(line);
		make_lowercase(line);

		/* This is a start of a new block */
		if ( line[0] == '[' ) {

			/* Number of geometries present in the file */
			if ( !strcmp(line, "[n_geo]") ) {
				free(line);
				line = readline(fd);
				stripline(line);
				n_geo = atoi(line);
				key = PyString_FromString("number_of_geometries");
				val = Py_BuildValue("i", n_geo);
				PyDict_SetItem(py_result, key, val);
				Py_DECREF(key);
				Py_DECREF(val);

				/* This will be used by arrays depending on number of geoms */
				dims[0] = n_geo;
				dims[1] = 3;
			}

			/* Energy of the subsequent geometries */
			else if ( !strcmp(line, "[geoconv]" ) ) {
				free(line);
				line = readline(fd);
				geoconv = (double*) malloc(n_geo * sizeof(double));
				if(geoconv == NULL) {
					PyErr_SetFromErrno(PyExc_MemoryError);
					return NULL; }
				for( i = 0; i < n_geo; i++) {
					free(line);
					line = readline(fd);
					stripline(line);
					geoconv[i] = atof(line);
				}
				py_geoconv = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, geoconv);
				key = PyString_FromString("energies");
				PyDict_SetItem(py_result, key, py_geoconv);
				Py_DECREF(key);
				Py_DECREF(py_geoconv);
			}

			/* Energy of the subsequent geometries */
			else if ( !strcmp(line, "[geometries] (xyz)" ) ) {
				py_geom = PyList_New(n_geo);
				for ( i = 0; i < n_geo; i++ )
					PyList_SetItem(py_geom, i, read_xyz(fd, 1.0));
				key = PyString_FromString("geometries");
				PyDict_SetItem(py_result, key, py_geom);
				Py_DECREF(key);
				Py_DECREF(py_geom);
			}

			/* Section 'atoms' present - this is one-geometry file */
			else if ( !strcmp(line, "[atoms] angs" ) ) {

				n_geo = 1;
				key = PyString_FromString("number_of_geometries");
				val = Py_BuildValue("i", n_geo);
				PyDict_SetItem(py_result, key, val);
				Py_DECREF(key);
				Py_DECREF(val);

				free(line);
				line = readline(fd);
				stripline(line);

				/* We don't know how many atoms are there, so we have to
                   store the lines. */
				nat = 0;
				line_store = (char**) malloc(10 * sizeof(char*));
				if ( line_store == NULL ) {
					PyErr_SetFromErrno(PyExc_MemoryError);
					return NULL; }
				while ( line[0] != '[' ) {
					line_store[nat++] = line;
					if ( nat % 10 == 0 ) {
						line_store = realloc(line_store, (nat + 10) * sizeof(char*));
						if( line_store == NULL ) {
							PyErr_SetFromErrno(PyExc_MemoryError);
						return NULL; }
					}
					line = readline(fd);
					stripline(line);
				}

				xyz = (float*) malloc(3 * nat * sizeof(float));
				if(xyz == NULL) {
					PyErr_SetFromErrno(PyExc_MemoryError);
					return NULL; }
				anum = (int*) malloc(nat * sizeof(int));
				if(anum == NULL) {
					PyErr_SetFromErrno(PyExc_MemoryError);
					return NULL; }

				py_syms = PyList_New(nat);

				/* Loop over atoms */
				for ( i = 0; i < nat; i++ ) {

					buffpos = line_store[i];
					token = strtok(buffpos, " ");
					val = Py_BuildValue("s", token);
					PyList_SetItem(py_syms, i, val);

					token = strtok(NULL, " ");
					/* not used */

					token = strtok(NULL, " ");
					anum[i] = atoi(token);

					token = strtok(NULL, " ");
					xyz[3*i+0] = atof(token);

					token = strtok(NULL, " ");
					xyz[3*i+1] = atof(token);

					token = strtok(NULL, " ");
					xyz[3*i+2] = atof(token);

					/* Get rid of the line. */
					free(line_store[i]);

				}

				/* Free the stored line pointers. */
				free(line_store);

				/* Add symbols to the dictionary */
				key = PyString_FromString("symbols");
				PyDict_SetItem(py_result, key, py_syms);
				Py_DECREF(key);
				Py_DECREF(py_syms);

				/* Add coordinates to the dictionary */
				dims[0] = nat;
				dims[1] = 3;
				py_geom = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, (float*) xyz);
				key = PyString_FromString("coordinates");
				PyDict_SetItem(py_result, key, py_geom);
				Py_DECREF(key);
				Py_DECREF(py_geom);

				/* Add atomic numbers to the dictionary */
				py_anum = PyArray_SimpleNewFromData(1, dims, NPY_INT, (int*) anum);
				key = PyString_FromString("atomic_numbers");
				PyDict_SetItem(py_result, key, py_anum);
				Py_DECREF(key);
				Py_DECREF(py_anum);

				/* This is to avoid reading the next line! */
				continue;
			}

		}

		line = readline(fd);

	}

	return py_result;
}



PyObject *read_fractional(FILE *fd) {

	int nofatoms, pos;
	char *buffer, *buffpos, *token;
	float *xyz;
	float *box;
	float transfm[3][3];
	float x, y, z;

	npy_intp dims[2];

	PyObject *key, *val, *py_result, *py_coord, *py_syms, *py_box;


	/* Create the dictionary that will be returned */
	py_result = PyDict_New();


	/* Read number of atoms */
	if( (buffer = readline(fd)) == NULL) {
		PyErr_SetFromErrno(PyExc_IOError);
		return NULL; }
	if( sscanf(buffer, "%d", &nofatoms) != 1 ) {
		PyErr_SetString(PyExc_IOError, "Incorrect atom number");
		return NULL; }

    val = Py_BuildValue("i", nofatoms);
	key = PyString_FromString("number_of_atoms");
	PyDict_SetItem(py_result, key, val);
	Py_DECREF(key);
	Py_DECREF(val);


	/* Instead of a comment line, there are six numbers:
       a[A] b[A] c[A] alpha[deg.] beta[deg.] gamma[deg.] */

	if ( (box = (float*) malloc(6 * sizeof(float))) == NULL ) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return NULL; }

	if((buffer = readline(fd)) == NULL) {
		PyErr_SetFromErrno(PyExc_IOError);
		return NULL; }
	buffer[strlen(buffer)-1] = '\0';
	buffpos = buffer;

	token = strtok(buffpos, " ");
	box[0] = atof(token);
	token = strtok(NULL, " ");
	box[1] = atof(token);
	token = strtok(NULL, " ");
	box[2] = atof(token);

	token = strtok(NULL, " ");
	box[3] = atof(token) / 180.0 * M_PI;
	token = strtok(NULL, " ");
	box[4] = atof(token) / 180.0 * M_PI;
	token = strtok(NULL, " ");
	box[5] = atof(token) / 180.0 * M_PI;


	dims[0] = 6;
	py_box = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT, box);
	key = PyString_FromString("lattice");
	PyDict_SetItem(py_result, key, py_box);
	Py_DECREF(key);
	Py_DECREF(py_box);

	/* Set the tranformation matrix */
	transfm[0][0] = box[0];
	transfm[0][1] = box[1] * cos(box[5]);
	transfm[0][2] = box[2] * cos(box[4]);
	transfm[1][0] = 0;
	transfm[1][1] = box[1] * sin(box[5]);
	transfm[1][2] = box[2] * (cos(box[3]) - cos(box[4])*cos(box[5])) / sin(box[5]);
	transfm[2][0] = 0;
	transfm[2][1] = 0;
	transfm[2][2] = box[2] * sqrt(1.0 - sq(cos(box[3])) - sq(cos(box[4])) - sq(cos(box[5]))
	                                  + 2*cos(box[3])*cos(box[4])*cos(box[5])) / sin(box[5]);

	/* Set-up the raw arrays for coordinates and charges */
	xyz = (float*) malloc(3 * nofatoms * sizeof(float));
	if(xyz == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return NULL; }

	py_syms = PyList_New(nofatoms);

	/* Atom loop */
	for(pos = 0; pos < nofatoms; pos++) {

		/* Get the whole line */
		if((buffer = readline(fd)) == NULL) {
			PyErr_SetFromErrno(PyExc_IOError);
			return NULL; }
		buffer[strlen(buffer)-1] = '\0';
		buffpos = buffer;

		/* Read symbol */
		token = strtok(buffpos, " ");
		val = Py_BuildValue("s", token);
		PyList_SetItem(py_syms, pos, val);

		/* Read coordinates */
		if ( (token = strtok(NULL, " ")) == NULL) {
			PyErr_SetString(PyExc_IOError, "Missing coordinate");
			return NULL; }
		x = atof(token);
		if ( (token = strtok(NULL, " ")) == NULL) {
			PyErr_SetString(PyExc_IOError, "Missing coordinate");
			return NULL; }
		y = atof(token);
		if ( (token = strtok(NULL, " ")) == NULL) {
			PyErr_SetString(PyExc_IOError, "Missing coordinate");
			return NULL; }
		z = atof(token);

		/* Convert from fractional to cartesian */
		xyz[3*pos + 0] = x * transfm[0][0] + y * transfm[0][1] + z * transfm[0][2];
		xyz[3*pos + 1] = y * transfm[1][1] + z * transfm[1][2];
		xyz[3*pos + 2] = z * transfm[2][2];

		/* Free the line buffer */
		free(buffer);
	}


	/* Add symbols to the dictionary */
	key = PyString_FromString("symbols");
	PyDict_SetItem(py_result, key, py_syms);
	Py_DECREF(key);
	Py_DECREF(py_syms);


	/* Add coordinates to the dictionary */
	dims[0] = nofatoms;
	dims[1] = 3;
	py_coord = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, (float*) xyz);
	/***************************************************************
	 * Do not free the raw array! It will be still used by Python! *

	free(*xyz);
	free(xyz);

     * when the ver. 1.7 arrives, use PyArray_SetBaseObject        *
     * to prevent memory leaks.                                    *
     ***************************************************************/

	key = PyString_FromString("coordinates");
	PyDict_SetItem(py_result, key, py_coord);
	Py_DECREF(key);
	Py_DECREF(py_coord);

	return py_result;
}



PyObject *read_gro(FILE *fd) {

	int nofatoms, pos;
	char *buffer;
	char symbuf[100];
	float *xyz, *vel;
	int *resid;
	unsigned short int velocities_present = 0;

	npy_intp dims[2];

	PyObject *key, *val, *py_result, *py_coord, *py_vel, *py_syms, *py_resn, *py_resid;


	/* Create the dictionary that will be returned */
	py_result = PyDict_New();


	/* Read the comment line */
	if((buffer = readline(fd)) == NULL) {
		PyErr_SetFromErrno(PyExc_IOError);
		return NULL; }
	buffer[strlen(buffer)-1] = '\0';
	stripline(buffer);

	val = Py_BuildValue("s", buffer);
	free(buffer);
	key = PyString_FromString("comment");
	PyDict_SetItem(py_result, key, val);
	Py_DECREF(key);
	Py_DECREF(val);

	/* Read number of atoms */
	if( (buffer = readline(fd)) == NULL) {
		PyErr_SetFromErrno(PyExc_IOError);
		return NULL; }
	if( sscanf(buffer, "%d", &nofatoms) != 1 ) {
		PyErr_SetString(PyExc_IOError, "Incorrect atom number");
		return NULL; }

    val = Py_BuildValue("i", nofatoms);
	key = PyString_FromString("number_of_atoms");
	PyDict_SetItem(py_result, key, val);
	Py_DECREF(key);
	Py_DECREF(val);


	/* Set-up the raw arrays for coordinates and charges */
	xyz = (float*) malloc(3 * nofatoms * sizeof(float));
	if(xyz == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return NULL; }
	vel = (float*) malloc(3 * nofatoms * sizeof(float));
	if(vel == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return NULL; }
	resid = (int*) malloc(nofatoms * sizeof(int));
	if(resid == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return NULL; }


	py_syms = PyList_New(nofatoms);
	py_resn = PyList_New(nofatoms);

	/* Atom loop */
	for(pos = 0; pos < nofatoms; pos++) {

		/* Get the whole line */
		if((buffer = readline(fd)) == NULL) {
			PyErr_SetFromErrno(PyExc_IOError);
			return NULL; }
		if( pos == 0 && strlen(buffer) > 50) velocities_present = 1;

		/* Read residue id */
		strncpy(symbuf, buffer, 5);
		symbuf[5] = '\0';
		stripline(symbuf);
		resid[pos] = atoi(symbuf);

		/* Read residue name */
		strncpy(symbuf, buffer+5, 5);
		symbuf[5] = '\0';
		stripline(symbuf);
		val = Py_BuildValue("s", symbuf);
		PyList_SetItem(py_resn, pos, val);

		/* Read atom name */
		strncpy(symbuf, buffer+10, 5);
		symbuf[5] = '\0';
		stripline(symbuf);
		val = Py_BuildValue("s", symbuf);
		PyList_SetItem(py_syms, pos, val);

		/* Read coordinates */
		strncpy(symbuf, buffer+20, 8);
		symbuf[8] = '\0';
		stripline(symbuf);
		xyz[3*pos + 0] = atof(symbuf);
		strncpy(symbuf, buffer+28, 8);
		symbuf[8] = '\0';
		stripline(symbuf);
		xyz[3*pos + 1] = atof(symbuf);
		strncpy(symbuf, buffer+36, 8);
		symbuf[8] = '\0';
		stripline(symbuf);
		xyz[3*pos + 2] = atof(symbuf);

		/* Read velocities */
		if(velocities_present) {
			strncpy(symbuf, buffer+44, 8);
			symbuf[8] = '\0';
			stripline(symbuf);
			vel[3*pos + 0] = atof(symbuf);
			strncpy(symbuf, buffer+52, 8);
			symbuf[8] = '\0';
			stripline(symbuf);
			vel[3*pos + 1] = atof(symbuf);
			strncpy(symbuf, buffer+60, 8);
			symbuf[8] = '\0';
			stripline(symbuf);
			vel[3*pos + 2] = atof(symbuf);
		}

		/* Free the line buffer */
		free(buffer);
	}


	/* Add symbols to the dictionary */
	key = PyString_FromString("atom_names");
	PyDict_SetItem(py_result, key, py_syms);
	Py_DECREF(key);
	Py_DECREF(py_syms);
	key = PyString_FromString("residue_names");
	PyDict_SetItem(py_result, key, py_resn);
	Py_DECREF(key);
	Py_DECREF(py_resn);


	/* Add coordinates to the dictionary */
	dims[0] = nofatoms;
	dims[1] = 3;
	py_coord = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, (float*) xyz);
	/***************************************************************
	 * Do not free the raw array! It will be still used by Python! *

	free(xyz);

     * when the ver. 1.7 arrives, use PyArray_SetBaseObject        *
     * to prevent memory leaks.                                    *
     ***************************************************************/

	key = PyString_FromString("coordinates");
	PyDict_SetItem(py_result, key, py_coord);
	Py_DECREF(key);
	Py_DECREF(py_coord);


	/* Add coordinates to the dictionary */
	if(velocities_present) {
		py_vel = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, (float*) vel);
		/***************************************************************
		 * Do not free the raw array! It will be still used by Python! *

		free(xyz);

    	 * when the ver. 1.7 arrives, use PyArray_SetBaseObject        *
	     * to prevent memory leaks.                                    *
    	 ***************************************************************/

		key = PyString_FromString("velocities");
		PyDict_SetItem(py_result, key, py_vel);
		Py_DECREF(key);
		Py_DECREF(py_vel);
	} else
		free(vel);


	return py_result;
}


