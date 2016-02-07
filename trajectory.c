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




static void Trajectory_dealloc(Trajectory* self)
{
	Py_XDECREF(self->symbols);
	Py_XDECREF(self->resids);
	Py_XDECREF(self->resnames);
	Py_XDECREF(self->atomicnumbers);
	self->ob_type->tp_free((PyObject*)self);
	free(self->filename);
	switch(self->type) {
		case XYZ:
		case MOLDEN:
		case GRO:
			if (self->fd != NULL) fclose(self->fd);
			break;
#ifdef HAVE_GROMACS
		case XTC:
			if (self->xd != NULL) close_xtc(self->xd);
#endif
	}
#ifdef HAVE_GROMACS
	sfree(self->xtcCoord);
#endif
}




static PyObject *Trajectory_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	Trajectory *self;

	self = (Trajectory *)type->tp_alloc(type, 0);
	if (self != NULL) {

	    self->filename = NULL;
		self->fd = NULL;
		self->xd = NULL;
		self->mode = 'r';
		self->type = GUESS;
		self->nofatoms = 0;
		self->nofframes = 0;
		self->xtcCoord = NULL;
		self->units = ANGS;
		self->lastFrame = -1;
		self->filePosition1 = -1;
		self->filePosition2 = -1;

		Py_INCREF(Py_None);
	    self->symbols = Py_None;
		
		Py_INCREF(Py_None);
	    self->atomicnumbers = Py_None;
		
		Py_INCREF(Py_None);
	    self->resids = Py_None;
		
		Py_INCREF(Py_None);
	    self->resnames = Py_None;
    }

    return (PyObject *)self;
}




/* This method should just guess the type, open the file and collect     *
 * general data like number of atoms and list of symbols. Actual reading *
 * of coordinates should go into read method                             */

static int Trajectory_init(Trajectory *self, PyObject *args, PyObject *kwds) {

	const char *filename;
	FILE *test;
	char *str_type = NULL;
	char ext[5];
	char *line;
	char *mode = NULL;
	char *units = NULL;
#ifdef HAVE_GROMACS
	int step;
	real time, prec;
	matrix box;
	gmx_bool bOK;
#endif

	static char *kwlist[] = {
		"filename", "format", "mode", "units", NULL };

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|sss", kwlist, &filename, &str_type, &mode, &units))
		return -1;

	self->filename = (char*) malloc(strlen(filename) * sizeof(char));
	strcpy(self->filename, filename);
	if (mode == NULL)
		self->mode = 'r';
	else
		self->mode = mode[0];

	/* Set the enum symbol of the file format */
	if ( str_type != NULL ) {
		if      ( !strcmp(str_type,    "XYZ") ) self->type = XYZ;
		else if ( !strcmp(str_type, "MOLDEN") ) self->type = MOLDEN;
		else if ( !strcmp(str_type,    "GRO") ) self->type = GRO;
		else if ( !strcmp(str_type,    "XTC") ) self->type = XTC;
	}

	/* Guess the file format, if not given explicitly */
	if ( self->type == GUESS ) {
		strcpy(ext, filename + strlen(filename) - 4);
		if      ( !strcmp(ext, ".xyz") ) self->type = XYZ;
		else if ( !strcmp(ext, ".gro") ) self->type = GRO;
		else if ( !strcmp(ext, ".xtc") ) self->type = XTC;
		else {
			/* Extract the first line */
			if ( (test = fopen(filename, "r")) == NULL ) {
				PyErr_SetFromErrno(PyExc_IOError);
				return -1; }
			if ( (line = readline(test)) == NULL ) {
				return -1; }
			make_lowercase(line);
			stripline(line);
			fclose(test);

			/* Perhaps it's Molden format? */
			if ( !strcmp(line, "[molden format]") ) self->type = MOLDEN;

			free(line);
		}
	}

	if (units == NULL) {
		switch(self->type) {
			case XYZ:
			case MOLDEN:
				self->units = ANGS;
				break;
			case GRO:
			case XTC:
				self->units = NM;
				break;
		}
	} else {
		if      (!strcmp(units, "angs")) self->units = ANGS;
		else if (!strcmp(units, "bohr")) self->units = BOHR;
		else if (!strcmp(units,   "nm")) self->units = NM;
	}

	/* Open the coordinate file */
	switch(self->type) {
		case XYZ:
		case GRO:
		case MOLDEN:
			if ( (self->fd = fopen(filename, "r")) == NULL ) {
				PyErr_SetFromErrno(PyExc_IOError);
				return -1; }
			break;
		case XTC:
#ifdef HAVE_GROMACS
			if( (self->xd = open_xtc(filename, "r")) == NULL) {
				PyErr_SetString(PyExc_IOError, "Error opening XTC file");
				return -1; }
#else
			PyErr_SetString(PyExc_SystemError, "The module has to be compiled with gromacs support to handle XTC files");
			return -1;
#endif
			break;
	}

	/* Router */
	switch(self->type) {
		case XYZ:
			if (read_topo_from_xyz(self) == -1) return -1;
			rewind(self->fd);
			self->filePosition1 = ftell(self->fd);
			self->filePosition2 = self->filePosition1;
			break;
		case MOLDEN:
			if (read_topo_from_molden(self) == -1) return -1;
			rewind(self->fd);
			/* read_topo_from_molden sets file position accordingly */
			//self->filePosition1 = ftell(self->fd);
			//self->filePosition2 = self->filePosition1;
			break;
		case GRO:
			if (read_topo_from_gro(self) == -1) return -1;
			rewind(self->fd);
			self->filePosition1 = ftell(self->fd);
			self->filePosition2 = self->filePosition1;
			break;
		case XTC:
#ifdef HAVE_GROMACS
			if (!read_first_xtc(self->xd, &(self->nofatoms), &step, &time,
                                box, &(self->xtcCoord), &prec, &bOK) && bOK) {
				PyErr_SetString(PyExc_IOError, "Error reading first frame");
				return -1; }
#endif
			break;
		/* If the file format is GUESS or different,
		   it means we've failed to guess :-(        */
		case GUESS:
		default:
			PyErr_SetString(PyExc_ValueError, "Unsupported file format");
			return -1;
	}

	return 0;
}



static PyObject *Trajectory_read(Trajectory *self) {

	PyObject *py_result = NULL;

	switch(self->type) {

		case XYZ:
			py_result = read_frame_from_xyz(self);
			if (py_result == Py_None) Py_RETURN_NONE;
			self->filePosition1 = ftell(self->fd);
			self->filePosition1 = self->filePosition2;
			self->lastFrame += 1;
			break;

		case GRO:
			py_result = read_frame_from_gro(self);
			if (py_result == Py_None) Py_RETURN_NONE;
			self->filePosition1 = ftell(self->fd);
			self->filePosition1 = self->filePosition2;
			self->lastFrame += 1;
			break;

		case MOLDEN:
			if (self->moldenStyle == MLATOMS)
				py_result = read_frame_from_molden_atoms(self);
			else
				py_result = read_frame_from_molden_geometries(self);
			if (py_result == Py_None) Py_RETURN_NONE;
			self->lastFrame += 1;
			break;

		case XTC:
			py_result = read_frame_from_xtc(self);
			break;

		default:
			break;
	}

	return py_result;

}


static PyObject* Trajectory_repr(Trajectory *self) {
	PyObject* str;
	char format[10];
	int len;

	switch(self->type) {
		case XYZ:
			strcpy(format,    "XYZ"); break;
		case MOLDEN:
			strcpy(format, "MOLDEN"); break;
		case GRO:
			strcpy(format,    "GRO"); break;
		case XTC:
			strcpy(format,    "XTC"); break;
		default:
			strcpy(format,       ""); break;
	}
	len = strlen(self->filename) + 50;
	str = PyString_FromFormat("Trajectory('%s', format='%s', mode='%c')",
								self->filename, format, self->mode);
	return str;
}



static PyMemberDef Trajectory_members[] = {
    {"symbols", T_OBJECT_EX, offsetof(Trajectory, symbols), 0,
     "Symbols of atoms"},
    {"atomicnumbers", T_OBJECT_EX, offsetof(Trajectory, atomicnumbers), 0,
     "Atomic numbers"},
    {"nofatoms", T_INT, offsetof(Trajectory, nofatoms), 0,
     "Number of atoms"},
    {"nofframes", T_INT, offsetof(Trajectory, nofframes), 0,
     "Number of frames"},
    {NULL}  /* Sentinel */
};


static PyMethodDef Trajectory_methods[] = {
    {"read", (PyCFunction)Trajectory_read, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"Trajectory.read()\n"
		"\n"
		"Read next frame from trajectory. Returns a dictionary:\n"
		"\n"
		"coordinates (ndarray)\n"
		"step (int)\n"
		"time (float)\n"
		"box (ndarray) [a, b, c, alpha, beta, gamma]\n"
		"\n" },
	{NULL}  /* Sentinel */
};


PyTypeObject TrajectoryType = {
	PyObject_HEAD_INIT(NULL)
	0,                         /*ob_size*/
	"moltools.Trajectory",       /*tp_name*/
	sizeof(Trajectory),          /*tp_basicsize*/
	0,                         /*tp_itemsize*/
	(destructor)Trajectory_dealloc, /*tp_dealloc*/
	0,                         /*tp_print*/
	0,                         /*tp_getattr*/
	0,                         /*tp_setattr*/
	0,                         /*tp_compare*/
	(reprfunc)Trajectory_repr, /*tp_repr*/
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
	"Trajectory class. Implements reading of trajectories from various "
	"file formats. Coordinates are stored in numpy array. "
	"Creating an instance:\n"
	"traj = Trajectory(filename, format='GUESS', mode='r', units='angs')\n"
	"Available formats include: XYZ, GRO, MOLDEN, XTC - guessed if not specified.\n"
	"Mode: 'r' (default), 'w', 'a'.\n"
	"Units: 'angs' (default), 'bohr', 'nm'.\n",           /* tp_doc */
	0,		               /* tp_traverse */
	0,		               /* tp_clear */
	0,		               /* tp_richcompare */
	0,		               /* tp_weaklistoffset */
	0,		               /* tp_iter */
	0,		               /* tp_iternext */
	Trajectory_methods,             /* tp_methods */
	Trajectory_members,             /* tp_members */
	0,                         /* tp_getset */
	0,                         /* tp_base */
	0,                         /* tp_dict */
	0,                         /* tp_descr_get */
	0,                         /* tp_descr_set */
	0,                         /* tp_dictoffset */
	(initproc)Trajectory_init,      /* tp_init */
	0,                         /* tp_alloc */
	Trajectory_new,                 /* tp_new */
};




int read_topo_from_xyz(Trajectory *self) {

	int nofatoms, pos, idx;
	char *buffer, *buffpos, *token;
	int *anum;
	extern Element element_table[];

	npy_intp dims[2];
	PyObject *val;

	/* Read number of atoms */
	if( (buffer = readline(self->fd)) == NULL) {
		return -1; }
	if( sscanf(buffer, "%d", &nofatoms) != 1 ) {
		PyErr_SetString(PyExc_IOError, "Incorrect atom number");
		return -1; }

	/* Read the comment line */
	if( (buffer = readline(self->fd)) == NULL) {
		return -1; }

	/* Get rid of Py_None in self->symbols */
	Py_DECREF(Py_None);
	self->symbols = PyList_New(nofatoms);

	anum = (int*) malloc(nofatoms * sizeof(int));
	if(anum == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return -1; }

	/* Atom loop */
	for(pos = 0; pos < nofatoms; pos++) {

		/* Get the whole line */
		if((buffer = readline(self->fd)) == NULL) {
			return -1; }
		buffer[strlen(buffer)-1] = '\0';
		buffpos = buffer;

		/* Read symbol */
		token = strtok(buffpos, " ");
		val = Py_BuildValue("s", token);
		PyList_SetItem(self->symbols, pos, val);

		idx = getElementIndexBySymbol(token);
		if (idx == -1)
			anum[pos] = -1;
		else
			anum[pos] = element_table[idx].number;

		/* Free the line buffer */
		free(buffer);
	}

	/* Add atomic numbers to the dictionary */
	dims[0] = nofatoms;
	dims[1] = 1;
	Py_DECREF(self->atomicnumbers);
	self->atomicnumbers = PyArray_SimpleNewFromData(1, dims, NPY_INT, (int*) anum);

	self->nofatoms = nofatoms;

	return 0;
}



/* Read molden file */
int read_topo_from_molden(Trajectory *self) {

	char *line = NULL;
	char *buffpos, *token;
	char **line_store;
	int n_geo, i, nat;
	int *anum;
	unsigned short int section_present = 0;

	npy_intp dims[2];
	PyObject *val, *py_anum;

	/* Loop over the lines */

	do {

		free(line);
		if ((line = readline(self->fd)) == NULL) return -1;

		stripline(line);
		make_lowercase(line);

		/* This is a start of a new block */
		if ( line[0] == '[' ) {

			/* Several geometries */
			if ( !strcmp(line, "[geometries] (xyz)" ) ) {

				/* Set filePosition1 for reading frames */
				self->filePosition1 = ftell(self->fd);
				self->moldenStyle = MLGEOMETRY;

				/*if (section_present) {
					PyErr_SetString(PyExc_RuntimeError, "Multiple geometry/atom sections");
					return -1; } */
				section_present = 1;
				read_topo_from_xyz(self);
				nat = self->nofatoms;
				break;
			}

			/* Section 'atoms' present - this is one-geometry file */
			else if ( !strcmp(line, "[atoms] angs") || !strcmp(line, "[atoms] au") ) {

				/* Set filePosition1 for reading frames */
				self->filePosition1 = ftell(self->fd);
				self->moldenStyle = MLATOMS;

				if ( !strcmp(line, "[atoms] angs") )
					self->units = ANGS;
				else if ( !strcmp(line, "[atoms] au") )
					self->units = BOHR;
				else {
					PyErr_SetString(PyExc_RuntimeError, "Unrecognized units");
					return -1;
				}

				anum = (int*) malloc(nat * sizeof(int));
				if(anum == NULL) {
					PyErr_SetFromErrno(PyExc_MemoryError);
					return -1; }

				/*if (section_present) {
					PyErr_SetString(PyExc_RuntimeError, "Multiple geometry/atom sections");
					return -1; }*/
				section_present = 1;
				free(line);
				if ( (line = readline(self->fd)) == NULL ) return -1;
				stripline(line);

				/* We don't know how many atoms are there, so we have to
                   store the lines. */
				nat = 0;
				line_store = (char**) malloc(10 * sizeof(char*));
				if ( line_store == NULL ) {
					PyErr_NoMemory();
					return -1; }
				while ( line[0] != '[' ) {
					line_store[nat++] = line;
					if ( nat % 10 == 0 ) {
						line_store = realloc(line_store, (nat + 10) * sizeof(char*));
						if( line_store == NULL ) {
							PyErr_NoMemory();
						return -1; }
					}
					if ( (line = readline(self->fd)) == NULL ) return -1;
					stripline(line);
				}

				/* Get rid of Py_None in self->symbols */
				Py_DECREF(Py_None);
				self->symbols = PyList_New(nat);

				/* Loop over atoms */
				for ( i = 0; i < nat; i++ ) {

					buffpos = line_store[i];
					token = strtok(buffpos, " ");
					val = Py_BuildValue("s", token);
					PyList_SetItem(self->symbols, i, val);

					token = strtok(NULL, " ");
					/* not used */

					token = strtok(NULL, " ");
					anum[i] = atoi(token);

					/* Get rid of the line. */
					free(line_store[i]);

				}

				/* Free the stored line pointers. */
				free(line_store);

				/* Add atomic numbers to the dictionary */
				dims[0] = nat;
				dims[1] = 1;
				Py_DECREF(self->atomicnumbers);
				self->atomicnumbers = PyArray_SimpleNewFromData(1, dims, NPY_INT, (int*) anum);

				break;
				/* This is to avoid reading the next line! */
				//continue;
			}

		}

	} while ( strlen(line) );
	free(line);


	rewind(self->fd);
	self->filePosition2 = -1;
	do {
		if ( (line = readline(self->fd)) == NULL ) return -1;
		stripline(line);
		make_lowercase(line);
		/* This is a start of a new block */
		if ( line[0] == '[' ) {
			/* Several geometries */
			if ( !strcmp(line, "[geoconv]" ) ) {
				free(line);
				/* consume the line 'energy' */
				if ( (line = readline(self->fd)) == NULL ) return -1;
				self->filePosition2 = ftell(self->fd);
			}
			/* Number of geometries present in the file */
			if ( !strcmp(line, "[n_geo]") ) {
				free(line);
				if ( (line = readline(self->fd)) == NULL ) return -1;
				stripline(line);
				self->nofframes = atoi(line);
			}
		}
	} while (strlen(line) != 0);
	free(line);

	if (!section_present) {
		PyErr_SetString(PyExc_RuntimeError, "geometry/atom section missing");
		return -1; }

	self->nofatoms = nat;
	return 0;
}




int read_topo_from_gro(Trajectory *self) {

	int nofatoms, pos;
	char *buffer;
	char symbuf[100];
	int *resid;

	PyObject *val;

	/* Read the comment line */
	if((buffer = readline(self->fd)) == NULL) {
		return -1; }
	buffer[strlen(buffer)-1] = '\0';
	stripline(buffer);

	/* Read number of atoms */
	if( (buffer = readline(self->fd)) == NULL) {
		return -1; }
	if( sscanf(buffer, "%d", &nofatoms) != 1 ) {
		PyErr_SetString(PyExc_IOError, "Incorrect atom number");
		return -1; }

	/* Get rid of Py_None in self->symbols etc. */
	Py_DECREF(Py_None);
	self->symbols = PyList_New(nofatoms);
	Py_DECREF(Py_None);
	self->resids = PyList_New(nofatoms);
	Py_DECREF(Py_None);
	self->resnames = PyList_New(nofatoms);

	resid = (int*) malloc(nofatoms * sizeof(int));
	if(resid == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return -1; }

	/* Atom loop */
	for(pos = 0; pos < nofatoms; pos++) {

		/* Get the whole line */
		if((buffer = readline(self->fd)) == NULL) {
			return -1; }

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
		PyList_SetItem(self->resnames, pos, val);

		/* Read atom name */
		strncpy(symbuf, buffer+10, 5);
		symbuf[5] = '\0';
		stripline(symbuf);
		val = Py_BuildValue("s", symbuf);
		PyList_SetItem(self->symbols, pos, val);

		/* Free the line buffer */
		free(buffer);
	}

	self->nofatoms = nofatoms;
	return 0;
}



/* This function is used by other readers, like    *
 * read_frame_from_molden_geometries for instance, *
 * so be careful with implementation.              */

PyObject *read_frame_from_xyz(Trajectory *self) {

	PyObject *py_result, *py_coord, *py_charges, *val, *key;
	char *buffer, *buffpos, *token;
	int pos, nat;
	float factor;
	float *xyz, *charges;
	unsigned short int charges_present;
	npy_intp dims[2];

	switch(self->units) {
		case ANGS: factor = 1.0; break;
		case NM: factor = 10.0; break;
		case BOHR: factor = BOHRTOANGS; break;
	}
	/* Create the dictionary that will be returned */
	py_result = PyDict_New();

	/* Read number of atoms */
	if( (buffer = readline(self->fd)) == NULL) {
		PyErr_SetFromErrno(PyExc_IOError);
		return NULL; }

	if (sscanf(buffer, "%d", &nat) != 1) {
		/* Possibly end of file */
		free(buffer);
		Py_DECREF(py_result);
	    Py_RETURN_NONE;
	}
	free(buffer);

	if (nat != self->nofatoms) {
		PyErr_SetString(PyExc_IOError, "Error reading number of atoms");
		return NULL; }

	/* Read the comment line */
	if((buffer = readline(self->fd)) == NULL) {
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
	xyz = (float*) malloc(3 * self->nofatoms * sizeof(float));
	if(xyz == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return NULL; }
	charges = (float*) malloc(self->nofatoms * sizeof(float));
	if(charges == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return NULL; }
	charges_present = 0;

	/* Atom loop */
	for(pos = 0; pos < self->nofatoms; pos++) {

		/* Get the whole line */
		if((buffer = readline(self->fd)) == NULL) {
			PyErr_SetFromErrno(PyExc_IOError);
			return NULL; }
		buffer[strlen(buffer)-1] = '\0';
		buffpos = buffer;

		/* Read symbol */
		token = strtok(buffpos, " ");

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

	/* Add coordinates to the dictionary */
	dims[0] = self->nofatoms;
	dims[1] = 3;
	py_coord = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, (float*) xyz);
	/***************************************************************
	 * Do not free the raw array! It will be still used by Python! *
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


PyObject *read_frame_from_molden_atoms(Trajectory *self) {

	char *line;
	char *token, *buffpos;
	int i;
	float *xyz;
	npy_intp dims[2];
	PyObject *py_result, *key, *val, *py_geom;

	/* Prepare dictionary */

	py_result = PyDict_New();

	xyz = (float*) malloc(3 * self->nofatoms * sizeof(float));
	if(xyz == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return NULL; }

	/* Loop over the lines */
	fseek(self->fd, self->filePosition1, SEEK_SET);
	for (i = 0; i < self->nofatoms; i++ ) {

		line = readline(self->fd);
		stripline(line);

		/* If not a letter, then we have run out of frames */
		if (!isalpha(line[0])) {
			Py_DECREF(py_result);
			free(xyz);
			free(line);
			Py_RETURN_NONE;
		}

		buffpos = line;
		token = strtok(buffpos, " ");

		token = strtok(NULL, " ");
		/* not used */

		token = strtok(NULL, " ");
		/* atomic number - not used here */

		token = strtok(NULL, " ");
		xyz[3*i+0] = atof(token);

		token = strtok(NULL, " ");
		xyz[3*i+1] = atof(token);

		token = strtok(NULL, " ");
		xyz[3*i+2] = atof(token);

		if (self->units == BOHR) {
			xyz[3*i+0] *= BOHRTOANGS;
			xyz[3*i+1] *= BOHRTOANGS;
			xyz[3*i+2] *= BOHRTOANGS;
		}

		free(line);

	}
	self->filePosition1 = ftell(self->fd);

	/* Add coordinates to the dictionary */
	dims[0] = self->nofatoms;
	dims[1] = 3;
	py_geom = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, (float*) xyz);
	key = PyString_FromString("coordinates");
	PyDict_SetItem(py_result, key, py_geom);
	Py_DECREF(key);
	Py_DECREF(py_geom);

	return py_result;

}




PyObject *read_frame_from_molden_geometries(Trajectory *self) {

	char *line;
	PyObject *py_result, *key, *val;

	fseek(self->fd, self->filePosition1, SEEK_SET);
	py_result = read_frame_from_xyz(self);
	if (py_result == Py_None) {
		Py_RETURN_NONE;
	}
	self->filePosition1 = ftell(self->fd);

	/* If filePosition2 == -1, the GEOCONV section was not found */
	if (self->filePosition2 >= 0) {
		fseek(self->fd, self->filePosition2, SEEK_SET);
		line = readline(self->fd);
		stripline(line);
		key = PyString_FromString("energy");
		val = Py_BuildValue("d", atof(line));
		PyDict_SetItem(py_result, key, val);
		Py_DECREF(key);
		Py_DECREF(val);
		free(line);
		self->filePosition2 = ftell(self->fd);
	}

	return py_result;

}


PyObject *read_frame_from_gro(Trajectory *self) {

	int nat, pos;
	char *buffer;
	float *xyz, *vel, *box;
	unsigned short int velocities_present = 0;

	npy_intp dims[2];

	PyObject *key, *val, *py_result, *py_coord, *py_vel, *py_box;


	/* Create the dictionary that will be returned */
	py_result = PyDict_New();

	/* Read the comment line */
	if((buffer = readline(self->fd)) == NULL) {
		PyErr_SetFromErrno(PyExc_IOError);
		return NULL; }

	if (buffer[0] == '\0') {
		/* EOF reached */
		free(buffer);
		Py_RETURN_NONE;
	}
	
	buffer[strlen(buffer)-1] = '\0';
	stripline(buffer);

	val = Py_BuildValue("s", buffer);
	free(buffer);
	key = PyString_FromString("comment");
	PyDict_SetItem(py_result, key, val);
	Py_DECREF(key);
	Py_DECREF(val);

	/* Read number of atoms */
	if( (buffer = readline(self->fd)) == NULL) {
		PyErr_SetFromErrno(PyExc_IOError);
		return NULL; }
	if( sscanf(buffer, "%d", &nat) != 1 || nat != self->nofatoms) {
		PyErr_SetString(PyExc_IOError, "Incorrect atom number");
		return NULL; }
	free(buffer);


	/* Set-up the raw arrays for coordinates and charges */
	xyz = (float*) malloc(3 * self->nofatoms * sizeof(float));
	if(xyz == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return NULL; }
	vel = (float*) malloc(3 * self->nofatoms * sizeof(float));
	if(vel == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return NULL; }
	box = (float*) malloc(9 * sizeof(float));
	if(box == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return NULL; }

	/* Atom loop */
	for(pos = 0; pos < self->nofatoms; pos++) {

		/* Get the whole line */
		if((buffer = readline(self->fd)) == NULL) {
			PyErr_SetFromErrno(PyExc_IOError);
			return NULL; }
		if(pos == 0 && strlen(buffer) > 50) velocities_present = 1;

		/* Read coordinates */
		xyz[3*pos + 0] = strPartFloat(buffer, 20, 8) * 10.0;
		xyz[3*pos + 1] = strPartFloat(buffer, 28, 8) * 10.0;
		xyz[3*pos + 2] = strPartFloat(buffer, 36, 8) * 10.0;

		/* Read velocities */
		if(velocities_present) {
			vel[3*pos + 0] = strPartFloat(buffer, 44, 8);
			vel[3*pos + 1] = strPartFloat(buffer, 52, 8);
			vel[3*pos + 2] = strPartFloat(buffer, 60, 8);
		}

		/* Free the line buffer */
		free(buffer);
	}

	/* Get the cell line */
	if((buffer = readline(self->fd)) == NULL) {
		PyErr_SetFromErrno(PyExc_IOError);
		return NULL; }
	box[3*0 + 0] = strPartFloat(buffer,  0, 10) * 10.0;
	box[3*1 + 1] = strPartFloat(buffer, 10, 10) * 10.0;
	box[3*2 + 2] = strPartFloat(buffer, 20, 10) * 10.0;
	if (strlen(buffer) > 31) box[3*0 + 1] = strPartFloat(buffer, 30, 10) * 10.0;
	else                     box[3*0 + 1] = 0.0;
	if (strlen(buffer) > 41) box[3*0 + 2] = strPartFloat(buffer, 40, 10) * 10.0;
	else                     box[3*0 + 2] = 0.0;
	if (strlen(buffer) > 51) box[3*1 + 0] = strPartFloat(buffer, 50, 10) * 10.0;
	else                     box[3*1 + 0] = 0.0;
	if (strlen(buffer) > 61) box[3*1 + 2] = strPartFloat(buffer, 60, 10) * 10.0;
	else                     box[3*1 + 2] = 0.0;
	if (strlen(buffer) > 71) box[3*2 + 0] = strPartFloat(buffer, 70, 10) * 10.0;
	else                     box[3*2 + 0] = 0.0;
	if (strlen(buffer) > 81) box[3*2 + 1] = strPartFloat(buffer, 80, 10) * 10.0;
	else                     box[3*2 + 1] = 0.0;
	free(buffer);

	/* Add coordinates to the dictionary */
	dims[0] = self->nofatoms;
	dims[1] = 3;
	py_coord = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, (float*) xyz);
	/***************************************************************
	 * Do not free the raw array! It will be still used by Python! *
     ***************************************************************/

	key = PyString_FromString("coordinates");
	PyDict_SetItem(py_result, key, py_coord);
	Py_DECREF(key);
	Py_DECREF(py_coord);


	/* Add velocities to the dictionary */
	if(velocities_present) {
		py_vel = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, (float*) vel);
		/***************************************************************
		 * Do not free the raw array! It will be still used by Python! *
    	 ***************************************************************/

		key = PyString_FromString("velocities");
		PyDict_SetItem(py_result, key, py_vel);
		Py_DECREF(key);
		Py_DECREF(py_vel);
	} else
		free(vel);

	dims[0] = 3;
	dims[1] = 3;
	py_box = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, (float*) box);
	key = PyString_FromString("box");
	PyDict_SetItem(py_result, key, py_box);
	Py_DECREF(key);
	Py_DECREF(py_box);

	return py_result;

}


PyObject *read_frame_from_xtc(Trajectory *self) {
	Py_RETURN_NONE;
}

