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


/* Methods of the Trajectory object */

static void Trajectory_dealloc(Trajectory* self)
{
	Py_XDECREF(self->symbols);
	Py_XDECREF(self->atomicnumbers);
	self->ob_type->tp_free((PyObject*)self);
	free(self->filename);
	switch(self->type) {
		case XYZ:
		case MOLDEN:
		case GRO:
		case FRAC:
			if (self->fd != NULL) fclose(self->fd);
			break;
#ifdef HAVE_GROMACS
		case XTC:
			if (self->xd != NULL) close_xtc(self->xd);
#endif
	}
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

		Py_INCREF(Py_None);
	    self->symbols = Py_None;
		
		Py_INCREF(Py_None);
	    self->atomicnumbers = Py_None;
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

	static char *kwlist[] = {
		"filename", "format", "mode", NULL };

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|ss", kwlist, &filename, &str_type, &mode))
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
		else if ( !strcmp(str_type,   "FRAC") ) self->type = FRAC;
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

	/* Open the coordinate file */
	switch(self->type) {
		case XYZ:
		case FRAC:
		case MOLDEN:
		case GRO:
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
		case FRAC:
			self->nofatoms = read_topo_from_xyz(self);
			rewind(self->fd);
			break;
		case MOLDEN:
			self->nofatoms = read_topo_from_molden(self);
			rewind(self->fd);
			break;
		case GRO:
			self->nofatoms = read_topo_from_gro(self);
			rewind(self->fd);
			break;
		case XTC:
			break;
		/* If the file format is GUESS or different,
		   it means we've failed to guess :-(        */
		case GUESS:
		default:
			PyErr_SetString(PyExc_ValueError, "Unsupported file format");
			return -1;
	}

	if ( self->nofatoms == -1 ) {
		return -1; }

	return 0;
}


//static PyObject *Trajectory_read(Trajectory *self, PyObject *args, PyObject *kwds) {
//
//	const char *filename;
//	char *line;
//	char *unit = NULL;
//	enum { GUESS, XYZ, MOLDEN, FRAC, GRO } type = GUESS;
//	char *str_type = NULL;
//	float factor;
//	FILE *fd, *test;
//	long int fpos;
//	struct stat fst;
//	char ext[5];
//	int nat;
//	char buffer[1000];
//
//	static char *kwlist[] = {"filename", "format", "unit", NULL};
//
//	/* Process the arguments */
//	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|ss",
//		kwlist, &filename, &str_type, &unit)) return NULL;
//
//	/* Get the unit measure for coordinates and set the factor accordingly */
//	if (unit == NULL) factor = 1.0;
//	else {
//		make_lowercase(unit);
//		if ( !strcmp(unit, "angs") ) factor = 1.0;
//		else if ( !strcmp(unit, "bohr") ) factor = BOHR;
//		else if ( !strcmp(unit, "nm") ) factor = 10.0;
//		else {
//			PyErr_SetString(PyExc_ValueError, "Unrecognized measure unit name");
//			return NULL; }
//	}
//
//	/* Set the enum symbol of the file format */
//	if ( str_type != NULL ) {
//		if      ( !strcmp(str_type,    "XYZ") ) type = XYZ;
//		else if ( !strcmp(str_type, "MOLDEN") ) type = MOLDEN;
//		else if ( !strcmp(str_type,   "FRAC") ) type = FRAC;
//		else if ( !strcmp(str_type,    "GRO") ) type = GRO;
//	}
//
//	/* Guess the file format, if not given explicitly */
//	if ( type == GUESS ) {
//		strcpy(ext, filename + strlen(filename) - 4);
//		if      ( !strcmp(ext, ".xyz") ) type = XYZ;
//		else if ( !strcmp(ext, ".gro") ) type = GRO;
//		else {
//			/* Extract the first line */
//			if ( (test = fopen(filename, "r")) == NULL ) {
//				PyErr_SetFromErrno(PyExc_IOError);
//				return NULL; }
//			if ( (line = readline(test)) == NULL ) {
//				PyErr_SetFromErrno(PyExc_IOError);
//				return NULL; }
//			make_lowercase(line);
//			stripline(line);
//			fclose(test);
//
//			/* Perhaps it's Molden format? */
//			if ( !strcmp(line, "[molden format]") ) type = MOLDEN;
//
//			free(line);
//		}
//	}
//
//	/* Open the coordinate file */
//	if ( (fd = fopen(filename, "r")) == NULL ) {
//		PyErr_SetFromErrno(PyExc_IOError);
//		return NULL; }
//
//	/* Router */
//	switch(type) {
//		case XYZ:
//			nat = read_frame_from_xyz(fd, factor, self);
//			if ( nat == -1 ) return NULL;
//			/* Support for multiframe XYZ: check the current position in
// 			 * the file; if we are near the end, finish; else, continue
// 			 * reading. The criteria is that there should be no more than
// 			 * 3 bytes left. */
//			fpos = ftell(fd);
//			if ( fpos == -1 ) {
//				PyErr_SetFromErrno(PyExc_IOError);
//				return NULL; }
//			if ( fstat(fileno(fd), &fst) == -1 ) {
//				PyErr_SetFromErrno(PyExc_IOError);
//				return NULL; }
//			if ( fst.st_size - fpos <= 3 ) {
//				break;
//			}
//			while ( fst.st_size - fpos > 3 ) {
//				nat = read_frame_from_xyz(fd, factor, self);
//				if ( nat == -1 ) return NULL;
//				fpos = ftell(fd);
//				if ( fpos == -1 ) {
//					PyErr_SetFromErrno(PyExc_IOError);
//					return NULL; }
//			}
//			break;
//		case MOLDEN:
//			//nat = read_frame_from_molden(fd, self);
//			if ( nat == -1 ) return NULL;
//			break;
//		case FRAC:
//			//nat = read_frame_from_fractional(fd, self);
//			if ( nat == -1 ) return NULL;
//			break;
//		case GRO:
//			//nat = read_frame_from_gro(fd, self);
//			if ( nat == -1 ) return NULL;
//			break;
//		/* If the file format is GUESS or different,
//		   it means we've failed to guess :-(        */
//		case GUESS:
//		default:
//			PyErr_SetString(PyExc_ValueError, "Unsupported file format");
//			return NULL;
//	}
//
//	fclose(fd);
//	Py_RETURN_NONE;
//}



static PyObject* Trajectory_repr(Trajectory *self) {
	PyObject* str;
	char format[10];
	int len;

	switch(self->type) {
		case XYZ:
			strcpy(format,    "XYZ"); break;
		case MOLDEN:
			strcpy(format, "MOLDEN"); break;
		case FRAC:
			strcpy(format,   "FRAC"); break;
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
     "Symbols of atoms."},
    {"atomicnumbers", T_OBJECT_EX, offsetof(Trajectory, atomicnumbers), 0,
     "Symbols of atoms."},
    {"nofatoms", T_INT, offsetof(Trajectory, nofatoms), 0,
     "Number of atoms."},
    {NULL}  /* Sentinel */
};


static PyMethodDef Trajectory_methods[] = {
//    {"read", (PyCFunction)Trajectory_read, METH_VARARGS | METH_KEYWORDS,
//		"\n"
//		"Trajectory.read(filename [, unit ] )\n"
//		"\n"
//		"Read frames from coordinate file and append to Trajectory.\n"
//		"Currently, supports only XYZ and Molden formats. In some file\n"
//		"formats the unit can be chosen between unit = angs, bohr, nm.\n"
//		"\n"
//		"XYZ:    supports extended files (charges in the fifth column)\n"
//		"        as well as standard files; coordinates are assumed to be\n"
//		"        in Angstroms. Also supports multiframe files.\n"
//		"\n"
//		"Molden: supports groups N_GEO, GEOCONV (energies only),\n"
//		"        GEOMETRIES (XYZ only), ATOMS (Angstroms only).\n"
//		"\n" },
    /*{"COM", (PyCFunction)Trajectory_COM, METH_VARARGS,
		"\n"
		"COM()\n"
		"\n"
		"Returns the center of mass as numpy array.\n"
		"\n" },*/
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
	"Trajectory class. Implements reading of trajectories from various"
	"file formats. Coordinates are stored in numpy array.",           /* tp_doc */
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

/* End of methods of the Trajectory object */


int read_topo_from_xyz(Trajectory *self) {

	int nofatoms, pos;
	char *buffer, *buffpos, *token;

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

		/* Free the line buffer */
		free(buffer);
	}

	return nofatoms;
}


/* Read molden file */
int read_topo_from_molden(Trajectory *self) {

	char *line;
	char *buffpos, *token;
	char **line_store;
	int n_geo, i, nat;
	int *anum;
	unsigned short int section_present = 0;

	PyObject *val;

	/* Loop over the lines */

	if ( (line = readline(self->fd)) == NULL) return -1;

	while ( strlen(line) ) {

		stripline(line);
		make_lowercase(line);

		/* This is a start of a new block */
		if ( line[0] == '[' ) {

			/* Several geometries */
			if ( !strcmp(line, "[geometries] (xyz)" ) ) {
				/*if (section_present) {
					PyErr_SetString(PyExc_RuntimeError, "Multiple geometry/atom sections");
					return -1; } */
				section_present = 1;
				nat = read_topo_from_xyz(self);
				break;
			}

			/* Section 'atoms' present - this is one-geometry file */
			else if ( !strcmp(line, "[atoms] angs" ) ) {

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

					/* Get rid of the line. */
					free(line_store[i]);

				}

				/* Free the stored line pointers. */
				free(line_store);

				break;
				/* This is to avoid reading the next line! */
				//continue;
			}

		}

		if ( (line = readline(self->fd)) == NULL ) return -1;

	}
	if (!section_present) {
		PyErr_SetString(PyExc_RuntimeError, "geometry/atom section missing");
		return -1; }

	return nat;
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

	/* Get rid of Py_None in self->symbols */
	Py_DECREF(Py_None);
	self->symbols = PyList_New(nofatoms);

	/* Atom loop */
	for(pos = 0; pos < nofatoms; pos++) {

		/* Get the whole line */
		if((buffer = readline(self->fd)) == NULL) {
			return -1; }

		/* Read atom name */
		strncpy(symbuf, buffer+10, 5);
		symbuf[5] = '\0';
		stripline(symbuf);
		val = Py_BuildValue("s", symbuf);
		PyList_SetItem(self->symbols, pos, val);

		/* Free the line buffer */
		free(buffer);
	}

	return nofatoms;
}



//int read_frame_from_xyz(FILE *fd, float factor, Trajectory *self) {
//
//	int nofatoms, pos, size;
//	char *buffer, *buffpos, *token;
//	float *ptr;
//	/* The following enum is used to decide if there are any charges to read.
//	 * If equal to NONE, it means that previous frames had no charges, so we
//	 * raise an exception if some are present.
//	 * If equal to REQUIRED, it means that previous frames had charges, so we
//	 * require them.
//	 * If equal to POSSIBLE, then we don't know yet. On the first line it will
//	 * be switched to NONE or REQUIRED, depending if we find the fourth number. */
//	enum { NONE, POSSIBLE, REQUIRED } charges_section;
//
//	PyObject *val;
//
//	/* Read number of atoms */
//	if( (buffer = readline(fd)) == NULL) {
//		PyErr_SetFromErrno(PyExc_IOError);
//		return -1; }
//	if( sscanf(buffer, "%d", &nofatoms) != 1 ) {
//		PyErr_SetString(PyExc_IOError, "Failed to read number of atoms");
//		return -1; }
//	if( nofatoms != self->nofatoms ) {
//		PyErr_SetString(PyExc_RuntimeError, "Incorrect number of atoms");
//		return -1; }
//
//	/* Read the comment line */
//	if((buffer = readline(fd)) == NULL) {
//		PyErr_SetFromErrno(PyExc_IOError);
//		return -1; }
//	buffer[strlen(buffer)-1] = '\0';
//
//	/* Determine wether we are expecting to find charges. */
//
//	if ( self->charges_dim[0] )
//		/* Charges must be present */
//		charges_section = REQUIRED;
//	else if ( ! self->charges_dim && self->frames_dim )
//		/* We have some frames but no charges, so there shouldn't be any */
//		charges_section = NONE;
//	else {
//		/* Possibly there will be charges */
//		charges_section = POSSIBLE;
//	}
//
//	/* Set-up the raw array and PyArray for coordinates */
//
//	if ( self->frames_raw == NULL ) {
//		self->frames_dim[0] = 1;
//		ptr = (float*) malloc(3 * self->frames_dim[1] * sizeof(float));
//	} else {
//		Py_DECREF(self->frames);
//		self->frames_dim[0] += 1;
//		ptr = realloc(self->frames_raw, self->frames_dim[0] * self->frames_dim[1] * 3 * sizeof(float));
//	}
//	if( ptr == NULL) {
//		PyErr_NoMemory();
//		return -1; }
//	self->frames_raw = ptr;
//	self->frames = PyArray_SimpleNewFromData(3, self->frames_dim, NPY_FLOAT, self->frames_raw);
//
//	/* Set-up the raw array and PyArray for charges.
//	 * If charges are not present, self->charges will be set back
//	 * to None at the end of this function. */
//
//	if ( self->charges_raw == NULL ) {
//		self->charges_dim[0] = 1;
//		ptr = (float*) malloc(self->charges_dim[1] * sizeof(float));
//	} else {
//		Py_DECREF(self->charges);
//		self->charges_dim[0] += 1;
//		ptr = realloc(self->charges_raw, self->charges_dim[0] * self->charges_dim[1] * sizeof(float));
//	}
//	if(ptr == NULL) {
//		PyErr_NoMemory();
//		return -1; }
//	self->charges_raw = ptr;
//	self->charges = PyArray_SimpleNewFromData(2, self->charges_dim, NPY_FLOAT, self->charges_raw);
//
//	/* Atom loop */
//	for(pos = 0; pos < nofatoms; pos++) {
//
//		/* Get the whole line */
//		if((buffer = readline(fd)) == NULL) {
//			PyErr_SetFromErrno(PyExc_IOError);
//			return -1; }
//		buffer[strlen(buffer)-1] = '\0';
//		buffpos = buffer;
//
//		/* Read symbol */
//		token = strtok(buffpos, " ");
//
//		/* Read coordinates */
//		if ( (token = strtok(NULL, " ")) == NULL) {
//			PyErr_SetString(PyExc_IOError, "Missing coordinate");
//			return -1; }
//		self->frames_raw[3*nofatoms*(self->frames_dim[0] - 1) + 3*pos + 0] = atof(token) * factor;
//		if ( (token = strtok(NULL, " ")) == NULL) {
//			PyErr_SetString(PyExc_IOError, "Missing coordinate");
//			return -1; }
//		self->frames_raw[3*nofatoms*(self->frames_dim[0] - 1) + 3*pos + 1] = atof(token) * factor;
//		if ( (token = strtok(NULL, " ")) == NULL) {
//			PyErr_SetString(PyExc_IOError, "Missing coordinate");
//			return -1; }
//		self->frames_raw[3*nofatoms*(self->frames_dim[0] - 1) + 3*pos + 2] = atof(token) * factor;
//
//		/* Read charge, if present */
//		token = strtok(NULL, " ");
//		if ( token != NULL ) {
//
//			/* This is bad: didn't expect the charges! */
//			if ( charges_section == NONE ) {
//				PyErr_SetString(PyExc_RuntimeError, "Found unexpected charge data");
//				return -1;
//			}
//
//			if ( charges_section == POSSIBLE ) charges_section = REQUIRED;
//
//			self->charges_raw[nofatoms*(self->charges_dim[0] - 1) + pos] = atof(token);
//
//		} else {
//
//			/* This is bad: we were expecting charges here and found nothing */
//			if ( charges_section == REQUIRED ) {
//				PyErr_SetString(PyExc_RuntimeError, "Missing charges");
//				return -1;
//			}
//
//			if ( charges_section == POSSIBLE ) charges_section = NONE;
//
//		}
//
//		/* Free the line buffer */
//		free(buffer);
//	}
//
//	self->nframes += 1;
//
//	if ( charges_section == NONE ) {
//
//		Py_DECREF(self->charges);
//		free(self->charges_raw);
//		self->charges_raw = NULL;
//		Py_INCREF(Py_None);
//		self->charges = Py_None;
//	
//	} else if ( charges_section == POSSIBLE ) {
//
//			PyErr_SetString(PyExc_AssertionError, "Internal error when reading charges");
//			return -1;
//	
//	}
//
//	return nofatoms;
//}

