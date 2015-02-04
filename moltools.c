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

//#include "periodic_table.c"



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
	    if (self->rscale == NULL) {
	        Py_DECREF(self);
	        return NULL; }
	    
	    self->rho = PyArray_EMPTY(1, dims, NPY_DOUBLE, 0);
	    if (self->rscale == NULL) {
	        Py_DECREF(self);
	        return NULL; }
	    
	    self->rho2 = PyArray_EMPTY(1, dims, NPY_DOUBLE, 0);
	    if (self->rscale == NULL) {
	        Py_DECREF(self);
	        return NULL; }
	    
	    self->Z = PyArray_EMPTY(1, dims, NPY_DOUBLE, 0);
	    if (self->rscale == NULL) {
	        Py_DECREF(self);
	        return NULL; }
	    
	    self->Z2 = PyArray_EMPTY(1, dims, NPY_DOUBLE, 0);
	    if (self->rscale == NULL) {
	        Py_DECREF(self);
	        return NULL; }
	    
	    self->F = PyArray_EMPTY(1, dims, NPY_DOUBLE, 0);
	    if (self->rscale == NULL) {
	        Py_DECREF(self);
	        return NULL; }
	    
	    self->F2 = PyArray_EMPTY(1, dims, NPY_DOUBLE, 0);
	    if (self->rscale == NULL) {
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
	//printf("%s\n", text);
	free(text);
	text = readline(eamfile);
	//printf("%s\n", text);
	free(text);
	fscanf(eamfile, "%d %lf %d %lf %*f", &nr, &dr, &nd, &dd);

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

static PyMethodDef EAMff_methods[] = {
    /*{"energy", (PyCFunction)EAMff_energy, METH_VARARGS,
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
		"\n" },*/
	{NULL}  /* Sentinel */
};

static PyTypeObject EAMffType = {
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



static PyObject *exposed_energy(PyObject *self, PyObject *args, PyObject *kwds) {
	FFType ff = FF_NONE;
	char *ff_type;
	float box[3] = { 0.0, 0.0, 0.0 };
	unsigned short int pbc = 0;
	double energy;
	int i, type;
	const char *ff_type_map[] = { "12-6", NULL };

	static char *kwlist[] = {
		"coordinates", "types", "ff_type", "ff", "box", NULL };

	PyObject *py_types, *py_ff;
	PyArrayObject *py_coords, *py_box = NULL;
	PyObject *py_result = NULL;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|sO!O!", kwlist,
			&PyArray_Type, &py_coords,
			&PyList_Type, &py_types,
			&ff_type,
			&PyList_Type, &py_ff,
			&PyArray_Type, &py_box))
		return NULL;

	if ( py_box != NULL ) {
		type = PyArray_TYPE(py_box);
		switch(type) {
			case NPY_FLOAT:
				box[0] = *( (float*) PyArray_GETPTR1(py_box, 0) );
				box[1] = *( (float*) PyArray_GETPTR1(py_box, 1) );
				box[2] = *( (float*) PyArray_GETPTR1(py_box, 2) );
				break;
			case NPY_DOUBLE:
				box[0] = *( (double*) PyArray_GETPTR1(py_box, 0) );
				box[1] = *( (double*) PyArray_GETPTR1(py_box, 1) );
				box[2] = *( (double*) PyArray_GETPTR1(py_box, 2) );
				break;
			default:
				PyErr_SetString(PyExc_ValueError, "Incorrect type in box vector");
				return NULL;
		}
		pbc = 1;
	}

	i = 0;
	while ( ff_type_map[i] != NULL ) {
		if (!strcmp(ff_type_map[i], ff_type)) ff = i;
		i += 1;
	}
	if (ff == FF_NONE) {
		PyErr_SetString(PyExc_ValueError, "Force field not implemented");
		return NULL; }

	energy = evaluate_energy(py_coords, py_types, ff, py_ff, box);
	py_result = PyFloat_FromDouble(energy);
	return py_result;
}



static PyObject *exposed_read(PyObject *self, PyObject *args, PyObject *kwds) {

	const char *filename;
	char *line;
	char *unit = NULL;
	enum { GUESS, XYZ, MOLDEN, FRAC, GRO } type = GUESS;
	char *str_type = NULL;
	float factor;
	FILE *fd, *test;
	long int fpos;
	struct stat fst;
	char ext[5];

	static char *kwlist[] = {"file", "format", "unit", NULL};

	PyObject *py_result, *py_list;

	py_result = NULL;

	/* Process the arguments */
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|ss",
		kwlist, &filename, &str_type, &unit)) return NULL;

	/* Get the unit measure for coordinates and set the factor accordingly */
	if (unit == NULL) factor = 1.0;
	else {
		make_lowercase(unit);
		if ( !strcmp(unit, "angs") ) factor = 1.0;
		else if ( !strcmp(unit, "bohr") ) factor = BOHR;
		else if ( !strcmp(unit, "nm") ) factor = 10.0;
		else {
			PyErr_SetString(PyExc_ValueError, "Unrecognized measure unit name");
			return NULL; }
	}

	/* Set the enum symbol of the file format */
	if ( str_type != NULL ) {
		if      ( !strcmp(str_type,    "XYZ") ) type = XYZ;
		else if ( !strcmp(str_type, "MOLDEN") ) type = MOLDEN;
		else if ( !strcmp(str_type,   "FRAC") ) type = FRAC;
		else if ( !strcmp(str_type,    "GRO") ) type = GRO;
	}

	/* Open the coordinate file */
	if ( (fd = fopen(filename, "r")) == NULL ) {
		PyErr_SetFromErrno(PyExc_IOError);
		return NULL; }

	/* Guess the file format, if not given explicitly */
	if ( type == GUESS ) {
		strcpy(ext, filename + strlen(filename) - 4);
		if      ( !strcmp(ext, ".xyz") ) type = XYZ;
		else if ( !strcmp(ext, ".gro") ) type = GRO;
		else {
			/* Extract the first line */
			if ( (test = fopen(filename, "r")) == NULL ) {
				PyErr_SetFromErrno(PyExc_IOError);
				return NULL; }
			if ( (line = readline(test)) == NULL ) {
				PyErr_SetFromErrno(PyExc_IOError);
				return NULL; }
			make_lowercase(line);
			stripline(line);
			fclose(test);

			/* Perhaps it's Molden format? */
			if ( !strcmp(line, "[molden format]") ) type = MOLDEN;

			free(line);
		}
	}

	/* Router */
	switch(type) {
		case XYZ:
			py_result = read_xyz(fd, factor);
			/* Support for multiframe XYZ: check the current position in
 			 * the file; if we are near the end, finish; else, continue
 			 * reading. The criteria is that there should be no more than
 			 * 3 bytes left. */
			fpos = ftell(fd);
			if ( fpos == -1 ) {
				PyErr_SetFromErrno(PyExc_IOError);
				return NULL; }
			if ( fstat(fileno(fd), &fst) == -1 ) {
				PyErr_SetFromErrno(PyExc_IOError);
				return NULL; }
			if ( fst.st_size - fpos > 3 ) {
				py_list = PyList_New(1);
				PyList_SetItem(py_list, 0, py_result);
			} else {
				break;
			}
			while ( fst.st_size - fpos > 3 ) {
				py_result = read_xyz(fd, factor);
				PyList_Append(py_list, py_result);
				fpos = ftell(fd);
				if ( fpos == -1 ) {
					PyErr_SetFromErrno(PyExc_IOError);
					return NULL; }
			}
			py_result = py_list;
			break;
		case MOLDEN:
			py_result = read_molden(fd);
			break;
		case FRAC:
			py_result = read_fractional(fd);
			break;
		case GRO:
			py_result = read_gro(fd);
			break;
		/* If the file format is GUESS or different,
		   it means we've failed to guess :-(        */
		case GUESS:
		default:
			PyErr_SetString(PyExc_ValueError, "Unsupported file format");
			return NULL;
	}

	fclose(fd);
	return py_result;
}


static PyObject *exposed_write(PyObject *self, PyObject *args, PyObject *kwds) {

	char *filename;
	char *comment = NULL;
	char *str_format = NULL;
	char *mode = NULL;
	enum { XYZ, GRO } format = XYZ;
	FILE *fd;
	int nat, i;

	static char *kwlist[] = {
		"file", "symbols", "coordinates", "comment", "residues",
		"residue_numbers", "box", "format", "mode", NULL };

	PyObject *py_symbols, *val;
	PyArrayObject *py_coords, *py_box = NULL;
	PyObject *py_resnam = NULL, *py_resid = NULL;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "sO!O!|sO!O!O!ss", kwlist,
			&filename,
			&PyList_Type, &py_symbols,
			&PyArray_Type, &py_coords,
			&comment,
			&PyList_Type, &py_resnam,
			&PyList_Type, &py_resid,
			&PyArray_Type, &py_box,
			&str_format,
			&mode))
		return NULL;

	if ( mode == NULL || !strcmp(mode, "w") ) {
		if( (fd = fopen(filename, "w")) == NULL ) {
			PyErr_SetFromErrno(PyExc_IOError);
			return NULL; }
	} else if ( !strcmp(mode, "a") ) {
		if( (fd = fopen(filename, "a")) == NULL ) {
			PyErr_SetFromErrno(PyExc_IOError);
			return NULL; }
	} else {
		PyErr_SetString(PyExc_ValueError, "Unsupported file mode");
		return NULL; }

	if ( str_format != NULL ) {
		if      ( !strcmp(str_format, "XYZ") ) format = XYZ;
		else if ( !strcmp(str_format, "GRO") ) format = GRO;
	}

	switch(format) {
		case XYZ:
			if ( write_xyz(fd, py_symbols, py_coords, comment) == -1 )
				return NULL;
			break;
		case GRO:
			if ( py_box == NULL ) {
				PyErr_SetString(PyExc_ValueError, "Box vector must be given, when GRO format is used");
				return NULL; }
			nat = PyList_Size(py_symbols);
			if ( py_resnam == NULL ) {
				py_resnam = PyList_New(nat);
				val = Py_BuildValue("s", "");
				for ( i = 0; i < nat; i++ ) {
					PyList_SetItem(py_resnam, i, val);
					Py_INCREF(val);
				}
				Py_DECREF(val);
			}
			if ( py_resid == NULL ) {
				py_resid = PyList_New(nat);
				val = Py_BuildValue("i", 1);
				for ( i = 0; i < nat; i++ ) {
					PyList_SetItem(py_resid, i, val);
					Py_INCREF(val);
				}
				Py_DECREF(val);
			}
			if ( write_gro(fd, py_symbols, py_coords, comment,
			               py_resnam, py_resid, py_box) == -1 )
				return NULL;
			break;
		default:
			PyErr_SetString(PyExc_RuntimeError, "This was unexpected...");
			return NULL;
	}

	fclose(fd);

	Py_RETURN_NONE;
}


static PyObject *find_bonds(PyObject *self, PyObject *args, PyObject *kwds) {
	extern Element element_table[];
	int i, j, nat, type, idx, start;
	float ax, ay, az, bx, by, bz;
	double ar, br, dist;
	npy_intp *numpyint;
	float factor = 1.3;
	char *format = NULL;
	enum Formats { FMT_LIST, FMT_DICT } fmt = FMT_LIST;

	static char *kwlist[] = {
		"coordinates", "types", "factor", "format", NULL };

	PyObject *val1, *val2, *tmp_list, *tmptup;
	PyObject *py_symbols;
	PyArrayObject *py_coords;
	PyObject *py_result = NULL;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|fs", kwlist,
			&PyList_Type, &py_symbols,
			&PyArray_Type, &py_coords,
			&factor,
			&format))
		return NULL;

	if(format != NULL && !strcmp(format, "list"))
		fmt = FMT_LIST;
	else if (format != NULL && !strcmp(format, "dict"))
		fmt = FMT_DICT;
	else if (format != NULL) {
		PyErr_SetString(PyExc_RuntimeError, "Unrecognized format.");
		return NULL;
	}

	nat = PyList_Size(py_symbols);
	numpyint = PyArray_DIMS(py_coords);
	if (numpyint[0] != nat) {
		PyErr_SetString(PyExc_RuntimeError, "Size of symbol and coordinate lists does not match.");
		return NULL;
	}

	type = PyArray_TYPE(py_coords);
	if (type != NPY_FLOAT && type != NPY_DOUBLE) {
		PyErr_SetString(PyExc_ValueError, "Coordinates must be of FLOAT or DOUBLE type.");
		return NULL;
	}

	if (fmt == FMT_DICT) {
		py_result = PyDict_New();
	} else {
		py_result = PyList_New(0);
	}

	for (i = 0; i < nat; i++) {
		if (type == NPY_FLOAT) {
			ax = *( (float*) PyArray_GETPTR2(py_coords, i, 0) );
			ay = *( (float*) PyArray_GETPTR2(py_coords, i, 1) );
			az = *( (float*) PyArray_GETPTR2(py_coords, i, 2) );
		} else {
			ax = *( (double*) PyArray_GETPTR2(py_coords, i, 0) );
			ay = *( (double*) PyArray_GETPTR2(py_coords, i, 1) );
			az = *( (double*) PyArray_GETPTR2(py_coords, i, 2) );
		}
		val1 = PyList_GetItem(py_symbols, i); // borrowed
		//val2 = PyDict_GetItem(py_types, val1); // borrowed
		//ar = PyFloat_AsDouble(val2);
		idx = getElementIndexBySymbol(PyString_AsString(val1));
		if(element_table[idx].number == -1) {
			PyErr_SetString(PyExc_RuntimeError, "Symbol unrecognized.");
			return NULL; }
		ar = element_table[idx].covalent_radius;
		if(ar < 0) {
			PyErr_SetString(PyExc_RuntimeError, "Covalent radius undefined.");
			return NULL; }
		tmp_list = PyList_New(0); // new
		val1 = PyInt_FromLong(i); // new
		if (fmt == FMT_DICT) start = 0;
		else start = i+1;
		for (j = start; j < nat; j++) {
			if (i == j) continue;
			if (type == NPY_FLOAT) {
				bx = *( (float*) PyArray_GETPTR2(py_coords, j, 0) );
				by = *( (float*) PyArray_GETPTR2(py_coords, j, 1) );
				bz = *( (float*) PyArray_GETPTR2(py_coords, j, 2) );
			} else {
				bx = *( (double*) PyArray_GETPTR2(py_coords, j, 0) );
				by = *( (double*) PyArray_GETPTR2(py_coords, j, 1) );
				bz = *( (double*) PyArray_GETPTR2(py_coords, j, 2) );
			}
			val2 = PyList_GetItem(py_symbols, j); // borrowed
			idx = getElementIndexBySymbol(PyString_AsString(val2));
			if(element_table[idx].number == -1) {
				PyErr_SetString(PyExc_RuntimeError, "Symbol unrecognized.");
				return NULL; }
			br = element_table[idx].covalent_radius;
			if(br < 0) {
				PyErr_SetString(PyExc_RuntimeError, "Covalent radius undefined.");
				return NULL; }
			dist = sq(bx-ax) + sq(by-ay) + sq(bz-az);
			//if (dist < sq((ar+br) * factor)) {
			if (sqrt(dist) < (ar+br) * factor) {
				val2 = PyInt_FromLong(j); // new
				if (fmt == FMT_DICT)
					PyList_Append(tmp_list, val2);
				else {
					tmptup = PyTuple_New(2);
					PyTuple_SetItem(tmptup, 0, val1);
					PyTuple_SetItem(tmptup, 1, val2);
					PyList_Append(py_result, tmptup);
					Py_DECREF(tmptup);
				}
				Py_DECREF(val2);
			}
		}
		if (fmt == FMT_DICT)
			if (PyList_Size(tmp_list))
				PyDict_SetItem(py_result, val1, tmp_list);
		Py_DECREF(val1);
		Py_DECREF(tmp_list);
	}

	return py_result;
}


static PyObject *find_molecules(PyObject *self, PyObject *args, PyObject *kwds) {
	int nbonds, i, j, k, molcount = 0, indexFound, natoms;
	unsigned int pyNAtoms;
	long int idx1, idx2, idx3;
	int *checkTable;

	static char *kwlist[] = {
		"natoms", "bonds", NULL };

	PyObject *bond, *mol, *atomIdx;
	PyObject *py_bonds;
	PyObject *py_result = NULL;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "IO!", kwlist,
			&pyNAtoms, &PyList_Type, &py_bonds))
		return NULL;

	checkTable = (int *)malloc(pyNAtoms * sizeof(int));
	for (i = 0; i < pyNAtoms; i++)
		checkTable[i] = 0;

	nbonds = PyList_Size(py_bonds);
	py_result = PyList_New(0);

	for (i = 0; i < nbonds; i++) {
		bond = PyList_GetItem(py_bonds, i); // borrowed
		atomIdx = PyTuple_GetItem(bond, 0);
		idx1 = PyInt_AsLong(atomIdx);
		atomIdx = PyTuple_GetItem(bond, 1);
		idx2 = PyInt_AsLong(atomIdx);
		if (!PyTuple_Check(bond)) {
			PyErr_SetString(PyExc_RuntimeError, "List of bonds should contain tuples.");
			return NULL; }
		if (PyTuple_Size(bond) != 2) {
			PyErr_SetString(PyExc_RuntimeError, "Bond tuple must have exactly two indices.");
			return NULL; }
		indexFound = 0;
		for (j = 0; j < molcount; j++) {
			mol = PyList_GetItem(py_result, j);
			natoms = PyList_Size(mol);
			for (k = 0; k < natoms; k++) {
				atomIdx = PyList_GetItem(mol, k);
				idx3 = PyInt_AsLong(atomIdx);
				if ( idx3 == idx1 ) {
					indexFound = 1;
					atomIdx = PyInt_FromLong(idx2);
					PyList_Append(mol, atomIdx);
					Py_DECREF(atomIdx);
				} else if ( idx3 == idx2 ) {
					indexFound = 1;
					atomIdx = PyInt_FromLong(idx1);
					PyList_Append(mol, atomIdx);
					Py_DECREF(atomIdx);
				}
			}
		}
		//printf("%d %ld %ld\n", indexFound, idx1, idx2);
		if (!indexFound) {
			mol = PyList_New(2);
			PyList_SetItem(mol, 0, PyInt_FromLong(idx1));
			PyList_SetItem(mol, 1, PyInt_FromLong(idx2));
			PyList_Append(py_result, mol);
			Py_DECREF(mol);
			molcount += 1;
		}
		checkTable[idx1] = 1;
		checkTable[idx2] = 1;
	}

	for (i = 0; i < pyNAtoms; i++) {
		if ( !checkTable[i] ) {
			mol = PyList_New(1);
			PyList_SetItem(mol, 0, PyInt_FromLong(i));
			PyList_Append(py_result, mol);
			Py_DECREF(mol);
			molcount += 1;
		}
	}

	free(checkTable);
	return py_result;
}


static PyObject *mass_list(PyObject *self, PyObject *args) {
	int nat, i, j;
	extern Element element_table[];
	char *symbol;
	double *masses;
	npy_intp dims[1];
	PyObject *py_symbols;
	PyObject *py_result = NULL;

	if(!PyArg_ParseTuple(args, "O!", &PyList_Type, &py_symbols))
		return NULL;

	nat = PyList_Size(py_symbols);
	masses = (double*) malloc( nat * sizeof(double) );

	for ( i = 0; i < nat; i++ ) {
		symbol = PyString_AsString(PyList_GetItem(py_symbols, i));
		j = 0;
		while ( strcmp(element_table[j].symbol, symbol) && element_table[j].number != -1 ) j++;
		if(element_table[j].number == -1) {
			PyErr_SetString(PyExc_RuntimeError, "Symbol unrecognized.");
			return NULL;
		}
		masses[i] = element_table[j].mass;
	}

	dims[0] = nat; 
	py_result = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, (double*) masses);

	return py_result;

}


static PyObject *centerofmass(PyObject *self, PyObject *args) {
	int i, nat;
	double totalmass = 0.0;
	float mx = 0.0, my = 0.0, mz = 0.0;
	float *com;
	npy_intp *numpyint;
	npy_intp dims[1];
	double mass;

	PyArrayObject *py_coords, *py_masses;
	PyObject *py_result = NULL;

	if(!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &py_coords, &PyArray_Type, &py_masses))
		return NULL;

	numpyint = PyArray_DIMS(py_masses);
	nat = numpyint[0];

	numpyint = PyArray_DIMS(py_coords);
	if ( nat != numpyint[0] ) {
		PyErr_SetString(PyExc_RuntimeError, "Coordinate- and mass arrays have different size.");
		return NULL;
	}

	for ( i = 0; i < nat; i++ ) {

		mass = *( (double*) PyArray_GETPTR1(py_masses, i) );

		totalmass += mass;
		
		mx += *( (float*) PyArray_GETPTR2(py_coords, i, 0) ) * mass;
		my += *( (float*) PyArray_GETPTR2(py_coords, i, 1) ) * mass;
		mz += *( (float*) PyArray_GETPTR2(py_coords, i, 2) ) * mass;
	}

	com = (float*) malloc( 3*sizeof(float) );
	com[0] = mx/totalmass;
	com[1] = my/totalmass;
	com[2] = mz/totalmass;

	dims[0] = 3;
	py_result = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT, (float*) com);

	return py_result;
}


static PyObject *rot_inertia(PyObject *self, PyObject *args, PyObject *kwds) {
	return NULL;
}


static PyObject *inertia(PyObject *self, PyObject *args, PyObject *kwds) {
	PyArrayObject *py_coords, *py_masses;
	PyObject *py_inertia;
	npy_intp *dim1, *dim2, dims[2];
	int nat, i, j, type;
	double *I;
	double mass;
	float x, y, z;

	static char *kwlist[] = {
		"coordinates", "masses", NULL };

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!", kwlist,
			&PyArray_Type, &py_coords,
			&PyArray_Type, &py_masses))
		return NULL;

	dim1 = PyArray_DIMS(py_coords);
	dim2 = PyArray_DIMS(py_masses);
	if (dim1[0] != dim2[0]) {
		PyErr_SetString(PyExc_RuntimeError, "Arrays not aligned.");
		return NULL;
	}

	nat = dim1[0];

	type = PyArray_TYPE(py_coords);
	if( type != NPY_FLOAT && type != NPY_DOUBLE) {
		PyErr_SetString(PyExc_ValueError, "Incorrect type of the coordinate set");
		return NULL;
	}

	I = (double*) malloc(9 * sizeof(double));
	if( I == NULL) {
		PyErr_SetString(PyExc_RuntimeError, "Could not allocate memory for the tensor.");
		return NULL;
	}

	for (i = 0; i < 9; i++)
			I[i] = 0.0;

	for (i = 0; i < nat; i++) {
		switch(type) {
			case NPY_FLOAT:
				x = *( (float*) PyArray_GETPTR2(py_coords, i, 0) );
				y = *( (float*) PyArray_GETPTR2(py_coords, i, 1) );
				z = *( (float*) PyArray_GETPTR2(py_coords, i, 2) );
				break;
			case NPY_DOUBLE:
				x = *( (double*) PyArray_GETPTR2(py_coords, i, 0) );
				y = *( (double*) PyArray_GETPTR2(py_coords, i, 1) );
				z = *( (double*) PyArray_GETPTR2(py_coords, i, 2) );
				break;
		}
		mass = *( (double*) PyArray_GETPTR1(py_masses, i) );
		I[0] += mass*(sq(y) + sq(z)); // (0,0)
		I[4] += mass*(sq(x) + sq(z)); // (1,1)
		I[8] += mass*(sq(x) + sq(y)); // (2,2)
		I[1] += -mass*x*y; // (0,1)
		I[2] += -mass*x*z; // (0,2)
		I[5] += -mass*y*z; // (1,2)
	}
	I[3] = I[1]; // (1,0)
	I[6] = I[2]; // (2,0)
	I[7] = I[5]; // (2,1)

	/* Add coordinates to the dictionary */
	dims[0] = 3;
	dims[1] = 3;
	py_inertia = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, (double*) I);

	return py_inertia;
}


static PyObject *mep_distance(PyObject *self, PyObject *args, PyObject *kwds) {

	PyArrayObject *py_coords1, *py_coords2;
	PyArrayObject *py_masses = NULL;
	PyObject *py_result;
	npy_intp *dim1, *dim2, *dim3;
	int nat, i, type1, type2;
	double mass;
	float ax, ay, az, bx, by, bz, delta, totalmass;

	static char *kwlist[] = {
		"coordinates1", "coordinates2", "masses", NULL };

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|O!", kwlist,
			&PyArray_Type, &py_coords1,
			&PyArray_Type, &py_coords2,
			&PyArray_Type, &py_masses))
		return NULL;

	dim1 = PyArray_DIMS(py_coords1);
	dim2 = PyArray_DIMS(py_coords2);
	if (dim1[0] != dim2[0]) {
		PyErr_SetString(PyExc_RuntimeError, "Arrays not aligned.");
		return NULL;
	}

	if ( py_masses != NULL ) {
		dim3 = PyArray_DIMS(py_masses);
		if (dim1[0] != dim3[0]) {
			PyErr_SetString(PyExc_RuntimeError, "Arrays not aligned.");
			return NULL;
		}
	}

	nat = dim1[0];

	delta = 0.0;
	totalmass = 0.0;

	type1 = PyArray_TYPE(py_coords1);
	type2 = PyArray_TYPE(py_coords2);
	if( type1 != NPY_FLOAT && type1 != NPY_DOUBLE) {
		PyErr_SetString(PyExc_ValueError, "Incorrect type of the first coordinate set");
		return NULL;
	}
	if( type2 != NPY_FLOAT && type2 != NPY_DOUBLE) {
		PyErr_SetString(PyExc_ValueError, "Incorrect type of the second coordinate set");
		return NULL;
	}

	for (i = 0; i < nat; i++) {
		if ( py_masses != NULL )
			mass = *( (double*) PyArray_GETPTR1(py_masses, i) );
		else
			mass = 1.0;

		switch(type1) {
			case NPY_FLOAT:
				ax = *( (float*) PyArray_GETPTR2(py_coords1, i, 0) );
				ay = *( (float*) PyArray_GETPTR2(py_coords1, i, 1) );
				az = *( (float*) PyArray_GETPTR2(py_coords1, i, 2) );
				break;
			case NPY_DOUBLE:
				ax = *( (double*) PyArray_GETPTR2(py_coords1, i, 0) );
				ay = *( (double*) PyArray_GETPTR2(py_coords1, i, 1) );
				az = *( (double*) PyArray_GETPTR2(py_coords1, i, 2) );
				break;
		}

		switch(type2) {
			case NPY_FLOAT:
				bx = *( (float*) PyArray_GETPTR2(py_coords2, i, 0) );
				by = *( (float*) PyArray_GETPTR2(py_coords2, i, 1) );
				bz = *( (float*) PyArray_GETPTR2(py_coords2, i, 2) );
				break;
			case NPY_DOUBLE:
				bx = *( (double*) PyArray_GETPTR2(py_coords2, i, 0) );
				by = *( (double*) PyArray_GETPTR2(py_coords2, i, 1) );
				bz = *( (double*) PyArray_GETPTR2(py_coords2, i, 2) );
				break;
		}

		//printf("%f %f %lf\n", az, bz, mass);
		delta += mass * (sq(ax - bx) + sq(ay - by) + sq(az - bz));
		totalmass += mass;
	}
	delta = sqrt(delta/totalmass);
	//delta /= sqrt(totalmass);

	return PyFloat_FromDouble((double)delta);
}


static PyMethodDef moltoolsMethods[] = {
    {"read", (PyCFunction)exposed_read, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"dict = read(filename [, unit ] )\n"
		"\n"
		"Reads a coordinate file and returns data as a dictionary.\n"
		"Currently, supports only XYZ and Molden formats. In some file\n"
		"formats the unit can be chosen between unit = angs, bohr, nm.\n"
		"\n"
		"XYZ:    supports extended files (charges in the fifth column)\n"
		"        as well as standard files; coordinates are assumed to be\n"
		"        in Angstroms. Also supports multiframe files. In that\n"
		"        case, returns a list of dictionaries.\n"
		"\n"
		"Molden: supports groups N_GEO, GEOCONV (energies only),\n"
		"        GEOMETRIES (XYZ only), ATOMS (Angstroms only).\n"
		"\n"
		"Output structure, dictionary including various keys, depending\n"
		"on what is found in the file:\n"
		"\n"
		"number_of_atoms - number of atoms (int)\n"
		"comment         - comment line extracted from second line (str)\n"
		"symbols         - atom symbols (list)\n"
		"coordinates     - array of [number_of_atoms,3] coordinates\n"
		"                  (numpy array)\n"
		"charges         - array of [number_of_atoms] point charges\n"
		"                  (numpy array)\n"
		"energies        - energies of the subsequent configurations\n"
		"                  (numpy array)\n"
		"geometries      - list of dictionaries, one for each configuration\n"
		"atomic_numbers  - atomic numbers (nuclear charge, numpy array, int)\n"
		"\n" },
    {"write", (PyCFunction)exposed_write, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"write(FILENAME, SYMBOLS, COORDS, [ comment=COMMENT,\n"
		"      residues=RESNAM, residue_numbers=RESNUM, box=PBCBOX,\n"
		"      format=FORMAT, mode=MODE ])\n"
		"\n"
		"Write structure in XYZ or GRO format. Following parameters are\n"
		"obligatory:\n"
		"FILENAME - file name (string)\n"
		"SYMBOLS  - atomic symbols (list)\n"
		"COORDS   - coordinates (numpy array)\n"
		"Following parameters are optional:\n"
		"FORMAT   - file format, 'XYZ' or 'GRO' (default 'XYZ')\n"
		"MODE     - writing mode, 'w' or 'a' (default 'w')\n"
		"COMMENT  - comment text (string), default is empty\n"
		"RESNAM   - residue names (list of strings), default are empty\n"
		"           strings\n"
		"RESNUM   - residue numbers (list if int), default is 1\n"
		"PBCBOX   - Dimentions of the box, used in PBC (numpy array);\n"
		"           required when GRO format is used\n"
		"\n" },
	{"findBonds", (PyCFunction)find_bonds, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"findBonds(symbols, coordinates, factor=1.3, format=list)\n"
		"\n"
		"Finds connections between atoms. 'symbols' is a list of atomic\n"
		"types, 'coordinates' is a numpy array of coordinates and factor\n"
		"is a factor used to multiply the covalent radii. Depending on the\n"
		"optional parameter 'format', returns either a dictionary of\n"
		"indices { <center> : (<center1>, <center2>, ... ) ... } for every\n"
		"atom or a list of unique bonds as tuples [ (center1, center2),\n"
		"(center1, center3), ... ]\n"
		"\n" },
	{"findMolecules", (PyCFunction)find_molecules, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"findMolecules(numberOfAtoms, bonds)\n"
		"\n"
		"Given a list of bonds as tuples, find sets of topologically\n"
		"connected atoms (molecules).\n"
		"\n" },
    {"energy", (PyCFunction)exposed_energy, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"energy(coordinates, types, ff_type, ff, box=None)\n"
		"\n"
		"Evaluate the energy of the configuration <coordinates>. <types> are assigned\n"
		"from <ff> and the energy is calculated according to <ff_type>, using PBC only\n"
		"if <box> dimensions are specified.\n"
		"\n" },
	{"mass_list", (PyCFunction)mass_list, METH_VARARGS,
		"\n"
		"mass_list(symbols)\n"
		"\n"
		"For a given list of atomic symbols, returns numpy array of masses.\n"
		"\n" },
	{"COM", (PyCFunction)centerofmass, METH_VARARGS,
		"\n"
		"COM(coordinates, masses)\n"
		"\n"
		"For given coordinates (numpy array) evaluate the center of mass\n"
		"(numpy array).\n"
		"\n" },
	{"MEP_distance", (PyCFunction)mep_distance, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"MEP_distance(coordinates1, coordinates2, masses=None)\n"
		"\n"
		"Computes mass-weighted distance of two structures, as used in\n"
		"IRC/MEP calculations. If masses are omitted, no weighting is done.\n"
		"\n" },
	{"inertia", (PyCFunction)inertia, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"inertia_tensor = inertia(coordinates, masses)\n"
		"\n"
		"Computes the moments of inertia tensor of a molecule.\n"
		"\n" },
	{"rot_inertia", (PyCFunction)rot_inertia, METH_VARARGS | METH_KEYWORDS,
		"\n"
		"coordinates_rot, principal_moments = rot_inertia(coordinates, inertia_tensor)\n"
		"\n"
		"Rotates the coordinates so that the principal moments of inertia\n"
		"are aligned with axes of the coordinate system.\n"
		"\n" },

    {NULL, NULL, 0, NULL}        /* Sentinel */
};


int build_tables(PyObject **list_symbols, PyObject **list_names,
                 PyObject **list_masses, PyObject **symbol2number) {

	extern Element element_table[];
	PyObject *key, *val;
	int elements_defined, idx;
	double mass;

	/* Check how many elements we know */
	elements_defined = 0;
	while ( element_table[elements_defined].number >= 0 )
		elements_defined += 1;

	/* Prepare lists of symbols, names and masses */
	*list_symbols  = PyList_New(elements_defined + 1);
	*list_names    = PyList_New(elements_defined + 1);
	*list_masses   = PyList_New(elements_defined + 1);
	*symbol2number = PyDict_New();

	/* Set the first element to None, so that indexing starts at 1 */
	val = Py_BuildValue("s", "");
	PyList_SetItem(*list_symbols, 0, val);
	PyList_SetItem(*list_names, 0, val);
	PyList_SetItem(*list_masses, 0, Py_None);

	/* Fill with data */
	for(idx = 0; idx < elements_defined; idx++) {

		/* Atom symbol */
		val = Py_BuildValue("s", element_table[idx].symbol);
		PyList_SetItem(*list_symbols, element_table[idx].number, val);

		/* Atom name */
		val = Py_BuildValue("s", element_table[idx].name);
		PyList_SetItem(*list_names, element_table[idx].number, val);

		/* Atom mass */
		mass = element_table[idx].mass;
		val = mass >= 0 ? Py_BuildValue("d", mass) : Py_None;
		PyList_SetItem(*list_masses, element_table[idx].number, val);

		/* Hash table symbol->number */
		key = Py_BuildValue("s", element_table[idx].symbol);
		val = Py_BuildValue("i", element_table[idx].number);
		PyDict_SetItem(*symbol2number, key, val);
		Py_DECREF(key);
		Py_DECREF(val);

	}

	return 0;
}


PyMODINIT_FUNC initmoltools(void)
{
    PyObject *md;
	PyObject *exposed_atom_symbols, *exposed_atom_names;
	PyObject *exposed_atom_masses, *exposed_symbol2number;

	/* Use system-wide locale, but make sure that decimal point is a point! */
	setlocale(LC_ALL, "");
	setlocale(LC_NUMERIC, "C");

	if (PyType_Ready(&EAMffType) < 0)
		return;
    md = Py_InitModule("moltools", moltoolsMethods);
    if(md == NULL) return;
	import_array();

	if( build_tables(&exposed_atom_symbols, &exposed_atom_names,
		&exposed_atom_masses, &exposed_symbol2number) ) return;

	/* Add the lists to the module */
	PyModule_AddObject(md, "AtomicSymbols", exposed_atom_symbols);
	PyModule_AddObject(md, "AtomicNames", exposed_atom_names);
	PyModule_AddObject(md, "AtomicMasses" , exposed_atom_masses);
	PyModule_AddObject(md, "AtomicNumbers" , exposed_symbol2number);

	Py_INCREF(&EAMffType);
	PyModule_AddObject(md, "EAMff", (PyObject *)&EAMffType);
}

