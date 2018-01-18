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
#include <numpy/halffloat.h>

#include "moltools.h"
#include "trajectory.h"





/* Method definitions */

static void Trajectory_dealloc(Trajectory* self)
{
    PyObject *tmp;

    tmp = self->symbols;
    self->symbols = NULL;
    Py_XDECREF(tmp);

    tmp = self->resids;
    self->resids = NULL;
    Py_XDECREF(tmp);

    tmp = self->resnames;
    self->resnames = NULL;
    Py_XDECREF(tmp);

    tmp = self->atomicnumbers;
    self->atomicnumbers = NULL;
    Py_XDECREF(tmp);

    tmp = self->moldenSections;
    self->moldenSections = NULL;
    Py_XDECREF(tmp);

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
            break;
#endif
        case GUESS:
        default:
            break;
    }
#ifdef HAVE_GROMACS
    sfree(self->xtcCoord);
#endif
    self->ob_type->tp_free((PyObject*)self);
}




static PyObject *Trajectory_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    Trajectory *self;

    self = (Trajectory *)type->tp_alloc(type, 0);
    if (self != NULL) {

        self->type = GUESS;
        self->units = ANGS;
        self->mode = 'r';
        self->filename = NULL;
        self->fd = NULL;
#ifdef HAVE_GROMACS
        self->xd = NULL;
        self->xtcCoord = NULL;
#endif
        self->filePosition1 = -1;
        self->filePosition2 = -1;
        self->moldenStyle = MLUNKNOWN;
        self->nofatoms = 0;
        self->nofframes = 0;
        self->lastFrame = -1;

        Py_INCREF(Py_None);
        self->symbols = Py_None;
        
        Py_INCREF(Py_None);
        self->atomicnumbers = Py_None;
        
        Py_INCREF(Py_None);
        self->resids = Py_None;
        
        Py_INCREF(Py_None);
        self->resnames = Py_None;
        
        Py_INCREF(Py_None);
        self->moldenSections = Py_None;
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

    PyObject *py_sym = NULL;
    PyObject *py_resid = NULL;
    PyObject *py_resn = NULL;;

    static char *kwlist[] = {
        "filename",
        "format", "mode", "units",
        "symbols", "resids", "resnames", NULL };

    if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|sssO!O!O!", kwlist,
            &filename,
            &str_type, &mode, &units,
            &PyList_Type, &py_sym,
            &PyArray_Type, &py_resid,
            &PyList_Type, &py_resn))
        return -1;

    if (py_sym != NULL) {
        Py_DECREF(self->symbols);
        self->symbols = py_sym;
        Py_INCREF(self->symbols);
    }

    if (py_resid != NULL) {
        Py_DECREF(self->resids);
        self->resids = py_resid;
        Py_INCREF(self->resids);
    }

    if (py_resn != NULL) {
        Py_DECREF(self->resnames);
        self->resnames = py_resn;
        Py_INCREF(self->resnames);
    }

    self->filename = (char*) malloc((strlen(filename)+1) * sizeof(char));
    strcpy(self->filename, filename);
    if (mode == NULL || mode[0] == 'r')
        self->mode = 'r';
    else if (mode[0] == 'w')
        self->mode = 'w';
    else if (mode[0] == 'a')
        self->mode = 'a';
    else {
        PyErr_SetString(PyExc_ValueError, "Incorrect mode");
        return -1;
    }

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
        else if (self->mode == 'r' || self->mode == 'a') {
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
        } else {
            PyErr_SetString(PyExc_ValueError, "Could not guess file format");
            return -1;
        }
    }

    /* Set correct units */
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
            case GUESS:
            default:
                PyErr_SetString(PyExc_ValueError, "Unsupported file format");
                return -1;
                break;
        }
    } else {
        if      (!strcmp(units, "angs")) self->units = ANGS;
        else if (!strcmp(units, "bohr")) self->units = BOHR;
        else if (!strcmp(units,   "nm")) self->units = NM;
    }

    if (self->mode == 'w' || self->mode == 'a') {

        if (self->symbols == Py_None) {
            PyErr_SetString(PyExc_ValueError, "Need atomic symbols");
            return -1; }

        self->nofatoms = PyList_Size(self->symbols);

        /* Open the coordinate file */
        switch(self->type) {
            case XYZ:
            case GRO:
                if ( (self->fd = fopen(filename, mode)) == NULL ) {
                    PyErr_SetFromErrno(PyExc_IOError);
                    return -1; }
                break;
            case MOLDEN:
            case XTC:
            default:
                PyErr_SetString(PyExc_NotImplementedError,
                                "Writing in this format is not implemented");
                return -1;
                break;
        }

    } else {

        /* Open the coordinate file */
        switch(self->type) {
            case XYZ:
            case GRO:
                if ( (self->fd = fopen(filename, "r")) == NULL ) {
                    PyErr_SetFromErrno(PyExc_IOError);
                    return -1; }
                break;
            case MOLDEN:
                if ( (self->fd = fopen(filename, "r")) == NULL ) {
                    PyErr_SetFromErrno(PyExc_IOError);
                    return -1; }
                Py_DECREF(self->moldenSections);
                self->moldenSections = read_molden_sections(self->fd);
                break;
            case XTC:
#ifdef HAVE_GROMACS
                if( (self->xd = open_xtc(filename, "r")) == NULL) {
                    PyErr_SetString(PyExc_IOError, "Error opening XTC file");
                    return -1; }
#else
                PyErr_SetString(PyExc_SystemError,
                    "The module has to be compiled with gromacs support to handle XTC files");
                return -1;
#endif
                break;
            case GUESS:
            default:
                PyErr_SetString(PyExc_ValueError, "Unsupported file format");
                return -1;
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
                close_xtc(self->xd);
                self->xd = open_xtc(filename, "r");
#endif
                break;
            /* If the file format is GUESS or different,
               it means we've failed to guess :-(        */
            case GUESS:
            default:
                PyErr_SetString(PyExc_ValueError, "Unsupported file format");
                return -1;
        }
    }

    return 0;
}



static PyObject *Trajectory_read(Trajectory *self) {

    PyObject *py_result = NULL;

    if (self->mode != 'r') {
        PyErr_SetString(PyExc_RuntimeError, "Trying to read in write mode");
        return NULL; }

    switch(self->type) {

        case XYZ:
            py_result = read_frame_from_xyz(self);
            if (py_result == Py_None) return py_result;
            self->filePosition1 = ftell(self->fd);
            self->filePosition1 = self->filePosition2;
            break;

        case GRO:
            py_result = read_frame_from_gro(self);
            if (py_result == Py_None) return py_result;
            self->filePosition1 = ftell(self->fd);
            self->filePosition1 = self->filePosition2;
            break;

        case MOLDEN:
            if (self->moldenStyle == MLATOMS)
                py_result = read_frame_from_molden_atoms(self);
            else
                py_result = read_frame_from_molden_geometries(self);
            if (py_result == Py_None) return py_result;
            break;

        case XTC:
            py_result = read_frame_from_xtc(self);
            break;

        default:
            break;
    }

    self->lastFrame += 1;
    return py_result;

}





static PyObject *Trajectory_write(Trajectory *self, PyObject *args, PyObject *kwds) {

    PyObject *py_coords = NULL;
    PyObject *py_vel = NULL;
    PyObject *py_box = NULL;
    char *comment = NULL;;
    int i, type;
    npy_half hx, hy, hz;
    long double lx, ly, lz;
    double dx, dy, dz;
    float x, y, z;
    char *s;

    static char *kwlist[] = {
        "coordinates", "velocities", "box", "comment", NULL };

    if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!|O!O!s", kwlist,
            &PyArray_Type, &py_coords,
            &PyArray_Type, &py_vel,
            &PyArray_Type, &py_box,
            &comment))
        return NULL;

    if (self->mode != 'a' && self->mode != 'w') {
        PyErr_SetString(PyExc_RuntimeError, "Trying to write in read mode");
        return NULL; }

    switch(self->type) {

        case XYZ:
            fprintf(self->fd, "%d\n", self->nofatoms);
            if( comment != NULL )
                fprintf(self->fd, "%s\n", comment);
            else
                fprintf(self->fd, "\n");

            type = PyArray_TYPE((PyArrayObject*)py_coords);
            for (i = 0; i < self->nofatoms; i++) {
                s = PyString_AsString(PyList_GetItem(self->symbols, i));
                switch(type) {
                    case NPY_HALF:
                        hx = *( (npy_half*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 0) );
                        hy = *( (npy_half*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 1) );
                        hz = *( (npy_half*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 2) );
                        x = npy_half_to_float(hx);
                        y = npy_half_to_float(hy);
                        z = npy_half_to_float(hz);
                        fprintf(self->fd, "%-3s  %8.3f  %8.3f  %8.3f\n", s, x, y, z);
                        break;
                    case NPY_FLOAT:
                        x = *( (float*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 0) );
                        y = *( (float*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 1) );
                        z = *( (float*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 2) );
                        fprintf(self->fd, "%-3s  %10.6f  %10.6f  %10.6f\n", s, x, y, z);
                        break;
                    case NPY_DOUBLE:
                        dx = *( (double*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 0) );
                        dy = *( (double*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 1) );
                        dz = *( (double*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 2) );
                        fprintf(self->fd, "%-3s  %12.8lf  %12.8lf  %12.8lf\n", s, dx, dy, dz);
                        break;
                    case NPY_LONGDOUBLE:
                        lx = *( (long double*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 0) );
                        ly = *( (long double*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 1) );
                        lz = *( (long double*) PyArray_GETPTR2((PyArrayObject*)py_coords, i, 2) );
                        fprintf(self->fd, "%-3s  %16.10Lf  %16.10Lf  %16.10Lf\n", s, lx, ly, lz);
                        break;
                    default:
                        PyErr_Format(PyExc_ValueError, "Incorrect type of coordinates array (%u)", type);
                        return NULL;
                }
            }
            break;

        case GRO:
            if (traj_write_gro(self, py_coords, py_vel, py_box, comment)) {
                PyErr_SetString(PyExc_RuntimeError, "Could not write");
                return NULL;
            }
            break;

        default:
            break;
    }

    self->lastFrame += 1;
    Py_RETURN_NONE;

}


static PyObject* Trajectory_repr(Trajectory *self) {
    PyObject* str;
    char format[10];

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
    str = PyString_FromFormat("Trajectory('%s', format='%s', mode='%c')",
                                self->filename, format, self->mode);
    return str;
}

/* End of method definitions */




/* Class definition */

static PyMemberDef Trajectory_members[] = {
    {"symbols", T_OBJECT_EX, offsetof(Trajectory, symbols), READONLY,
     "A list of atomic symbols"},
    {"atomicNumbers", T_OBJECT_EX, offsetof(Trajectory, atomicnumbers), READONLY,
     "An ndarray with atomic numbers"},
    {"resIDs", T_OBJECT_EX, offsetof(Trajectory, resids), READONLY,
     "An ndarray with residue numbers - one number per atom"},
    {"resNames", T_OBJECT_EX, offsetof(Trajectory, resnames), READONLY,
     "A list of residue names"},
    {"nOfAtoms", T_INT, offsetof(Trajectory, nofatoms), READONLY,
     "Number of atoms (int)"},
    {"nOfFrames", T_INT, offsetof(Trajectory, nofframes), READONLY,
     "Number of frames (int)"},
    {"lastFrame", T_INT, offsetof(Trajectory, lastFrame), READONLY,
     "The number of the Last frame read or written"},
    {"moldenSections", T_OBJECT_EX, offsetof(Trajectory, moldenSections), READONLY,
     "Dictionary containing byte offsets to sections in Molden file"},
    {"filename", T_STRING, offsetof(Trajectory, filename), READONLY,
     "File name (str)"},
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
        "box (ndarray) shape=3,3\n"
        "\n" },

    {NULL}  /* Sentinel */
};




PyTypeObject TrajectoryType = {

    PyObject_HEAD_INIT(NULL)
    0,                           /*ob_size*/
    "moltools.Trajectory",       /*tp_name*/
    sizeof(Trajectory),          /*tp_basicsize*/
    0,                           /*tp_itemsize*/

    /* Methods to implement standard operations */
    (destructor)Trajectory_dealloc, /*tp_dealloc*/
    0,                              /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    (reprfunc)Trajectory_repr, /*tp_repr*/

    /* Method suites for standard classes */
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/

    /* More standard operations (here for binary compatibility) */
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/

    /* Functions to access object as input/output buffer */
    0,                         /*tp_as_buffer*/

    /* Flags to define presence of optional/expanded features */
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/

    /* Documentation string */
    "Trajectory class. Implements reading of trajectories from XYZ. Molden, GRO\n"
    "and XTC. Writing is implemented for XYZ and GRO. The process is two-step;\n"
    "first, the object must be created, by specifying filename (for reading) or\n"
    "topology information (for writing). Second, frames can be read/saved\n"
    "repeteadly. Reading examples:\n"
    "  traj = Trajectory('my.xyz')\n"
    "  frame1 = traj.read()\n"
    "  frame2 = traj.read()\n"
    "Object of the class Trajectory contains such fields as:\n"
    "symbols, atomicnumbers, resids, resnames, nofatoms, nofframes, lastframe,\n"
    "moldenSections, filename. Method read() returns a dictionary with items\n"
    "depending on the file format, but at least 'coordinates' are present.\n"
    "Writing example:\n"
    "  traj = Trajectory('my.xyz', symbols_list)\n"
    "  traj.write(coordinates1)\n"
    "  traj.write(coordinates2)\n"
    "When writing a trajectory, at least the file name and the list of symbols\n"
    "must be specified.\n"
    "Creating an instance for reading:\n"
    "  traj = Trajectory(filename, format='GUESS', mode='r', units='angs')\n"
    "Available formats include: XYZ, GRO, MOLDEN, XTC - guessed if not specified.\n"
    "Mode: 'r' (default), 'w', 'a'.\n"
    "Units: 'angs' (default), 'bohr', 'nm'.\n"
    "Creating an instance for writing:\n"
    "  traj = Trajectory(filename, format='GUESS', mode='w', symbols=, resids=,\n"
    "                    resnames=)\n"
    "symbols and resnames are lists, while resids, coodinates and velocities are\n"
    "ndarrays.\n",           /* tp_doc */

    /* Assigned meaning in release 2.0 */
    /* call function for all accessible objects */
    0,                       /* tp_traverse */

    /* delete references to contained objects */
    0,                       /* tp_clear */

    /* Assigned meaning in release 2.1 */
    /* rich comparisons */
    0,                       /* tp_richcompare */

    /* weak reference enabler */
    0,                       /* tp_weaklistoffset */

    /* Added in release 2.2 */
    /* Iterators */
    0,                       /* tp_iter */
    0,                       /* tp_iternext */

    /* Attribute descriptor and subclassing stuff */
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
    0, /* tp_free; Low-level free-memory routine */
    0, /* tp_is_gc; For PyObject_IS_GC */
    0, /* tp_bases; */
    0, /* tp_mro; method resolution order */
    0, /* tp_cache; */
    0, /* tp_subclasses; */
    0, /* tp_weaklist; */
};

/* End of class definition */
