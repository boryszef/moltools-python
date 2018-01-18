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


/* Local helper functions */

static int read_topo_from_xyz(Trajectory *self) {

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
    Py_DECREF(self->symbols);
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




/* Read MOLDEN file and get the sections present. Store offset to the    *
 * particular section too. Returns a python dictionary or NULL on error. */

PyObject* read_molden_sections(FILE *fd) {

    long filepos;
    char *strptr, *line;
    char buffer[256];
    int len;
    PyObject *dict, *key, *val;

    dict = PyDict_New();

    rewind(fd);

    filepos = ftell(fd);
    if ((line = readline(fd)) == NULL) return NULL;
    stripline(line);
    make_lowercase(line);

    do {

        /* Start of a section */
        if(line[0] == '[') {

            /* Get the name */
            strptr = strchr(line, ']');
            len = (int)(strptr - line) - 1;
            strncpy(buffer, line+1, len);
            buffer[len] = '\0';

            /* Put the name and position in list */
            key = PyString_FromString(buffer);
            val = Py_BuildValue("i", filepos);
            PyDict_SetItem(dict, key, val);
            Py_DECREF(key);
            Py_DECREF(val);

        }
        
        free(line);
        filepos = ftell(fd);
        if ((line = readline(fd)) == NULL) return NULL;
        stripline(line);
        make_lowercase(line);

    } while (line[0] != '\0');
    free(line);

    rewind(fd);

    return dict;
}



static int read_topo_from_molden(Trajectory *self) {

    char *line = NULL;
    char *buffpos, *token;
    char **line_store;
    int i, nat;
    int *anum;
    long filepos;

    npy_intp dims[2];
    PyObject *val;
    PyObject *keyGeo, *keyAtom, *keyN, *keyConv;

    keyGeo = PyString_FromString("geometries");
    keyAtom = PyString_FromString("atoms");
    keyN = PyString_FromString("n_geo");
    keyConv = PyString_FromString("geoconv");

    /* Make sure that the sections are done */
    if (self->moldenSections == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Molden sections have not been read");
        return -1;
    }

    /* Determine the style of this file */
    if (PyDict_Contains(self->moldenSections, keyGeo)) {

        filepos = PyInt_AsLong(PyDict_GetItem(self->moldenSections, keyGeo));
        fseek(self->fd, filepos, SEEK_SET);
        if ((line = readline(self->fd)) == NULL) return -1;
        free(line);
        self->filePosition1 = ftell(self->fd);
        self->moldenStyle = MLGEOMETRY;
        read_topo_from_xyz(self);

    } else if (PyDict_Contains(self->moldenSections, keyAtom)) {
        filepos = PyInt_AsLong(PyDict_GetItem(self->moldenSections, keyAtom));
        fseek(self->fd, filepos, SEEK_SET);
        if ((line = readline(self->fd)) == NULL) return -1;
        stripline(line);
        make_lowercase(line);
        if ( !strcmp(line, "[atoms] angs") )
            self->units = ANGS;
        else if ( !strcmp(line, "[atoms] au") )
            self->units = BOHR;
        else {
            PyErr_SetString(PyExc_RuntimeError, "Unrecognized units");
            return -1; }
        free(line);
        self->filePosition1 = ftell(self->fd);
        self->moldenStyle = MLATOMS;
        /* We don't know how many atoms are there, so we have to
           store the lines. */
        nat = 0;
        line_store = (char**) malloc(10 * sizeof(char*));
        if ( line_store == NULL ) {
            PyErr_NoMemory();
            return -1; }
            if ((line = readline(self->fd)) == NULL) return -1;
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
            free(line);

        self->nofatoms = nat;

        /* Get rid of Py_None in self->symbols */
        Py_DECREF(self->symbols);
        self->symbols = PyList_New(nat);

        anum = (int*) malloc(nat * sizeof(int));
        if(anum == NULL) {
            PyErr_SetFromErrno(PyExc_MemoryError);
            return -1; }

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

    } else {

        PyErr_SetString(PyExc_RuntimeError, "geometry/atom section missing");
        return -1; }
    
    if (PyDict_Contains(self->moldenSections, keyN)) {
        filepos = PyInt_AsLong(PyDict_GetItem(self->moldenSections, keyN));
        fseek(self->fd, filepos, SEEK_SET);
        if ((line = readline(self->fd)) == NULL) return -1;
        free(line);
        if ((line = readline(self->fd)) == NULL) return -1;
        stripline(line);
        self->nofframes = atoi(line);
        free(line);
    }

    if (PyDict_Contains(self->moldenSections, keyConv)) {
        filepos = PyInt_AsLong(PyDict_GetItem(self->moldenSections, keyConv));
        fseek(self->fd, filepos, SEEK_SET);
        if ((line = readline(self->fd)) == NULL) return -1;
        free(line);
        if ((line = readline(self->fd)) == NULL) return -1;
        free(line);
        self->filePosition2 = ftell(self->fd);
    }
    Py_DECREF(keyGeo);
    Py_DECREF(keyAtom);
    Py_DECREF(keyN);
    Py_DECREF(keyConv);

    return 0;
}




static int read_topo_from_gro(Trajectory *self) {

    Py_ssize_t pos;
    int nofatoms;
    char *buffer;
    char symbuf[100];
    int *resid;
    npy_intp dims[2];
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

    self->nofatoms = nofatoms;

    /* Get rid of Py_None in self->symbols etc. */
    Py_DECREF(self->symbols);
    self->symbols = PyList_New(nofatoms);
    Py_DECREF(self->resnames);
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

    /* Add residue IDs to the dictionary */
    dims[0] = nofatoms;
    dims[1] = 1;
    Py_DECREF(self->resids);
    self->resids = PyArray_SimpleNewFromData(1, dims, NPY_INT, (int*) resid);
    /***************************************************************
     * Do not free the raw array! It will be still used by Python! *
     ***************************************************************/

    return 0;
}




/* This function is used by other readers, like    *
 * read_frame_from_molden_geometries for instance, *
 * so be careful with implementation.              */

static PyObject *read_frame_from_xyz(Trajectory *self) {

    PyObject *py_result, *py_coord, *py_charges, *val, *key;
    char *buffer, *buffpos, *token;
    int pos, nat;
    float factor;
    ARRAY_REAL *xyz, *charges;
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
        Py_DECREF(py_result);
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
        Py_DECREF(py_result);
        return NULL; }

    /* Read the comment line */
    if((buffer = readline(self->fd)) == NULL) {
        PyErr_SetFromErrno(PyExc_IOError);
        Py_DECREF(py_result);
        return NULL; }
    buffer[strlen(buffer)-1] = '\0';

    val = Py_BuildValue("s", buffer);
    free(buffer);
    key = PyString_FromString("comment");
    PyDict_SetItem(py_result, key, val);
    Py_DECREF(key);
    Py_DECREF(val);

    /* Set-up the raw arrays for coordinates and charges */
    xyz = (ARRAY_REAL*) malloc(3 * self->nofatoms * sizeof(ARRAY_REAL));
    if(xyz == NULL) {
        PyErr_SetFromErrno(PyExc_MemoryError);
        Py_DECREF(py_result);
        return NULL; }
    charges = (ARRAY_REAL*) malloc(self->nofatoms * sizeof(ARRAY_REAL));
    if(charges == NULL) {
        PyErr_SetFromErrno(PyExc_MemoryError);
        Py_DECREF(py_result);
        return NULL; }
    charges_present = 0;

    /* Atom loop */
    for(pos = 0; pos < self->nofatoms; pos++) {

        /* Get the whole line */
        if((buffer = readline(self->fd)) == NULL) {
            PyErr_SetFromErrno(PyExc_IOError);
        Py_DECREF(py_result);
            return NULL; }
        buffer[strlen(buffer)-1] = '\0';
        buffpos = buffer;

        /* Read symbol */
        token = strtok(buffpos, " ");

        /* Read coordinates */
        if ( (token = strtok(NULL, " ")) == NULL) {
            PyErr_SetString(PyExc_IOError, "Missing coordinate");
            Py_DECREF(py_result);
            return NULL; }
        xyz[3*pos + 0] = atof(token) * factor;
        if ( (token = strtok(NULL, " ")) == NULL) {
            PyErr_SetString(PyExc_IOError, "Missing coordinate");
            Py_DECREF(py_result);
            return NULL; }
        xyz[3*pos + 1] = atof(token) * factor;
        if ( (token = strtok(NULL, " ")) == NULL) {
            PyErr_SetString(PyExc_IOError, "Missing coordinate");
            Py_DECREF(py_result);
            return NULL; }
        xyz[3*pos + 2] = atof(token) * factor;

        /* Read charge, if present */
        token = strtok(NULL, " ");
        if ( token != NULL ) {

            /* This is bad: until now, there were no charges */
            if ( pos > 0 && !charges_present ) {
                PyErr_SetString(PyExc_IOError, "Unexpected charges found");
                Py_DECREF(py_result);
                return NULL;
            }

            charges_present = 1;
            charges[pos] = atof(token);

        } else {

            /* This is bad: we were expecting charges here and found nothing */
            if ( pos > 0 && charges_present ) {
                PyErr_SetString(PyExc_IOError, "Missing charges");
                Py_DECREF(py_result);
                return NULL;
            }
        }

        /* Free the line buffer */
        free(buffer);
    }

    /* Add coordinates to the dictionary */
    dims[0] = self->nofatoms;
    dims[1] = 3;
    py_coord = PyArray_SimpleNewFromData(2, dims, NPY_ARRAY_REAL, (ARRAY_REAL*) xyz);
    /***************************************************************
     * Do not free the raw array! It will be still used by Python! *
     ***************************************************************/

    key = PyString_FromString("coordinates");
    PyDict_SetItem(py_result, key, py_coord);
    Py_DECREF(key);
    Py_DECREF(py_coord);


    /* Add charges, if present */
    if ( charges_present ) {
        py_charges = PyArray_SimpleNewFromData(1, dims, NPY_ARRAY_REAL, (ARRAY_REAL*) charges);
        key = PyString_FromString("charges");
        PyDict_SetItem(py_result, key, py_charges);
        Py_DECREF(key);
        Py_DECREF(py_charges);
    } else
        /* Free the charges ONLY if the Python object was not created! */
        free(charges);

    return py_result;

}




static PyObject *read_frame_from_molden_atoms(Trajectory *self) {

    char *line;
    char *token, *buffpos;
    int i;
    ARRAY_REAL *xyz;
    npy_intp dims[2];
    PyObject *py_result, *key, *py_geom;

    /* Prepare dictionary */

    py_result = PyDict_New();

    xyz = (ARRAY_REAL*) malloc(3 * self->nofatoms * sizeof(ARRAY_REAL));
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
    py_geom = PyArray_SimpleNewFromData(2, dims, NPY_ARRAY_REAL, (ARRAY_REAL*) xyz);
    key = PyString_FromString("coordinates");
    PyDict_SetItem(py_result, key, py_geom);
    Py_DECREF(key);
    Py_DECREF(py_geom);

    return py_result;

}




static PyObject *read_frame_from_molden_geometries(Trajectory *self) {

    char *line;
    PyObject *py_result, *key, *val;

    fseek(self->fd, self->filePosition1, SEEK_SET);
    py_result = read_frame_from_xyz(self);
    if (py_result == Py_None) {
        return py_result;
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




static PyObject *read_frame_from_gro(Trajectory *self) {

    int nat, pos;
    char *buffer;
    ARRAY_REAL *xyz, *vel, *box;
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
    xyz = (ARRAY_REAL*) malloc(3 * self->nofatoms * sizeof(ARRAY_REAL));
    if(xyz == NULL) {
        PyErr_SetFromErrno(PyExc_MemoryError);
        return NULL; }
    vel = (ARRAY_REAL*) malloc(3 * self->nofatoms * sizeof(ARRAY_REAL));
    if(vel == NULL) {
        PyErr_SetFromErrno(PyExc_MemoryError);
        return NULL; }
    box = (ARRAY_REAL*) malloc(9 * sizeof(ARRAY_REAL));
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
    py_coord = PyArray_SimpleNewFromData(2, dims, NPY_ARRAY_REAL, (ARRAY_REAL*) xyz);
    /***************************************************************
     * Do not free the raw array! It will be still used by Python! *
     ***************************************************************/

    key = PyString_FromString("coordinates");
    PyDict_SetItem(py_result, key, py_coord);
    Py_DECREF(key);
    Py_DECREF(py_coord);


    /* Add velocities to the dictionary */
    if(velocities_present) {
        py_vel = PyArray_SimpleNewFromData(2, dims, NPY_ARRAY_REAL, (ARRAY_REAL*) vel);
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
    py_box = PyArray_SimpleNewFromData(2, dims, NPY_ARRAY_REAL, (ARRAY_REAL*) box);
    key = PyString_FromString("box");
    PyDict_SetItem(py_result, key, py_box);
    Py_DECREF(key);
    Py_DECREF(py_box);

    return py_result;

}




static PyObject *read_frame_from_xtc(Trajectory *self) {

    PyObject *py_dict, *val, *key, *py_coord, *py_box;
    matrix mbox;
    float *box, *xyz;
    float time, prec;
    npy_intp dims[2];
    gmx_bool bOK;
    int i, step;

    /* Create the dictionary that will be returned */
    py_dict = PyDict_New();

    if (!read_next_xtc(self->xd, self->nofatoms, &step, &time, mbox, self->xtcCoord, &prec, &bOK)) {
        Py_DECREF(py_dict);
        Py_RETURN_NONE;
    }

    if (!bOK) {
        PyErr_SetString(PyExc_IOError, "Corrupted frame");
        return NULL;
    }

    val = Py_BuildValue("i", step);
    key = PyString_FromString("step");
    PyDict_SetItem(py_dict, key, val);
    Py_DECREF(key);
    Py_DECREF(val);

    val = Py_BuildValue("f", time);
    key = PyString_FromString("time");
    PyDict_SetItem(py_dict, key, val);
    Py_DECREF(key);
    Py_DECREF(val);

    box = (float*) malloc(9 * sizeof(float));
    if(box == NULL) {
        PyErr_SetFromErrno(PyExc_MemoryError);
        return NULL; }

    /* Only orthogonal boxes; implement other later */
    box[0] = (float)mbox[0][0];
    box[1] = (float)mbox[0][1];
    box[2] = (float)mbox[0][2];
    box[3] = (float)mbox[1][0];
    box[4] = (float)mbox[1][1];
    box[5] = (float)mbox[1][2];
    box[6] = (float)mbox[2][0];
    box[7] = (float)mbox[2][1];
    box[8] = (float)mbox[2][2];

    dims[0] = 3;
    dims[1] = 3;
    py_box = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, (float*)box);
    key = PyString_FromString("box");
    PyDict_SetItem(py_dict, key, py_box);
    Py_DECREF(key);
    Py_DECREF(py_box);

    /* Set-up the raw arrays for coordinates and charges */
    xyz = (float*) malloc(3 * self->nofatoms * sizeof(float));
    if(xyz == NULL) {
        PyErr_SetFromErrno(PyExc_MemoryError);
        return NULL; }

    for (i = 0; i < self->nofatoms; i++) {
        /* Times 10, because converting from nm */
        xyz[i*3    ] = (float)(self->xtcCoord[i][0] * 10.0);
        xyz[i*3 + 1] = (float)(self->xtcCoord[i][1] * 10.0);
        xyz[i*3 + 2] = (float)(self->xtcCoord[i][2] * 10.0);
    }

    /* Add coordinates to the dictionary */
    dims[0] = self->nofatoms;
    dims[1] = 3;
    py_coord = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, (float*) xyz);
    /***************************************************************
     * Do not free the raw array! It will be still used by Python! *
     ***************************************************************/

    key = PyString_FromString("coordinates");
    PyDict_SetItem(py_dict, key, py_coord);
    Py_DECREF(key);
    Py_DECREF(py_coord);

    return py_dict;
}

/* End of helper functions */


