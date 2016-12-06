#include "moltools.h"


int lookupStringInList(char *needle, char **stack, int len) {
	int i;

	for (i = 0; i < len; i++)
		if (!strcmp(needle, stack[i])) return 1;
	return 0;
}


PyObject *findHBonds(PyObject *self, PyObject *args, PyObject *kwds) {
	PyArrayObject *py_coords, *py_box;
	PyObject *py_syms, *py_accs, *out, *tuple;
	char *atomSymbol1, *atomSymbol2;
	npy_intp *numpyint;
    int use_pbc, type, naccept, natoms;
	int h, i, j, nearest;
	double box[3], half[3];
	double H[3], I[3], J[3], other[3], nearHI[3], nearHJ[3];
	double *A, *B, *C;
	double d2hi, d2hj, d2ij, cosval, cos_cutoff, test;
	double cutoff = 3.5, carbon_cutoff = -1.0;
	double acceptor_cutoff = 2.6, angle_cutoff = 30.;
	char **acceptors;

	static char *kwlist[] = {
		"symbols", "coordinates", "acceptors", "box", "cutoff",
		"carbon_cutoff", "acceptor_cutoff", "angle_cutoff", NULL };

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!O!|O!dddd", kwlist,
			&PyList_Type, &py_syms,
			&PyArray_Type, &py_coords,
			&PyList_Type, &py_accs,
			&PyArray_Type, &py_box,
			&cutoff, &carbon_cutoff, &acceptor_cutoff, &angle_cutoff))
		return NULL;

	if (carbon_cutoff < 0) carbon_cutoff = cutoff;

	cutoff *= cutoff;
	carbon_cutoff *= carbon_cutoff;
	acceptor_cutoff *= acceptor_cutoff;

	cos_cutoff = cos(angle_cutoff / 180. * M_PI);

	if (py_box == NULL) use_pbc = 0;
	else use_pbc = 1;

    if (use_pbc) {
		vectorToDouble(box, py_box);
        half[0] = box[0] / 2.0;
        half[1] = box[1] / 2.0;
        half[2] = box[2] / 2.0;
    }

	type = PyArray_TYPE(py_coords);
	if( type != NPY_FLOAT && type != NPY_DOUBLE) {
		PyErr_SetString(PyExc_ValueError, "Incorrect type of the coordinate set");
		return NULL; }

	naccept = PyList_Size(py_accs);
	acceptors = (char**) malloc(naccept * sizeof(char*));
	for (i = 0; i < naccept; i++)
		acceptors[i] = PyString_AsString(PyList_GetItem(py_accs, i));

	numpyint = PyArray_DIMS(py_coords);
	natoms = numpyint[0];

	out = PyList_New(0);

	for ( h = 0; h < natoms; h++ ) {

		atomSymbol1 = PyString_AsString(PyList_GetItem(py_syms, h));

		if(strcmp(atomSymbol1, "H")) continue;

		switch(type) {
			case NPY_FLOAT:
				H[0] = *( (float*) PyArray_GETPTR2(py_coords, h, 0) );
				H[1] = *( (float*) PyArray_GETPTR2(py_coords, h, 1) );
				H[2] = *( (float*) PyArray_GETPTR2(py_coords, h, 2) );
				break;
			case NPY_DOUBLE:
				H[0] = *( (double*) PyArray_GETPTR2(py_coords, h, 0) );
				H[1] = *( (double*) PyArray_GETPTR2(py_coords, h, 1) );
				H[2] = *( (double*) PyArray_GETPTR2(py_coords, h, 2) );
				break;
		}
		if (use_pbc) wrapCartesian(H, box);

		// Find the nearest heavy atom

		d2hi = -1.0;

		for (i = 0; i < natoms; i++) {

			if (h == i) continue;

			atomSymbol1 = PyString_AsString(PyList_GetItem(py_syms, i));
			//if (!lookupStringInList(atomSymbol1, acceptors, naccept)) continue;

			if(!strcmp(atomSymbol1, "H")) continue;

			switch(type) {
				case NPY_FLOAT:
					I[0] = *( (float*) PyArray_GETPTR2(py_coords, i, 0) );
					I[1] = *( (float*) PyArray_GETPTR2(py_coords, i, 1) );
					I[2] = *( (float*) PyArray_GETPTR2(py_coords, i, 2) );
					break;
				case NPY_DOUBLE:
					I[0] = *( (double*) PyArray_GETPTR2(py_coords, i, 0) );
					I[1] = *( (double*) PyArray_GETPTR2(py_coords, i, 1) );
					I[2] = *( (double*) PyArray_GETPTR2(py_coords, i, 2) );
					break;
			}
			if (use_pbc) {
				wrapCartesian(I, box);
				copyPoint(nearHI, I);
				nearestImage(H, nearHI, half);
				test = distanceSquare(H, nearHI);
			} else
				test = distanceSquare(H, I);

			//if (d2hi > cutoff) continue;
			if (d2hi < 0 || d2hi > test) {
				d2hi = test;
				nearest = i;
			}
		}

		// Skip if the nearest is not in the list
		atomSymbol1 = PyString_AsString(PyList_GetItem(py_syms, nearest));
		if (!lookupStringInList(atomSymbol1, acceptors, naccept)) continue;

		switch(type) {
			case NPY_FLOAT:
				I[0] = *( (float*) PyArray_GETPTR2(py_coords, nearest, 0) );
				I[1] = *( (float*) PyArray_GETPTR2(py_coords, nearest, 1) );
				I[2] = *( (float*) PyArray_GETPTR2(py_coords, nearest, 2) );
				break;
			case NPY_DOUBLE:
				I[0] = *( (double*) PyArray_GETPTR2(py_coords, nearest, 0) );
				I[1] = *( (double*) PyArray_GETPTR2(py_coords, nearest, 1) );
				I[2] = *( (double*) PyArray_GETPTR2(py_coords, nearest, 2) );
				break;
		}

		//if (!strcmp(atomSymbol1, "C") && d2hi > carbon_cutoff) continue;

		for (j = 0; j < natoms; j++) {

			if (j == nearest || j == h) continue;

			atomSymbol2 = PyString_AsString(PyList_GetItem(py_syms, j));
			if (!lookupStringInList(atomSymbol2, acceptors, naccept)) continue;

			if(!strcmp(atomSymbol1, "C")) continue;

			switch(type) {
				case NPY_FLOAT:
					J[0] = *( (float*) PyArray_GETPTR2(py_coords, j, 0) );
					J[1] = *( (float*) PyArray_GETPTR2(py_coords, j, 1) );
					J[2] = *( (float*) PyArray_GETPTR2(py_coords, j, 2) );
					break;
				case NPY_DOUBLE:
					J[0] = *( (double*) PyArray_GETPTR2(py_coords, j, 0) );
					J[1] = *( (double*) PyArray_GETPTR2(py_coords, j, 1) );
					J[2] = *( (double*) PyArray_GETPTR2(py_coords, j, 2) );
					break;
			}

			if (use_pbc) {
				wrapCartesian(J, box);
				copyPoint(nearHJ, J);
				nearestImage(H, nearHJ, half);
				d2hj = distanceSquare(H, nearHJ);
			} else
				d2hj = distanceSquare(H, J);

			if (d2hj > acceptor_cutoff) continue;

			//if (!strcmp(atomSymbol2, "C") && d2hj > carbon_cutoff) continue;

			if (use_pbc) {
				copyPoint(other, J);
				nearestImage(I, other, half);
				d2ij = distanceSquare(I, other);
			} else
				d2ij = distanceSquare(I, J);

			if (d2ij > cutoff) continue;

			if (!strcmp(atomSymbol1, "C") && d2ij > carbon_cutoff) continue;

			A = H;
			if (d2hi < d2hj) {
				B = I;
				C = J;
			} else {
				B = J;
				C = I;
			}

			if (use_pbc) {
				nearestImage(B, A, half);
				nearestImage(B, C, half);
			}

			cosval = threePointAngleCosine(A, B, C);

			if (cosval > cos_cutoff) {
				tuple = Py_BuildValue("(iii)", h, i, j);
				PyList_Append(out, tuple);
				Py_DECREF(tuple);
			}
		}
	}

	free(acceptors);

	return out;
}


PyObject *measureAngleCosine(PyObject *self, PyObject *args, PyObject *kwds) {
	PyArrayObject *py_atom1, *py_atom2, *py_atom3, *py_box = NULL;
	PyObject *out;
	double box[3], half[3];
	double A[3], B[3], C[3], p[3], q[3];
	double lp, lq, cos;
	int use_pbc, type;
	npy_intp *numpyint;

	static char *kwlist[] = {
		"atom1", "atom2", "atom3", "box", NULL };

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!O!|O!", kwlist,
			&PyArray_Type, &py_atom1,
			&PyArray_Type, &py_atom2,
			&PyArray_Type, &py_atom3,
			&PyArray_Type, &py_box))
		return NULL;

	if (py_box == NULL) use_pbc = 0;
	else use_pbc = 1;

    if (use_pbc) {
		vectorToDouble(box, py_box);
        half[0] = box[0] / 2.0;
        half[1] = box[1] / 2.0;
        half[2] = box[2] / 2.0;
    }

	vectorToDouble(A, py_atom1);
	vectorToDouble(B, py_atom2);
	vectorToDouble(C, py_atom3);

	if (use_pbc) {
		wrapCartesian(A, box);
		wrapCartesian(B, box);
		wrapCartesian(C, box);
		nearestImage(B, A, half);
		nearestImage(B, C, half);
	}

	cos = threePointAngleCosine(A, B, C);

	return PyFloat_FromDouble(cos);
}


PyObject *distanceMatrix(PyObject *self, PyObject *args, PyObject *kwds) {
	PyArrayObject *py_coords, *py_box = NULL, *py_squared = NULL;
	PyObject *py_dist;
	int i, j, natoms, type;
	npy_intp *numpyint;
	double ax, ay, az, bx, by, bz, dx, dy, dz, dist;
	double box[3], half[3];
	npy_intp dims[2];
	float *distances;
	int use_pbc, use_square;

	static char *kwlist[] = {
		"coordinates", "box", "squared", NULL };

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!|O!O!", kwlist,
			&PyArray_Type, &py_coords,
			&PyArray_Type, &py_box,
			&PyBool_Type, &py_squared))
		return NULL;

	if (py_box == NULL) use_pbc = 0;
	else use_pbc = 1;

	if(py_squared == NULL || !PyBool_Check(py_squared)) use_square = 0;
	else use_square = 1;

	numpyint = PyArray_DIMS(py_coords);
	natoms = numpyint[0];

	distances = (float*) malloc(natoms * natoms * sizeof(float));
	if(distances == NULL) {
		PyErr_SetFromErrno(PyExc_MemoryError);
		return NULL; }

	if (use_pbc) {
		vectorToDouble(box, py_box);
		half[0] = box[0] / 2.0;
		half[1] = box[1] / 2.0;
		half[2] = box[2] / 2.0;
	}

	type = PyArray_TYPE(py_coords);
	if( type != NPY_FLOAT && type != NPY_DOUBLE) {
		PyErr_SetString(PyExc_ValueError, "Incorrect type of the coordinate set");
		return NULL; }

	for ( i = 0; i < natoms; i++ ) {
		switch(type) {
			case NPY_FLOAT:
				ax = *( (float*) PyArray_GETPTR2(py_coords, i, 0) );
				ay = *( (float*) PyArray_GETPTR2(py_coords, i, 1) );
				az = *( (float*) PyArray_GETPTR2(py_coords, i, 2) );
				break;
			case NPY_DOUBLE:
				ax = *( (double*) PyArray_GETPTR2(py_coords, i, 0) );
				ay = *( (double*) PyArray_GETPTR2(py_coords, i, 1) );
				az = *( (double*) PyArray_GETPTR2(py_coords, i, 2) );
				break;
		}
		if (use_pbc) {
			ax = fmod(ax, box[0]);
			ay = fmod(ay, box[1]);
			az = fmod(az, box[2]);
			if (ax < 0) ax += box[0];
			if (ay < 0) ay += box[1];
			if (az < 0) az += box[2];
		}
		distances[i*natoms + i] = 0.0;
		for ( j = i + 1; j < natoms; j++ ) {
			switch(type) {
				case NPY_FLOAT:
					bx = *( (float*) PyArray_GETPTR2(py_coords, j, 0) );
					by = *( (float*) PyArray_GETPTR2(py_coords, j, 1) );
					bz = *( (float*) PyArray_GETPTR2(py_coords, j, 2) );
					break;
				case NPY_DOUBLE:
					bx = *( (double*) PyArray_GETPTR2(py_coords, j, 0) );
					by = *( (double*) PyArray_GETPTR2(py_coords, j, 1) );
					bz = *( (double*) PyArray_GETPTR2(py_coords, j, 2) );
					break;
			}
			if (use_pbc) {
				bx = fmod(bx, box[0]);
				by = fmod(by, box[1]);
				bz = fmod(bz, box[2]);
				if (bx < 0) bx += box[0];
				if (by < 0) by += box[1];
				if (bz < 0) bz += box[2];
			}
			dx = fabs(ax-bx);
			dy = fabs(ay-by);
			dz = fabs(az-bz);
			if (use_pbc) {
				if (dx > half[0]) dx = box[0] - dx;
				if (dy > half[1]) dy = box[1] - dy;
				if (dz > half[2]) dz = box[2] - dz;
			}
			dist = sq(dx) + sq(dy) + sq(dz);
			if(!use_square) dist = sqrt(dist);
			distances[i*natoms + j] = dist;
			distances[j*natoms + i] = dist;
		}
	}
	dims[0] = natoms;
	dims[1] = natoms;
	py_dist = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, (float*) distances);
	PyArray_ENABLEFLAGS((PyArrayObject*)py_dist, NPY_ARRAY_OWNDATA);

	return py_dist;
}


PyObject *centerofmass(PyObject *self, PyObject *args) {
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


PyObject *inertia(PyObject *self, PyObject *args, PyObject *kwds) {
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
		return NULL; }

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


PyObject *mep_distance(PyObject *self, PyObject *args, PyObject *kwds) {

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

		delta += mass * (sq(ax - bx) + sq(ay - by) + sq(az - bz));
		totalmass += mass;
	}
	delta = sqrt(delta/totalmass);

	return PyFloat_FromDouble((double)delta);
}


