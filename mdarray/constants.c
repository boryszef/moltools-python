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

#include "periodic_table.h"


int build_tables(PyObject **list_symbols, PyObject **list_names,
                 PyObject **list_masses, PyObject **symbol2number,
                 PyObject **dict_radii) {

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
	*dict_radii    = PyDict_New();

	/* Set the first element to None, so that indexing starts at 1 */
	val = Py_BuildValue("s", "");
	PyList_SetItem(*list_symbols, 0, val);
	Py_INCREF(val);
	PyList_SetItem(*list_names, 0, val);
	Py_INCREF(Py_None);
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
		if (mass >= 0) {
			val = Py_BuildValue("d", mass);
		} else {
			val = Py_None;
			Py_INCREF(Py_None);
		}
		PyList_SetItem(*list_masses, element_table[idx].number, val);

		/* Hash tables symbol->number and symbol->covalent radius */
		key = Py_BuildValue("s", element_table[idx].symbol);
		val = Py_BuildValue("i", element_table[idx].number);
		PyDict_SetItem(*symbol2number, key, val);
		Py_DECREF(val);
		val = Py_BuildValue("f", element_table[idx].covalent_radius);
		PyDict_SetItem(*dict_radii, key, val);
		Py_DECREF(val);
		Py_DECREF(key);

	}

	return 0;
}


