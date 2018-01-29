
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



# Introduction

**mdarray** is a Python module for (mostly) reading trajectories from
Molecular Dynamics. The central part of the module is Trajectory class that
reads coordinates from the trajectory (plus some other data), which can be
further processed in Python scripts. The module has been developed with two
assumptions in mind:
* it should be fast and easy to use
* the data should be stored using standard types, so that the user can
  easily build own programs on top of it.
Because of that, each frame read is returned as dictionary with items such as
coordinates, velocities, box, etc. Arrays of numbers are stored using ndarray
type provided by the Numpy package. The choice of Numpy is dictated by the
fact, that the data stored in array can be built directly in C code as a
contiguous C array and is not further transformed by Numpy, but only
interfaced in Python. As a result, reading and processing the data is fast.
Also, Numpy is a versatile package with dozens of functions and operators that
allow fast manipulation of arrays.

The `Trajectory` class currently supports reading XYZ, GRO and XTC formats, the
latter one only if gromacs library is available in the system. Writting is
supported in XYZ and GRO formats.

# Usage

Simple reading examples:

```Python
>>> import mdarray
>>> traj = mdarray.Trajectory('meoh.xyz')
```

This will determine the type of the file and read all topology information
- whatever is available in the provided format. In case of xyz file that means
atomic symbols only:

```Python
>>> print(traj.nAtoms)
6
>>> print(traj.symbols)
['C', 'H', 'H', 'H', 'O', 'H']
```

The instance contains also standard atomic masses:
```Python
>>> print(traj.masses)
[ 12.011   1.008   1.008   1.008  15.999   1.008]
```

To read single frame from file, call `read()` method:
```Python
>>> frame = traj.read()
>>> print(frame)
{'coordinates': array([[ 0.   ,  0.   ,  0.   ],
       [ 0.   ,  0.   ,  1.089],
       [ 1.027,  0.   , -0.363],
       [-0.513, -0.889, -0.363],
       [-0.66 ,  1.143, -0.467],
       [-0.66 ,  1.143, -1.414]]), 'comment': 'Methanol molecule'}
```

In case of formats that support multiple frames, such as XYZ and XTC, reading
is done sequentially, since this is more economic and faster. Each call to
`read()` reads just a single frame. The method will return None if there are
no more frames to read:
```Python
>>> frame = traj.read()
>>> print(frame)
None
```

Reading GRO file is similar:
```Python
>>> import mdarray
>>> traj = mdarray.Trajectory('gas.gro')
```

GRO files contain additional topology information:
```Python
>>> traj.resids
array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=int32)
>>> traj.resNames
['CHL', 'CHL', 'CHL', 'CHL', 'CHL', 'CHL', 'CHL', 'CHL', 'CHL', 'CHL']
>>> traj.symbols
['CA', 'HA1', 'HA2', 'OA', 'HA3', 'CB', 'HB1', 'HB2', 'OB', 'HB3']
```

Reading frame from GRO is exactly like from XYZ:
```Python
>>> frame = traj.read()
>>> print(frame)
{'comment': 'Gas phase', 'box': array([[ 35.0909996 ,   0.        ,   0.        ],
       [  0.        ,  34.57499981,   0.        ],
       [  0.        ,   0.        ,  37.08400011]]), 'coordinates': array([[  9.99000013,  13.83000016,  15.03999949],
       [  9.63999987,  13.05999994,  15.49999952],
       [  9.59999979,  14.70000029,  15.39999962],
       [ 11.39999986,  13.83000016,  15.21999955],
       [ 11.74000025,  14.63000059,  15.04999995],
       [  9.59999979,  13.68000031,  13.5800004 ],
       [  8.70999992,  13.60999942,  13.53999972],
       [ 10.        ,  12.92999983,  13.20000052],
       [ 10.06000042,  14.8300004 ,  12.84999967],
       [ 10.34999967,  14.50999975,  12.13999987]])}
```
Note that GRO format has information about PBC, hence the dictionary has the
`box` key with (3,3) array.

**mdarray** always converts coordinates to Angstroms. It is assumed that XYZ
files are in Angstroms, while GRO and XTC formats are in nm. However, it is
possible to set input units (angs, nm, bohr) like this:
```Python
>>> traj = mdarray.Trajectory('meoh.xyz', units="bohr")
```

This will convert atomic units into Angstroms while reading.

Writing will be illustrated with the following example: let's take coordinates
from GRO file, shift all atoms by a vector (10, -10, 0) and save to XYZ.

Import mdarray and numpy as well - for vector operations:
```Python
>>> import mdarray
>>> import numpy
```

Load input trajectory and read a frame:
```Python
>>> traj = mdarray.Trajectory('gas.gro')
>>> frame = traj.read()
```

Translate coordinates by a vector:
```Python
>>> crd = frame['coordinates']
>>> vec = numpy.array([10, -10, 0])
>>> crd += vec
```

Open new file for writing ('w' mode); symbols must be specified if XYZ format
is used - we'll take them from GRO trajectory:
```Python
>>> out = mdarray.Trajectory('out.xyz', 'w', symbols=traj.symbols)
```

Write out the modified coordinates with original comment line from GRO file:
```Python
>>> out.write(crd, comment=frame['comment'])
```

Since the `write()` method can be called several times, in order to save
subsequent frames, the file is not being closed at once. The class keeps the
file open until the instance is destroyed. Because of that, when you write
to a file, better delete the instance as soon as you don't need it anymore:
```Python
>>> del out
```

# Installation

To install the package, simply type
```
python setup.py install
```
To install as a regular user, add `--user` flag.

Note: the module supports only Python 3 and requires Numpy package. To use XTC
files, the module uses libgromacs shared library, which is automatically
detected using `pkg-config`.
