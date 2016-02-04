from numpy.distutils.core import setup, Extension

moltools = Extension('moltools',
	sources = [ 'periodic_table.c', 'utils.c', 'readers.c',
	            'measure.c', 'constants.c', 'writers.c', 'topology.c',
	            'molecule.c', 'moltools.c', ],
)

setup (name = 'moltools-python', version = '0.1',
       description = 'Python module for manipulation of atomic coordinates',
       ext_modules = [ moltools ] )
