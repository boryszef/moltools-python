from numpy.distutils.core import setup, Extension

moltools = Extension('moltools',
	sources = [ 'eam.c', 'periodic_table.c', 'utils.c', 'readers.c',
				'measure.c', 'constants.c', 'writers.c', 'energy.c',
	            'molecule.c', 'moltools.c', ],
	extra_compile_args=["-g"],
	extra_link_args=["-g"]
)

setup (name = 'moltools-python', version = '0.1',
       description = 'Python module for manipulation of atomic coordinates',
       ext_modules = [ moltools ] )
