from numpy.distutils.core import setup, Extension

# Check if gromacs is present
if extraPackagePresent(libgromacs):
    cflags =
    libs =

moltools = Extension('moltools',
	sources = [ 'eam.c', 'periodic_table.c', 'utils.c', 'readers.c',
	            'measure.c', 'constants.c', 'writers.c', 'topology.c',
	            'molecule.c', 'moltools.c', ],
	#extra_compile_args=["-O0 -g"],
	#extra_link_args=["-O0 -g"]
)

setup (name = 'moltools-python', version = '0.1.1',
       description = 'Python module for manipulation of atomic coordinates',
       author = 'Borys Szefczyk', author_email = 'boryszef@gmail.com',
       url = 'https://github.com/boryszef/moltools-python',
       ext_modules = [ moltools ] )
