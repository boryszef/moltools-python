import sys
if sys.version_info.major < 3:
    print("This package works with Python3 only.")
    exit(1)

from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import get_info

extraCFlags = [ ]#"-O0 -g -Wall -Wextra" ]
extraLFlags = [ ]#"-O0 -g" ]
extraInfo = get_info('npymath')

moltools = Extension('moltools',

	sources = [ 'trajectory.c', 'trajectory.h',
                    'utils.c', 'utils.h',
                    'readers.c', 'readers.h',
                    'writers.c', 'writers.h',
                    'moltools.c', 'moltools.h' ],
	extra_compile_args=extraCFlags,	extra_link_args=extraLFlags,

)

setup (name = 'moltools-python', version = '0.1.2',
       description = 'Python module for manipulation of atomic coordinates',
       author = 'Borys Szefczyk', author_email = 'boryszef@gmail.com',
       url = 'https://github.com/boryszef/moltools-python',
       ext_modules = [ moltools ] )
