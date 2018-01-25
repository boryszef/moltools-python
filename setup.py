import sys
if sys.version_info.major < 3:
    print("This package works with Python3 only.")
    exit(1)

from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import get_info
import subprocess

def extraPackagePresent(pkg):
    if subprocess.call(['pkg-config', '--exists', pkg]):
        return False
    else:
        return True

def getPackageFlags(pkg):
    cflags = subprocess.check_output(['pkg-config', '--cflags', pkg])
    libs = subprocess.check_output(['pkg-config', '--libs', pkg])
    return cflags.strip(), libs.strip()

extraCFlags = [ "-O0 -g -Wall -Wextra" ]
extraLFlags = [ "-O0 -g" ]
extraInfo = get_info('npymath')

for libdir in extraInfo['library_dirs']:
    extraLFlags.append("-L"+libdir)
for lib in extraInfo['libraries']:
    extraLFlags.append("-l"+lib)
for hdir in extraInfo['include_dirs']:
    extraLFlags.append("-I"+hdir)

# Check if gromacs is present
if extraPackagePresent('libgromacs'):
    cflags, libs = getPackageFlags('libgromacs')
    extraCFlags.append('-DHAVE_GROMACS')
    extraCFlags.append(cflags)
    extraLFlags.append(libs)

moltools = Extension('moltools',
	sources = [ 'trajectory.c', 'trajectory.h',
                    'utils.c', 'utils.h',
                    'periodic_table.c', 'periodic_table.h',
                    'moltools.c', 'moltools.h' ],
	extra_compile_args=extraCFlags,	extra_link_args=extraLFlags,

)

setup (name = 'moltools-python', version = '0.1.2',
       description = 'Python module for manipulation of atomic coordinates',
       author = 'Borys Szefczyk', author_email = 'boryszef@gmail.com',
       url = 'https://github.com/boryszef/moltools-python',
       ext_modules = [ moltools ] )
