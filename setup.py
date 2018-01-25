import sys
if sys.version_info.major < 3:
    print("This package works with Python3 only.")
    exit(1)

from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import get_info
import subprocess
from distutils.ccompiler import new_compiler
from distutils.sysconfig import get_python_inc, get_python_lib
import glob

def extraPackagePresent(pkg):
    if subprocess.call(['pkg-config', '--exists', pkg]):
        return False
    else:
        return True

def getPackageFlags(pkg):
    cflags = subprocess.check_output(['pkg-config', '--cflags', pkg])
    libs = subprocess.check_output(['pkg-config', '--libs', pkg])
    return cflags.strip(), libs.strip()

extraCFlags = [ "-O3 -march=native" ]
#extraCFlags = [ "-O0 -g -Wall -Wextra" ]
extraLFlags = []
#extraLFlags = [ "-O0 -g" ]
extraInfo = get_info('npymath')

# This is ugly and needs better approach
cc = new_compiler(verbose=1)
cc.add_include_dir(".")
cc.add_include_dir(get_python_inc())
cc.add_library_dir(get_python_lib())
cc.add_library('cunit')
for libdir in extraInfo['library_dirs']:
    extraLFlags.append("-L"+libdir)
    cc.add_library_dir(libdir)
for lib in extraInfo['libraries']:
    extraLFlags.append("-l"+lib)
    cc.add_library(lib)
for hdir in extraInfo['include_dirs']:
    extraLFlags.append("-I"+hdir)
    cc.add_include_dir(hdir)

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
	extra_compile_args=extraCFlags,	extra_link_args=extraLFlags)

x=setup(name = 'moltools-python', version = '0.1.2',
      description = 'Python module for manipulation of atomic coordinates',
      author = 'Borys Szefczyk', author_email = 'boryszef@gmail.com',
      url = 'https://github.com/boryszef/moltools-python',
      ext_modules = [ moltools ] )
cc.add_include_dir(x.include_dirs)
cc.add_library(x.include_dirs[0].split("/")[-1])

def find_object(obj):
    for obj in glob.iglob('build/**/'+obj, recursive=True):
        return obj

# Make tests
obj = cc.compile(["tests/test_utils.c"])
obj.append(find_object("utils.o"))
obj.append(find_object("periodic_table.o"))
cc.link_executable(obj, "tests/test_utils")
