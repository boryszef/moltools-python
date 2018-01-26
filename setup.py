import sys
if sys.version_info.major < 3:
    print("This package works with Python3 only.")
    exit(1)

from setuptools import setup, Extension
from numpy.distutils.misc_util import get_info
import glob
import unittest
#import subprocess

#def extraPackagePresent(pkg):
#    if subprocess.call(['pkg-config', '--exists', pkg]):
#        return False
#    else:
#        return True

#def getPackageFlags(pkg):
#    cflags = subprocess.check_output(['pkg-config', '--cflags', pkg])
#    libs = subprocess.check_output(['pkg-config', '--libs', pkg])
#    return cflags.strip(), libs.strip()

#extraCFlags = [ "-O3 -march=native" ]
extraCFlags = [ "-O0 -g" ]
extraLFlags = []

# Check if gromacs is present
#if extraPackagePresent('libgromacs'):
#    cflags, libs = getPackageFlags('libgromacs')
#    extraCFlags.append('-DHAVE_GROMACS')
#    extraCFlags.append(cflags)
#    extraLFlags.append(libs)

npymathInfo = get_info('npymath')
mdarray = Extension(
    'mdarray',
    sources = glob.glob('mdarray/*.c'),
    include_dirs = npymathInfo['include_dirs'],
    library_dirs = npymathInfo['library_dirs'],
    libraries = npymathInfo['libraries'],
    extra_compile_args = extraCFlags,
    extra_link_args=extraLFlags)

def test_suite():
    test_loader = unittest.TestLoader()
    suite = test_loader.discover('tests', pattern='test_*.py')
    return suite

setup(
    name = 'mdarray',
    version = '0.1.2',
    description = 'Python module for manipulation of atomic coordinates and trajectories',
    keywords = 'molecular dynamics modeling modelling trajectory xyz',
    author = 'Borys Szefczyk',
    author_email = 'boryszef@gmail.com',
    url = 'https://github.com/boryszef/mdarray',
    license = 'GPL-3.0',
    install_requires = [ 'numpy>=1.10' ],
    ext_modules = [ mdarray ],
    test_suite = "setup.test_suite",
    zip_safe = False
    )
