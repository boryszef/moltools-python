import sys
if sys.version_info.major < 3:
    print("This package works with Python3 only.")
    exit(1)

from setuptools import setup, Extension
from numpy.distutils.misc_util import get_info
import glob
import unittest
import subprocess

def extraPackagePresent(pkg):
    if subprocess.call(['pkg-config', '--exists', pkg]): return False
    else: return True

def getPackageFlags(pkg):
    flagMap = { '-I' : 'inc_dirs', '-L' : 'lib_dirs', '-l' : 'libs', '-D' : 'define' }
    out = { k : [] for k in flagMap.values() }
    out['extra'] = []
    flags = subprocess.check_output(['pkg-config', '--cflags', '--libs', pkg])
    for token in flags.decode('ASCII').split():
        if token[:2] in flagMap: key = flagMap[token[:2]]
        else: key = 'extra'
        out[key].append(token[2:])
    # Convert '=' in define's into tuples
    tmp = []
    for token in out['define']:
        if '=' in token: tmp.append(token.split('='))
        else: tmp.append((token, None))
    out['define'] = tmp
    return out

inc_dirs = []
lib_dirs = []
libs = []
define = []

extraCFlags = [ "-O3", "-march=native" ]
#extraCFlags = [ "-O0", "-g" ]
extraLFlags = []
#extraLFlags = [ "-O0", "-g" ]

npymathInfo = get_info('npymath')
inc_dirs.extend(npymathInfo['include_dirs'])
lib_dirs.extend(npymathInfo['library_dirs'])
libs.extend(npymathInfo['libraries'])

# Check if gromacs is present
if extraPackagePresent('libgromacs'):
    flags = getPackageFlags('libgromacs')
    inc_dirs.extend(flags['inc_dirs'])
    lib_dirs.extend(flags['lib_dirs'])
    libs.extend(flags['libs'])
    define.extend(flags['define'])
    define.append(('HAVE_GROMACS', None))
    extraCFlags.extend(flags['extra'])
    extraLFlags.extend(flags['extra'])

print("INC_DIRS:", inc_dirs)
print("LIB_DIRS:", lib_dirs)
print("LIBS:", libs)
print("DEFINE:", define)
print("EXTRA C FLAGS:", extraCFlags)
print("EXTRA LINK FLAGS:", extraLFlags)

mdarray = Extension(
    'mdarray',
    glob.glob('mdarray/*.c'),
    include_dirs = inc_dirs,
    library_dirs = lib_dirs,
    libraries = libs,
    define_macros = define,
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
    long_description = open('README.md').read(),
    author = 'Borys Szefczyk',
    author_email = 'boryszef@gmail.com',
    url = 'https://github.com/boryszef/mdarray',
    license = 'GPL-3.0',
    install_requires = [ 'numpy>=1.10' ],
    ext_modules = [ mdarray ],
    test_suite = "setup.test_suite",
    zip_safe = False
    )
