from numpy.distutils.core import setup, Extension
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

extraCFlags = [ "-O0 -g" ]
extraLFlags = [ "-O0 -g" ]

# Check if gromacs is present
if extraPackagePresent('libgromacs'):
    cflags, libs = getPackageFlags('libgromacs')
    extraCFlags.append('-DHAVE_GROMACS')
    extraCFlags.append(cflags)
    extraLFlags.append(libs)

moltools = Extension('moltools',
	sources = [ 'eam.c', 'periodic_table.c', 'utils.c', 'readers.c',
	            'measure.c', 'constants.c', 'writers.c', 'topology.c',
	            'molecule.c', 'moltools.c', ],
	extra_compile_args=extraCFlags,	extra_link_args=extraLFlags,
)

setup (name = 'moltools-python', version = '0.1.1',
       description = 'Python module for manipulation of atomic coordinates',
       author = 'Borys Szefczyk', author_email = 'boryszef@gmail.com',
       url = 'https://github.com/boryszef/moltools-python',
       ext_modules = [ moltools ] )
