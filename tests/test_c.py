import os
import subprocess
from sys import abiflags, version_info
from distutils.ccompiler import new_compiler
from distutils.sysconfig import get_python_inc, get_python_lib
from numpy.distutils.misc_util import get_info

#testDir = os.path.dirname(os.path.realpath(__file__))

cc = new_compiler(verbose=1)

cc.add_include_dir("mdarray")
cc.add_include_dir(get_python_inc())
cc.add_library_dir(get_python_lib())
cc.add_library('cunit')
cc.add_library('python%d.%d%s' % (version_info.major, version_info.minor, abiflags))

npymathInfo = get_info('npymath')
for inc in npymathInfo['include_dirs']:
    cc.add_include_dir(inc)
for lib in npymathInfo['library_dirs']:
    cc.add_library_dir(lib)
for lib in npymathInfo['libraries']:
    cc.add_library(lib)

# Make tests
obj = cc.compile(["tests/test_utils.c", "mdarray/utils.c", "mdarray/periodic_table.c" ])
cc.link_executable(obj, "tests/test_utils.x")

#subprocess.run([testDir+"/test_utils"], check=True)
subprocess.run(["./tests/test_utils.x"], check=True)
