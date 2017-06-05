from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext as build_pyx
import numpy
import os

os.environ['CC'] = 'gcc-7'

compile_args = ['-O3', '-fopenmp']
link_args = ['-O3', '-fopenmp']

setup(name='test',
      ext_modules=[Extension('fmmpy', ['fmmpy.pyx', 'fmm.c'],
                             include_dirs=[numpy.get_include()],
                             extra_compile_args=compile_args,
                             extra_link_args=link_args)],
      cmdclass={'build_ext': build_pyx},
      gdb_debug=True)
