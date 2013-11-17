from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as NP

ext_modules = [Extension("c_burstsearch", ["c_burstsearch.pyx"])]

setup(
  name = 'Burst search',
  cmdclass = {'build_ext': build_ext},
  include_dirs = [NP.get_include()],
  ext_modules = ext_modules
)
