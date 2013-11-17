from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("c_fuse", ["c_fuse.pyx"])]

setup(
  name = 'Burst fusing',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
