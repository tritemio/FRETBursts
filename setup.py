from setuptools import setup
from setuptools.extension import Extension
import numpy as np
import versioneer

## Metadata
project_name = 'fretbursts'
long_description = """
FRETBursts
==========

**FRETBursts** is a software toolkit for burst analysis of confocal
single-molecule FRET (smFRET) measurements. It can analyze both single-spot
and multi-spot smFRET data with or without alternating laser excitation (ALEX).

Quick links:

- `FRETBursts Homepage <http://tritemio.github.io/FRETBursts/>`_
- `FRETBursts Reference Documentation <http://fretbursts.readthedocs.org>`_
- `FRETBursts Tutorials <https://github.com/tritemio/FRETBursts_notebooks>`_
"""

## Configure versioneer
versioneer.VCS = 'git'
versioneer.versionfile_source = project_name + '/_version.py'
versioneer.versionfile_build = project_name + '/_version.py'
versioneer.tag_prefix = '' # tags are like 1.2.0
versioneer.parentdir_prefix = project_name + '-'

## Configuration to build Cython extensions
try:
    from Cython.Distutils import build_ext
except ImportError:
    # cython is not installed: do not build extensions
    has_cython = False
    ext_modules = []
else:
    # cython is installed: register the extensions to be built
    has_cython = True
    ext_modules = [Extension("burstsearchlib_c",
                             [project_name + \
                             "/burstsearch/burstsearchlib_c.pyx"])]

## Configure setup.py commands
cmdclass = versioneer.get_cmdclass()
if has_cython:
    cmdclass.update(build_ext=build_ext)


setup(name = project_name,
      version = versioneer.get_version(),
      cmdclass = cmdclass,
      include_dirs = [np.get_include()],
      ext_modules = ext_modules,
      author = 'Antonino Ingargiola',
      author_email = 'tritemio@gmail.com',
      url          = 'http://github.com/tritemio/FRETBursts/',
      download_url = 'http://github.com/tritemio/FRETBursts/',
      install_requires = ['numpy', 'scipy', 'matplotlib', 'lmfit', 'seaborn',
                          'phconvert', 'future'],
      license = 'GPLv2',
      description = ("Burst analysis toolkit for single and multi-spot "
                     "smFRET data."),
      long_description = long_description,
      platforms = ('Windows', 'Linux', 'Mac OS X'),
      classifiers=['Intended Audience :: Science/Research',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3.4',
                   'Topic :: Scientific/Engineering',
                   ],
      packages = ['fretbursts', 'fretbursts.utils', 'fretbursts.fit',
                  'fretbursts.burstsearch', 'fretbursts.dataload',
                  'fretbursts.tests'],
      keywords = 'single-molecule FRET smFRET burst-analysis biophysics',
      )

