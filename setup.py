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

For more info please refer to:

- **FRETBursts: An Open Source Toolkit for Analysis of Freely-Diffusing Single-Molecule FRET**
  *Ingargiola et. al.* (2016). PLoS ONE doi: `10.1371/journal.pone.0160716 <10.1371/journal.pone.0160716>`__.


Quick links:

- `FRETBursts Homepage <http://tritemio.github.io/FRETBursts/>`_
- `FRETBursts Reference Documentation <http://fretbursts.readthedocs.org>`_
- `FRETBursts Tutorials <https://github.com/tritemio/FRETBursts_notebooks>`_

See also `Release Notes <http://fretbursts.readthedocs.io/en/latest/releasenotes.html>`__.
"""

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
    ext_modules = [Extension("burstsearch_c",
                             [project_name + "/phtools/burstsearch_c.pyx"]),
                   Extension("phrates_c",
                             [project_name + "/phtools/phrates_cy.pyx"],
                             include_dirs = ["."],)]

## Configure setup.py commands
cmdclass = versioneer.get_cmdclass()
if has_cython:
    cmdclass.update(build_ext=build_ext)
else:
    print('WARNING: No cython found. Fast routines will not be installed.')


setup(name = project_name,
      version = versioneer.get_version(),
      cmdclass = cmdclass,
      include_dirs = [np.get_include()],
      ext_modules = ext_modules,
      author = 'Antonino Ingargiola',
      author_email = 'tritemio@gmail.com',
      url          = 'http://tritemio.github.io/FRETBursts/',
      download_url = 'http://tritemio.github.io/FRETBursts/',
      install_requires = ['numpy', 'scipy', 'matplotlib', 'lmfit', 'seaborn',
                          'phconvert', 'future'],
      include_package_data = True,
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
                   'Programming Language :: Python :: 3.5',
                   'Programming Language :: Python :: 3.6',
                   'Topic :: Scientific/Engineering',
                   ],
      packages = ['fretbursts', 'fretbursts.utils', 'fretbursts.fit',
                  'fretbursts.phtools', 'fretbursts.dataload',
                  'fretbursts.tests'],
      keywords = 'single-molecule FRET smFRET burst-analysis biophysics',
      )
