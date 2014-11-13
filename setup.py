from setuptools import setup
from setuptools.extension import Extension
from Cython.Distutils import build_ext
import numpy as np
import versioneer

ext_modules = [Extension("burstsearchlib_c",
                         ["fretbursts/burstsearch/burstsearchlib_c.pyx"])]

project_name = 'fretbursts'
versioneer.VCS = 'git'
versioneer.versionfile_source = project_name + '/_version.py'
versioneer.versionfile_build = project_name + '/_version.py'
versioneer.tag_prefix = '' # tags are like 1.2.0
versioneer.parentdir_prefix = project_name + '-'

long_description = """
FRETBursts
==========

**FRETBursts** is a software toolkit for burst analysis of confocal 
single-molecule FRET (smFRET) measurements. It can analyze both single-spot
and multi-spot smFRET data with or without alternating laser excitation (ALEX).

Quick links: 

- `FRETBursts Homepage <https://github.com/tritemio/FRETBursts>`_
- `Reference documentation <http://fretbursts.readthedocs.org/index.html>`_
- `FRETBursts tutorials <https://github.com/tritemio/FRETBursts_notebooks>`_
"""

cmdclass = versioneer.get_cmdclass()
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
      requires = ('numpy', 'scipy', 'matplotlib', 'ipython'),
      license = 'GPLv2',
      description = ("Burst analysis toolkit for single and multi-spot "
                     "smFRET data."),
      long_description = long_description,
      platforms = ('Windows', 'Linux', 'Mac OS X'),
      classifiers=['Intended Audience :: Science/Research',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python',
                   'Programming Language :: Python :: 2.7',
                   'Topic :: Scientific/Engineering',
                   ],
      packages = ['fretbursts', 'fretbursts.utils', 'fretbursts.fit',
                  'fretbursts.burstsearch', 'fretbursts.dataload'],
      keywords = 'single-molecule FRET smFRET burst-analysis biophysics',
      )

