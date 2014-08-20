
from setuptools import setup

import codecs
import os
import re

here = os.path.abspath(os.path.dirname(__file__))

# Read the version number from a source file.
# Why read it, and not import?
# see https://groups.google.com/d/topic/pypa-dev/0PkjVpcxTzQ/discussion
def find_version(*file_paths):
    # Open in Latin-1 so that we avoid encoding errors.
    # Use codecs.open for Python 2 compatibility
    with codecs.open(os.path.join(here, *file_paths), 'r', 'latin1') as f:
        version_file = f.read()

    # The version line must have the form
    # __version__ = 'ver'
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

# Try to use pandoc to convert markdown to rst
try:
    from pypandoc import convert
    read_md = lambda f: convert(f, to='rst', format='markdown')
except ImportError:
    print("warning: pypandoc module not found, could not convert Markdown to RST")
    read_md = lambda f: open(f, 'r').read()


setup(name = 'fretbursts',
      version = find_version('fretbursts', '__init__.py'),
      author='Antonino Ingargiola',
      author_email='tritemio@gmail.com',
      url          = 'http://github.com/tritemio/FRETBursts/',
      download_url = 'http://github.com/tritemio/FRETBursts/',
      requires = ('numpy', 'scipy', 'matplotlib', 'ipython'),
      license = 'GPLv2',
      description = "Burst analysis toolkit for single and multi-spot smFRET data.",
      long_description = read_md('README.md'),
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

