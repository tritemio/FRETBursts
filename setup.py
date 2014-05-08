
from setuptools import setup
import fretbursts

try:
    from pypandoc import convert
    read_md = lambda f: convert(f, 'rst')
except ImportError:
    print("warning: pypandoc module not found, could not convert Markdown to RST")
    read_md = lambda f: open(f, 'r').read()


setup(name = 'fretbursts',
      version = fretbursts.__version__,
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
      packages = ['fretbursts'],
      keywords = 'single-molecule FRET smFRET burst-analysis biophysics',
      )

