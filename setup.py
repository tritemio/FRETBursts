from setuptools import setup
import versioneer


project_name = 'fretbursts'
versioneer.VCS = 'git'
versioneer.versionfile_source = project_name + '/_version.py'
versioneer.versionfile_build = project_name + '/_version.py'
versioneer.tag_prefix = '' # tags are like 1.2.0
versioneer.parentdir_prefix = project_name + '-' 

# Try to use pandoc to convert markdown to rst
try:
    from pypandoc import convert
    read_md = lambda f: convert(f, to='rst', format='markdown')
except ImportError:
    print("warning: pypandoc module not found, could not convert "
          "Markdown to RST")
    read_md = lambda f: open(f, 'r').read()


setup(name = 'fretbursts',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
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

