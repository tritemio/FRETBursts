FRETBursts Dependencies
=======================

For documentation purposes, this is the list of dependencies to run FRETBursts:

- Python 2.7
- Numpy 1.6 or newer
- Scipy 0.13 or newer
- Matplotlib 1.3 or newer, with QT4 backend (either PyQT4 or PySide).
- IPython 1.x (2.x recommended)
- PyTables 3.x. To load/save the :ref:`Photon-HDF5 <hdf5-format>`.
- lmfit 0.8 or newer, used for flexible histogram fitting.
- Pandas, currently used just for nice table representations.

If you want to compile the cython extensions (optional) you also need:

- cython 0.20 or newer.
- a C compiler

For developing FRETBursts you should also install

- sphinx 1.2.2 with the napoleon extension (sphinxcontrib-napoleon) to
  build this documentation.
- pytest to execute the unit tests.

Note that, unless you know what you are doing, you should never install these
dependencies manually. Use a scientific python distribution like
`Continuum Anaconda <https://store.continuum.io/cshop/anaconda/>`__
instead.
