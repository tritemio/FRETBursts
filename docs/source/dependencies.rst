FRETBursts Dependencies
=======================

For documentation purposes, this is the list of dependencies to run FRETBursts:

- Python 3.5+ or 2.7 (deprecated)
- Numpy 1.6+
- Scipy 0.17+
- Matplotlib 1.5+ or 2+, with QT4 backend (either PyQT4 or PySide) or QT5.
- PyTables 3.x. To load/save the :ref:`Photon-HDF5 <hdf5-format>`.
- lmfit 0.9.3+, used for flexible histogram fitting.
- Jupyter environment: notebook, ipython, ipywidgets.
- Pandas, for nice table representation and exporting data.

If you want to compile the cython extensions (optional) you also need:

- cython 0.20 or newer.
- a C compiler

For developing FRETBursts you should also install

- sphinx 1.3+ (we use napoleon extension) to build this documentation.
- pytest to execute the unit tests.

Note that, unless you know what you are doing, you should never install these
dependencies manually. Use a scientific python distribution like
`Continuum Anaconda <https://store.continuum.io/cshop/anaconda/>`__
instead.
