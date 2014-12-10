.. _hdf5-format:

HDF5-based smFRET file format
=============================

.. module:: fretbursts

We developed an HDF5-based format called **Photon-HDF5** for smFRET
and other measurements involving series of photon timestamps.
The specifications of the Photon-HDF5 format can be found in
`Photon-HDF5 format 0.2 (draft) <http://photon-hdf5.readthedocs.org/>`_.

For a general overview on the importance of a standard file format
for smFRET see also :ref:`hdf5-advantages`.


Read and write HDF5 smFRET files
--------------------------------

To load a smFRET data contained in HDF5-Ph-Data use the
function :func:`loader.hdf5`.

Any measurements data loaded in a :class:`burstlib.Data` object can be saved
in HDF5-Ph-Data format by using the :func:`hdf5.store` function.

FRETBursts `hdf5` module
------------------------

The module :mod:`fretbursts.hdf5` provides the :func:`hdf5.store` function
and other utility functions to quickly print structure (hierarchy)
and attributes (metadata) of HDF5 files.

.. automodule:: fretbursts.hdf5
    :members:
