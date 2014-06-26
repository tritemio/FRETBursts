HDF5-based smFRET file format
=============================

.. module:: fretbursts

Read and write HDF5 smFRET files
--------------------------------

.. warning::

    At the moment the HDF5 smFRET format is tentative and open to discussion. 
    However we already reserved a version attribute to be able to upgrade 
    the format in a backward compatible manner.
    
    If you have suggestions or want to contribute developing the format
    please contact us.

In FRETBursts we support reading/writing from/to HDF5 format.

To load a smFRET dataset store in HDF5 format use the loader function
:func:`loader.hdf5`.

A dataset loaded in a :class:`burstlib.Data` object can be saved in HDF5
with :func:`hdf5.store`.

.. module:: fretbursts.hdf5
.. autofunction:: store

