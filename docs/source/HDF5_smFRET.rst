.. _hdf5-format:

HDF5-based smFRET file format
=============================

.. module:: fretbursts

We developed an HDF5-based format called **Photon-HDF5** for smFRET
and other measurements involving series of photon timestamps.
The specifications of the Photon-HDF5 format can be found in
`Photon-HDF5 format <http://photon-hdf5.readthedocs.org/>`_.

For a general overview on the importance of a standard file format
for smFRET see also :ref:`hdf5-advantages`.


Read and write HDF5 smFRET files
--------------------------------

To load a smFRET data contained in HDF5-Ph-Data use the
function :func:`loader.photon_hdf5`.

You can convert files from any format to Photon-HDF5 by using
`phconvert <https://github.com/tritemio/phconvert>
(already pre-installed with FRETBursts).

