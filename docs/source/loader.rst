.. currentmodule:: fretbursts.loader

Loader functions
================

While FRETBursts can load data files from different file formats, we
advocate using `Photon-HDF5 <http://photon-hdf5.readthedocs.org/>`_,
a file format specifically designed for freely-diffusing single-molecule
spectroscopy data.

Photon-HDF5 files can be loaded with the function :func:`photon_hdf5`,
regardless of the type of excitation or number of spots.

Single-spot μs-ALEX measurement stored in SM files can be loaded via
:func:`usalex` and single-spot ns-ALEX measurement stored in SPC files
(Beckr & Hickl) can be loaded via :func:`nsalex`.

Note that regardless of the format, for alternated excitation data,
you need to apply the alternation parameters using
:func:`alex_apply_period`.

.. contents::


List of loader functions
------------------------

.. automodule:: fretbursts.loader
    :members:


Load data manually
------------------

In case the data is available in a format not directly supported by
FRETBursts it is possible to manually create a :class:`fretbursts.burstslib.Data` object.
For example, for non-ALEX smFRET data, two arrays of same length are
needed: the timestamps and the acceptor-mask. The timestamps need to be
an int64 numpy array containing the recorded photon timestamps in arbitrary
units (usually dictated by the acquisition hardware clock period).
The acceptor-mask needs to be a numpy boolean array that is `True`
when the corresponding timestamps comes from the acceptor channel and
`False` when it comes from the donor channel. Having these arrays a
`Data` object can be manually created with::


    d = Data(ph_times_m=[timestamps], A_em=[acceptor_mask],
             clk_p=10e-9, ALEX=False, nch=1, fname='file_name')

In the previous example, we set the timestamp unit (`clk_p`) to 10~ns
and we specify that the data is not from an ALEX measurement. Creating
`Data` objects for ALEX and ns-ALEX measurements follows the same lines.
