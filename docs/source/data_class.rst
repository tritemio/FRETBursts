.. currentmodule:: fretbursts.burstlib

The "Data()" class
==================

The :class:`Data` class is the main container for smFRET measurements.
It contains timestamps, detectors and all the results of data processing
such as background estimation, burst data, fitted FRET and so on.

The reference documentation of the class follows.

.. contents::


"Data()" class: description and attributes
------------------------------------------

A description of the :class:`Data` class and its main attributes.

.. module:: fretbursts.burstlib
.. autoclass:: Data


Analysis methods
----------------

List of :class:`Data` methods used to perform different analysis.

.. class:: Data

    .. automethod:: calc_bg

    .. automethod:: burst_search_t

    .. automethod:: calc_fret

    .. automethod:: calc_ph_num

    .. automethod:: fuse_bursts


Basic info methods
------------------

List of :class:`Data` methods that output basic information.

.. class:: Data

    .. automethod:: time_max

    .. automethod:: num_bursts

    .. automethod:: burst_sizes_ich

    .. automethod:: status

    .. automethod:: name

    .. automethod:: Name


Burst corrections methods
-------------------------

List of :class:`Data` methods used to apply burst corrections.

.. class:: Data

    .. automethod:: update_leakage

    .. automethod:: update_dir_ex

    .. automethod:: update_gamma

    .. automethod:: background_correction_t

    .. automethod:: leakage_correction

    .. automethod:: dither


Other burst methods
-------------------

List of :class:`Data` methods not falling in previous categories.

.. class:: Data

    .. automethod:: calc_sbr

    .. automethod:: calc_max_rate

    .. automethod:: fit_E_generic

    .. automethod:: fit_E_m

    .. automethod:: fit_E_ML_poiss

    .. automethod:: fit_E_minimize

    .. automethod:: fit_E_two_gauss_EM


Utility methods
---------------

List of :class:`Data` methods used to get (or iterate over) the different
arrays of timestamps or burst data.

.. class:: Data

    .. automethod:: get_ph_times

    .. automethod:: iter_ph_times

    .. automethod:: get_ph_mask

    .. automethod:: iter_ph_masks

    .. automethod:: expand

    .. automethod:: copy

    .. automethod:: slice_ph


Photon selection
----------------

.. module:: fretbursts.ph_sel

The class :class:`Ph_sel` is used to specify which sub-set of
photons/timestamps are "selected" (i.e. all-photons, Donor-excitation-period
photons, etc...).

.. autoclass:: Ph_sel

