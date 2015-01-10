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


Basic info methods
------------------

List of :class:`Data` methods that output basic information.

.. class:: Data

    .. autoattribute:: time_max

    .. autoattribute:: time_min

    .. autoattribute:: ph_data_sizes

    .. autoattribute:: num_bursts

    .. automethod:: burst_sizes

    .. automethod:: burst_sizes_ich

    .. autoattribute:: burst_widths

    .. automethod:: ph_in_bursts_ich

    .. automethod:: ph_in_bursts_mask_ich

    .. automethod:: status

    .. automethod:: name

    .. automethod:: Name


Analysis methods
----------------

Methods for background estimation, burst search and burst-data calculations.

.. class:: Data

    .. automethod:: calc_bg

    .. automethod:: burst_search

    .. automethod:: calc_fret

    .. automethod:: calc_ph_num

    .. automethod:: fuse_bursts

    .. automethod:: calc_sbr

    .. automethod:: calc_max_rate


Burst correction methods
------------------------

List of :class:`Data` methods used to apply burst corrections.

.. class:: Data

    .. automethod:: update_leakage

    .. automethod:: update_dir_ex

    .. automethod:: update_gamma

    .. automethod:: background_correction

    .. automethod:: leakage_correction

    .. automethod:: dither


Fitting methods
-------------------

Some fitting methods for burst data. Note that E and S histogram fitting
with generic models is now handled with the new
:ref:`fitting framework <fit-section>`.

.. class:: Data

    .. automethod:: fit_E_generic

    .. automethod:: fit_E_m

    .. automethod:: fit_E_ML_poiss

    .. automethod:: fit_E_minimize

    .. automethod:: fit_E_two_gauss_EM


Data access methods
-------------------

These methods are used to access (or iterate over) the different
arrays of timestamps or burst data.

.. class:: Data

    .. automethod:: get_ph_times

    .. automethod:: iter_ph_times

    .. automethod:: get_ph_mask

    .. automethod:: iter_ph_masks

    .. automethod:: iter_bursts_ph

    .. automethod:: expand

    .. automethod:: copy

    .. automethod:: slice_ph
