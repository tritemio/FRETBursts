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


Summary information
-------------------

List of :class:`Data` attributes and
methods providing summary information on the measurement:

.. class:: Data

    .. autoattribute:: time_max

    .. autoattribute:: time_min

    .. autoattribute:: ph_data_sizes

    .. autoattribute:: num_bursts

    .. automethod:: burst_sizes

    .. automethod:: burst_sizes_pax_ich

    .. automethod:: burst_sizes_ich

    .. automethod:: get_naa_corrected

    .. autoattribute:: burst_widths

    .. automethod:: ph_in_bursts_ich

    .. automethod:: ph_in_bursts_mask_ich

    .. automethod:: status

    .. autoattribute:: name

    .. automethod:: Name


Analysis methods
----------------

The following methods perform background estimation, burst search and
burst-data calculations:

- :meth:`Data.calc_bg`
- :meth:`Data.burst_search`
- :meth:`Data.calc_fret`
- :meth:`Data.calc_ph_num`
- :meth:`Data.fuse_bursts`
- :meth:`Data.calc_sbr`
- :meth:`Data.calc_max_rate`

The methods documentation follows:

.. class:: Data

    .. automethod:: calc_bg

    .. automethod:: burst_search

    .. automethod:: calc_fret

    .. automethod:: calc_ph_num

    .. automethod:: fuse_bursts

    .. automethod:: calc_sbr

    .. automethod:: calc_max_rate


Burst corrections
-----------------

Correction factors
^^^^^^^^^^^^^^^^^^

The following are the various burst correction factors. They are `Data`
properties, so setting their value automatically updates all the burst
quantities (including `E` and `S`).

.. class:: Data

    .. autoattribute:: gamma

    .. autoattribute:: leakage

    .. autoattribute:: dir_ex

    .. autoattribute:: chi_ch



Correction methods
^^^^^^^^^^^^^^^^^^

List of :class:`Data` methods used to apply burst corrections.

.. class:: Data

    .. automethod:: background_correction

    .. automethod:: leakage_correction

    .. automethod:: dither



Burst selection methods
-----------------------

:class:`Data` methods that allow to filter bursts according to different rules.
See also :ref:`burst_selection`.

.. class:: Data

    .. automethod:: select_bursts

    .. automethod:: select_bursts_mask

    .. automethod:: select_bursts_mask_apply


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

The following methods are used to access (or iterate over) the arrays of timestamps
(for different photon streams), timestamps masks and burst data.

- :meth:`Data.get_ph_times`
- :meth:`Data.iter_ph_times`
- :meth:`Data.get_ph_mask`
- :meth:`Data.iter_ph_masks`
- :meth:`Data.iter_bursts_ph`
- :meth:`Data.expand`
- :meth:`Data.copy`
- :meth:`Data.slice_ph`

The methods documentation follows:

.. class:: Data

    .. automethod:: get_ph_times

    .. automethod:: iter_ph_times

    .. automethod:: get_ph_mask

    .. automethod:: iter_ph_masks

    .. automethod:: iter_bursts_ph

    .. automethod:: expand

    .. automethod:: copy

    .. automethod:: slice_ph
