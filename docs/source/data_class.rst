The "Data()" class
==================

The :class:`Data` class is the main container for smFRET measurements.
It contains timestamps, detectors and all the results of data processing
such as background estimation, burst data, fitted FRET and so on.

The reference documentation of the class follows.

Description and attributes
--------------------------

.. module:: fretbursts.burstlib
.. autoclass:: Data


Analysis methods
----------------

.. currentmodule:: fretbursts.burstlib

.. class:: Data

    .. automethod:: calc_bg

    .. automethod:: burst_search_t

    .. automethod:: calc_fret

    .. automethod:: calc_ph_num

    .. automethod:: fuse_bursts


Basic info methods
------------------

.. class:: Data

    .. automethod:: time_max

    .. automethod:: num_bu

    .. automethod:: status

    .. automethod:: name

    .. automethod:: Name


Burst corrections methods
-------------------------

.. class:: Data

    .. automethod:: update_bt

    .. automethod:: update_gamma

    .. automethod:: background_correction_t

    .. automethod:: bleed_through_correction

    .. automethod:: dither


Other burst methods
-------------------

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

.. class:: Data

    .. automethod:: get_ph_times

    .. automethod:: iter_ph_times

    .. automethod:: get_ph_mask

    .. automethod:: iter_ph_masks

    .. automethod:: expand

    .. automethod:: copy

    .. automethod:: slice_ph

