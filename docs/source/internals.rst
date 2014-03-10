
The ``Data()`` object
=====================

This is the object that contains all the information of a dataset (name,
timestamps, bursts, correction factors) and provides several methods to
perform analisys (background estimation, burst search, FRET fitting,
etc...).

It's a dictionary in which the items are also attributes. To add or
delete attributes use method ``.add()`` or ``delete()``.

How to copy ``Data()``
----------------------

To create a copy of a variable of type ``Data()`` you can create
(istantiate) a new ``Data()`` variable (``d_new``) using the old
variable (``d_old``) for intialization:

::

    d_new = Data(**d_old)

In this way each field of the new variable will point to the old
variable data. You can reassign attributes in one variable without
influencing the other. However if you modify the content of an attribute
(change some elements of an array) the modification will be reflected in
both variables.

To copy to a new variable duplicating all the data (so that will double
the RAM usage but the data is completely separated) use the method
``.copy()``:

::

    d_new = d_old.copy()


List of ``Data()`` attributes
=============================

***NOTE:*** In the following attributes marked with *(list)* contain one
element per channel.

Measurement attributes
----------------------

-  ``nch``: number of channels
-  ``clk_p``: clock period in seconds (for ``ph_times``)
-  ``ph_times_m``, ``A_em``: (list) ph times and relative bool mask for
   acceptor ph
-  ``BT``: bleed-through or leakage percentage
-  ``gamma``: gamma factor may be scalar or same size as nch

ALEX Specific:
^^^^^^^^^^^^^^

-  ``D_em``, ``A_em``, ``D_ex``, ``A_ex``, they are lists (1 per ch) and
   each element is a boolean mask for ``ph_times_m[i]``.
-  ``D_ON``, ``A_ON``: tuples of int (start-end values) for donor ex.
   and acceptor ex. selection.
-  ``switch_window``: *(int)* lenght of the alternation window in clk
   cycles.


BG Attributes
-------------

-  ``bg``, ``bg_dd``, ``bg_ad``, ``bg_aa``: *(list)* bg for each channel
   calculated every X sec
-  ``nperiods``: number of periods in which ph are splitted for bg
   calculation
-  ``bg_fun``: function used for bg calculation
-  ``Lim``: *(list)* for each ch, is a list of pairs of index of
   ``.ph_times_m[i]`` that identify first and last photon in each period
-  ``Ph_p``: *(list)* for each ch, is a list of pairs of arrival time
   for the first and the last photon in each period

Old attributes (now just the per-ch mean of ``bg``, ``bg_dd``, ``bg_ad``
and ``bg_aa``):

::

    `rate_m`: array of bg rates for D+A channel pairs (ex. 4 for 4 spots)
    `rate_dd`: array of bg rates for D em (and D ex if ALEX)
    `rate_da`: array of bg rates for A em (and D ex if ALEX)
    `rate_aa`: array of bg rates for A em and A ex (only for ALEX)


Burst search parameters (user input)
------------------------------------

Burst search introduction
~~~~~~~~~~~~~~~~~~~~~~~~~

A burst-start is detected if :math:`m` consecutive photons/timestamps
lie within an time interval :math:`T`. Basically the minimum rate for
burst start is:

.. math:: Min\, Rate = \frac{m}{T}

:math:`T` (and thus the minimum rate) is computed from the BG rate
according to different criteria. When the :math:`m`-photon rate falls
below the minimmu rate the burst-end is marked. Then, if the total
number of photons in burst is :math:`>L` the burst is saved, otherwise
is rejected. Setting :math:`L = m` means that no burst is rejected.

Parameters
~~~~~~~~~~

-  ``ph_sel``: *(string)* type of ph selection for burst search: 'D',
   'A' or 'DA' (default)
-  ``m``, ``L`` : parameters for burt search
-  ``P``: 1-prob. to have a burst start due to BG (assuming Poisson
   distrib.).
-  ``F``: multiplying factor for BG used to compute ``TT`` and ``T``
   (see next section)

If ``P`` is ``None`` then values in ``TT`` are computed as
:math:`T = \frac{m}{F \cdot BG}` otherwise they are computed from the
`Erlang
distribution <http://en.wikipedia.org/wiki/Erlang_distribution>`__ with
probability ``P`` to have a burst start from a Poisson process with rate
:math:`F\cdot BG` and considering :math:`m` consecutive timestamps (in
other words, the shape parameter :math:`k` of the Erlang distribution is
:math:`m`).

*See also:*
~~~~~~~~~~~

-  ``burst_search_t()``: method of ``Data()`` to perform the burst
   search
-  ``poisson_threshold.py``: file where the function to compute ``TT``
   from the Erlang distribution is defined.


Burst search data (available after burst search)
------------------------------------------------

-  ``mburst``: *(list)* array containing burst data:
   ``[tstart, width, #ph, istart]``
-  ``TT``: (same size as .bg) T values (in sec.) for burst search
-  ``T``: *(array)* per-channel mean of TT parameter

-  ``nd``,\ ``na``,\ ``nt`` : *(list)* number of donor, acceptor and
   total ph in each burst,
-  ``bp``: *(list)* index of the time period in which the burst happens.
   Same length as nd. Needed to identify which bg value to use.
-  ``bg_bs`` *(list)*: BG used for threshold in burst search (points to
   ``bg``, ``bg_dd`` or ``bg_ad``)

-  ``fuse``: is not None, contains the parameters used for fusing bursts

-  ``E``: *(list)* FRET efficiency value for each burst
   :math:`E = n_a/(n_a+n_d)`
-  ``S``: *[ALEX only] (list)* stochiometry value for each burst
   :math:`S = n_t/(n_t+n_{aa})`
    these are eventually BG and Lk corrected.
-  ``naa``: *[ALEX only] (list)* number of ph in the A (acceptor) ch
   during A excitation in each burst

**NOTE**: Burst data such as ``nd``, ``na``, ``nt``, ``E``, ``S``,
etc..., are list (one element per channel) of arrays. Each array
contains an element for each burst.

``ph_sel`` attribute
====================


Question:
~~~~~~~~~

The **``ph_sel``** attribute in ``Data()`` is something related to the
burst search. Why is set also on ``calc_bg_cache()`` when an old bg
cache file is opened? And why is saved by ``calc_bg_cache()``?

Answer
~~~~~~

When BG is computed we also compute **``Ph_p``** and **``Lim``** (resp.
timestamps and index): they define how ``ph_times`` is sliced to compute
the (time-dependent) BG.

So, if **``ph_sel == 'DA'``** then ``Lim`` and ``Ph_p`` are index and
elements of ``ph_times_m[ich]``, otherwise they are index (and elements)
of ``ph_times_m[ich][a_em]`` (or ``d_em``).
