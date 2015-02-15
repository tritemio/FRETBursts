.. _burst_selection:

Burst selection
===============

.. module:: fretbursts.burstlib

After performing a burst search is common to select bursts according to
different criteria (burst size, FRET efficiency, etc...).

In FRETBursts this can be easily accomplished using the method
:meth:`Data.select_bursts`. This method takes
a :ref:`selection function <selection_functions>`
as parameters. :meth:`Data.select_bursts` returns a new :class:`Data` object
containing only the new sub-set of bursts. A new selection can be applied
to this new object as well. In this way, different
selection criteria can be freely combined in order to obtain a
burst population satisfying arbitrary constrains.

FRETBursts provides a large number of
:ref:`selection functions <selection_functions>`. Moreover, creating a new
selection function is extremely simple, requiring (usually) 2-3 lines of code.
You can take the functions in `select_bursts.py` as examples to create your
own selection rule.

In the next section we list all the selection functions. You may also want
to check the :class:`Data` methods that deal with burst selection:

- :meth:`Data.select_bursts`
- :meth:`Data.select_bursts_mask`
- :meth:`Data.select_bursts_mask_apply`


.. _selection_functions:

Selection functions
-------------------

.. automodule:: fretbursts.select_bursts
    :members:
