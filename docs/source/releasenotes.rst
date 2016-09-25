FRETBursts Release Notes
========================

Version 0.5.7 (Sep. 2016)
-------------------------

Refactoring and expansion of gamma and beta corrections.
Briefly, in all the places where corrected burst sizes are being computed,
we removed the `gamma1` argument and added a flag `donor_ref`.
Additionally, the values `Data.S` are now beta corrected.

These changes affected
several components as described below.

Data Class
~~~~~~~~~~

- Data methods `Data.burst_sizes_ich` and `Data.burst_sizes`) now accept
  arguments ``gamma``, ``beta`` and ``donor_ref``. The argument ``gamma1``
  was removed.
  The two conventions of corrected burst sizes are chosen with the boolean
  flag ``donor_ref``.
  See the `burst_sizes_ich docs <http://fretbursts.readthedocs.io/en/latest/data_class.html?highlight=get_naa#fretbursts.burstlib.Data.burst_sizes_ich>`__
  for details.

- New method `get_naa_corrected` returns the array of `naa` burst counts
  corrected with the passed ``gamma`` and ``beta`` values. Like for the burst
  size, the argument ``donor_ref`` selects the convention for the correction.
  See the `get_naa_corrected docs <http://fretbursts.readthedocs.io/en/latest/data_class.html?highlight=get_naa#fretbursts.burstlib.Data.get_naa_corrected>`__
  for details.

- A new `Data` attribute ``beta`` (default: 1) stores a beta value that is used
  to compute the corrected S. This value is never implicitly used to compute
  corrected burst sizes or naa (for these a `beta` arguments needs to be
  passed explicitly).


Plot functions
~~~~~~~~~~~~~~

Plot functions `hist_size` and `hist_brightness` accept the new arguments
for corrected burst size (``gamma``, ``beta`` and ``donor_ref``).

Burst selection
~~~~~~~~~~~~~~~

Burst selection by `size` and `naa` accept the new arguments
for corrected burst size (``gamma``, ``beta`` and ``donor_ref``).

Burst Weights
~~~~~~~~~~~~~

Functions that accept weights don't accept the gamma1 argument anymore,
but they don't (yet) support the arguments `donor_ref` and `beta`.
As a result, for the purpose of weighting, there is only one expression
for corrected burst size (``na + gamma*nd``), with the option to add ``naa``
but without beta correction.


All these changes are covered by unit tests.

Installation via conda-forge
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since version 0.5.6 we started distributing conda packages for FRETBursts
through the `conda-forge <https://conda-forge.github.io/>`__ channel
(a community supported repository, as opposed to a private channel we were using before).
To install or update FRETBursts you should now use::

    conda install fretbursts -c conda-forge

Using the conda-forge channel simplifies our release process since
their infrastructure automatically builds packages for multiple
platforms and python versions. Please report any issues in installing
or upgrading FRETBursts on the
`GitHub Issues <https://github.com/tritemio/FRETBursts/issues>`__ page.

For more detailed installation instructions see the
`Getting Started <http://fretbursts.readthedocs.io/en/latest/getting_started.html>`__
documentation.


Version 0.5.6
-------------

For older release notes see  `GitHub Releases Page <https://github.com/tritemio/FRETBursts/releases/>`__.
