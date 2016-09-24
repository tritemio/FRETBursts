FRETBursts Release Notes
========================

Version 0.5.7 (Sep. 2016) RC1
-----------------------------

Refactoring and expansion of gamma and beta corrections.

- Data methods that return burst sizes (i.e. `Data.burst_sizes_ich` and
  `Data.burst_sizes`) now accept the arguments ``gamma``, ``beta``
  and ``donor_ref``, while the argument ``gamma1`` was removed.
  The two conventions of corrected burst size are chosen with the boolean
  flag ``donor_ref``.
  See the `burst_sizes_ich docs <http://fretbursts.readthedocs.io/en/latest/data_class.html?highlight=get_naa#fretbursts.burstlib.Data.burst_sizes_ich>`__
  for details.

- New method `get_naa_corrected` returns the array of `naa` burst counts
  corrected with the passed ``gamma`` and ``beta`` values. Like for the burst
  size, the argument ``donor_ref`` selects the convention for the correction.
  See the `get_naa_corrected docs <http://fretbursts.readthedocs.io/en/latest/data_class.html?highlight=get_naa#fretbursts.burstlib.Data.get_naa_corrected>`__
  for details.

- Plot functions `hist_size` and `hist_brightness` accept the new arguments
  for corrected burst size (``gamma``, ``beta`` and ``donor_ref``).

- Burst selection by `size` and `naa` accept the new arguments
  for corrected burst size (``gamma``, ``beta`` and ``donor_ref``).

- A new `Data` attribute ``beta`` (default: 1) stores a beta value that is used
  to compute the corrected S. This value is never implicitly used to compute
  corrected burst sizes or naa (for these a `beta` arguments needs to be
  passed explicitly).

- Unit test for new features.

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
