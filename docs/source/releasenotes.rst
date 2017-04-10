FRETBursts Release Notes
========================

Version 0.6.2 (Apr. 2017)
--------------------------

This is a technical release that removes the hard dependency on QT
and solves some installation issues due to QT pinning on conda-forge.


Version 0.6.1 (Apr. 2017)
--------------------------

For this version of FRETBursts, conda packages are distributed for
python 2.7, 3.5, 3.6 and numpy 1.11 and 1.12. FRETBursts still works
with python 3.4 but conda packages are not provided anymore.
Python 2.7 is now deprecated. Support for python 2.7 will be removed
in a future version.

The current release includes the following changes:

- SangYoon Chung (@chungjjang80) found that the `L` argument in
  burst search was ignored and submitted a fix to the problem in
  `PR #57 <https://github.com/tritemio/FRETBursts/pull/57>`__.
  Tests were added to avoid future regressions.
- Fix access to the deprecated background attributes (introduced in 0.6).
  See `b850a5 <https://github.com/tritemio/FRETBursts/commit/b850a595033c27cc66f8f4a748b1d0bf68366750>`__.
- Add plot wrapper for 16-ch data.
- Improved example notebook showing how to export burst data.
  See `Exporting Burst Data <https://github.com/tritemio/FRETBursts/blob/49a45dd815b40602c5e754a162c66a837bbd2477/notebooks/Example%20-%20Exporting%20Burst%20Data%20Including%20Timestamps.ipynb>`__.
- Re-enable background rate caching.
  See `PR #53 <https://github.com/tritemio/FRETBursts/pull/53>`__.
- Support Path objects as filename in `loader.photon_hdf5()`.
  See `201b5c <https://github.com/tritemio/FRETBursts/commit/201b5c089eca0f0867ceb453c3c111c54a21704d>`__.
- Improve `Ph_sel` string representation, added factory method `Ph_sel.from_str`
  and added new tests.
  See `3dc5f0 <https://github.com/tritemio/FRETBursts/commit/3dc5f078c678ca3c806f49b27223a2e1cd6df64a>`__.


Version 0.6 (Jan. 2017)
-----------------------

.. module:: fretbursts.burstlib

- Improvements to the layout of 48-spot plots.
- Simplify background computation avoiding useless recomputations.
  This results in 3x speed increase in background computation
  for measurement loaded with `ondisk=True` and 30% speed increase
  when using `ondisk=False`.
  Now all background rates are stored in the dictionary :attr:`Data.bg`,
  while the mean background rate in the dictionary :attr:`Data.bg_mean`.
  The old attributes `Data.bg_*` and `Data.rate_*` have been deprecated
  and will be removed in a future release (see below).
- Fix loading files with `ondisk=True`. With this option timestamps are not
  kept in RAM but loaded spot-by-spot when needed. This option has no effect
  on single-spot measurements but will save RAM in multi-spot measurements.
- Add new plot functions
  `hist_interphoton <http://fretbursts.readthedocs.io/en/latest/plots.html#fretbursts.burst_plot.hist_interphoton>`__
  and `hist_interphoton_single <http://fretbursts.readthedocs.io/en/latest/plots.html#fretbursts.burst_plot.hist_interphoton_single>`__
  to plot the interphoton delay distribution. In previous versions the
  function `hist_bg` (and `hist_bg_single`) did the same plot but required
  the background to be fitted. `hist_interphoton*` do not require any prior
  background fit and also have a cleaner and improved API.
- Detect and handle smFRET files (no ALEX) with counts not only in D or A channels
  (`f0e33d <https://github.com/tritemio/FRETBursts/commit/f0e33d855d6dfb31c89f282b249f80d845472124>`__).
- Better error message when a burst filtering function fails
  (`c7826d <https://github.com/tritemio/FRETBursts/commit/c7826d5190a034578b1fdb9c4325f8fbfe2c01d4>`__).

Backward-incompatible changes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Effect on burst search
""""""""""""""""""""""
Version 0.6 introduced a small change in how the auto-threshold
for background estimation is computed. This results in slightly different
background rates. As a consequence, burst searches setting a threshold
as function of the background, will set a slightly different threshold and
will find different number of bursts. The difference is not dramatic,
but can result in slight numeric changes in estimated parameters.

Details of auto-threshold changes
"""""""""""""""""""""""""""""""""
The refactor included a change in how the background is computed when using
`tail_min_us='auto'`. As before, with this setting, the background is
estimated iteratively in two steps. A first raw estimation with a fixed
threshold (250us), and second estimation with a threshold function of the
rate computed in the first step. Before version 0.6, the first step estimated
a single rate for the whole measurement. Now the first-step estimation is
performed in each background period separately. As before, the second step
computes the background separately in each background period.
This change was motivated by the need to simplify the internal logic
of background estimation, and to increase the computation efficiency
and accuracy.

Background attributes
"""""""""""""""""""""
The background refactor resulted in an incompatible change in the
:attr:`Data.bg` attribute. Users upgrading to version 0.6, may need to replace
`Data.bg` with `Data.bg[Ph_sel('all')]` in their notebooks. Note that
no official FRETBursts notebook was using `Data.bg`, so most users will not be
affected.

Compatibility layer
"""""""""""""""""""
All the old background-related attributes (bg_dd, bg_ad, bg_da, bg_aa,
rate_dd, rate_ad, rate_da, rate_aa, rate_m) are still present but deprecated.
The same data is now contained in the dictionaries
:attr:`Data.bg` and :attr:`Data.bg_mean`.
When using the deprecated attributes, a message will indicate the new syntax.
If you see the deprecation warning, please update the notebook
to avoid future errors.

Details of changed attributes
"""""""""""""""""""""""""""""

Before version 0.6, `Data.bg` contained background rates
fitted for **all-photons** stream. `Data.bg` was a list of arrays:
one array per spot, one array element per background period.
In version 0.6+, `Data.bg` contains the background rates for **all** the fitted
photon streams. `Data.bg` is now a dict using `Ph_sel` objects as keys.
Each dict entry is a list of array, one array per spot and one array element
per background period. For more details please refer to the following
documentation :attr:`Data.bg` and :attr:`Data.bg_mean`.


Version 0.5.9 (Sep. 2016)
-------------------------

- Added support for pyqt and qt 5+.
- Fix burst selection with multispot data.
  See `this commit <https://github.com/tritemio/FRETBursts/commit/f05e807cbd032e748580af9cc310585bcde97e40>`__.

There may still be some glitches when using
the QT5 GUIs from the notebook, but installing (and importing) FRETBursts
does not require QT4 anymore (QT5 is the current default in anaconda).
Please report any issue.


Version 0.5.7 (Sep. 2016)
-------------------------

Refactoring and expansion of gamma and beta corrections.
Briefly, in all the places where corrected burst sizes are being computed,
we removed the `gamma1` argument and added a flag `donor_ref`.
Additionally, the values `Data.S` are now beta corrected.

These changes affected
several components as described below.

Data Class
^^^^^^^^^^

- Methods `Data.burst_sizes_ich` and `Data.burst_sizes` now accept the
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
^^^^^^^^^^^^^^

Plot functions `hist_size` and `hist_brightness` accept the new arguments
for corrected burst size (``gamma``, ``beta`` and ``donor_ref``).

Burst selection
^^^^^^^^^^^^^^^

Burst selection by `size` and `naa` accept the new arguments
for corrected burst size (``gamma``, ``beta`` and ``donor_ref``).

Burst Weights
^^^^^^^^^^^^^

Functions that accept weights don't accept the gamma1 argument anymore,
but they don't (yet) support the arguments `donor_ref` and `beta`.
As a result, for the purpose of weighting, there is only one expression
for corrected burst size (``na + gamma*nd``), with the option to add ``naa``
but without beta correction.


All these changes are covered by unit tests.

Installation via conda-forge
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
