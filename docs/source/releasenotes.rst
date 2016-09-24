FRETBursts Release Notes
========================

Version 0.5.7 (Sep. 2016)
-------------------------

Refactoring and expansion of gamma and beta corrections.

- Methods that return burst size now accept the arguments ``gamma``, ``beta``
  and ``donor_ref``, while the argument ``gamma1`` was removed.
  The two conventions of corrected burst size are chosen with the boolean
  flag ``donor_ref``. See the docs for details.

- New method `get_naa_corrected` returns the array of `naa` burst counts
  corrected with the passed ``gamma`` and ``beta`` values. Like for the burst
  size, the argument ``donor_ref`` selects the convention for the correction.
  See docs for details.

- Plot functions `hist_size` and `hist_brightness` accept the new arguments
  for corrected burst size (``gamma``, ``beta`` and ``donor_ref``).

- Burst selection by `size` and `naa` accept the new arguments
  for corrected burst size (``gamma``, ``beta`` and ``donor_ref``).

- A new `Data` attribute ``beta`` (default: 1) stores a beta value that is used
  to compute the corrected S. This value is never implicitly used to compute
  corrected burst sizes or naa (for these a `beta` arguments needs to be
  passed explicitly).

- Unit test for new features.


Version 0.5.6
-------------

See `Release Note on GitHub <https://github.com/tritemio/FRETBursts/releases/tag/0.5.6>`__.
