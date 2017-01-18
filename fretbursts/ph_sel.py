#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014-2016 The Regents of the University of California,
#               Antonino Ingargiola <tritemio@gmail.com>
#
"""
In this module we define the class :class:`Ph_sel` used to specify a
"selection" of a sub-set of photons/timestamps  (i.e. all-photons,
Donor-excitation-period photons, etc...).

A photon selection is one of the base *photon streams* or a combination of
them. Base *photon streams* are photon from the donor (or acceptor) emission
channel detected during the donor (or acceptor) excitation period. For
non-ALEX data there is only the donor excitation period.

The following table shows base *photon streams* for smFRET data (non-ALEX):

=================   =====================
Photon selection    Syntax
=================   =====================
D-emission          `Ph_sel(Dex='Dem')`
A-emission          `Ph_sel(Dex='Aem')`
=================   =====================

and for ALEX data:

==============================   =====================
Photon selection                 Syntax
==============================   =====================
D-emission during D-excitation   `Ph_sel(Dex='Dem')`
A-emission during D-excitation   `Ph_sel(Dex='Aem')`
D-emission during A-excitation   `Ph_sel(Aex='Dem')`
A-emission during A-excitation   `Ph_sel(Aex='Aem')`
==============================   =====================

Additionally, all the photons can be selected with `Ph_sel('all')` (that is a
shortcut for  `Ph_sel(Dex='DAem', Aex='DAem')`.

Examples:

    - `Ph_sel(Dex='DAem', Aex='DAem')` or `Ph_sel('all')` select all photons.
    - `Ph_sel(Dex='DAem')` selects only donor and acceptor photons
      emitted during donor excitation. These are all the photons for
      non-ALEX data.
    - `Ph_sel(Dex='Aem', Aex='Aem')` selects all the photons detected from
      the acceptor-emission channel.

The documentation for the :class:`Ph_sel` class follows.

"""

from collections import namedtuple


# Implementation Rationale:
#
#   This class is implemented as a named tuple that is an immutable object.
#   This means the Ph_sel objects can be used as dictionary keys. This would
#   not be possible if Ph_sel were a dict. The __new__ method is just a nicety
#   to allow only valid arguments values and to throw meaningful exceptions
#   otherwise.
#
class Ph_sel(namedtuple('Ph_sel', ['Dex', 'Aex'])):
    """Class that describes a selection of photons.

    This class takes two arguments `Dex` and `Aex`.
    Valid values for the arguments are the strings 'DAem', 'Dem', 'Aem' or
    `None`. These values select, respectively, donor+acceptor, donor-only,
    acceptor-only or no photons during an excitation period (`Dex` or `Aex`).

    The class must be called with at least one keyword argument or using
    the string 'all' as the only argument. Calling `Ph_sel('all')` is
    equivalent to `Ph_sel(Dex='DAem', Aex='DAem')`.
    Not specifying a keyword argument is equivalent to setting it to None.

    """
    valid_values = ('DAem', 'Dem', 'Aem', None)

    def __new__(cls, Dex=None, Aex=None):
        if Dex is None and Aex is None:
            raise ValueError("You need to specify at least one argument "
                             "(Dex, Aex or 'all').")
        if Dex == 'all':
            return super(Ph_sel, cls).__new__(cls, 'DAem', 'DAem')
        if Dex not in cls.valid_values or Aex not in cls.valid_values:
            raise ValueError("Invalid value %s. Valid values are "
                             "'DAem', 'Dem' or 'Aem' (or None)." %
                             str((Dex, Aex)))
        return super(Ph_sel, cls).__new__(cls, Dex, Aex)

    @classmethod
    def _get_str_mapping(cls, invert=False):
        mapping = {cls('all'): 'all',
                   cls(Dex='Dem'): 'DexDem', cls(Dex='Aem'): 'DexAem',
                   cls(Aex='Aem'): 'AexAem', cls(Aex='Dem'): 'AexDem',
                   cls(Dex='DAem'): 'Dex', cls(Aex='DAem'): 'Aex',
                   cls(Dex='Dem', Aex='Dem'): 'Dem',
                   cls(Dex='Aem', Aex='Aem'): 'Aem',
                   cls(Dex='DAem', Aex='Aem'): 'DexDAem_AexAem'}
        if invert:
            mapping = {v: k for k, v in mapping.items()}
        return mapping

    @classmethod
    def from_str(cls, ph_sel_str):
        str_to_ph_sel = cls._get_str_mapping(invert=True)
        return str_to_ph_sel[ph_sel_str]

    def __str__(self):
        labels = self._get_str_mapping()
        return labels.get(self, repr(self))
