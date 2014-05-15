#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
A class to describe the photon selection.
"""

from collections import namedtuple


class Ph_sel(namedtuple('Ph_sel', ['Dex', 'Aex'])):
    """Class that describes a selection of photons.

    This class takes two arguments 'Dex' and 'Aex'.
    Valid values for the arguments are 'DAem', 'Dem', 'Aem' or None. These
    values select, respectively, donor+acceptor, donor-only, acceptor-only
    or no photons during an excitation period (Dex or Aex).

    The class must be called with at least one keyword argument or using
    the string 'all' as the only argument. Calling `Ph_sel('all')` is
    equivalent to `Ph_sel(Dex='DAem', Aex='DAem')`.
    Not specifying a keyword argument is equivalet to setting it to None.

    Example:
        `Ph_sel(Dex='DAem', Aex='DAem')` or `Ph_sel('all')` select all photons.
        `Ph_sel(Dex='DAem')` selects only donor and acceptor photons
        emitted during donor excitation.
        `Ph_sel(Dex='Aem', Aex='Aem')` selects all the acceptor photons.
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
                             "'DAem', 'Dem' or 'Aem' (or None)." %\
                             str((Dex, Aex)))
        return super(Ph_sel, cls).__new__(cls, Dex, Aex)
