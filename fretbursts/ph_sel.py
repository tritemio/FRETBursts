#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
Photon selection classes.
"""

class ph_sel:
    DAex_DAem = 99
    all = DAex_DAem
    Dex, Aex, Dem, Aem, Dex_Dem, Dex_Aem, Aex_Dem, Aex_Aem = range(8)
    Dex_DAem_Aex_Aem = 98
    Dex_DAem_Aex_Dem = 97
    Dex_Dem_Aex_DAem = 96
    Dex_Aem_Aex_DAem = 95
    Dex_Dem_Aex_Aem = 94
    Dex_Aem_Aex_Dem = 93


class Ph_sel(dict):
    def __init__(self, *args, **kwargs):
        """Object used to specify the photon selection.

        Example:
            To select all photons use `Ph_sel(Dex='DAem', Aex='DAem')` or
            `Ph_sel('all')`. To select only donor and acceptor photons
            emitted during donor excitiation use `Ph_sel(Dex='DAem')`.

        The object is a custom dict in which the only valid keys are 'Dex'
        or 'Aex', for the donor and acceptor excitation periods. Furthermore
        values for a key can only be 'DAem', 'Dem', 'Aem' or None. These
        values select donor+acceptor, donor-only, acceptor-only or none
        photons during an excitation period (defined by the key).

        The object must be instatiated specifying at least one keyword
        argument or using the string 'all' as the only argument. Calling
        `Ph_sel('all')` is equivalent to `Ph_sel(Dex='DAem', Aex='DAem')`.
        Not specifying an keyword argument is equivalet to setting it to None.
        The returned object is guaranteed to always have both the keys.

        For consistency the object should not be modified after intialization.
        """

        super(Ph_sel, self).__init__(Dex=None, Aex=None)
        if len(args) > 0:
            if len(args) > 1 or args[0].lower() != 'all':
                print args
                raise ValueError('Invalid arguments %s. The only valid '
                                  'non-keyword argument is "all".' % str(args))
            else:
                self['Dex'] = 'DAem'
                self['Aex'] = 'DAem'
        else:
            if len(kwargs) == 0:
                raise ValueError('You need to specify at least one keyword '
                                 'argument (Dex or Aex).')
            for key, value in kwargs.items():
                if key not in ['Dex', 'Aex']:
                    raise ValueError("Invalid keyword name '%s'. Valid names"
                                     " are 'Dex' or 'Aex'." % key)
                if value not in ['DAem', 'Dem', 'Aem']:
                    raise ValueError("Invalid keyword value '%s'. Valid "
                                     "values are 'DAem', 'Dem' or 'Aem'." %\
                                     value)
            self.update(**kwargs)
