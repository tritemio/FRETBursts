#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
This module provides functions to store and load background fit data
to and from an HDF5 file.

The functions here assume to find an open pyTables file reference in
the :class:`Data` attribute `.bg_data_file`.

Implementation details
----------------------

Background caching only works when `bg_fun = bg.exp_fit` (MLE tail fit) and
assumes that `bg_ph_sel == Ph_sel('all')`.

Background estimation results are indentified by the `Data.calc_bg` arguments::

    time_s, tail_min_us, F_bg, error_metrics, fit_allph

These are serialized and used as group name in `/background`.
The HDF5 file stores::

    bg

the arrays::

    Lim, Ph_p

and optionally the field `bg_auto_th_us0` used when `tail_min_us == 'auto'`.
The following `Data` attributes are computed (and not stored):

- bg_th_us_user: is tail_min_us when `tail_min_us != 'auto'`
- bg_auto_F_bg: is F_bg when `tail_min_us == 'auto'`

Finally the following `Data` attributes are fixed:

- bg_fun: fixed to bg.exp_fit, not saved but implied
- bg_ph_sel: fixed to Ph_sel('all'), not saved but implied

The following properties are not stored since they are compute on-fly every time:

- nperiods
- bg_mean

Attributes not saved and not restored (so far):

- Th_us: actual threshold to select the tail of the interphoton delays
  distribution. Dict of lists just like Data.bg.
- BG_err: metric for fit error estimation. Dict of lists just like Data.bg.

"""

from __future__ import absolute_import
from builtins import range, zip

from pathlib import Path
import json
import numpy as np
import tables

from .utils.misc import pprint
from .ph_sel import Ph_sel
from .background import exp_fit


def to_signature(time_s, tail_min_us, F_bg, error_metrics, fit_allph):
    return json.dumps(dict(time_s=time_s, tail_min_us=tail_min_us, F_bg=F_bg,
                           error_metrics=error_metrics, fit_allph=fit_allph))


def from_signature(string):
    return {k: tuple(v) if isinstance(v, list) else v
            for k, v in json.loads(string).items()}


def remove_cache(dx):
    """Remove all the saved background data."""
    assert 'bg_data_file' in dx
    pprint(' * Removing all the cached background data ... ')
    if '/background' in dx.bg_data_file:
        dx.bg_data_file.remove_node('/background', recursive=True)
        dx.bg_data_file.flush()
    pprint('[DONE]\n')


def _save_bg_data(bg, Lim, Ph_p, bg_calc_kwargs, h5file, bg_auto_th_us0=None):
    """Save background data to HDF5 file."""
    # Save the bg data
    group_name = to_signature(**bg_calc_kwargs)
    if _bg_is_cached(h5file, bg_calc_kwargs):
        h5file.remove_node('/background', group_name, recursive=True)

    bg_group = h5file.create_group('/background', group_name,
                                   createparents=True)

    for ph_sel in bg:
        h5file.create_array(bg_group, 'BG_%s' % str(ph_sel), obj=bg[ph_sel])
    h5file.create_array(bg_group, 'Lim', obj=Lim)
    h5file.create_array(bg_group, 'Ph_p', obj=Ph_p)
    if bg_calc_kwargs['tail_min_us'] == 'auto':
        assert bg_auto_th_us0 is not None
        h5file.create_array(bg_group, 'bg_auto_th_us0', obj=bg_auto_th_us0)


def _load_bg_data(bg_calc_kwargs, h5file):
    """Load background data from a HDF5 file."""
    group_name = to_signature(**bg_calc_kwargs)
    if group_name not in h5file.root.background:
        msg = 'Group "%s" not found in the HDF5 file.' % group_name
        raise ValueError(msg)
    bg_auto_th_us0 = None
    bg_group = h5file.get_node('/background/', group_name)

    pprint('\n - Loading bakground data: ')
    bg = {}
    for node in bg_group._f_iter_nodes():
        if node._v_name.startswith('BG_'):
            ph_sel = Ph_sel.from_str(node._v_name[len('BG_'):])
            bg[ph_sel] = [np.asfarray(b) for b in node.read()]

    Lim = bg_group.Lim.read()
    Ph_p = bg_group.Ph_p.read()
    if 'bg_auto_th_us0' in bg_group:
        bg_auto_th_us0 = bg_group.bg_auto_th_us0.read()
    return bg, Lim, Ph_p, bg_auto_th_us0


def _bg_is_cached(h5file, bg_calc_kwargs):
    """Returns signature matches a group in /backgroung.
    """
    group_name = to_signature(**bg_calc_kwargs)
    return ('background' in h5file.root and
            group_name in h5file.root.background)


def get_h5file(dx):
    datafile = Path(dx.fname)
    cachefile = datafile.with_name(datafile.stem + '_cache.hdf5')
    return tables.open_file(str(cachefile), mode='a')


def calc_bg_cache(dx, fun, time_s, tail_min_us, F_bg, error_metrics, fit_allph,
                  recompute=False):
    """Cached version of `.calc_bg()` method."""
    assert fun == exp_fit, 'Cache only supports the bg.exp_fit function.'
    assert error_metrics is None, 'Cache only support `error_metrics=None`.'
    bg_calc_kwargs = dict(time_s=time_s, tail_min_us=tail_min_us,
                          F_bg=F_bg, error_metrics=error_metrics,
                          fit_allph=fit_allph)
    h5file = get_h5file(dx)
    if _bg_is_cached(h5file, bg_calc_kwargs) and not recompute:
        # Background found in cache. Load it.
        pprint(' * Loading BG rates from cache ... ')
        bg, Lim, Ph_p, bg_auto_th_us0 = _load_bg_data(bg_calc_kwargs, h5file)

        bg_dict = dict(bg_fun=exp_fit, bg_ph_sel=Ph_sel('all'))     # fixed
        bg_dict.update(bg=bg, Lim=Lim, Ph_p=Ph_p, bg_time_s=time_s)
        if bg_calc_kwargs['tail_min_us'] == 'auto':
            bg_dict['bg_auto_F_bg'] = bg_calc_kwargs['F_bg']
            assert bg_auto_th_us0 is not None
            bg_dict['bg_auto_th_us0'] = bg_auto_th_us0
        else:
            bg_dict['bg_th_us_user'] = tail_min_us
        dx._clean_bg_data()
        dx.add(**bg_dict)
        pprint(' [DONE]\n')
    else:
        pprint(' * Computing BG rates:\n')
        dx.calc_bg(fun=fun, **bg_calc_kwargs)
        bg_auto_th_us0 = dx.get('bg_auto_th_us0', None)
        _save_bg_data(dx.bg, dx.Lim, dx.Ph_p, bg_calc_kwargs, h5file,
                      bg_auto_th_us0=bg_auto_th_us0)
        pprint(' [DONE]\n')
    h5file.close()
