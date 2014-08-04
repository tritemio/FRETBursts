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
"""

import os
import numpy as np
from fretbursts.utils.misc import pprint


def remove_cache(dx):
    """Remove all the saved background data."""
    assert 'bg_data_file' in dx
    pprint(' * Removing all the cached background data ... ')
    if '/background' in dx.bg_data_file:
        dx.bg_data_file.remove_node('/background', recursive=True)
        dx.bg_data_file.flush()
    pprint('[DONE]\n')

def _get_bg_arrays_info(dx):
    """Return a dict of background numeric data description.

    The keys of the returned dict are names and the values are a string
    description of numeric data related to the background fit.
    This info is used to set the title attribute in the HDF5 arrays.
    """
    bg_arrays_info = dict(
        bg = 'BG rate in each CH vs time (Total)',
        bg_dd = 'BG rate in each CH vs time (D_em during D_ex)',
        bg_ad = 'BG rate in each CH vs time (A_em during D_ex)',
        bg_da = 'BG rate in each CH vs time (D_em during A_ex)',
        bg_aa = 'BG rate in each CH vs time (A_em during A_ex)',
        Lim = 'Index of first and last timestamp in each period',
        Ph_p = 'First and last timestamp in each period',
        nperiods = 'Number of time periods in which BG is computed',
        bg_time_s = 'Time duration of the period (windows in which computing BG)',
        bg_auto_th = '1 if the bg threshold was computed automatically, else 0'
    )

    bg_arrays_info_auto = dict(
        bg_auto_th_us0 = 'Threshold used for the initial bg rate estimation.',
        bg_auto_F_bg = 'Factor that multiplies the initial rate estimation to '
                       'get the auto threshold',
    )
    bg_arrays_info_noauto = dict(
        bg_th_us_user = ('Waiting time thresholds for BG fit. '
                         'This array contains 5 thresholds for different '
                         'photon selections. In the order: all photons, '
                         ' D_em-D_ex, A_em-D_ex, D_em-A_ex, A_em-A_ex.')
    )

    if dx.bg_auto_th:
        bg_arrays_info.update(**bg_arrays_info_auto)
    else:
        bg_arrays_info.update(**bg_arrays_info_noauto)
    return bg_arrays_info

def bg_save_hdf5(dx):
    """Save background data to HDF5 file."""
    assert 'bg_data_file' in dx
    assert 'bg' in dx

    bg_arrays_info = _get_bg_arrays_info(dx)
    bg_attr_names = ['bg_fun', 'bg_fun_name']
    if '/background' not in dx.bg_data_file:
        dx.bg_data_file.create_group('/', 'background',
                                  title='Background estimation data')
    ## Save the bg data
    group_name = _get_bg_groupname(dx)
    if group_name in dx.bg_data_file:
        dx.bg_data_file.remove_node(group_name, recursive=True)

    bg_group = dx.bg_data_file.create_group(os.path.dirname(group_name),
                                         os.path.basename(group_name),
                                         createparents=True)
    # Save arrays and scalars
    pprint('\n - Saving arrays/scalars: ')
    for name, info in bg_arrays_info.items():
        pprint(name + ', ')
        arr = np.array(dx[name])
        dx.bg_data_file.create_array(bg_group, name, obj=arr, title=info)

    # Save the attributes
    pprint('\n - Saving HDF5 attributes: ')
    for attr in bg_attr_names:
        pprint(attr + ', ')
        bg_group._v_attrs[attr] = dx[attr]
    pprint('\n')
    dx.bg_data_file.flush()

def bg_load_hdf5(dx, group_name):
    """Load background data from a HDF5 file."""
    assert 'bg_data_file' in dx
    if group_name not in dx.bg_data_file:
        print 'Group "%s" not found in the HDF5 file.' % group_name
        return

    ## Load the bg data
    bg_arrays = dict()
    bg_attrs = dict()

    bg_group = dx.bg_data_file.get_node(group_name)

    # Load arrays and scalars
    pprint('\n - Loading arrays/scalars: ')
    for node in bg_group._f_list_nodes():
        name = node.name
        pprint(name + ', ')
        #title = node.title
        arr = bg_group._f_get_child(name)
        bg_arrays[name] = arr.read()
    dx.add(**bg_arrays)

    # Load the attributes
    pprint('\n - Loading HDF5 attributes: ')
    for attr in bg_group._v_attrs._f_list():
        pprint(attr + ', ')
        bg_attrs[attr] = bg_group._v_attrs[attr]
    dx.add(**bg_attrs)

    pprint('\n - Generating additional fields: ')
    in_map = ['', '_dd', '_ad', '_da', '_aa']
    out_map = ['_m', '_dd', '_ad', '_da', '_aa']
    new_attrs = {}
    for in_s, out_s in zip(in_map, out_map):
        assert 'bg' + in_s in dx
        pprint('bg' + in_s + ', ')
        new_attrs['rate' + out_s] = [bg.mean() for bg in dx['bg' + in_s]]
    dx.add(**new_attrs)
    pprint('\n')

def _get_bg_groupname(dx, time_s=None):
    """Get the HDF5 group name for the background data.

    Arguments:
        dx (Data): object containing the measurement data
        time_s (None or float): if None the value is taken from dx.bg_time_s
    Returns:
        A string for the background group name.
    """
    if time_s is None and 'bg' not in dx:
        print 'You need to compute the background or provide time_s.'
        return
    time_slice = time_s if time_s is not None else dx.bg_time_s
    return '/background/time_%ds' % time_slice

def _bg_is_cached(dx, signature):
    """Returns wheter background with given `signature` is in disk the cache.
    """
    if 'bg_data_file' in dx:
        bg_groupname = _get_bg_groupname(dx, time_s=signature['time_s'])
        if bg_groupname in dx.bg_data_file:
            # At least we have a group with the right time_s
            # Check whether its signature matches the current one
            group = dx.bg_data_file.get_node(bg_groupname)
            if signature == group._v_attrs.signature:
                return True
    return False

def _bg_add_signature(dx, signature):
    """Add the signature attr to current background group.
    """
    assert 'bg_data_file' in dx
    bg_groupname = _get_bg_groupname(dx, time_s=signature['time_s'])
    assert bg_groupname in dx.bg_data_file

    group = dx.bg_data_file.get_node(bg_groupname)
    group._v_attrs.signature = signature
    dx.bg_data_file.flush()

def calc_bg_cache(dx, fun, time_s=60, tail_min_us=500, F_bg=2, recompute=False):
    """Cached version of `.calc_bg()` method."""

    curr_call_signature = dict(fun_name=fun.__name__, time_s=time_s,
                               tail_min_us=tail_min_us, F_bg=F_bg)
    if _bg_is_cached(dx, curr_call_signature) and not recompute:
        # Background found in cache. Load it.
        pprint(' * Loading BG rates from cache ... ')
        bg_groupname = _get_bg_groupname(dx, time_s=time_s)
        dx._clean_bg_data()
        bg_load_hdf5(dx, bg_groupname)
        pprint(' [DONE]\n')
    else:
        # Background not found in cache. Compute it.
        pprint(' * No cached BG rates, recomputing:\n')
        dx.calc_bg(fun=fun, time_s=time_s, tail_min_us=tail_min_us, F_bg=F_bg)
        if 'bg_data_file' in dx:
            pprint(' * Storing BG  to disk ... ')
            # And now store it to disk
            bg_save_hdf5(dx)
            _bg_add_signature(dx, curr_call_signature)
            pprint(' [DONE]\n')


def test_calc_bg_cache(dx):
    pass

