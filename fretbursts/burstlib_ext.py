# -*- coding: utf-8 -*-
#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
This module contains extensions to `burstslib.py`.

Functions here defined operate on :class:`Data()` objects extending the
functionality beyond the core functions and methods defined in
:module:`burstlib`.
"""
from __future__ import division
import os

import numpy as np
from scipy.stats import erlang
from scipy.optimize import leastsq
from fretbursts.utils.misc import pprint


def get_bg_distrib_erlang(d, ich=0, m=10, ph_sel='DA', bp=(0, -1)):
    """Return a frozen erlang distrib. with rate equal to the bg rate.
    """
    assert ph_sel in ['DA', 'D', 'A']
    if np.size(bp) == 1: bp = (bp, bp)
    periods = slice(d.Lim[ich][bp[0]][0], d.Lim[ich][bp[1]][1] + 1)
    # Compute the BG distribution
    if ph_sel == 'DA':
        bg_ph = d.bg_dd[ich] + d.bg_ad[ich]
    elif ph_sel == 'D':
        bg_ph = d.bg_dd[ich]
    elif ph_sel == 'A':
        bg_ph = d.bg_ad[ich]

    rate_ch_kcps = bg_ph[periods].mean()/1e3   # bg rate in kcps
    bg_dist = erlang(a=m, scale=1./rate_ch_kcps)
    return bg_dist

def calc_mdelays_hist(d, ich=0, m=10, bp=(0, -1), bins_s=(0, 10, 0.02),
                      ph_sel='DA', bursts=False, bg_fit=True, bg_F=0.8):
    """Compute histogram of m-photons delays (or waiting times).

    Arguments:
        m (int): number of photons used to compute each delay.

        bp (int or 2-element tuple): index of the period to use. If tuple,
            the period range between bp[0] and bp[1] (included) is used.
        bins_s (3-element tuple): start, stop and step for the bins
        ph_sel (string): Photons to use. 'DA' donor + acceptor, 'D' donor,
            'A' acceptor.
    Returns:
        bin_x: array of bins centers
        histograms_y: arrays of histograms, contains 1 or 2 histograms
            (when `bursts` is False or True)
        bg_dist: erlang distribution with same rate as background (kcps)
        a, rate_kcps (floats, optional): amplitude and rate for an Erlang
            distribution fitted to the histogram for bin_x > bg_mean*bg_F.
            Returned only if `bg_fit` is True.
    """
    assert ph_sel in ['DA', 'D', 'A']
    if np.size(bp) == 1: bp = (bp, bp)
    periods = slice(d.Lim[ich][bp[0]][0], d.Lim[ich][bp[1]][1] + 1)
    bins = np.arange(*bins_s)

    if ph_sel == 'DA':
        ph = d.ph_times_m[ich][periods]
        if bursts:
            phb = ph[d.ph_in_burst[ich][periods]]
    elif ph_sel == 'D':
        donor_ph_period = -d.A_em[ich][periods]
        ph = d.ph_times_m[ich][periods][donor_ph_period]
        if bursts:
            phb = ph[d.ph_in_burst[ich][periods][donor_ph_period]]
    elif ph_sel == 'A':
        accept_ph_period = d.A_em[ich][periods]
        ph = d.ph_times_m[ich][periods][accept_ph_period]
        if bursts:
            phb = ph[d.ph_in_burst[ich][periods][accept_ph_period]]

    ph_mdelays = np.diff(ph[::m])*d.clk_p*1e3        # millisec
    if bursts:
        phb_mdelays = np.diff(phb[::m])*d.clk_p*1e3  # millisec
        phb_mdelays = phb_mdelays[phb_mdelays < 5]

    # Compute the PDF through histograming
    hist_kwargs = dict(bins=bins, normed=True)
    mdelays_hist_y, _ = np.histogram(ph_mdelays, **hist_kwargs)
    bin_x = bins[:-1] + 0.5*(bins[1] - bins[0])
    if bursts:
        mdelays_b_hist_y, _ = np.histogram(phb_mdelays, **hist_kwargs)
        mdelays_b_hist_y *= phb_mdelays.size/ph_mdelays.size
    if bursts:
        histograms_y = np.vstack([mdelays_hist_y, mdelays_b_hist_y])
    else:
        histograms_y = mdelays_hist_y

    results = [bin_x, histograms_y]

    # Compute the BG distribution
    bg_dist = get_bg_distrib_erlang(d, ich=ich, m=m, bp=bp, ph_sel=ph_sel)
    bg_mean = bg_dist.mean()
    results.append(bg_dist)

    if bg_fit:
        ## Fitting the BG portion of the PDF to an Erlang
        _x = bin_x[bin_x > bg_mean*bg_F]
        _y = mdelays_hist_y[bin_x > bg_mean*bg_F]
        fit_func = lambda x, a, rate_kcps: a*erlang.pdf(x, a=m,
                                                        scale=1./rate_kcps)
        err_func = lambda p, x, y: fit_func(x, p[0], p[1]) - y
        p, flag = leastsq(err_func, x0=[0.9, 3.], args=(_x,_y))
        print p, flag
        a, rate_kcps = p
        results.extend([a, rate_kcps])

    return results

def burst_data_period_mean(dx, burst_data):
    """Compute mean `bursts_data` in each period.

    Arguments:
        burst_data (list of arrays): one array per channel,
            each array hase one element per burst.

    Returns:
        2D of arrays with shape (nch, nperiods).

    Example:
        burst_period_mean(dx, dx.nt)
    """
    mean_burst_data = np.zeros((dx.nch, dx.nperiods))
    for ich, (b_data_ch, period) in enumerate(zip(burst_data, dx.bp)):
        for iperiod in xrange(dx.nperiods):
            mean_burst_data[ich, iperiod] = b_data_ch[period == iperiod].mean()
    return mean_burst_data


##
# Store and load background to and from an HDF5 file
#

def _get_bg_arrays_info(dx):
    bg_arrays_info = dict(
        bg = 'BG rate in each CH vs time (Total)',
        bg_dd = 'BG rate in each CH vs time (D_em during D_ex)',
        bg_ad = 'BG rate in each CH vs time (A_em during D_ex)',
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
        bg_th_us_all = 'Waiting time threshold for BG fit of all timestamps',
        bg_th_us_DD = 'Waiting time threshold for BG fit of D_em D_ex timestamps',
        bg_th_us_AD = 'Waiting time threshold for BG fit of A_em D_ex timestamps',
        bg_th_us_AA = 'Waiting time threshold for BG fit of A_em A_ex timestamps',
    )

    if dx.bg_auto_th:
        bg_arrays_info.update(**bg_arrays_info_auto)
    else:
        bg_arrays_info.update(**bg_arrays_info_noauto)
    return bg_arrays_info


def bg_save_hdf5(dx):
    """Save background to HDF5 file (experimental, manta-only)."""
    if 'bg' not in dx:
        print "Compute background first!"
        return
    if 'data_file' not in dx:
        print "No open HDF5 file found."
        return

    bg_arrays_info = _get_bg_arrays_info(dx)
    bg_attr_names = ['bg_fun', 'bg_fun_name']
    if '/background' not in dx.data_file:
        dx.data_file.create_group('/', 'background',
                                  title='Background estimation data')
    ## Save the bg data
    group_name = _get_bg_groupname(dx)
    if group_name in dx.data_file:
        dx.data_file.remove_node(group_name, recursive=True)

    bg_group = dx.data_file.create_group(os.path.dirname(group_name),
                                         os.path.basename(group_name),
                                         createparents=True)
    # Save arrays and scalars
    for name, info in bg_arrays_info.items():
        arr = np.array(dx[name])
        print name
        dx.data_file.create_array(bg_group, name, obj=arr, title=info)

    # Save the attributes
    for attr in bg_attr_names:
        bg_group._v_attrs[attr] = dx[attr]

    dx.data_file.flush()

def bg_load_hdf5(dx, group_name):
    """Load background from a HDF5 file (experimental, manta-only)."""
    if group_name not in dx.data_file:
        print 'Group not found in the HDF5 file.'
        return

    ## Load the bg data
    bg_arrays = dict()
    bg_attrs = dict()

    bg_group = dx.data_file.get_node(group_name)

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
    in_map = ['', '_dd', '_ad', '_aa']
    out_map = ['_m', '_dd', '_ad', '_aa']
    new_attrs = {}
    for in_s, out_s in zip(in_map, out_map):
        assert 'bg' + in_s in dx
        pprint('bg' + in_s + ', ')
        new_attrs['rate' + out_s] = [bg.mean() for bg in dx['bg' + in_s]]
    dx.add(**new_attrs)
    pprint('\n')

def _get_bg_groupname(dx, time_s=None):
    if time_s is None and 'bg' not in dx:
        print 'You need to compute the background or provide time_s.'
        return
    time_slice = time_s if time_s is not None else dx.bg_time_s
    return '/background/time_%ds' % time_slice

def _bg_is_cached(dx, signature):
    if 'data_file' in dx:
        bg_groupname = _get_bg_groupname(dx, time_s=signature['time_s'])
        if bg_groupname in dx.data_file:
            # At least we have a group with the right time_s
            # Check whether its signature matches the current one
            group = dx.data_file.get_node(bg_groupname)
            if signature == group._v_attrs.signature:
                return True
    return False

def _bg_add_signature(dx, signature):
    assert 'data_file' in dx
    bg_groupname = _get_bg_groupname(dx, time_s=signature['time_s'])
    assert bg_groupname in dx.data_file

    group = dx.data_file.get_node(bg_groupname)
    group._v_attrs.signature = signature
    dx.data_file.flush()

def calc_bg_cache(dx, fun, time_s=60, tail_min_us=500, F_bg=2, recompute=False):
    """Cached version of `.calc_bg()` method."""
    if not 'data_file' in dx:
        print 'No open HDF5 file found.'
        return

    curr_call_signature = dict(fun_name=fun.__name__,
                             time_s=time_s, tail_min_us=tail_min_us, F_bg=F_bg)
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

        pprint(' * Storing BG  to disk ... ')
        # And now store it to disk
        bg_save_hdf5(dx)
        _bg_add_signature(dx, curr_call_signature)
        pprint(' [DONE]\n')


def test_calc_bg_cache(dx):
    pass

