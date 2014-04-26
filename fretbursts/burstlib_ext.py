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

import numpy as np
from scipy.stats import erlang
from scipy.optimize import leastsq


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
            each array has one element per burst.

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

def join_data(d_list, gap=1):
    """Joins burst data of different measurements in a single Data() object.

    This function requires that all the passed data objects use the same
    period (bg_time_s). For each measurement, the time of burst start is
    offset by the duration of the previous measurement + an additional `gap`.

    The index of the first/last photon in the burst (returned by `b_istart()`
    and `b_iend()`) are keept unmodified and refer to the original timestamp
    array. The timestamp arrays are not copied: the new `Data` object will
    have empty lists instead of timestamp arrays. This may cause error if
    calling functions that require the timestamps arrays.

    Burst widths and sizes are kept unchanged.

    A new attribute (`i_origin`), containing, for each burst, the index of
    the original data object in the list, is saved in the returned object.

    Returns:
        A `Data` object containing bursts from the all the objects in `d_list`.
    """
    from fretbursts.burstlib import Data, itstart

    nch = d_list[0].nch
    bg_time_s = d_list[0].bg_time_s
    for d in d_list:
        assert d.nch == nch
        assert d.bg_time_s == bg_time_s

    new_d = Data(**d_list[0])
    new_d.add(ph_times_m = [[]]*nch)

    # Set all the bursts fields by simple concatenation along axis = 0
    for name in Data.burst_fields:
        if name in new_d:
            new_d.add(**{ name: [np.array([])]*nch })
            for ich in xrange(nch):
                new_size = np.sum([d[name][ich].shape[0] for d in d_list])
                if new_size == 0:
                    continue  # -> No bursts in this ch
                value = np.concatenate([d[name][ich] for d in d_list])
                new_d[name][ich] = value
                assert new_d[name][ich].shape[0] == new_size

    # Set the i_origin attribute
    new_d.add(i_origin = [])
    for ich in xrange(nch):
        i_origin_ch = np.concatenate([i_d*np.ones(d.num_bu()[ich])
                        for i_d, d in enumerate(d_list)])
        new_d.i_origin.append(i_origin_ch)

    # Modify the new mburst so the time of burst start is monotonic
    offset_clk = 0
    for i_orig, d_orig in enumerate(d_list):
        for ich in xrange(nch):
            mask = new_d.i_origin[ich] == i_orig
            new_d.mburst[ich][mask, itstart] += offset_clk
        offset_clk += (d_orig.time_max() + gap)/d_orig.clk_p

    return new_d
