# -*- coding: utf-8 -*-
#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
The module `burtlib_ext.py` contains extensions to `burstslib.py`. It can be
though as a simple plugin system for FRETBursts.

Functions here defined operate on :class:`fretbursts.burstlib.Data()` objects
extending the functionality beyond the core functions and methods defined in
`burstlib`.
"""
from __future__ import division

import numpy as np
from scipy.stats import erlang
from scipy.optimize import leastsq

from ph_sel import Ph_sel
import burstsearch.burstsearchlib as bslib
from utils.misc import pprint


def _get_bg_distrib_erlang(d, ich=0, m=10, ph_sel=Ph_sel('all'), bp=(0, -1)):
    """Return a frozen (scipy) erlang distrib. with rate equal to the bg rate.
    """
    assert ph_sel in [Ph_sel('all'), Ph_sel(Dex='Dem'), Ph_sel(Dex='Aem')]
    if np.size(bp) == 1: bp = (bp, bp)
    periods = slice(d.Lim[ich][bp[0]][0], d.Lim[ich][bp[1]][1] + 1)
    # Compute the BG distribution
    if ph_sel == Ph_sel('all'):
        bg_ph = d.bg_dd[ich] + d.bg_ad[ich]
    elif ph_sel == Ph_sel(Dex='Dem'):
        bg_ph = d.bg_dd[ich]
    elif ph_sel == Ph_sel(Dex='Aem'):
        bg_ph = d.bg_ad[ich]

    rate_ch_kcps = bg_ph[periods].mean()/1e3   # bg rate in kcps
    bg_dist = erlang(a=m, scale=1./rate_ch_kcps)
    return bg_dist

def calc_mdelays_hist(d, ich=0, m=10, bp=(0, -1), bins_s=(0, 10, 0.02),
                      ph_sel=Ph_sel('all'), bursts=False, bg_fit=True,
                      bg_F=0.8):
    """Compute histogram of m-photons delays (or waiting times).

    Arguments:
        dx (Data object): contains the burst data to process.
        ich (int): the channel number. Default 0.
        m (int): number of photons used to compute each delay.
        bp (int or 2-element tuple): index of the period to use. If tuple,
            the period range between bp[0] and bp[1] (included) is used.
        bins_s (3-element tuple): start, stop and step for the bins
        ph_sel (Ph_sel object): photon selection to use.
    Returns:
        Tuple of values:

            * bin_x (array): array of bins centers
            * histograms_y (array): arrays of histograms, contains 1 or 2
              histograms (when `bursts` is False or True)
            * bg_dist (random distribution): erlang distribution with same
              rate as background (kcps)
            * a, rate_kcps (floats, optional): amplitude and rate for an
              Erlang distribution fitted to the histogram for
              bin_x > bg_mean*bg_F. Returned only if `bg_fit` is True.
    """
    assert ph_sel in [Ph_sel('all'), Ph_sel(Dex='Dem'), Ph_sel(Dex='Aem')]
    if np.size(bp) == 1: bp = (bp, bp)
    periods = slice(d.Lim[ich][bp[0]][0], d.Lim[ich][bp[1]][1] + 1)
    bins = np.arange(*bins_s)

    if ph_sel == Ph_sel('all'):
        ph = d.ph_times_m[ich][periods]
        if bursts:
            phb = ph[d.ph_in_burst[ich][periods]]
    elif ph_sel == Ph_sel(Dex='Dem'):
        donor_ph_period = -d.A_em[ich][periods]
        ph = d.ph_times_m[ich][periods][donor_ph_period]
        if bursts:
            phb = ph[d.ph_in_burst[ich][periods][donor_ph_period]]
    elif ph_sel == Ph_sel(Dex='Aem'):
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
    bg_dist = _get_bg_distrib_erlang(d, ich=ich, m=m, bp=bp, ph_sel=ph_sel)
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
        dx (Data object): contains the burst data to process
        burst_data (list of arrays): one array per channel, each array
            has one element of "burst data" per burst.

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
    """Joins burst data of different measurements in a single `Data` object.

    This function requires that all the passed data objects use the same
    period (bg_time_s). For each measurement, the time of burst start is
    offset by the duration of the previous measurement + an additional `gap`.

    The index of the first/last photon in the burst (returned by `b_istart()`
    and `b_iend()`) are kept unmodified and refer to the original timestamp
    array. The timestamp arrays are not copied: the new `Data` object will
    not contain any timestamp arrays (ph_times_m). This may cause error when
    calling functions that require the timestamps data.

    The background arrays (bg, bg_dd, etc...) are concatenated. The burst
    attribute `bp` is updated to refer to these new concatenated arrays.
    The attributes `Lim` and `Ph_p` are concatenated and left unchanged.
    Therefore different sections will refer to different original timestamp
    arrays.

    Burst widths and sizes are kept unchanged.

    A new attribute (`i_origin`), containing for each burst the index of the
    original data object in the list, is saved in the returned object.

    Returns:
        A `Data` object containing bursts from the all the objects in `d_list`.
    """

    from fretbursts.burstlib import Data, itstart, itend

    nch = d_list[0].nch
    bg_time_s = d_list[0].bg_time_s
    for d in d_list:
        assert d.nch == nch
        assert d.bg_time_s == bg_time_s

    new_d = Data(**d_list[0])
    new_d.delete('ph_times_m')

    # Set the bursts fields by concatenation along axis = 0
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

    # Set the background fields by concatenation along axis = 0
    new_nperiods = np.sum([d.nperiods for d in d_list])
    for name in Data.bg_fields:
        if name in new_d:
            new_d.add(**{ name: [] })
            for ich in xrange(nch):
                value = np.concatenate([d[name][ich] for d in d_list])
                new_d[name].append(value)
                assert new_d[name][ich].shape[0] == new_nperiods

    # Set the i_origin burst attribute
    new_d.add(i_origin = [])
    for ich in xrange(nch):
        i_origin_ch = np.concatenate([i_d*np.ones(d.num_bursts()[ich])
                        for i_d, d in enumerate(d_list)])
        new_d.i_origin.append(i_origin_ch)

    # Update the `bp` attribute to refer to the background period in
    # the new concatenated background arrays.
    sum_nperiods = np.cumsum([d.nperiods for d in d_list])
    for i_d, d in zip(xrange(1, len(d_list)), d_list[1:]):
        for ich in xrange(nch):
            # Burst "slice" in new_d coming from current d
            b_mask = new_d.i_origin[ich] == i_d
            # Add the nperiods of all the previous measurements
            new_d.bp[ich][b_mask] = new_d.bp[ich][b_mask] + sum_nperiods[i_d-1]

    # Modify the new mburst so the time of burst start/end is monotonic
    offset_clk = 0
    for i_orig, d_orig in enumerate(d_list):
        for ich in xrange(nch):
            if np.size(new_d.mburst[ich]) == 0: continue
            mask = new_d.i_origin[ich] == i_orig
            new_d.mburst[ich][mask, itstart] += offset_clk
            new_d.mburst[ich][mask, itend] += offset_clk
        offset_clk += (d_orig.time_max() + gap)/d_orig.clk_p

    return new_d

def burst_search_and_gate(dx, F=6, m=10, ph_sel1=Ph_sel(Dex='DAem'),
                          ph_sel2=Ph_sel(Aex='Aem'), mute=False):
    """Return a Data object containing bursts obtained by and-gate burst-search.

    The and-gate burst search is a composition of 2 burst searches performed
    on different photon selections. The bursts in the and-gate burst search
    are the overlapping bursts in the 2 initial burst searches, and their
    duration is the intersection of the two overlapping bursts.

    By default the 2 photon selections are D+A photons during D excitation
    (`Ph_sel(Dex='DAem')`) and A photons during A excitation
    (`Ph_sel(Aex='Aex')`).

    Arguments:
        dx (Data object): contains the data on which to perform the burst
            search. Background estimation must be performed before the search.
        F (float): Burst search parameter F.
        m (int): Burst search parameter m.
        ph_sel1 (Ph_sel object): photon selections used for bursts search 1.
        ph_sel2 (Ph_sel object): photon selections used for bursts search 2.
        mute (bool): if True nothing is printed. Default: False.

    Return:
        A new `Data` object containing bursts from the and-gate search.

    See also :meth:`fretbursts.burstlib.Data.burst_search_t`.
    """
    dx_d = dx
    dx_a = dx.copy()
    dx_and = dx.copy()

    dx_d.burst_search_t(L=m, m=m, F=F, ph_sel=ph_sel1)
    dx_a.burst_search_t(L=m, m=m, F=F, ph_sel=ph_sel2)

    mburst_and = []
    for mburst_d, mburst_a in zip(dx_d.mburst, dx_a.mburst):
        mburst_and.append(bslib.burst_and(mburst_d, mburst_a))

    dx_and.add(mburst=mburst_and)

    pprint(" - Calculating burst periods ...", mute)
    dx_and._calc_burst_period()                       # writes bp
    pprint("[DONE]\n", mute)

    # Note: dx_and.bg_bs will not be meaningful
    dx_and.add(m=m, L=m, F=F, P=None, ph_sel=Ph_sel(Dex='DAem'))
    dx_and.add(bg_corrected=False, leakage_corrected=False,
               dir_ex_corrected=False, dithering=False)

    pprint(" - Counting D and A ph and calculating FRET ... \n", mute)
    dx_and.calc_fret(count_ph=True, corrections=True, mute=mute)
    pprint("   [DONE Counting D/A]\n", mute)

    return dx_and

