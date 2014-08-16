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
import select_bursts
import burstlib as bl
import burstsearch.burstsearchlib as bslib
from utils.misc import pprint

import fret_fit
from fit.weighted_kde import gaussian_kde_w


def fit_E_kde_peak(dx, E_range=(-0.1, 1.1), bandwidth=0.03, E_ax=None,
                   weights='size', return_pdf=False):
    """Fit E by finding the KDE maximum on all the channels.

    Parameters
        dx (Data): `Data` object containing the FRET data
        E_range (tuple of floats): min-max range where to search for the peak.
            Used to select a single peak in a multi-peaks distribution.
        bandwidth (float): bandwidth for the Kernel Densisty Estimation
        E_ax (array or None): FRET efficiency axis (i.e. x-axis) used to
            evaluate the Kernel Density
        weights (string or None): kind of burst weights.
            See :func:`fretbursts.fret_fit.get_weights`.
        return_pdf (bool): if True returns also the (x, y) values of the
            PDF computed by the KDE. Default False.

    Returns
        An array of E values (one per ch) and, if return_pdf is True,
        the array of E values (size M) and the array of PDF (size nch x M).
    """
    if E_ax is None:
        E_ax = np.arange(-0.2, 1.2, 0.0002)

    E_fit_mch = np.zeros(dx.nch)
    if return_pdf:
        E_pdf_mch = np.zeros((dx.nch, E_ax.size))

    for ich in range(dx.nch):

        res = fit_E_kde_peak_single_ch(dx, ich=ich,
                              E_range=E_range, bandwidth=bandwidth, E_ax=E_ax,
                              weights=weights, return_pdf=return_pdf)

        if return_pdf:
            E_fit_mch[ich] = res[0]
            E_pdf_mch[ich] = res[2]
        else:
             E_fit_mch[ich] = res

    if return_pdf:
        return E_fit_mch, E_ax, E_pdf_mch
    else:
        return E_fit_mch

def fit_E_kde_peak_single_ch(dx, ich=0, E_range=(-0.1, 1.1), bandwidth=0.03,
                   E_ax=None, weights='size', return_pdf=False):
    """Fit E by finding the KDE maximum on channel `ich`.

    Parameters
        dx (Data): `Data` object containing the FRET data
        ich (int): channel number for multi-spot data. Default 0.
        E_range (tuple of floats): min-max range where to search for the peak.
            Used to select a single peak in a multi-peaks distribution.
        bandwidth (float): bandwidth for the Kernel Densisty Estimation
        E_ax (array or None): FRET efficiency axis (i.e. x-axis) used to
            evaluate the Kernel Density
        weights (string or None): kind of burst weights.
            See :func:`fretbursts.fret_fit.get_weights`.
        return_pdf (bool): if True returns also the (x, y) values of the
            PDF computed by the KDE. Default False.

    Returns
        Fitted E value and (optionally) the KDE of the PDF.
    """
    E_ax, E_pdf = compute_E_kde(dx, ich=ich, bandwidth=bandwidth,
                                E_ax=E_ax, weights=weights)
    E_fit = find_max(E_ax, E_pdf, xmin=E_range[0], xmax=E_range[1])

    if return_pdf:
        return E_fit, E_ax, E_pdf
    else:
        return E_fit

def find_max(x, y, xmin=None, xmax=None):
    """Find peak position of a curve (x, y) between `xmin` and `xmax`.
    """
    if xmin is None:
        xmin = x.min()
    if xmax is None:
        xmax = x.max()

    mask = np.where((x >= xmin)*(x <= xmax))
    return x[mask][y[mask].argmax()]

def compute_E_kde(dx, ich=0, bandwidth=0.03, E_ax=None, weights='size',
                  E_range=None):
    """Compute the KDE for E values in `dx`, channel `ich`.

    Parameters
        dx (Data): `Data` object containing the FRET data
        ich (int): channel number for multi-spot data. Default 0.
        bandwidth (float): bandwidth for the Kernel Densisty Estimation
        E_ax (array or None): FRET efficiency axis (i.e. x-axis) used to
            evaluate the Kernel Density
        weights (string or None): kind of burst weights.
            See :func:`fretbursts.fret_fit.get_weights`.
        E_range (None or tuple): if not None, the KDE is computed on the
            bursts between min/max E values in `E_range`.

    Returns
        Arrays of E (x-axis) and KDE PDF (y-axis)
    """

    if E_ax is None:
        E_ax = np.arange(-0.2, 1.2, 0.0002)

    if E_range is not None:
        dx = bl.Sel(dx, select_bursts.E, E1=E_range[0], E2=E_range[1])

    w = fret_fit.get_weights(dx.nd[ich], dx.na[ich], weights=weights)
    kde = gaussian_kde_w(dx.E[ich], bw_method=bandwidth, weights=w)
    E_pdf = kde.evaluate(E_ax)

    return E_ax, E_pdf


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
    dx_and.add(m=m, L=m, F=F, P=None, ph_sel='AND-gate')
    dx_and.add(bg_corrected=False, leakage_corrected=False,
               dir_ex_corrected=False, dithering=False)

    pprint(" - Counting D and A ph and calculating FRET ... \n", mute)
    dx_and.calc_fret(count_ph=True, corrections=True, mute=mute)
    pprint("   [DONE Counting D/A]\n", mute)

    return dx_and

