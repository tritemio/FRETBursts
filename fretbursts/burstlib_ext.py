# -*- coding: utf-8 -*-
#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
The module `burtlib_ext.py` (usually imported as `bext`) contains extensions
to `burstslib.py`. It can be though as a simple plugin system for FRETBursts.

Functions here defined operate on :class:`fretbursts.burstlib.Data()` objects
extending the functionality beyond the core functions and methods defined in
`burstlib`. This modularization allows to implement new functions without
overloading the :class:`fretbursts.burstlib.Data` with an high number
of non-core methods.

The type of functions here implemented are quite diverse. A short summary
follows.

* :func:`bursts_fitter` and :func:`fit_bursts_kde_peak` help to build and
  fit histograms and KDEs for E or S.

* :func:`burst_search_and_gate` performs the AND-gate burst search taking
  intersection of the bursts detected in two photons streams.

* :func:`calc_mdelays_hist` computes the histogram of the m-delays
  distribution of photon intervals.

* :func:`burst_data_period_mean` computes a mean of any "burst data"
  for each background period.

* :func:`join_data` joins different measuremets to create a single
  "virtual" measurement from a series of measurements.

Finally a few functions deal with burst timestamps:

* :func:`get_burst_photons` returns a list of timestamps for each burst.
* :func:`ph_burst_stats` compute any statistics (for example mean or median)
  on the timestamps of each burst.
* :func:`asymmetry` returns a burst "asymmetry index" based on the difference
  between Donor and Acceptor timestamps.

"""
from __future__ import division

import numpy as np
from scipy.stats import erlang
from scipy.optimize import leastsq
import pandas as pd
import tables

from ph_sel import Ph_sel
import burstsearch.burstsearchlib as bslib
import background as bg
from utils.misc import pprint

import burstlib
import fret_fit
import mfit


def calc_mean_lifetime(dx, t1=0, t2=np.inf, ph_sel=Ph_sel('all')):
    """Compute the mean lifetime in each burst.

    Arguments:
        t1, t2 (floats): min and max value (in TCSPC bin units) for the
            nanotime to be included in the mean
        ph_sel (Ph_sel object): object defining the photon selection.
            See :class:`fretbursts.ph_sel.Ph_sel` for details.

    Returns:
        List of arrays of per-burst mean lifetime. One array per channel.
    """
    mean_lifetimes = []

    for bursts, nanot, mask in zip(dx.mburst, dx.nanotimes,
                                   dx.iter_ph_masks(ph_sel)):
        selection = (nanot > t1)*(nanot < t2)
        # Select photons in ph_sel AND with nanotime in [t1, t2]
        if mask != slice(None):
            selection *= mask
        mean_lifetimes.append(
            burstlib.burst_ph_stats(nanot, bursts, mask=selection,
                                    func=np.mean))

    return mean_lifetimes


def _store_bg_data(store, base_name, min_ph_delays_us, best_bg, best_th,
                   BG_data, BG_data_e):
    if not base_name.endswith('/'):
        base_name = base_name + '/'
    store_name = store.filename
    group_name = '/' + base_name[:-1]
    store.create_carray(group_name, 'min_ph_delays_us', obj=min_ph_delays_us,
                        createparents=True)
    for ph_sel, values in BG_data.items():
        store.create_carray(group_name, str(ph_sel), obj=values)
    for ph_sel, values in BG_data_e.items():
        store.create_carray(group_name, str(ph_sel) + '_err', obj=values)
    store.close()
    store = pd.HDFStore(store_name)
    store[base_name + 'best_bg'] = best_bg
    store[base_name + 'best_th'] = best_th
    store.close()

def _load_bg_data(store, base_name, ph_streams):
    if not base_name.endswith('/'):
        base_name = base_name + '/'
    store_name = store.filename
    group_name = '/' + base_name[:-1]
    min_ph_delays = store.get_node(group_name, 'min_ph_delays_us')[:]
    BG_data = {}
    for ph_sel in ph_streams:
        BG_data[ph_sel] = store.get_node(group_name, str(ph_sel))[:]
    BG_data_e = {}
    for ph_sel in ph_streams:
        BG_data_e[ph_sel] = store.get_node(group_name, str(ph_sel) + '_err')[:]
    store.close()
    store = pd.HDFStore(store_name)
    best_bg = store[base_name + 'best_bg']
    best_th = store[base_name + 'best_th']
    store.close()
    return best_th, best_bg, BG_data, BG_data_e, min_ph_delays

def calc_bg_brute_cache(dx, min_ph_delay_list=None, return_all=False,
                        error_metrics='KS', force_recompute=False):
    """Compute background for all the ch, ph_sel and periods caching results.

    This function performs a brute-force search of the min ph delay
    threshold. The best threshold is the one the minimizes the error
    function. The best background fit is the rate fitted using the
    best threshold.

    Results are cached to disk and loaded trasparentlty when needed.
    The cache file is an HDF5 file named `dx.fname[:-5] + '_BKG.hdf5'`.

    Arguments:
        min_ph_delay_list (sequence): sequence of values used for the
            brute-force search. Background and error will be computed
            for each value in `min_ph_delay_list`.
        return_all (bool): if True return all the fitted backgrounds and
            error functions. Default False.
        error_metrics (string): Specifies the error metric to use.
            See :func:`fretbursts.background.exp_fit` for more details.
        force_recompute (bool): if True, recompute results even if a cache
            is found.

    Returns:
        Two arrays with best threshold (us) and best background. If
        `return_all = True` also returns the dictionaries containing all the
        fitted backgrounds and errors.
    """
    if min_ph_delay_list is None:
        min_ph_delay_list = np.arange(100, 8500, 100)
    else:
        min_ph_delay_list = np.asfarray(min_ph_delay_list)

    base_name = 'bg_%s_%ds_%s/' % (dx.bg_fun_name, dx.bg_time_s, error_metrics)
    store_fname = dx.fname[:-5] + '_BKG.hdf5'
    comp_filter = tables.Filters(complevel=6, complib='zlib')
    store = tables.open_file(store_fname, mode='a', filters=comp_filter)
    loaded = False
    if '/' + base_name + 'min_ph_delays_us' in store:
        Th = store.get_node('/', base_name + 'min_ph_delays_us').read()
        if np.all(Th == min_ph_delay_list) and not force_recompute:
            print ' - Loading BG from cache'
            res = _load_bg_data(store, base_name, dx.ph_streams)
            loaded = True
    if not loaded:
        print ' - Computing BG'
        res = calc_bg_brute(dx, min_ph_delay_list=min_ph_delay_list,
                            return_all=True, error_metrics=error_metrics)
        best_th, best_bg, BG_data, BG_data_e, min_ph_delays = res
        _store_bg_data(store, base_name, min_ph_delays, best_bg, best_th,
                       BG_data, BG_data_e)

    if return_all:
        return res
    else:
        return res[:-3]

def calc_bg_brute(dx, min_ph_delay_list=None, return_all=False,
                  error_metrics='KS'):
    """Compute background for all the ch, ph_sel and periods.

    This function performs a brute-force search of the min ph delay
    threshold. The best threshold is the one the minimizes the error
    function. The best background fit is the rate fitted using the
    best threshold.

    Parameters
        min_ph_delay_list (sequence): sequence of values used for the
            brute-force search. Background and error will be computed
            for each value in `min_ph_delay_list`.
        return_all (bool): if True return all the fitted backgrounds and
            error functions. Default False.
        error_metrics (string): Specifies the error metric to use.
            See :func:`fretbursts.background.exp_fit` for more details.

    Returns
        Two arrays with best threshold (us) and best background. If
        `return_all = True` also returns the dictionaries containing all the
        fitted backgrounds and errors.
    """
    if min_ph_delay_list is None:
        min_ph_delay_list = np.arange(100, 8500, 100)
    else:
        min_ph_delay_list = np.asfarray(min_ph_delay_list)

    ph_sel_labels = [str(p) for p in dx.ph_streams]

    BG_data, BG_data_e = {}, {}

    index = pd.MultiIndex.from_product([range(dx.nch), ph_sel_labels],
                                       names=['CH', 'ph_sel'])
    best_th = pd.DataFrame(index=index, columns=np.arange(dx.nperiods))
    best_th.columns.name = 'period'
    best_th.sortlevel(inplace=True)
    best_bg = best_th.copy()

    for ph_sel in dx.ph_streams:
        # Compute BG and error for all ch, periods and thresholds
        # Shape: (nch, nperiods, len(thresholds))
        BG_data[ph_sel], BG_data_e[ph_sel] = bg.fit_varying_min_delta_ph(
                dx, min_ph_delay_list, bg_fit_fun=bg.exp_fit, ph_sel=ph_sel,
                error_metrics=error_metrics)

        # Compute the best Th and BG estimate for all ch and periods
        for ich in range(dx.nch):
            for period in range(dx.nperiods):
                b = BG_data_e[ph_sel][ich, period, :]
                i_best_th = b[-np.isnan(b)].argmin()
                best_th.loc[(ich, str(ph_sel)), period] = \
                        min_ph_delay_list[i_best_th]
                best_bg.loc[(ich, str(ph_sel)), period] = \
                        BG_data[ph_sel][ich, period, i_best_th]
    if return_all:
        return best_th, best_bg, BG_data, BG_data_e, min_ph_delay_list
    else:
        return best_th, best_bg


def burst_data(dx, ich=0, include_bg=False, include_ph_index=False):
    """Return a pandas Dataframe (one row per bursts) with all the burst data.
    """
    if dx.ALEX:
        nd, na, naa, bg_d, bg_a, bg_aa, wid = dx.expand(ich=ich, alex_naa=True,
                                                        width=True)
        nt = nd + na + naa + dx.nda[ich]
    else:
        nd, na, bg_d, bg_a, wid = dx.expand(ich=ich, width=True)
        nt = nd + na

    size_raw = burstlib.b_size(dx.mburst[ich])
    t_start = burstlib.b_start(dx.mburst[ich])*dx.clk_p
    t_end = burstlib.b_end(dx.mburst[ich])*dx.clk_p
    i_start = burstlib.b_istart(dx.mburst[ich])
    i_end = burstlib.b_iend(dx.mburst[ich])
    asym = asymmetry(dx, dropnan=False)

    data_dict = dict(size_raw=size_raw, nt=nt, width_ms=wid*1e3,
                     t_start=t_start, t_end=t_end, asymmetry=asym)

    if include_ph_index:
        data_dict.update(i_start=i_start, i_end=i_end)

    if include_bg:
        data_dict.update(bg_d=bg_d, bg_a=bg_a)
        if dx.ALEX:
            data_dict.update(bg_aa=bg_aa)

    burst_fields = dx.burst_fields[:]
    burst_fields.remove('mburst')
    for field in burst_fields:
        if field in dx:
            data_dict[field] = dx[field][ich]

    return pd.DataFrame.from_dict(data_dict)


def fit_bursts_kde_peak(dx, burst_data='E', bandwidth=0.03, weights=None,
                        gamma=1, add_naa=False, x_range=(-0.1, 1.1),
                        x_ax=None, save_fitter=True):
    """Fit burst data (typ. E or S) by finding the KDE max on all the channels.

    Parameters
        dx (Data): `Data` object containing the FRET data
        burst_data (string): name of burst-data attribute (i.e 'E' or 'S').
        bandwidth (float): bandwidth for the Kernel Densisty Estimation
        weights (string or None): kind of burst-size weights.
            See :func:`fretbursts.fret_fit.get_weights`.
        gamma (float): gamma factor passed to `get_weights()`.
        add_naa (bool): if True adds `naa` to the burst size.
        save_fitter (bool): if True save the `MultiFitter` object in the
            `dx` object with name: burst_data + '_fitter'.
        x_range (tuple of floats): min-max range where to search for the peak.
            Used to select a single peak in a multi-peaks distribution.
        x_ax (array or None): x-axis used to evaluate the Kernel Density

    Returns
        An array of max peak positions (one per ch). If the number of
        channels is 1 returns a scalar.
    """
    if x_ax is None:
        x_ax = np.arange(-0.2, 1.2, 0.0002)

    fitter = bursts_fitter(dx, burst_data=burst_data, save_fitter=save_fitter,
                           weights=weights, gamma=gamma, add_naa=add_naa)
    fitter.calc_kde(bandwidth=bandwidth)
    fitter.find_kde_max(x_ax, xmin=x_range[0], xmax=x_range[1])
    KDE_max_mch = fitter.kde_max_pos
    if dx.nch == 1:
        KDE_max_mch = KDE_max_mch[0]
    return KDE_max_mch

def bursts_fitter(dx, burst_data='E', save_fitter=True,
                  weights=None, gamma=1, add_naa=False):
    """Create a mfit.MultiFitter object (for E or S) add it to `dx`.

    A MultiFitter object allows to fit multi-channel data with the same
    model.

    Parameters
        dx (Data): `Data` object containing the FRET data
        save_fitter (bool): if True save the `MultiFitter` object in the
            `dx` object with name: burst_data + '_fitter'.
        burst_data (string): name of burst-data attribute (i.e 'E' or 'S').
        weights (string or None): kind of burst-size weights.
            See :func:`fretbursts.fret_fit.get_weights`.
        gamma (float): gamma factor passed to `get_weights()`.
        add_naa (bool): if True adds `naa` to the burst size.

    Returns
        The `mfit.MultiFitter` object with the specified burst-size weights.
    """
    assert burst_data in dx
    fitter = mfit.MultiFitter(dx[burst_data])
    if weights is not None:
        weight_kwargs = dict(weights=weights, gamma=gamma, nd=dx.nd, na=dx.na)
        if add_naa:
            weight_kwargs.update(naa = dx.naa)
        fitter.set_weights_func(weight_func = fret_fit.get_weights,
                                weight_kwargs = dict(weights=weights,
                                                     gamma=gamma,
                                                     nd=dx.nd, na=dx.na))
    if save_fitter:
        dx.add(**{burst_data + '_fitter': fitter,
                  'burst_weights': (weights, float(gamma), add_naa)})
    return fitter

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

def _get_bg_erlang(d, ich=0, m=10, ph_sel=Ph_sel('all'), period=0):
    """Return a frozen (scipy) erlang distrib. with rate equal to the bg rate.
    """
    bg_rate = d.bg_from(ph_sel=ph_sel)[ich][period]
    #bg_rate_kcps = bg_rate*1e-3
    bg_dist = erlang(a=m, scale=1./bg_rate)
    return bg_dist

def histogram_mdelays(d, ich=0, m=10, period=None, ph_sel=Ph_sel('all'),
                      binwidth=1e-3, dt_max= 10e-3, bins=None, bursts=False,
                      pdf=True):
    """Compute histogram of m-photons delays (or waiting times).

    Arguments
        dx (Data object): contains the burst data to process.
        ich (int): the channel number. Default 0.
        m (int): number of photons used to compute each delay.
        period (int): index of the period to use.
        ph_sel (Ph_sel object): photon selection to use.

    Returns
        Two arrays for the histogram counts and the bins centers. If `bursts`
        is True the counts arrays contains a second row that is the histogram
        for only the photons inside bursts.
    """
    if bins is None:
        bins = np.arange(0, dt_max, binwidth)

    if period is None:
        ph = d.get_ph_times(ich=ich, ph_sel=ph_sel)
    else:
        ph = d.get_ph_times_period(period=period, ich=ich, ph_sel=ph_sel)

    ph_mdelays = np.diff(ph[::m])*d.clk_p
    if bursts:
        if period is not None:
            print "WARNING: the burst-ph histogram is built from all periods"
        ph_in_bursts = d.ph_in_bursts_ich(ich=ich, ph_sel=ph_sel)
        phb_mdelays = np.diff(ph_in_bursts[::m])*d.clk_p

    # Compute the histogram
    hist_kwargs = dict(bins=bins, density=False)
    counts_tot, _ = np.histogram(ph_mdelays, **hist_kwargs)
    bin_center = bins[:-1] + 0.5*(bins[1] - bins[0])
    pdf_tot = counts_tot / (counts_tot.sum()*binwidth)
    histograms_y = pdf_tot if pdf else counts_tot
    if bursts:
        counts_bursts, _ = np.histogram(phb_mdelays, **hist_kwargs)
        pdf_bursts_normalized = counts_bursts / (counts_tot.sum()*binwidth)
        if pdf:
            histograms_y = np.vstack([pdf_tot, pdf_bursts_normalized])
        else:
            histograms_y = np.vstack([counts_tot, counts_bursts])
    return histograms_y, bin_center

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
        i_origin_ch = np.concatenate([i_d*np.ones(d.num_bursts[ich])
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
    dx_a = dx.copy(mute=mute)
    dx_and = dx.copy(mute=mute)

    dx_d.burst_search_t(L=m, m=m, F=F, ph_sel=ph_sel1, mute=mute)
    dx_a.burst_search_t(L=m, m=m, F=F, ph_sel=ph_sel2, mute=mute)

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


##
#  Burst asymmetry
#

def get_burst_photons(d, ich=0, ph_sel=Ph_sel('all')):
    """Return a list of arrays of photon timestamps in each burst.

    Arguments:
        d (Data): Data() object
        ich (int): channel index
        ph_sel (Ph_sel): photon selection. It allows to select timestamps
            from a specific photon selection. Example ph_sel=Ph_sel(Dex='Dem').
            See :class:`fretbursts.ph_sel.Ph_sel` for details.

    Returns:
        A list of arrays of photon timestamps (one array per burst).
    """
    bursts = d.mburst[ich]
    i_start, i_end = burstlib.b_istart(bursts), burstlib.b_iend(bursts)

    ph_times = d.get_ph_times(ich)
    burst_slices = [slice(i1, i2 + 1) for i1, i2 in zip(i_start, i_end)]
    burst_photons = [ph_times[slice_i] for slice_i in burst_slices]

    if ph_sel != Ph_sel('all'):
        ph_times_mask = d.get_ph_mask(ich, ph_sel=ph_sel)
        photon_masks = [ph_times_mask[slice_i] for slice_i in burst_slices]
        burst_photons = [ph[mask] for ph, mask in zip(burst_photons,
                                                      photon_masks)]
    return burst_photons

def ph_burst_stats(d, ich=0, func=np.mean, ph_sel=Ph_sel('all')):
    """Applies function `func` to the timestamps of each burst.

    Arguments:
        d (Data): Data() object
        ich (int): channel index
        func (function): a function that take an array of burst-timestamps
            and return a scalar. Default `numpy.mean`.
        ph_sel (Ph_sel): photon selection. It allows to select timestamps
            from a specific photon selection. Default Ph_sel('all').
            See :class:`fretbursts.ph_sel.Ph_sel` for details.

    Returns:
        An array containing per-burst timestamp statistics.
    """
    burst_photons = get_burst_photons(d, ich, ph_sel=ph_sel)
    stats = [func(times) for times in burst_photons]
    return np.array(stats)

def asymmetry(dx, ich=0, func=np.mean, dropnan=True):
    """Compute an asymmetry index for each burst in channel `ich`.

    It computes each burst the difference func({t_D}) - func({t_A})
    where `func` is a function (default `mean`) that computes some statistics
    on the timestamp and {t_D} and {t_A} are the sets of D or A timestamps
    in a bursts (during D excitation).

    Arguments:
        d (Data): Data() object
        ich (int): channel index
        func (function): the function to be used to extract D and A photon
            statistics in each bursts.

    Returns:
        An arrays of photon timestamps (one array per burst).
    """
    stats_d = ph_burst_stats(dx, ich=ich, func=func, ph_sel=Ph_sel(Dex='Dem'))
    stats_a = ph_burst_stats(dx, ich=ich, func=func, ph_sel=Ph_sel(Dex='Aem'))

    #b_size = d.burst_sizes(ich, add_naa=False)
    #b_width = burstlib.b_width(d.mburst[ich])
    burst_asym = (stats_d - stats_a)*dx.clk_p*1e3
    if dropnan:
        burst_asym = burst_asym[-np.isnan(burst_asym)]
    return burst_asym
