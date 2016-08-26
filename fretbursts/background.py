#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2013-2016 The Regents of the University of California,
#               Antonino Ingargiola <tritemio@gmail.com>
#
"""
Routines to compute the background from an array of timestamps. This module
is normally imported as `bg` when fretbursts is imported.

The important functions are :func:`exp_fit` and :func:`exp_cdf_fit` that
provide two (fast) algorithms to estimate the background without binning.
These functions are not usually called directly but passed to
:meth:`Data.calc_bg` to compute the background of a measurement.

See also :func:`exp_hist_fit` for background estimation using an histogram fit.
"""

from __future__ import absolute_import
from builtins import zip, range

import numpy as np
from .ph_sel import Ph_sel
from .utils.misc import pprint
from .fit import exp_fitting
from .fit.gaussian_fitting import gaussian_fit_hist


def raw_fit(ph, clk_p=12.5e-9, residuals=False, tail_min_us=None):
    """Compute the "raw" rate: (number of ph / duration).

    Note:
        `tail_min_us` argument is there only to have the same signature
        of the other background fit functions but is ignored.
    """
    Lambda = ph.size/((ph[-1]-ph[0])*clk_p)
    if residuals:
        resid = exp_fitting.get_residuals(np.diff(ph), 1./(Lambda*clk_p))
        return Lambda, np.abs(resid).max()*100
    else:
        return Lambda, 0


def _compute_error(residuals, x_residuals, error_metrics):
    assert error_metrics in ['KS', 'CM', None]
    error = None
    if error_metrics == 'KS':
        error = np.abs(residuals).max()*100
    elif error_metrics == 'CM':
        error = np.trapz(residuals**2, x=x_residuals)
    return error

def _exp_fit_generic(ph, fit_fun, tail_min_us=None, tail_min_p=0.1,
                     clk_p=12.5e-9, error_metrics=None):
    """Computes BG rates on timestamp delays above a min. value.

    Compute a background rate, selecting waiting-times (delays) larger than a
    minimum threshold.

    You need to pass the specific fitting function as `fit_fun`.

    Arguments:
        ph (array): timestamps from which estimete the background.
        fit_fun (function): function to use for background fitting.
        tail_min_us (float): minimum photon separation threshold in us for
            photons to be considered as background.
        tail_min_p (int): min threshold in percentace. ONly used when
            `tail_min_us` is None. Deprecated.
        clk_p (float): unit of timestamps in `ph` (seconds).
        error_metrics (string or None): Valid values are 'KS' or 'CM'.
            'KS' (Kolmogorov-Smirnov statistics) computes the error as the
            max of deviation of the empirical CDF from the fitted CDF.
            'CM' (Crames-von Mises) uses the L^2 distance.
            If None, no error metric is computed (returns None).

    Returns:
        2-Tuple: Estimated background rate in cps, and a "quality of fit"
        index (the lower the better) according to the chosen metric.
        If error_metrics==None, the returned "quality of fit" is None.
    """
    dph = np.diff(ph)
    if tail_min_us is None:
        tail_min = dph.max()*tail_min_p
    else:
        tail_min = tail_min_us*1e-6/clk_p

    res = fit_fun(dph, s_min=tail_min, calc_residuals=error_metrics is not None)
    Lambda, residuals, x_residuals, s_size = res

    error = _compute_error(residuals, x_residuals, error_metrics)
    Lambda /= clk_p
    return Lambda, error


def exp_fit(ph, tail_min_us=None, clk_p=12.5e-9, error_metrics=None):
    """Return a background rate using the MLE of mean waiting-times.

    Compute the background rate, selecting waiting-times (delays) larger
    than a minimum threshold.

    This function performs a Maximum Likelihood (ML) fit. For
    exponentially-distributed waiting-times this is the empirical mean.

    Arguments:
        ph (array): timestamps array from which to extract the background
        tail_min_us (float): minimum waiting-time in micro-secs
        clk_p (float): clock period for timestamps in `ph`
        error_metrics (string or None): Valid values are 'KS' or 'CM'.
            'KS' (Kolmogorov-Smirnov statistics) computes the error as the
            max of deviation of the empirical CDF from the fitted CDF.
            'CM' (Crames-von Mises) uses the L^2 distance.
            If None, no error metric is computed (returns None).

    Returns:
        2-Tuple: Estimated background rate in cps, and a "quality of fit"
        index (the lower the better) according to the chosen metric.
        If error_metrics==None, the returned "quality of fit" is None.

    See also:
        :func:`exp_cdf_fit`, :func:`exp_hist_fit`
    """
    return _exp_fit_generic(ph, fit_fun=exp_fitting.expon_fit,
                            tail_min_us=tail_min_us, clk_p=clk_p,
                            error_metrics=error_metrics)

def exp_cdf_fit(ph, tail_min_us=None, clk_p=12.5e-9, error_metrics=None):
    """Return a background rate fitting the empirical CDF of waiting-times.

    Compute the background rate, selecting waiting-times (delays) larger
    than a minimum threshold.

    This function performs a least square fit of an exponential Cumulative
    Distribution Function (CDF) to the empirical CDF of waiting-times.

    Arguments:
        ph (array): timestamps array from which to extract the background
        tail_min_us (float): minimum waiting-time in micro-secs
        clk_p (float): clock period for timestamps in `ph`
        error_metrics (string or None): Valid values are 'KS' or 'CM'.
            'KS' (Kolmogorov-Smirnov statistics) computes the error as the
            max of deviation of the empirical CDF from the fitted CDF.
            'CM' (Crames-von Mises) uses the L^2 distance.
            If None, no error metric is computed (returns None).

    Returns:
        2-Tuple: Estimated background rate in cps, and a "quality of fit"
        index (the lower the better) according to the chosen metric.
        If error_metrics==None, the returned "quality of fit" is None.

    See also:
        :func:`exp_fit`, :func:`exp_hist_fit`
    """
    return _exp_fit_generic(ph, fit_fun=exp_fitting.expon_fit_cdf,
                            tail_min_us=tail_min_us, clk_p=clk_p,
                            error_metrics=error_metrics)


def exp_hist_fit(ph, tail_min_us, binw=50e-6, clk_p=12.5e-9,
                 weights='hist_counts', error_metrics=None):
    """Compute background rate with WLS histogram fit of waiting-times.

    Compute the background rate, selecting waiting-times (delays) larger
    than a minimum threshold.

    This function performs a Weighed Least Squares (WLS) fit of the
    histogram of waiting times to an exponential decay.

    Arguments:
        ph (array): timestamps array from which to extract the background
        tail_min_us (float): minimum waiting-time in micro-secs
        binw (float): bin width for waiting times, in seconds.
        clk_p (float): clock period for timestamps in `ph`
        weights (None or string): if None no weights is applied.
            if is 'hist_counts', each bin has a weight equal to its counts
            if is 'inv_hist_counts', the weight is the inverse of the counts.
        error_metrics (string or None): Valid values are 'KS' or 'CM'.
            'KS' (Kolmogorov-Smirnov statistics) computes the error as the
            max of deviation of the empirical CDF from the fitted CDF.
            'CM' (Crames-von Mises) uses the L^2 distance.
            If None, no error metric is computed (returns None).

    Returns:
        2-Tuple: Estimated background rate in cps, and a "quality of fit"
        index (the lower the better) according to the chosen metric.
        If error_metrics==None, the returned "quality of fit" is None.

    See also:
        :func:`exp_fit`, :func:`exp_cdf_fit`
    """
    assert np.size(ph) > 0
    dph = np.diff(ph)
    tail_min = tail_min_us*1e-6/clk_p
    binw_clk = binw/clk_p
    bins = np.arange(0, dph.max() - tail_min + 1, binw_clk)

    res = exp_fitting.expon_fit_hist(dph, bins=bins, s_min=tail_min,
                                     weights=weights,
                                     calc_residuals=error_metrics is not None)

    Lambda, residuals, x_residuals, s_size = res
    error = _compute_error(residuals, x_residuals, error_metrics)
    Lambda /= clk_p
    return Lambda, error

##
# Fit background as function of th
#
def fit_varying_min_delta_ph(d, min_delta_ph_list, bg_fit_fun=exp_fit,
                             ph_sel=Ph_sel('all'), **kwargs):
    """
    Fit the background as a function of the min photon interval threshold.

    The background is fitted for all the channels, all the periods and all the
    values in `min_ph_iterval_list`.

    Parameters
        d (Data object): the Data object containing the timestamps. You need to
            call d.calc_bg() before calling this function in order to create
            the background periods.
        min_delta_ph_list (list or array): list of minimum photon separation
            values above which photons are considered background. Unit: us.
        bg_fit_fun (function): function used to fit the background.
        ph_sel (Ph_sel object): photon selection on which the background is
            computed. See :class:`fretbursts.ph_sel.Ph_sel` for details.

    Returns
        Two arrays for background rate and fit-error of shape
        (nch, nperiods, len(min_delta_ph_list)).
    """

    BG = np.zeros((d.nch, d.nperiods, np.size(min_delta_ph_list)))
    BG_err = np.zeros_like(BG)
    BG[:], BG_err[:] = None, None

    for ich in range(d.nch):
        for period, ph in enumerate(
                d.iter_ph_times_period(ich=ich, ph_sel=ph_sel)):
            for i_min, min_delta_ph in enumerate(min_delta_ph_list):
                try:
                    BG[ich, period, i_min], BG_err[ich, period, i_min] = \
                            bg_fit_fun(ph, tail_min_us=min_delta_ph,
                                       clk_p=d.clk_p, **kwargs)
                except AssertionError:
                    # There are not enough delays with current threshold
                    break   # Skip remaining values in min_delta_ph_list
    return BG, BG_err


def fit_var_tail_us(d, Tail_min_us_list, t_max_s,
                    bg_fit_fun=exp_fit, ph_sel=Ph_sel('all'), **kwargs):
    """
    Fit BG of a ph_sel on all CH for all the values in `Tail_min_us_list`.

    The BG is fitted from the first `t_max_s` seconds of the measurement.

    Returns
        Two arrays for background rate and fit-error of shape
        (nch, len(min_delta_ph_list)).
    """
    assert ph_sel in [Ph_sel('all'), Ph_sel(Dex='Dem'), Ph_sel(Dex='Aem')]
    BG = np.zeros((d.nch, np.size(Tail_min_us_list)))
    BG_err = np.zeros((d.nch, np.size(Tail_min_us_list)))
    t_max_clk = t_max_s/d.clk_p

    ph_times_m_slice, A_em_slice = [], []
    for ph, a_em in zip(d.ph_times_m, d.A_em):
        ph_times_m_slice.append(ph[ph < t_max_clk])
        A_em_slice.append(a_em[ph < t_max_clk])

    Ph_times = []
    for ph, a_em in zip(ph_times_m_slice, A_em_slice):
        if ph_sel == Ph_sel(Dex='Dem'):
            Ph_times.append(ph[-a_em])
        elif ph_sel == Ph_sel(Dex='Aem'):
            Ph_times.append(ph[a_em])
        else:
            Ph_times.append(ph)

    for ch, ph_t in enumerate(Ph_times):
        for it, t in enumerate(Tail_min_us_list):
            try:
                BG[ch, it], BG_err[ch, it] = bg_fit_fun(
                    ph_t, tail_min_us=t, clk_p=d.clk_p, **kwargs)
            except:
                break # Skip remaining Tail_min
    return BG, BG_err


##
# Other functions
#
def histo(ph, bin_ms=10., t_max_s=None, clk_p=12.5e-9):
    """Returns an histogram and bins-centers of ph (ph arrival times)."""
    if t_max_s is not None:
        ph = ph[ph <= t_max_s/clk_p]
    bins = np.arange(ph[0], ph[-1]+1, (bin_ms*1e-3)/clk_p)
    H = np.histogram(ph, bins=bins)
    tt = H[0]
    ti = H[1][:-1]+0.5*(H[1][1]-H[1][0])
    return tt, ti

def gauss_fit(ph, bin_ms=10, clk_p=12.5e-9):
    """Returns the BG rate of ph calculated from the hist (PDF) of timetrace.
    """
    assert np.size(ph) > 0
    tt, ti = histo(ph=ph, bin_ms=bin_ms, clk_p=clk_p)
    #mu, sig = gaussian_fit(tt, mu_sigma_guess=[tt.mean(), tt.std()])
    mu, sig = gaussian_fit_hist(tt, mu0=tt.mean(), sigma0=tt.std())
    mask = (tt < (mu+3*sig))*(tt > (mu-3*sig))
    tt2 = tt[mask]
    #mu2, sig2 = gaussian_fit(tt2, mu_sigma_guess=[mu,sig])
    mu2, sig2 = gaussian_fit_hist(tt2, mu0=mu, sigma0=sig)
    return mu2/(bin_ms*1e-3)#, sig2/(bin_ms*1e-3)


##
# Experimental functions
#
def smart_bg(d, ich=0, bin_=50e-3, step=1):
    """BG calculation through binning (WARNING: very slow!)."""
    bg = []
    t = d.ph_times_m[ich]*d.clk_p
    t_max = np.floor(t.max())
    pprint(" Calculation started:")
    for s in np.arange(step, t_max, step):
        #if (s % (t_max/50) == 0): pprint(" %d %%" % (s/t_max*100))
        h = np.histogram(t[(t < s)*(t > (s-step))],
                         bins=np.arange(s-step, s+1e-3, bin_))
        #print(h[0])
        bg.append(h[0].min())
    pprint('\n')
    return np.array(bg)/bin_
