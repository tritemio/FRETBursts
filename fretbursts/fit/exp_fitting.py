#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014-2016 Antonino Ingargiola <tritemio@gmail.com>
#
"""
Generic functions to fit exponential populations.

These functions can be used directly, or, in a typical FRETBursts workflow
they are passed to higher level methods.

*See also:*

* :doc:`background`
"""

from builtins import range, zip

import numpy as np
from scipy.stats import linregress, expon
from scipy.optimize import leastsq


def get_ecdf(s, offset=0.5):
    """Return arrays (x, y) for the empirical CDF curve of sample `s`.

    See the code for more info (is a one-liner!).

    Arguments:
        s (array of floats): sample
        offset (float, default 0.5): Offset to add to the y values of the CDF

    Returns:
        (x, y) (tuple of arrays): the x and y values of the empirical CDF
    """
    return np.sort(s), np.arange(offset, s.size+offset)*1./s.size

def get_residuals(s, tau_fit, offset=0.5):
    """Returns residuals of sample `s` CDF vs an exponential CDF.

    Arguments:
        s (array of floats): sample
        tau_fit (float): mean waiting-time of the exponential distribution
            to use as reference
        offset (float): Default 0.5. Offset to add to the empirical CDF.
            See :func:`get_ecdf` for details.

    Returns:
        residuals (array): residuals of empirical CDF compared with analytical
        CDF with time constant `tau_fit`.
    """
    x, y = get_ecdf(s, offset=offset)
    ye = expon.cdf(x, scale=tau_fit)
    residuals = y - ye
    return x, residuals

def expon_fit(s, s_min=0, offset=0.5, calc_residuals=True):
    """Fit sample `s` to an exponential distribution using the ML estimator.

    This function computes the rate (Lambda) using the maximum likelihood (ML)
    estimator of the mean waiting-time (Tau), that for an exponentially
    distributed sample is the sample-mean.

    Arguments:
        s (array): array of exponetially-distributed samples
        s_min (float): all samples < `s_min` are discarded
            (`s_min` must be >= 0).
        offset (float): offset for computing the CDF. See :func:`get_ecdf`.
        calc_residuals (bool): if True compute the residuals of the fitted
            exponential versus the empirical CDF.

    Returns:
        A 4-tuple of the fitted rate (1/life-time), residuals array,
        residuals x-axis array, sample size after threshold.
    """
    if s_min > 0: s = s[s >= s_min] - s_min
    assert s.size > 10

    # Maximum likelihood estimator of the waiting-time
    Lambda = 1./s.mean()

    x_residuals, residuals = None, None
    if calc_residuals:
        x_residuals, residuals = get_residuals(s, tau_fit=1./Lambda,
                                               offset=offset)
    return Lambda, residuals, x_residuals, s.size

def expon_fit_cdf(s, s_min=0, offset=0.5, calc_residuals=True):
    """Fit of an exponential model to the empirical CDF of `s`.

    This function computes the rate (Lambda) fitting a line (linear
    regression) to the log of the empirical CDF.

    Arguments:
        s (array): array of exponetially-distributed samples
        s_min (float): all samples < `s_min` are discarded
            (`s_min` must be >= 0).
        offset (float): offset for computing the CDF. See :func:`get_ecdf`.
        calc_residuals (bool): if True compute the residuals of the fitted
            exponential versus the empirical CDF.

    Returns:
        A 4-tuple of the fitted rate (1/life-time), residuals array,
        residuals x-axis array, sample size after threshold.
    """
    if s_min > 0: s = s[s >= s_min] - s_min
    assert s.size > 10

    # Line fit the log of the eCDF to compute the rate (Lambda)
    ecdf = get_ecdf(s, offset=offset)
    decr_line = np.log(1-ecdf[1])
    L = linregress(ecdf[0], decr_line)
    Lambda = -L[0]

    x_residuals, residuals = None, None
    if calc_residuals:
        x_residuals, residuals = get_residuals(s, tau_fit=1./Lambda,
                                               offset=offset)
    return Lambda, residuals, x_residuals, s.size


def expon_fit_hist(s, bins, s_min=0, weights=None, offset=0.5,
                   calc_residuals=True):
    """Fit of an exponential model to the histogram of `s` using least squares.

    Arguments:
        s (array): array of exponetially-distributed samples
        bins (float or array): if float is the bin width, otherwise is the
            array of bin edges (passed to `numpy.histogram`)
        s_min (float): all samples < `s_min` are discarded
            (`s_min` must be >= 0).
        weights (None or string): if None no weights is applied.
            if is 'hist_counts', each bin has a weight equal to its counts
            if is 'inv_hist_counts', the weight is the inverse of the counts.
        offset (float): offset for computing the CDF. See :func:`get_ecdf`.
        calc_residuals (bool): if True compute the residuals of the fitted
            exponential versus the empirical CDF.

    Returns:
        A 4-tuple of the fitted rate (1/life-time), residuals array,
        residuals x-axis array, sample size after threshold.
    """
    if s_min > 0: s = s[s >= s_min] - s_min
    assert s.size > 10

    counts, bins = np.histogram(s, bins=bins, density=True)
    x = bins[:-1] + 0.5*(bins[1] - bins[0])  # bin center position
    y = counts
    x = x[y > 0]
    y = y[y > 0]

    if weights is None:
        w = np.ones(y.size)
    elif weights == 'hist_counts':
        w = np.sqrt(y*s.size*(bins[1]-bins[0]))
    elif weights == 'inv_hist_counts':
        w = np.sqrt(1./y*s.size*(bins[1]-bins[0]))
    else:
        raise ValueError('Weighting scheme not valid (use: None, or '
                         '"hist_counts")')

    exp_fun = lambda x, rate: rate*np.exp(-x*rate)
    err_fun = lambda rate, x, y, w: (exp_fun(x, rate) - y)*w

    res = leastsq(err_fun, x0=1./(s.mean()), args=(x, y, w))
    Lambda = res[0]

    x_residuals, residuals = None, None
    if calc_residuals:
        x_residuals, residuals = get_residuals(s, tau_fit=1./Lambda,
                                               offset=offset)
    return Lambda, residuals, x_residuals, s.size
