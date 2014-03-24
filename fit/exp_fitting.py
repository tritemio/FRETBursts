#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
Generic functions to fit exponential populations.
"""

import numpy as np
from scipy.stats import linregress
from scipy.optimize import leastsq


def expon_fit(s, s_min=0):
    """Exponential fit of samples s using MLE.

    Arguments:
        s (array): array of exponetially-distributed samples
        s_min (float): all samples < `s_min` are discarded 
            (`s_min` must be >= 0).

    Returns:
        The lambda parameter (1/life-time) of the exponential.
    """
    if s_min > 0: s = s[s >= s_min]
    assert s.size > 0    
    Tau = s.mean() - s_min
    return 1./Tau

def expon_fit_cdf(s, s_min=0):
    """Exponential fit of samples s using a curve fit to the empirical CDF.

    Arguments:
        s (array): array of exponetially-distributed samples
        s_min (float): all samples < `s_min` are discarded 
            (`s_min` must be >= 0).

    Returns:
        The lambda parameter (1/life-time) of the exponential.
    """
    if s_min > 0: s = s[s >= s_min]
    assert s.size > 0    
    ## Empirical CDF
    ecdf = [np.sort(s), np.arange(0.5,s.size+0.5)*1./s.size]
    decr_line = np.log(1-ecdf[1])
    L = linregress(ecdf[0], decr_line)
    Lambda = -L[0]
    return Lambda

def expon_fit_hist(s, bins, s_min=0, weights=None):
    """Exponential fit of samples s using a curve fit of the histogram.

    Arguments:
        s (array): array of exponetially-distributed samples
        bins (float or array): if float is the bin width, otherwise is the
            array of bin edges (passed to `numpy.histogram`)
        s_min (float): all samples < `s_min` are discarded 
            (`s_min` must be >= 0).
        weights (None or string): if None no weights is applied.
            if is 'hist_counts', each bin has a weight equal to its counts
            in the minimization.
    
    Returns:
        The lambda parameter (1/life-time) of the exponential.
    """
    if s_min > 0: s = s[s >= s_min]
    assert s.size > 0
    
    counts, bins = np.histogram(s, bins=bins, density=True)
    x = bins[:-1] + 0.5*(bins[1] - bins[0])        
    y = counts
    x = x[y > 0]
    y = y[y > 0]
    
    if weights is None:
        w = np.ones(y.size)
    elif weights == 'hist_counts':
        w = np.sqrt(y*s.size*(bins[1]-bins[0]))
    else:
        raise ValueError('Weighting scheme not valid (use: None, or '
                         '"hist_counts")')
    
    exp_fun = lambda x, rate, x_min: rate*np.exp(-(x - x_min)*rate)
    err_fun = lambda rate, x, y, x_min, w: (exp_fun(x, rate, x_min) - y)*w
    
    res = leastsq(err_fun, x0=1./(s.mean() - s_min), 
                  args=(x, y, s_min, w))
    rate = res[0]
    return rate