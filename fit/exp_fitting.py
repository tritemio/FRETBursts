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

def expon_fit_hist(s, bins, s_min=0):
    """Exponential fit of samples s using a curve fit of the histogram.

    Arguments:
        s (array): array of exponetially-distributed samples
        bins (float or array): if float is the bin width, otherwise is the
            array of bin edges (passed to `numpy.histogram`)
        s_min (float): all samples < `s_min` are discarded 
            (`s_min` must be >= 0).
    
    Returns:
        The lambda parameter (1/life-time) of the exponential.
    """
    if s_min > 0: s = s[s >= s_min]
    assert s.size > 0
    
    H = np.histogram(s, bins=bins, density=True)
    x = H[1][:-1] + 0.5*(H[1][1] - H[1][0])        
    x = x[H[0] > 0]
    y = np.log(H[0][H[0] > 0])
    
    exp_fun = lambda x, rate, x_min: (-(x - x_min)*rate)
    err_fun = lambda rate, x, y, x_min: exp_fun(x, rate, x_min) - y
    
    res = leastsq(err_fun, x0=1./(s.mean() - s_min), args=(x, y, s_min))
    rate = res[0]
    return rate

def expon_fit_histw(s, bins, s_min=0, weights=None):
    """Exponential fit of samples s using a curve fit of the histogram.

    Arguments:
        s (array): array of exponetially-distributed samples
        bins (float or array): if float is the bin width, otherwise is the
            array of bin edges (passed to `numpy.histogram`)
        s_min (float): all samples < `s_min` are discarded 
            (`s_min` must be >= 0).
    
    Returns:
        The lambda parameter (1/life-time) of the exponential.
    """
    assert weights is None or weights in ['y', 'inv_y']
    if s_min > 0: s = s[s >= s_min]
    assert s.size > 0
    
    H = np.histogram(s, bins=bins, density=True)
    x = H[1][:-1] + 0.5*(H[1][1] - H[1][0])
    y = H[0]
    if weights is None:
        w = np.ones(y.size)
    elif weights == 'y':
        w = y
    elif weights == 'inv_y':
        w = 1./y
    else:
        raise ValueError('Weighting scheme not valid (use: None, "y" or '
                         '"inv_y")')

    exp_fun = lambda x, rate, x_min: np.exp(-(x - x_min)*rate)
    err_fun = lambda rate, x, y, x_min, w: (exp_fun(x, rate, x_min) - y)*w
    
    res = leastsq(err_fun, x0=1./(s.mean() - s_min), 
                  args=(x, y, s_min, w))
    rate = res[0]
    return rate
