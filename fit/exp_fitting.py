#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
Generic functions to fit an exponential populations.
"""

import numpy as np
from scipy.stats import linregress
from scipy.optimize import leastsq


def expon_fit(s, s_min=0):
    """Eponential fit of samples s using MLE.
    All samples < s_min are discarded (s_min must be >= 0).
    Returns the lambda parameter (1/life-time) of the exponential.
    """
    if s_min > 0: s = s[s >= s_min]
    assert s.size > 0    
    Tau = s.mean() - s_min
    return 1./Tau

def expon_fit_cdf(s, s_min=0):
    """Eponential fit of samples s using a curve fit to the empirical CDF.
    All samples < s_min are discarded (s_min must be >= 0).
    Returns the lambda parameter (1/life-time) of the exponential.
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
    """Eponential fit of samples s using a curve fit of the histogram.
    All samples < s_min are discarded (s_min must be >= 0).
    Returns the lambda parameter (1/life-time) of the exponential.
    """
    if s_min > 0: s = s[s >= s_min]
    assert s.size > 0
    
    H = np.histogram(s, bins=bins, density=True)
    x = H[1][:-1] + 0.5*(H[1][1] - H[1][0])
    y = H[0]
    
    exp_fun = lambda x, rate, x_min: np.exp(-(x - x_min)*rate)
    err_fun = lambda rate, x, y, x_min: exp_fun(x, rate, x_min) - y
    
    res = leastsq(err_fun, x0=1./(s.mean() - s_min), args=(x, y, s_min))
    rate = res[0]
    return rate

