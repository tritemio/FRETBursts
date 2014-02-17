"""
This module contains routines for background rate estimation.
"""

import numpy as np
from utils.misc import pprint
from fit.exp_fitting import expon_fit, expon_fit_cdf, expon_fit_hist
from fit.gaussian_fitting import gaussian_fit_hist


def exp_fit(ph, fit_fun=expon_fit, tail_min_p=0.1, tail_min_us=None,
                     clk_p=12.5e-9):
    """Return BG rate for ph computed as mean of delays (above a min value).
    """
    dph = np.diff(ph)
    if tail_min_us is None:
        tail_min = dph.max()*tail_min_p
    else:
        tail_min = tail_min_us*1e-6/clk_p
    Lambda = fit_fun(dph, s_min=tail_min)/clk_p
    return Lambda

def exp_cdf_fit(ph, tail_min_p=0.1, tail_min_us=None, clk_p=12.5e-9):
    """Return BG rate for ph computed fitting the exponential CDF of delays
    """
    return exp_fit(ph, fit_fun=expon_fit_cdf, tail_min_p=tail_min_p,
                       tail_min_us=tail_min_us, clk_p=clk_p)

def raw_fit(ph, clk_p=12.5e-9):
    """Compute the "raw" rate: (number of ph / duration). """
    return ph.size/((ph[-1]-ph[0])*clk_p)

def hist_fit(ph, tail_min_us, binw=0.1e-3, clk_p=12.5e-9):
    """Returns the BG rate of ph calculated from the hist (PDF) of timetrace.
    """
    assert np.size(ph) > 0
    dph = np.diff(ph)
    tail_min = tail_min_us*1e-6/clk_p
    binw_clk = binw/clk_p
    bins = np.arange(0, dph.max(), binw_clk)
    Lambda = expon_fit_hist(dph, bins=bins, s_min=tail_min)/clk_p
    return Lambda

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
    mask = (tt<(mu+3*sig))*(tt>(mu-3*sig))
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
        h = np.histogram(t[(t<s)*(t>(s-step))],
                bins=np.arange(s-step, s+1e-3, bin_))
        print h[0]
        bg.append(h[0].min())
    pprint('\n')
    return np.array(bg)/bin_

