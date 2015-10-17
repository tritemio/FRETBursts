#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""

The module :mod:`select_bursts` defines functions to select
bursts according to different criteria.

These functions are usually passed to
:meth:`Data.select_bursts() <fretbursts.burstlib.Data.select_bursts>`.
For example::

    ds = d.select_bursts(select_bursts.E, th1=0.2, th2=0.6)

returns a new object `ds` containing only the bursts of `d` that pass the
specified selection criterium (`E` between 0.2 and 0.6 in this case).

"""

from __future__ import print_function, absolute_import
from builtins import range, zip

import numpy as np
from scipy import stats

from .utils.misc import clk_to_s as _clk_to_s


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  BURSTS SELECTION FUNCTIONS
#

def str_G(gamma, gamma1):
    """Return a string to indicate if and how gamma or gamma1 were used."""
    if gamma1 is None: s = "G%.1f" % gamma
    else: s = "G1_%.1f" % gamma1
    return s

## Selection on E or S values
def E(d, ich=0, E1=-np.inf, E2=np.inf):
    """Select bursts with E between E1 and E2."""
    assert E1 <= E2, 'Threshold E1 (%.2f) must be <= of E2 (%.2f)' % (E1, E2)
    burst_mask = (d.E[ich] >= E1)*(d.E[ich] <= E2)
    return burst_mask, ''

def S(d, ich=0, S1=-np.inf, S2=np.inf):
    assert S1 <= S2, 'Threshold S1 (%.2f) must be <= of S2 (%.2f)' % (S1, S2)
    """Select bursts with S between S1 and S2."""
    burst_mask = (d.S[ich] >= S1)*(d.S[ich] <= S2)
    return burst_mask, ''

def ES(d, ich=0, E1=-np.inf, E2=np.inf, S1=-np.inf, S2=np.inf, rect=True):
    """Select bursts with E between E1 and E2 and S between S1 and S2.

    When `rect` is True the selection is rectangular otherwise is elliptical.
    """
    assert E1 <= E2, 'Threshold E1 (%.2f) must be <= of E2 (%.2f)' % (E1, E2)
    assert S1 <= S2, 'Threshold S1 (%.2f) must be <= of S2 (%.2f)' % (S1, S2)
    if rect:
        return ES_rect(d, ich, E1=E1, E2=E2, S1=S1, S2=S2)
    else:
        return ES_ellips(d, ich, E1=E1, E2=E2, S1=S1, S2=S2)

def ES_rect(d, ich=0, E1=-np.inf, E2=np.inf, S1=-np.inf, S2=np.inf):
    """Select bursts inside the rectangle defined by E1, E2, S1, S2.
    """
    burst_mask = (d.S[ich] >= S1)*(d.S[ich] <= S2) * \
            (d.E[ich] >= E1)*(d.E[ich] <= E2)
    return burst_mask, ''

def ES_ellips(d, ich=0, E1=-1e3, E2=1e3, S1=-1e3, S2=1e3):
    """Select bursts with E-S inside an ellipsis inscribed in E1, E2, S1, S2.
    """
    def ellips(x, y, x1, x2, y1, y2):
        rx, ry = 0.5*abs(x2-x1), 0.5*abs(y2-y1)
        return ((x - np.mean([x1,x2]))/rx)**2 + ((y - np.mean([y1,y2]))/ry)**2

    burst_mask = (ellips(d.E[ich], d.S[ich], E1, E2, S1, S2) <= 1)
    return burst_mask, ''

## Selection on static burst size, width or period
def period(d, ich=0, bp1=0, bp2=None):
    """Select bursts from period bp1 to period bp2 (included)."""
    if bp2 is None: bp2 = d.bp[ich].max()
    burst_mask = (d.bp[ich] >= bp1)*(d.bp[ich] <= bp2)
    return burst_mask, ''

def time(d, ich=0, time_s1=0, time_s2=None):
    """Select the burst starting from time_s1 to time_s2 (in seconds)."""
    burst_start = d.mburst[ich].start * d.clk_p
    if time_s2 is None: time_s2 = burst_start.max()
    burst_mask = (burst_start >= time_s1)*(burst_start <= time_s2)
    return burst_mask, ''

def nd(d, ich=0, th1=20, th2=np.inf):
    """Select bursts with (nd >= th1) and (nd <= th2)."""
    assert th1 <= th2, 'th1 (%.2f) must be <= of th2 (%.2f)' % (th1, th2)
    bursts_mask = (d.nd[ich] >= th1)*(d.nd[ich] <= th2)
    return bursts_mask, ''

def na(d, ich=0, th1=20, th2=np.inf):
    """Select bursts with (na >= th1) and (na <= th2)."""
    assert th1 <= th2, 'th1 (%.2f) must be <= of th2 (%.2f)' % (th1, th2)
    bursts_mask = (d.na[ich] >= th1)*(d.na[ich] <= th2)
    return bursts_mask, ''

def naa(d, ich=0, th1=20, th2=np.inf):
    """Select bursts with (naa >= th1) and (naa <= th2)."""
    assert th1 <= th2, 'th1 (%.2f) must be <= of th2 (%.2f)' % (th1, th2)
    bursts_mask = (d.naa[ich] >= th1)*(d.naa[ich] <= th2)
    return bursts_mask, ''

def size(d, ich=0, th1=20, th2=np.inf, gamma=1., gamma1=None,
         add_naa=False, beta=None):
    """Select bursts with burst sizes between th1 and th2.

    The parameters `gamma`, `gamma1`, `add_naa` and `beta` are passed to
    :meth:`fretbursts.burstlib.Data.burst_sizes_ich` to compute the burst size.
    """
    assert th1 <= th2, 'th1 (%.2f) must be <= of th2 (%.2f)' % (th1, th2)
    burst_size = d.burst_sizes_ich(ich, gamma, gamma1, add_naa, beta)
    if d.nch > 1 and (np.size(th1) == d.nch): th1 = th1[ich]
    if d.nch > 1 and (np.size(th2) == d.nch): th2 = th2[ich]
    bursts_mask = (burst_size >= th1)*(burst_size <= th2)
    s = "size_th%d" % th1
    if th2 < 1000: s +="_th2_%d" % th2
    return bursts_mask, s+str_G(gamma, gamma1)

def width(d, ich=0, th1=0.5, th2=np.inf):
    """Select bursts with width between th1 and th2 (ms)."""
    assert th1 <= th2, 'th1 (%.2f) must be <= of th2 (%.2f)' % (th1, th2)
    th1, th2 = th1*1e-3/d.clk_p, th2*1e-3/d.clk_p
    burst_width = d.mburst[ich].width
    bursts_mask = (burst_width >= th1)*(burst_width <= th2)
    return bursts_mask, ''

def sbr(d, ich=0, th1=0, th2=np.inf):
    """Select bursts with SBR between `th1` and `th2`."""
    assert th1 <= th2, 'th1 (%.2f) must be <= of th2 (%.2f)' % (th1, th2)
    if 'sbr' not in d:
        d.calc_sbr()
    sbr_ich = d.sbr[ich]
    bursts_mask = (sbr_ich >= th1)*(sbr_ich <= th2)
    return bursts_mask, ''

def peak_phrate(d, ich=0, th1=0, th2=np.inf):
    """Select bursts with peak phtotons rate between th1 and th2 (cps).

    Note that this function requires to compute the peak photon rate
    first using :meth:`fretbursts.burstlib.Data.calc_max_rate`.
    """
    rate = d.max_rate[ich]
    mask = (rate >= th1)*(rate <= th2)
    return mask, ''

def brightness(d, ich=0, th1=0, th2=np.inf, add_naa=False, gamma=1):
    """Select bursts with size/width between th1 and th2 (cps).
    """
    sizes = d.burst_sizes(gamma=gamma, add_naa=add_naa)[ich]
    brightness = sizes/d.burst_widths[ich]
    mask = (brightness >= th1)*(brightness <= th2)
    return mask, ''


def nda_percentile(d, ich=0, q=50, low=False, gamma=1., gamma1=None,
                   add_naa=False):
    """Select bursts with SIZE >= q-percentile (or <= if `low` is True)

    `gamma`, `gamma1` and `add_naa` are passed to
    :meth:`fretbursts.burstlib.Data.burst_sizes_ich` to compute the burst size.
    """
    burst_size = d.burst_sizes_ich(ich, gamma, gamma1, add_naa)
    q_percentile = np.percentile(burst_size, q=q)
    if low: bursts_mask = (burst_size <= q_percentile)
    else: bursts_mask = (burst_size >= q_percentile)
    return bursts_mask, 'perc%d' % q

def topN_nda(d, ich=0, N=500, gamma=1., gamma1=None, add_naa=False):
    """Select the N biggest bursts in the channel.

    `gamma`, `gamma1` and `add_naa` are passed to
    :meth:`fretbursts.burstlib.Data.burst_sizes_ich` to compute the burst size.
    """
    burst_size = d.burst_sizes_ich(ich, gamma, gamma1, add_naa)
    index_sorted = burst_size.argsort()
    burst_mask = np.zeros(burst_size.size, dtype=bool)
    burst_mask[index_sorted[-N:]] = True
    return burst_mask, 'topN%d%s' % (N, str_G(gamma,gamma1))

def topN_max_rate(d, ich=0, N=500):
    """Select `N` bursts with the highest max burst rate.
    """
    max_burst_rate = d.max_rate[ich]
    index_sorted = max_burst_rate.argsort()
    burst_mask = np.zeros(max_burst_rate.size, dtype=bool)
    burst_mask[index_sorted[-N:]] = True
    return burst_mask, 'topN_MaxRate%d' % N

def topN_sbr(d, ich=0, N=200):
    """Select the top `N` bursts with hightest SBR."""
    if 'sbr' not in d:
        d.calc_sbr()
    sbr_ich = d.sbr[ich]
    index_sorted = sbr_ich.argsort()
    burst_mask = np.zeros(sbr_ich.size, dtype=bool)
    burst_mask[index_sorted[-N:]] = True
    return burst_mask, 'top%d_sbr' % N



## Selection on burst time (nearby, overlapping or isolated bursts)
def single(d, ich=0, th=1):
    """Select bursts that are at least th millisec apart from the others."""
    th = th*1e-3/d.clk_p
    burst_start = d.mburst[ich].start
    burst_stop = d.mburst[ich].stop
    gap_mask = (burst_start[1:] - burst_stop[:-1]) >= th
    bursts_mask = np.hstack([gap_mask,False])*np.hstack([False,gap_mask])
    return bursts_mask, ''

def consecutive(d, ich=0, th1=0, th2=np.inf, kind='both'):
    """Select consecutive bursts with th1 <= separation <= th2 (in sec.).

    Arguments:
        kind (string): valid values are 'first' to select the first burst
            of each pair, 'second' to select the second burst of each pair
            and 'both' to select both bursts in each pair.
    """
    assert th1 <= th2, 'th1 (%.2f) must be <= of th2 (%.2f)' % (th1, th2)
    kind_valids = ['first', 'second', 'both']
    assert kind in kind_valids, \
        "Invalid value for 'kind'. Valid values are %s" % str(kind_valids)

    bseparation = d.mburst[ich].separation*d.clk_p
    pair_mask = (bseparation >= th1)*(bseparation <= th2)

    bursts_mask = np.zeros(d.num_bursts[ich], dtype=bool)
    if kind in ['first', 'both']:
        bursts_mask += np.hstack([pair_mask, (False,)])
    if kind in ['second', 'both']:
        bursts_mask += np.hstack([(False,), pair_mask])
    return bursts_mask, ''

## Selection on burst size vs BG
def nd_bg(d, ich=0, F=5):
    """Select bursts with (nd >= bg_dd*F)."""
    bg_burst = d.bg_dd[ich][d.bp[ich]] * d.mburst[ich].width * d.clk_p
    bursts_mask = (d.nd[ich] >= F*bg_burst)
    return bursts_mask, ''

def na_bg(d, ich=0, F=5):
    """Select bursts with (na >= bg_ad*F)."""
    bg_burst = d.bg_ad[ich][d.bp[ich]] * d.mburst[ich].width * d.clk_p
    bursts_mask = (d.na[ich] >= F*bg_burst)
    return bursts_mask, ''

def naa_bg(d, ich=0, F=5):
    """Select bursts with (naa >= bg_aa*F)."""
    bg_burst = d.bg_aa[ich][d.bp[ich]] * d.mburst[ich].width * d.clk_p
    bursts_mask = (d.naa[ich] >= F*bg_burst)
    return bursts_mask, ''

def nt_bg(d, ich=0, F=5):
    """Select bursts with (nt >= bg*F)."""
    bg_burst = d.bg[ich][d.bp[ich]] * d.mburst[ich].width * d.clk_p
    bursts_mask = (d.nt[ich] > F*bg_burst)
    return bursts_mask, ''

## Selection on burst size vs BG (probabilistic)
def na_bg_p(d, ich=0, P=0.05, F=1.):
    """Select bursts w/ AD signal using P{F*BG>=na} < P."""
    accept_ch_bg_rate = d.rate_ad[ich]
    bursts_width = _clk_to_s(d.mburst[ich].width)
    max_num_bg_ph = stats.poisson(F*accept_ch_bg_rate*bursts_width).isf(P)
    #print("Min num. ph = ",  max_num_bg_ph)
    bursts_mask = (d.na[ich] >= max_num_bg_ph)
    return bursts_mask, ''

def nd_bg_p(d, ich=0, P=0.05, F=1.):
    """Select bursts w/ DD signal using P{F*BG>=nd} < P."""
    donor_ch_bg_rate = d.rate_dd[ich]
    bursts_width = _clk_to_s(d.mburst[ich].width)
    max_num_bg_ph = stats.poisson(F*donor_ch_bg_rate*bursts_width).isf(P)
    #print("Min num. ph = ", max_num_bg_ph)
    bursts_mask = (d.nd[ich] >= max_num_bg_ph)
    return bursts_mask, ''

def naa_bg_p(d, ich=0, P=0.05, F=1.):
    """Select bursts w/ AA signal using P{F*BG>=naa} < P."""
    A_em_ex_bg_rate = d.rate_aa[ich]
    bursts_width = _clk_to_s(d.mburst[ich].width)
    max_num_bg_ph = stats.poisson(F*A_em_ex_bg_rate*bursts_width).isf(P)
    #print("Min num. ph = ", max_num_bg_ph)
    bursts_mask = (d.naa[ich] >= max_num_bg_ph)
    return bursts_mask, ''

def nt_bg_p(d, ich=0, P=0.05, F=1.):
    """Select bursts w/ signal using P{F*BG>=nt} < P."""
    bg_rate = d.rate_m[ich]
    bursts_width = _clk_to_s(d.mburst[ich].width)
    max_num_bg_ph = stats.poisson(F*bg_rate*bursts_width).isf(P)
    #print("Min num. ph = ", max_num_bg_ph)
    #print("burst width (ms) = ", bursts_width*1e3)
    #print("Poisson rate = ", bg_rate*bursts_width)
    #print("rate = ", bg_rate)
    bursts_mask = (d.nt[ich] >= max_num_bg_ph)
    return bursts_mask, ''


## Old selection functions

### Selection on burst skeness
### this uses bleaching() from burstlib_misc.py
#def centered(d, ich=0, th=0.1):
#    """Select bursts with absolute value of skewness index less than th."""
#    skew_index,_,_ = bleaching(d, ich=ich, exclude_nan=False)
#    bursts_mask = (skew_index <= th)*(skew_index >= -th)
#    return bursts_mask
#
### this uses bleaching() from burstlib_misc.py
#def skewness(d, ich=0, th1=0.2, th2=1, **kwargs):
#    """Select bursts with skeness between th1 and th2."""
#    skew_index,_,_ = bleaching(d, ich=ich, **kwargs)
#    bursts_mask = (skew_index <= th2)*(skew_index >= th1)
#    return bursts_mask

def size_noise(d, ich=0, th=2):
    """Select bursts w/ size th times above the noise on both D and A ch."""
    burst_width = d.mburst[ich].width
    noise_d, noise_a = burst_width*d.rate_dd[ich], burst_width*d.rate_ad[ich]
    bursts_mask = (d.nd[ich] >= th*noise_d)*(d.na[ich] >= th*noise_a)
    return bursts_mask, ''

def size_noise_or(d, ich=0, th=2):
    """Select bursts w/ size th times above the noise on D or A ch."""
    burst_width = d.mburst[ich].width
    noise_d, noise_a = burst_width*d.rate_dd[ich], burst_width*d.rate_ad[ich]
    bursts_mask = (d.nd[ich] >= th*noise_d)+(d.na[ich] >= th*noise_a)
    return bursts_mask, ''

## this uses prob_to_be_bg from burstlib_misc.py
#def no_bg(d, ich=0, P=0.005, NF=1.):
#    """Select bursts with prob. to be from BG < P."""
#    burst_prob = prob_to_be_bg(d, ich=ich, NF=NF)
#    bursts_mask = (burst_prob < P)
#    return bursts_mask

def fret_value(d, ich=0, E=0.5, P_th=0.01):
    """Select bursts with prob. > P_th to have FRET equal to `E`."""
    bsizes = np.around(d.nd[ich]+d.na[ich]).astype(np.uint16)
    bursts_mask = np.zeros(bsizes.size, dtype=bool)
    for burst_size in range(bsizes.min(), bsizes.max()+1):
        indexes = np.where(bsizes == burst_size)
        RV = stats.binom(burst_size, E)
        #accept_num = arange(burst_size+1)
        #y = RV.cdf(accept_num)
        #min_accept_num = interp(th, y,accept_num)
        min_accept_num = RV.ppf(P_th) # ppf: percent point function (cdf^-1)
        bursts_mask[indexes] = (d.na[ich][indexes] > min_accept_num)
    return bursts_mask, ''
