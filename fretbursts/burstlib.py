#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
This is the main library to import (or run) in order to load all the analysis
functions.

`burstslib.py` defines the fundamental object `Data()` that contains both the
experimental data (attributes) and the high-level analysis routines (methods).

Furthermore it loads all the remaining **FRETBursts** modules (except for
`loaders.py`).

For usage example see the IPython Notebooks in sub-folder "notebooks".
"""

import os
import hashlib
import numpy as np
import copy
from numpy import zeros, size, r_
import scipy.stats as SS

from utils.misc import pprint, clk_to_s, deprecate
from poisson_threshold import find_optimal_T_bga
import fret_fit
import bg_cache
from ph_sel import Ph_sel
from fretmath import gamma_correct_E, gamma_uncorrect_E

from burstsearch import burstsearchlib as bslib
from burstsearch.burstsearchlib import (
        itstart, iwidth, inum_ph, iistart, iiend, itend,
        # Burst search function
        bsearch,
        # Burst data functions
        b_start, b_end, b_width, b_istart, b_iend, b_size,
        b_rate, b_separation,
        # Photon counting function,
        mch_count_ph_in_bursts
        )

import background as bg
import select_bursts
import fit
from fit.gaussian_fitting import (gaussian_fit_hist,
                                  gaussian_fit_cdf,
                                  two_gaussian_fit_hist,
                                  two_gaussian_fit_hist_min,
                                  two_gaussian_fit_hist_min_ab,
                                  two_gaussian_fit_EM,
                                  two_gauss_mix_pdf,
                                  two_gauss_mix_ab,)


# Redefine some old functions that have been renamed so old scripts will not
# break but will print a warning
bg_calc_exp = deprecate(bg.exp_fit, 'bg_calc_exp', 'bg.exp_fit')
bg_calc_exp_cdf = deprecate(bg.exp_fit, 'bg_calc_exp_cdf', 'bg.exp_cdf_fit')


def _get_bsearch_func(pure_python=False):
    if pure_python:
        # return the python version
        return bslib.bsearch_py
    else:
        # or what is available
        return bsearch

def _get_mch_count_ph_in_bursts_func(pure_python=False):
    if pure_python:
        # return the python version
        return bslib.mch_count_ph_in_bursts_py
    else:
        # or what is available
        return mch_count_ph_in_bursts
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##  GLOBAL VARIABLES
##
#
# itstart, iwidth, inum_ph, iistart and others defined in burstsearch/bs.py
#

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  BURST SELECTION FUNCTIONS
#

def Select_bursts(d_orig, filter_fun, negate=False, nofret=False, **kwargs):
    """Uses `filter_fun` to select a sub-set of bursts from `d_orig`.

    Arguments:
        d_orig (Data object): the original Data() object to filter
        filter_fun (fuction): function used for burst selection
        negate (boolean): If True, negates (i.e. take the complementary)
            of the selection returned by `filter_fun`. Default `False`.
        nofret (boolean): If True do not recompute burst correction and FRET
            efficiencies in the new retuerned object. Default `False`.

    Kwargs:
        Any additional keyword argument is passed to `filter_fun()`.

    Returns:
        d_new: a new Data() object containing the selected bursts.

    Warning:
        This function is also avaliable with the shortcut name `Sel`.

    Note:
        In order to save RAM `ph_times_m` is shared between `d_orig` and
        `d_new`. However, all the bursts data (`mburst`, `nd`, `na`, etc...)
        are new objects.
    """
    Masks, str_sel = Sel_mask(d_orig, filter_fun, negate=negate,
            return_str=True, **kwargs)
    d_sel = Sel_mask_apply(d_orig, Masks, nofret=nofret, str_sel=str_sel)
    return d_sel

# Alias for :func:`Select_bursts`
Sel = Select_bursts

def Sel_mask(d_orig, filter_fun, negate=False, return_str=False, **kwargs):
    """Returns a list of nch masks to select bursts according to filter_fun.
    The passed 'filter_fun' is used to compute the mask for each channel.

    Use this function only if you want to apply a selection from one object
    to a second object. Otherwise use :func:`Select_bursts`.

    See also:
        :func:`Select_bursts`, :func:`Sel_mask_apply`
    """
    ## Create the list of bool masks for the bursts selection
    M = [filter_fun(d_orig,i,**kwargs) for i in range(d_orig.nch)]
    Masks = [-m[0] if negate else m[0] for m in M]
    str_sel = M[0][1]
    if return_str: return Masks, str_sel
    else: return Masks

def Sel_mask_apply(d_orig, Masks, nofret=False, str_sel=''):
    """Returns a new Data object with bursts select according to Masks.
    Note that 'ph_times_m' is shared to save RAM, but `mburst`, `nd`, `na`,
    `nt`, `bp` (and `naa` if ALEX) are new objects.

    Use this function only if you want to apply a selection from one object
    to a second object. Otherwise use :func:`Select_bursts`.

    See also:
        :func:`Select_bursts`, :func:`Sel_mask`,
    """
    ## Attributes of ds point to the same objects of d_orig
    ds = Data(**d_orig)

    ## Copy the per-burst fields that must be filtered
    used_fields = [field for field in Data.burst_fields if field in d_orig]
    for name in used_fields:

        # Recreate the current attribute as a new list to avoid modifying
        # the old list that is also in the original object.
        # The list is initialized with empty arrays because this is the valid
        # value when a ch has no bursts.
        ds.add(**{ name: [np.array([])]*d_orig.nch })

        # Assign the new data
        for ich, mask in enumerate(Masks):
            if d_orig[name][ich].size == 0: continue # -> no bursts in ch
            # Note that boolean masking implies numpy array copy
            # On the contrary slicing only makes a new view of the array
            ds[name][ich] = d_orig[name][ich][mask]

    # Recompute E and S
    if not nofret:
        ds.calc_fret(count_ph=False)
    # Add the annotation about the filter function
    ds.s = list(d_orig.s+[str_sel]) # using append would modify also d_orig
    return ds


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Bursts and Timestamps utilities
#
def get_alex_fraction(on_range, alex_period):
    """Get the fraction of period beween two numbers indicating a range.
    """
    assert len(on_range) == 2
    if on_range[0] < on_range[1]:
        fraction = 1.*(on_range[1] - on_range[0])/alex_period
    else:
        fraction = 1.*(alex_period + on_range[1] - on_range[0])/alex_period
    return fraction

def top_tail(nx, a=0.1):
    """Return for each ch the mean size of the top `a` fraction.
    nx is one of nd, na, nt from Data() (list of burst size in each ch).
    """
    assert a>0 and a<1
    return np.r_[[n[n>n.max()*(1-a)].mean() for n in nx]]

# Quick functions to calculate rate-trace from ph_times
def ph_rate(m, ph):
    """Return an array of m-photons rates of size ph.size - m + 1."""
    return 1.*m/(ph[m-1:] - ph[:ph.size-m+1])   # rate

def ph_rate_t(m, ph):
    """Return the mean time for each rate compute by `ph_rate`."""
    return 0.5*(ph[m-1:] + ph[:ph.size-m+1])  # time for rate

def ph_select(ph, mburst):
    """Return bool mask to select all ph inside any burst"""
    mask = zeros(ph.size, dtype=bool)
    iBurstStart, iBurstEnd = b_istart(mburst), b_iend(mburst)
    for istart, iend in zip(iBurstStart, iBurstEnd):
        mask[istart:iend+1] = True
    return mask

def mch_ph_select(PH, MBurst):
    """Multi-ch version of ph_select."""
    Mask = [ph_select(ph, mb) for ph, mb in zip(PH, MBurst)]
    return Mask


def b_ph_times1(b, ph_times, pad=0):
    """Returns a slice of ph_times inside one burst."""
    return ph_times[b[iistart]-pad:b[iiend]+pad+1]
def b_ph_times_v(bursts, ph_times, pad=0):
    """Returns a list of arrays containing ph_times inside each burst."""
    PH = [ph_times[b[iistart]-pad:b[iiend]+pad+1] for b in bursts]
    return PH

def b_rate_max(ph, m, mburst, mask=None):
    """Returns the max m-photons rate reached inside each burst.

    Arguments
        ph (1D array): array of photons timestamps
        m (int): number of timestamps to use to compute the rate
        mburst (2D array): array of burst data as returned by the burst search
            function
        mask (boolean mask or None): if not None, is a boolean mask
            to select photons in `ph` (for example Donor-ch photons).

    Return
        Array of max photon rate reached inside each burst.
    """
    PHB = []
    for burst in mburst:
        burst_slice = slice(burst[iistart], burst[iiend] + 1)
        burst_ph = ph[burst_slice]
        if mask is not None:
            burst_ph = burst_ph[mask[burst_slice]]
            if burst_ph.size < m:
                PHB.append(None)
                continue
        PHB.append(burst_ph)
    rates_max = np.array(
        [ph_rate(m=m, ph=phb).max() if phb is not None else 0 for phb in PHB]
        )
    return rates_max

def b_irange(bursts, b_index, pad=0):
    """Returns range of indices of ph_times inside one burst"""
    pad = np.array(pad).astype(bursts.dtype) # to avoid unwanted conversions
    _i_start = bursts[b_index, iistart]
    _i_end = bursts[b_index, iiend]
    return np.arange(_i_start-pad, _i_end+pad+1)

def b_ph_times(bursts, b_index, ph_times, pad=0):
    """Returns ph_times inside one burst, with "pad" ph before and after."""
    return ph_times[b_irange(bursts, b_index, pad=pad)]

def b_rates_inside(ph, b, bi, m=3, pad=0):
    """Returns all the m-ph-rates of burst #bi."""
    return ph_rate(m, b_ph_times(b, bi, ph, pad=pad))


def find_burst(bursts, size, width_ms, clk_p=12.5e-9):
    """Find b_index of burst(s) of given size AND width."""
    width = (width_ms*1e-3)/clk_p
    th = 0.01e-3/clk_p # 800clk or 10us @ clk_p=12.5e-9s
    return np.where((b_size(bursts) == size)*(abs(b_width(bursts)-width) < th))

def fuse_bursts_direct(mburst, ms=0, clk_p=12.5e-9, verbose=True):
    """Fuse bursts separated by less than `ms` (milli-secs).

    This function is a direct implementation using a single loop.
    For a faster implementation see :func:`fuse_bursts_iter`.

    Parameters:
        bursts (2D array): Nx6 array of burst data, one row per burst
            See `burstseach.burstseachlib.py` for details.
        ms (float):
            minimum waiting time between bursts (in millisec). Burst closer
            than that will be fuse in asingle burst.
        clk_p (float): clock period or timestamp units in seconds.
        verbose (bool): if True print a summary of fused bursts.

    Returns:
        new_bursts (2D array): new array of burst data
    """
    max_delay_clk = (ms*1e-3)/clk_p

    fused_bursts = []
    fused_burst = None
    for burst1, burst2 in zip(mburst[:-1], mburst[1:]):
        if fused_burst is not None:
            burst1c = fused_burst
        else:
            burst1c = burst1

        separation = burst2[itstart] - burst1c[itend]
        if separation <= max_delay_clk:
            fb_tstart = burst1c[itstart]
            fb_istart = burst1c[iistart]
            fb_tend = burst2[itend]
            fb_iend = burst2[iiend]
            fb_width = burst1c[iwidth] + burst2[iwidth]
            fb_num_ph = burst1c[inum_ph] + burst2[inum_ph]
            if burst1c[iiend] >= burst2[iistart]:
                n_overlap_ph = burst1c[iiend] - burst2[iistart] + 1
                fb_num_ph -= n_overlap_ph
                t_overlap = burst1c[itend] - burst2[itstart]
                fb_width -= t_overlap
            fused_burst = [fb_tstart, fb_width, fb_num_ph,
                           fb_istart, fb_iend, fb_tend]
        else:
            if fused_burst is not None:
                fused_bursts.append(fused_burst)
                fused_burst = None
            else:
                fused_bursts.append(burst1c)

    # Append the last bursts (either a fused one or isolated)
    if fused_burst is not None:
        fused_bursts.append(fused_burst)
    else:
        fused_bursts.append(burst2)

    init_nburst = mburst.shape[0]
    delta_b = init_nburst - len(fused_bursts)
    pprint(" --> END Fused %d bursts (%.1f%%)\n\n" %\
            (delta_b, 100.*delta_b/init_nburst), mute=-verbose)
    return np.array(fused_bursts)

def fuse_bursts_iter(bursts, ms=0, clk_p=12.5e-9, verbose=True):
    """Fuse bursts separated by less than `ms` (milli-secs).

    This function calls iteratively :func:`b_fuse` until there are no more
    bursts to fuse. See also :func:`fuse_bursts_direct`.

    Parameters:
        bursts (2D array): Nx6 array of burst data, one row per burst
            See `burstseach.burstseachlib.py` for details.
        ms (float):
            minimum waiting time between bursts (in millisec). Burst closer
            than that will be fused in a single burst.
        clk_p (float): clock period or timestamp units in seconds.
        verbose (bool): if True print a summary of fused bursts.

    Returns:
        new_bursts (2D array): new array of burst data
    """

    z = 0
    init_nburst = bursts.shape[0]
    new_nburst, nburst = 0, 1  # starting condition
    while (new_nburst < nburst):
        z += 1
        nburst = bursts.shape[0]
        bursts = b_fuse(bursts, ms=ms, clk_p=clk_p)
        new_nburst = bursts.shape[0]
    delta_b = init_nburst-nburst
    pprint(" --> END Fused %d bursts (%.1f%%, %d iter)\n\n" %\
            (delta_b, 100.*delta_b/init_nburst, z), mute=-verbose)
    return bursts

def b_fuse(mburst, ms=0, clk_p=12.5e-9):
    """Fuse bursts separated by less than `ms` (milli-secs).

    This is a low-level function that only fuses 2 consecutive bursts
    separated by less than `ms` millisec. If there are 3 or consecutive
    bursts separated by less than `ms` only the first 2 are fused.
    See :func:`fuse_bursts_iter` or :func:`fuse_bursts_direct` for
    higher level functions.

    Parameters:
        mburst (2D array): Nx6 array of burst data, one row per burst
            See `burstseach.burstseachlib.py` for details.
        ms (float):
            minimum waiting time between bursts (in millisec). Burst closer
            than that will be fused in a single burst.
        clk_p (float): clock period or timestamp units in seconds.

    Returns:
        new_mburst (2D array): new array of burst data
    """
    max_delay_clk = (ms*1e-3)/clk_p
    # Nearby bursts masks
    delays = (b_separation(mburst) <= max_delay_clk)
    first_burst = np.hstack([delays, (False,)])
    second_burst = np.hstack([(False,), delays])
    # Maintain just the 1st in case there were more than 2 consecutive bursts
    first_burst -= (second_burst*first_burst)
    second_burst = np.hstack([(False,), first_burst[:-1]])
    both_burst = first_burst + second_burst

    # istart is from the first bursts, iend is from the second burst
    fused_burst1 = mburst[first_burst, :] # slicing makes a copy
    fused_burst2 = mburst[second_burst, :]

    # pure bool mask, no copy, b will be changed
    #fused_burst1 = mburst[first_burst]
    #fused_burst2 = mburst[second_burst]

    num_ph = b_size(fused_burst1) + b_size(fused_burst2)
    overlap = b_iend(fused_burst1) - b_istart(fused_burst2) + 1
    # NOTE: overlap == 0 means touching but not overlapping bursts, if ph[i] is
    #       the last ph of burst1 ph[i+1] is the first ph of burst2
    overlap[overlap < 0] = 0
    num_ph -= overlap
    #[2] overlap_non_neg = overlap >= 0
    #[2] num_ph[overlap_non_neg] -= overlap[overlap_non_neg]
    #assert (num_ph <= (b_size(fused_burst1) + b_size(fused_burst2))).all()

    width = b_width(fused_burst1) + b_width(fused_burst2)
    t_overlap = b_end(fused_burst1) - b_start(fused_burst2)
    #assert (t_overlap[overlap > 0] >= 0).all()
    #assert (t_overlap[overlap == 1] == 0).all()
    width[overlap > 0] -= t_overlap[overlap > 0]
    #[2] width[overlap_non_neg] -= t_overlap[overlap_non_neg]
    # NOTE for [2]: overlap_non_neg includes also cases of overlap==0 for which
    #       t_overlap is negative, it's an arbitrary choice if in this case we
    #       should add (s2-e1) = -t_overlap to the new width. See paper notes.
    #assert (width <= (b_width(fused_burst1) + b_width(fused_burst2))).all()

    # Assign the new burst data
    # fused_burst1 has alredy the right tstart and istart
    fused_burst1[:, inum_ph] = num_ph
    fused_burst1[:, iwidth] = width
    fused_burst1[:, iiend] = b_iend(fused_burst2)
    fused_burst1[:, itend] = b_end(fused_burst2)

    new_burst = np.vstack([fused_burst1, mburst[-both_burst, :]])
    reorder = new_burst[:, itstart].argsort()
    return new_burst[reorder, :]

def mch_fuse_bursts(MBurst, ms=0, clk_p=12.5e-9, verbose=True):
    """Multi-ch version of `fuse_bursts`. `MBurst` is a list of arrays.
    """
    mburst = [b.copy() for b in MBurst] # safety copy
    new_mburst = []
    ch = 0
    for mb in mburst:
        ch += 1
        pprint(" - - - - - CHANNEL %2d - - - - \n" % ch, -verbose)
        if mb.size == 0:
            continue
        new_bursts = fuse_bursts_iter(mb, ms=ms, clk_p=clk_p, verbose=verbose)
        new_mburst.append(new_bursts)
    return new_mburst

def stat_burst(d, ich=0, fun=np.mean):
    """Compute a per-ch statistics (`fun`) for bursts.
    """
    tstart, width, num_ph, istart = 0, 1, 2, 3
    statg, statr = zeros(d.mburst[ich].shape[0]), zeros(d.mburst[ich].shape[0])
    statt = zeros(d.mburst[ich].shape[0])
    for i, b in enumerate(d.mburst[ich]):
        burst_slice = slice(b[istart], b[istart]+b[num_ph])
        # Arrival times of (all) ph in burst b
        ph = d.ph_times_m[ich][burst_slice]
        # Select only one color ch if requested
        r_mask = d.A_em[ich][burst_slice]
        g_mask = -r_mask
        #if r_mask.sum() == 0 or g_mask.sum() == 0:
        #    pprint("WARNING: Zero-size bursts.\n")
        statg[i] = fun(ph[g_mask].astype(float))
        statr[i] = fun(ph[r_mask].astype(float))
        statt[i] = fun(ph.astype(float))
    return statg, statr, statt

def burst_stats(mburst, clk_p=12.5*1e9):
    """Compute average duration, size and burst-delay for bursts in mburst.
    """
    width_stats = np.array([[b[:, 1].mean(), b[:, 1].std()] for b in mburst
        if len(b) > 0]).T
    height_stats = np.array([[b[:, 2].mean(), b[:, 2].std()] for b in mburst
        if len(b) > 0]).T
    mean_burst_delay = np.array([np.diff(b[:, 0]).mean() for b in mburst
        if len(b) > 0])
    return (clk_to_s(width_stats, clk_p)*1e3, height_stats,
            clk_to_s(mean_burst_delay, clk_p))

def print_burst_stats(d):
    """Print some bursts statistics."""
    nch = len(d.mburst)
    width_ms, height, delays = burst_stats(d.mburst, d.clk_p)
    s = "\nNUMBER OF BURSTS: m = %d, L = %d" % (d.m, d.L)
    s += "\nPixel:          "+"%7d "*nch % tuple(range(1, nch+1))
    s += "\n#:              "+"%7d "*nch % tuple([b.shape[0] for b in d.mburst])
    s += "\nT (us) [BS par] "+"%7d "*nch % tuple(np.array(d.T)*1e6)
    s += "\nBG Rat T (cps): "+"%7d "*nch % tuple(d.rate_m)
    s += "\nBG Rat D (cps): "+"%7d "*nch % tuple(d.rate_dd)
    s += "\nBG Rat A (cps): "+"%7d "*nch % tuple(d.rate_ad)
    s += "\n\nBURST WIDTH STATS"
    s += "\nPixel:          "+"%7d "*nch % tuple(range(1, nch+1))
    s += "\nMean (ms):      "+"%7.3f "*nch % tuple(width_ms[0, :])
    s += "\nStd.dev (ms):   "+"%7.3f "*nch % tuple(width_ms[1, :])
    s += "\n\nBURST SIZE STATS"
    s += "\nPixel:          "+"%7d "*nch % tuple(range(1, nch+1))
    s += "\nMean (# ph):    "+"%7.2f "*nch % tuple(height[0, :])
    s += "\nStd.dev (# ph): "+"%7.2f "*nch % tuple(height[1, :])
    s += "\n\nBURST MEAN DELAY"
    s += "\nPixel:          "+"%7d "*nch % tuple(range(1, nch+1))
    s += "\nDelay (s):      "+"%7.3f "*nch % tuple(delays)
    return s

def ES_histog(E, S, bin_step=0.05, E_bins=None, S_bins=None):
    """Returns 2D (ALEX) histogram and bins of bursts (E,S).
    """
    if E_bins is None:
        E_bins = np.arange(-0.6, 1.6+1e-4, bin_step)
    if S_bins is None:
        S_bins = np.arange(-0.6, 1.6+1e-4, bin_step)
    H, E_bins, S_bins = np.histogram2d(E, S, bins=[E_bins, S_bins])
    return H, E_bins, S_bins

def delta(x):
    """Return x.max() - x.min()"""
    return x.max() - x.min()

def mask_empty(mask):
    """Returns True if `mask` is empty, otherwise False.

    `mask` can be a boolean array or a slice object.
    """
    if type(mask) is slice:
        is_slice_empty = (mask.stop == 0)
        return is_slice_empty
    else:
        # Bolean array
        return not mask.any()

class DataContainer(dict):
    """
    Generic class for storing data.

    It's a dictionary in which each key is also an attribute d['nt'] or d.nt.
    """
    def __init__(self, **kwargs):
        dict.__init__(self, **kwargs)
        for k in self:
            dict.__setattr__(self, k, self[k])

    def add(self, **kwargs):
        """Adds or updates elements (attributes and/or dict entries). """
        self.update(**kwargs)
        for k, v in kwargs.items():
            setattr(self, k, v)
    def delete(self, *args):
        """Delete an element (attribute and/or dict entry). """
        for name in args:
            try:
                self.pop(name)
            except KeyError:
                print ' WARNING: Name %s not found (dict).' % name
            try:
                delattr(self, name)
            except AttributeError:
                print ' WARNING: Name %s not found (attr).' % name


class Data(DataContainer):
    """
    Container for all the information (timestamps, bursts) of a dataset.

    Data() contains all the information of a dataset (name, timestamps, bursts,
    correction factors) and provides several methods to perform analysis
    (background estimation, burst search, FRET fitting, etc...).

    When loading a measurement file a Data() object is created by one
    of the loader functions in `loaders.py`. Data() objects can be also
    created by some methods (`.copy()`, `.fuse_bursts()`, etc...) or by a
    "burst selection" with `Sel()` function.

    To add or delete data-attributes use `.add()` or `.delete()` methods.
    All the standard data-attributes are listed below.

    Note:
        Attributes of type "*list*" contain one element per channel.
        Each element, in turn, can be an array. For example `.ph_times_m[i]`
        is the array of timestamps for channel `i`; or `.nd[i]` is the array
        of donor counts in each burst for channel `i`.

    **Measurement attributes**

    Attributes:
        fname (string): measurements file name
        nch (int): number of channels
        clk_p (float): clock period in seconds for timestamps in `ph_times_m`
        ph_times_m (list): list of timestamp arrays (int64). Each array
            contains all the timestamps (donor+acceptor) in one channel.
        A_em (list): list of boolean arrays marking acceptor timestamps. Each
            array is a boolean mask for the corresponding ph_times_m array.
        leakage (float or array of floats): leakage (or bleed-through) fraction.
            May be scalar or same size as nch.
        gamma (float or array of floats): gamma factor.
            May be scalar or same size as nch.
        D_em (list of boolean arrays):  **[ALEX-only]**
            boolean mask for `.ph_times_m[i]` for donor emission
        D_ex, A_ex (list of boolean arrays):  **[ALEX-only]**
            boolean mask for `.ph_times_m[i]` during donor or acceptor
            excitation
        D_ON, A_ON (2-element tuples of int ): **[ALEX-only]**
            start-end values for donor and acceptor excitation selection.
        alex_period (int): **[ALEX-only]**
            duration of the alternation period in clock cycles.

    **Background Attributes**

    These attributes contain the estimated background rate. Each attribute is
    a list (one element per channel) of arrays. Each array contains several
    rates computed every `time_s` seconds of measurements. `time_s` is an
    argument passed to `.calc_bg()` method. Each `time_s` measurement slice
    is here called background **`period`**.

    Attributes:
        bg (list of arrays):  total background (donor + acceptor) during donor
            excitation
        bg_dd (list of arrays): background of donor emission during donor
            excitation
        bg_ad (list of arrays): background of acceptor emission during donor
            excitation
        bg_aa (list of arrays): background of acceptor emission during
            acceptor excitation
        nperiods (int): number of periods in which timestamps are split for
            background calculation
        bg_fun (function): function used to compute the background rates
        Lim (list): each element of this list is a list of index pairs for
            `.ph_times_m[i]` for **first** and **last** photon in each period.
        Ph_p (list): each element in this list is a list of timestamps pairs
            for **first** and **last** photon of each period.
        bg_ph_sel (Ph_sel object): photon selection used by Lim and Ph_p.
            See :class:`fretbursts.ph_sel.Ph_sel` for details.

    Other attributes (per-ch mean of `bg`, `bg_dd`, `bg_ad` and `bg_aa`)::

        rate_m: array of bg rates for D+A channel pairs (ex. 4 for 4 spots)
        rate_dd: array of bg rates for D em (and D ex if ALEX)
        rate_da: array of bg rates for A em (and D ex if ALEX)
        rate_aa: array of bg rates for A em and A ex (only for ALEX)

    **Burst search parameters (user input)**

    These are the parameters used to perform the burst search
    (see :meth:`burst_search_t`).

    Attributes:
        ph_sel (Ph_sel object): photon selection used for burst search.
            See :class:`fretbursts.ph_sel.Ph_sel` for details.
        m (int): number of consecutive timestamps used to compute the
            local rate during burst search
        L (int): min. number of photons for a burst to be identified and saved
        P (float, probability): valid values [0..1].
            Probability that a burst-start is due to a Poisson background.
            The employed Poisson rate is the one computed by `.calc_bg()`.
        F (float): `(F * background_rate)` is the minimum rate for burst-start

    **Burst search data (available after burst search)**

    When not specified, parameters marked as (list of arrays) contains arrays
    with one element per bursts. `mburst` arrays contain one "row" per burst.
    `TT` arrays contain one element per `period` (see above: background
    attributes).

    Attributes:
        mburst (list of arrays): list of 2-D arrays containing burst data.
            Each row contains data for 1 burst. Each column contains burst
            fields like burst start, end, duration, size and indexes.
            The proper way to access these fields is through the
            function b_* (ex: b_start, b_end, b_width, b_size, etc...).
            For more details see :mod:`fretbursts.burstsearch.burstsearchlib`.

        TT (list of arrays): list of arrays of *T* values (in sec.). A *T*
            value is the maximum delay between `m` photons to have a
            burst-start. Each channels has an array of *T* values, one for
            each background "period" (see above).
        T (array): per-channel mean of `TT`

        nd, na (list of arrays): number of donor or acceptor photons during
            donor excitation in each burst
        nt (list of arrays): total number photons (nd+na+naa)
        naa (list of arrays): number of acceptor photons in each bursts
            during acceptor excitation **[ALEX only]**
        bp (list of arrays): time period for each burst. Same shape as `nd`.
            This is needed to identify the background rate for each burst.
        bg_bs (list): background rates used for threshold computation in burst
            search (is a reference to `bg`, `bg_dd` or `bg_ad`).

        fuse (None or float): if not None, the burst separation in ms below
            which bursts have been fused (see `.fuse_bursts()`).

        E (list): FRET efficiency value for each burst:
                    E = na/(na + gamma*nd).
        S (list): stoichiometry value for each burst:
                    S = (gamma*nd + na) /(gamma*nd + na + naa)
    """

    # Attribute names containing per-photon data.
    # Each attribute is a list (1 element per ch) of arrays (1 element
    # per photon).
    ph_fields = ['ph_times_m', 'A_em', 'D_em', 'A_ex', 'D_ex']

    # Attribute names containing background data.
    # Each attribute is a list (1 element per ch) of sequences (1 element per
    # background period). For example `.bg` is a list of arrays, while `.Lim`
    # and `.Ph_p` are lists of lists-of-tuples (one tuple per background
    # period). These attributes do not exist before computing the background.
    bg_fields = ['bg', 'bg_dd', 'bg_ad', 'bg_da', 'bg_aa', 'Lim', 'Ph_p']

    # Attribute names containing per-burst data.
    # Each attribute is a list (1 element per ch) of arrays (1 element
    # per burst).
    # They do not necessarly exist. For example 'naa' exists only for ALEX
    # data. Also none of them exist before performing a burst search.
    burst_fields = ['E', 'S', 'mburst', 'nd', 'na', 'nt', 'bp', 'nda', 'naa',
                    'max_rate', 'sbr']

    # List of photon selections on which the background is computed
    _ph_streams = [Ph_sel('all'), Ph_sel(Dex='Dem'), Ph_sel(Dex='Aem'),
                   Ph_sel(Aex='Dem'), Ph_sel(Aex='Aem')]

    @property
    def ph_streams(self):
        if self.ALEX:
            return self._ph_streams
        else:
            return [Ph_sel('all'), Ph_sel(Dex='Dem'), Ph_sel(Dex='Aem')]

    def __init__(self, **kwargs):
        # Default values
        init_kw = dict(ALEX=False, leakage=0., gamma=1., chi_ch=1., dir_ex=0.,
                       s=[])
        # Override with user data
        init_kw.update(**kwargs)
        DataContainer.__init__(self, **init_kw)

    ##
    # Infrastructure methods
    #
    def copy(self, mute=False):
        """Copy data in a new object. All arrays copied except for ph_times_m
        """
        pprint('Deep copy executed.\n', mute)
        new_d = Data(**self) # this make a shallow copy (like a pointer)

        ## Deep copy (not just reference) or array data
        for field in self.burst_fields + self.bg_fields:
            # Making sure k is defined
            if field in self:

                # Make a deepcopy of the per-channel lists
                new_d[field] = copy.deepcopy(self[field])

                # Set the attribute: new_d.k = new_d[k]
                setattr(new_d, field, new_d[field])
        return new_d

    def ph_times_hash(self, hash_name='md5', hexdigest=True):
        """Return an hash for the timestamps arrays.
        """
        m = hashlib.new(hash_name)
        for ph in self.iter_ph_times():
            if type(ph) is np.ndarray:
                m.update(ph.data)
            else:
                # TODO Handle ph_times in PyTables files
                raise NotImplementedError
        if hexdigest:
            return m.hexdigest()
        else:
            return m

    def iter_ph_masks(self, ph_sel=Ph_sel('all')):
        """Iterator returning masks for `ph_sel` photons.

        Arguments:
            ph_sel (Ph_sel object): object defining the photon selection.
                See :class:`fretbursts.ph_sel.Ph_sel` for details.
        """
        for ich in xrange(self.nch):
            yield self.get_ph_mask(ich, ph_sel=ph_sel)

    def _check_ph_sel(self, ph_sel):
        """Check consistency of `ph_sel` with current data (ALEX vs not ALEX).
        """
        assert type(ph_sel) is Ph_sel

        if not self.ALEX and ph_sel == Ph_sel(Dex='DAem'):
            ph_sel = Ph_sel('all')

        # If Aex != None and not ALEX print a warning and set Aex to None
        not_ph_sel_all = ph_sel != Ph_sel('all')
        if not_ph_sel_all and (not self.ALEX) and (ph_sel.Aex is not None):
            print 'WARNING: Use of acceptor excitation with non-ALEX data.'
            ph_sel = Ph_sel(Dex=ph_sel.Dex)

        return ph_sel

    def get_ph_mask(self, ich=0, ph_sel=Ph_sel('all')):
        """Returns a mask for `ph_sel` photons in channel `ich`.

        The masks are either boolean arrays or slices (full or empty). In
        both cases they can be used to index the timestamps of the
        corresponding channel.

        Arguments:
            ph_sel (Ph_sel object): object defining the photon selection.
                See :class:`fretbursts.ph_sel.Ph_sel` for details.
        """
        ph_sel = self._check_ph_sel(ph_sel)

        # This is the only case in which Aex='DAem' for non-ALEX data is OK
        if ph_sel == Ph_sel('all'):
            # Note that slice(None) is equivalent to [:].
            # Also, numpy arrays are not copies when sliced.
            # So getting all photons with this mask is efficient
            # Note: the drawback is that the slice cannot be indexed
            #       (where a normal boolean array would)
            return slice(None)

        # Base selections
        elif ph_sel == Ph_sel(Dex='Dem'):
            return self.get_D_em_D_ex(ich)
        elif ph_sel == Ph_sel(Dex='Aem'):
            return self.get_A_em_D_ex(ich)
        elif ph_sel == Ph_sel(Aex='Dem'):
            return self.get_D_em(ich)*self.get_A_ex(ich)
        elif ph_sel == Ph_sel(Aex='Aem'):
            return self.get_A_em(ich)*self.get_A_ex(ich)

        # Selection of all photon in one emission ch
        elif ph_sel == Ph_sel(Dex='Dem', Aex='Dem'):
            return self.get_D_em(ich)
        elif ph_sel == Ph_sel(Dex='Aem', Aex='Aem'):
            return self.get_A_em(ich)

        # Selection of all photon in one excitation period
        elif ph_sel == Ph_sel(Dex='DAem'):
            return self.get_D_ex(ich)
        elif ph_sel == Ph_sel(Aex='DAem'):
            return self.get_A_ex(ich)

        # Selection of all photons except for Dem during Aex
        elif ph_sel == Ph_sel(Dex='DAem', Aex='Aem'):
            return self.get_D_ex(ich) + self.get_A_em(ich)*self.get_A_ex(ich)

        else:
            raise ValueError('Selection not implemented.')

    def iter_ph_times(self, ph_sel=Ph_sel('all')):
        """Iterator that returns the arrays of timestamps in `.ph_times_m`.

        Arguments:
            ph_sel (Ph_sel object): object defining the photon selection.
                See :class:`fretbursts.ph_sel.Ph_sel` for details.
        """
        for ich in xrange(self.nch):
            yield self.get_ph_times(ich, ph_sel=ph_sel)

    def get_ph_times(self, ich=0, ph_sel=Ph_sel('all')):
        """Returns the timestamps array for channel `ich`.

        This method always returns in-memory arrays, even when ph_times_m
        is a disk-backed list of arrays.

        Arguments:
            ph_sel (Ph_sel object): object defining the photon selection.
                See :class:`fretbursts.ph_sel.Ph_sel` for details.
        """
        ph = self.ph_times_m[ich]

        # If not a list is an on-disk array, and needs to be loaded (with [:])
        if type(self.ph_times_m) is not list:
            ph = ph[:]

        return ph[self.get_ph_mask(ich, ph_sel=ph_sel)]

    def _get_ph_mask_single(self, ich, mask_name, negate=False):
        """Get the bool array `mask_name` for channel `ich`.
        If the internal "bool array" is a scalar return a slice (full or empty)
        """
        mask = np.asarray(getattr(self, mask_name)[ich])
        if negate:
            mask = np.logical_not(mask)
        if len(mask.shape) == 0:
            # If mask is a boolean scalar, select all or nothing
            mask = slice(None) if mask else slice(0)
        return mask

    def get_A_em(self, ich=0):
        """Returns a mask to select photons detected in the acceptor ch."""
        return self._get_ph_mask_single(ich, 'A_em')

    def get_D_em(self, ich=0):
        """Returns a mask to select photons detected in the donor ch."""
        return self._get_ph_mask_single(ich, 'A_em', negate=True)

    def get_A_ex(self, ich=0):
        """Returns a mask to select photons in acceptor-excitation periods."""
        return self._get_ph_mask_single(ich, 'A_ex')

    def get_D_ex(self, ich=0):
        """Returns a mask to select photons in donor-excitation periods."""
        if self.ALEX:
            return self._get_ph_mask_single(ich, 'D_ex')
        else:
            return slice(None)

    def get_D_em_D_ex(self, ich=0):
        """Returns a mask of donor photons during donor-excitation."""
        if self.ALEX:
            return self.get_D_em(ich)*self.get_D_ex(ich)
        else:
            return self.get_D_em(ich)

    def get_A_em_D_ex(self, ich=0):
        """Returns a mask of acceptor photons during donor-excitation."""
        if self.ALEX:
            return self.get_A_em(ich)*self.get_D_ex(ich)
        else:
            return self.get_A_em(ich)

    def iter_ph_times_period(self, ich=0, ph_sel=Ph_sel('all')):
        """Iterate through arrays of ph timestamps in each background period.
        """
        mask = self.get_ph_mask(ich=ich, ph_sel=ph_sel)
        for period in range(self.nperiods):
            yield self.get_ph_times_period(period, ich=ich, mask=mask)

    def get_ph_times_period(self, period, ich=0, ph_sel=Ph_sel('all'),
                            mask=None):
        """Return the array of ph_times in `period`, `ich` and `ph_sel`.
        """
        istart, iend = self.Lim[ich][period]
        period_slice = slice(istart, iend + 1)

        ph_times = self.get_ph_times(ich=ich)
        if mask is None:
            mask = self.get_ph_mask(ich=ich, ph_sel=ph_sel)

        if type(mask) is slice and mask == slice(None):
            ph_times_period = ph_times[period_slice]
        else:
            ph_times_period = ph_times[period_slice][mask[period_slice]]
        return ph_times_period


    def slice_ph(self, time_s1=0, time_s2=None, s='slice'):
        """Return a new Data object with ph in [`time_s1`,`time_s2`] (seconds)
        """
        if time_s2 is None: time_s2 = self.time_max()
        if time_s2 >= self.time_max() and time_s1 <= 0:
                return self.copy()

        t1_clk, t2_clk = time_s1/self.clk_p, time_s2/self.clk_p
        assert np.array([t1_clk < ph.max() for ph in self.ph_times_m]).all()

        masks = [(ph >= t1_clk)*(ph <= t2_clk) for ph in self.iter_ph_times()]

        new_d = Data(**self)
        for name in self.ph_fields:
            if name in self:
                #if name == 'A_em':
                #    raise ValueError
                new_d[name] = [a[mask] for a, mask in zip(self[name], masks)]
                setattr(new_d, name, new_d[name])
        new_d.delete_burst_data()

        # Shift timestamps to start from 0 to avoid problems with BG calc
        for ich in range(self.nch):
            ph_i = new_d.get_ph_times(ich)
            ph_i -= t1_clk
        new_d.s.append(s)
        return new_d

    def burst_sizes_ich(self, ich=0, gamma=1., gamma1=None, add_naa=False):
        """Return gamma corrected burst sizes for channel `ich`.

        The gamma corrected size is computed as::

            1) nd + na/gamma  (so th1 is the min. burst size for donly bursts)
            2) nd*gamma1 + na (so th1 is the min. burst size for high FRET
            bursts)

        If `gamma1` is not None, the definition (2) is used.
        If `d.ALEX` and `add_naa` are True, `naa` is added to the burst size.

        Returns
            Array of burst sizes for channel `ich`.
        """
        if gamma1 is not None:
            burst_size = 1.*self.nd[ich]*gamma1 + self.na[ich]
        else:
            burst_size = self.nd[ich] + 1.*self.na[ich]/gamma
        if self.ALEX and add_naa:
            burst_size += self.naa[ich]
        return burst_size

    def bursts_slice(self, N1=0, N2=-1):
        """Return new Data object with bursts between `N1` and `N2`
        `N1` and `N2` can be scalars or lists (one per ch).
        """
        if np.isscalar(N1): N1 = [N1]*self.nch
        if np.isscalar(N2): N2 = [N2]*self.nch
        assert len(N1) == len(N2) == self.nch
        d = Data(**self)
        d.add(mburst=[b[n1:n2, :] for b, n1, n2 in zip(d.mburst, N1, N2)])
        d.add(nt=[nt[n1:n2] for nt, n1, n2 in zip(d.nt, N1, N2)])
        d.add(nd=[nd[n1:n2] for nd, n1, n2 in zip(d.nd, N1, N2)])
        d.add(na=[na[n1:n2] for na, n1, n2 in zip(d.na, N1, N2)])
        if self.ALEX:
            d.add(naa=[aa[n1:n2] for aa, n1, n2 in zip(d.naa, N1, N2)])
        d.calc_fret() # recalc fret efficiency
        return d

    def collapse(self, update_gamma=True):
        """Returns an object with 1-ch data joining the multi-ch data.

        The argument `update_gamma` (bool, default True) allows to avoid
        recomputing gamma as the average of the original gamma. This flag
        should be always True. Set False only for testing/debugging.
        """
        dc = Data(**self)

        mburst = np.vstack(self.mburst)
        burst_start = b_start(mburst)
        sort_index = burst_start.argsort()

        ich_burst = [i*np.ones(nb) for i, nb in enumerate(self.num_bursts())]
        dc.add(ich_burst=np.hstack(ich_burst)[sort_index])

        for name in self.burst_fields:
            if name in self:
                # Concatenate arrays along axis = 0
                value = [ np.concatenate(self[name])[sort_index] ]
                dc.add(**{name: value})
        dc.add(nch=1)
        dc.add(chi_ch=1.)
        # NOTE: Updating gamma has the side effect of recomputing E
        #       (and S if ALEX). We need to update gamma because, in general,
        #       gamma can be an array with a value for each ch.
        if update_gamma:
            dc.update_gamma(np.mean(self.get_gamma_array()))
        return dc

    ##
    # Utility methods
    #
    def get_params(self):
        """Returns a plain dict containing only parameters and no arrays.
        This can be used as a summary of data analysis parameters.
        Additional keys `name' and `Names` are added with values
        from `.name()` and `.Name()`.
        """
        p_names = ['fname', 'clk_p', 'nch', 'ph_sel', 'L', 'm', 'F', 'P',
                   'leakage', 'dir_ex', 'gamma', 'bg_time_s', 'nperiods',
                   'rate_dd', 'rate_ad', 'rate_aa', 'rate_m', 'T', 'rate_th',
                   'bg_corrected', 'leakage_corrected', 'dir_ex_corrected',
                   'dithering', 'chi_ch', 's', 'ALEX']
        p_dict = dict(self)
        for name in p_dict.keys():
            if name not in p_names:
                p_dict.pop(name)
        p_dict.update(name=self.name(), Name=self.Name())
        return p_dict

    def expand(self, ich=0, alex_naa=False, width=False):
        """Return per-burst D and A sizes (nd, na) and background (bg_d, bg_a).

        This method returns for each bursts the corrected signal counts and
        background counts in donor and acceptor channels. Optionally, the
        burst width is also returned.

        Arguments:
            ich (int): channel for the bursts (can be not 0 only in multi-spot)
            alex_naa (bool): if True and self.ALEX, returns burst sizes and
                background also for acceptor photons during accept. excitation
            width (bool): whether return the burst duration (in seconds).

        Returns:
            List of arrays: nd, na, donor bg, acceptor bg.
            If `alex_naa` is True returns: nd, na, naa, bg_d, bg_a, bg_aa.
            If `width` is True returns the bursts duration (in sec.) as last
            element.
        """
        period = self.bp[ich]
        w = b_width(self.mburst[ich])*self.clk_p
        bg_a = self.bg_ad[ich][period]*w
        bg_d = self.bg_dd[ich][period]*w
        res = [self.nd[ich], self.na[ich]]
        if self.ALEX and alex_naa:
            bg_aa = self.bg_aa[ich][period]*w
            res.extend([self.naa[ich], bg_d, bg_a, bg_aa])
        else:
            res.extend([bg_d, bg_a])
        if width:
            res.append(w)
        return res

    def time_max(self):
        """Return the measurement time (last photon) in seconds."""
        if 'ph_times_m' in self:
            return max([t[-1]*self.clk_p for t in self.ph_times_m])
        elif 'ph_times_t' in self:
            return self.ph_times_t.max()*self.clk_p
        elif 'mburst' in self:
            return max([b_end(mb)[-1]*self.clk_p for mb in self.mburst])
        else:
            raise ValueError("No timestamps or bursts found.")

    def num_bu(self):
        """Shortcut for :meth:`num_bursts`."""
        return self.num_bursts()

    def num_bursts(self):
        """Return an array with the number of bursts in each channel."""
        return np.r_[[mb.shape[0] for mb in self.mburst]]

    def ph_select(self):
        """Return masks of ph inside bursts for all the channels."""
        return mch_ph_select(self.ph_times_m, self.mburst)

    def ph_in_bursts(self, ich=0, ph_sel=Ph_sel('all')):
        """Return photons inside bursts.

        Returns
            Array photon timestamps in channel `ich` and photon
            selection `ph_sel` that are inside bursts.
        """
        ph_all = self.get_ph_times(ich=ich)
        bursts_mask = ph_select(ph_all, self.mburst[ich])
        if ph_sel == Ph_sel('all'):
            return ph_all[bursts_mask]
        else:
            ph_sel_mask = self.get_ph_mask(ich=ich, ph_sel=ph_sel)
            return ph_all[ph_sel_mask*bursts_mask]

    def calc_max_rate(self, m, ph_sel=Ph_sel('all')):
        """Compute the max m-photon rate reached in each burst.

        Arguments:
            m (int): number of timestamps to use to compute the rate
            ph_sel (Ph_sel object): object defining the photon selection.
                See :class:`fretbursts.ph_sel.Ph_sel` for details.
        """
        if ph_sel == Ph_sel('all'):
            Max_Rate = [b_rate_max(ph=ph, m=m, mburst=mb)
                    for ph, mb in zip(self.iter_ph_times(), self.mburst)]
        else:
            Max_Rate = [b_rate_max(ph=ph, m=m, mburst=mb, mask=mask)
                    for ph, mask, mb in zip(self.iter_ph_times(),
                                            self.iter_ph_masks(ph_sel=ph_sel),
                                            self.mburst)]

        Max_Rate = [mr/self.clk_p - bg[bp] for bp, bg, mr in
                    zip(self.bp, self.bg_from(ph_sel), Max_Rate)]
        self.add(max_rate=Max_Rate)

    def delete_burst_data(self):
        """Erase all the burst data"""
        for k in self.burst_fields + ['fuse', 'lsb']:
            if k in self:
                self.delete(k)
    ##
    # Background analysis methods
    #
    def calc_bg_cache(self, fun, time_s=60, tail_min_us=500, F_bg=2,
                      recompute=False):
        """Compute time-dependent background rates for all the channels.

        This version is the cached version of :meth:`calc_bg`.
        This method tries to load the background data from the HDF5 file in
        self.bg_data_file. If a saved background data is not found, it computes
        the background and stores the data to the HDF5 file.

        The arguments are the same as :meth:`calc_bg` with the only addition
        of `recompute` (bool) to force a background recomputation even if
        a cached version is found.

        Form more details on the other arguments see :meth:`calc_bg`.
        """
        bg_cache.calc_bg_cache(self, fun, time_s=time_s,
                               tail_min_us=tail_min_us, F_bg=F_bg,
                               recompute=recompute)

    def _get_auto_bg_th_arrays(self, F_bg=2, tail_min_us0=250):
        """Return a dict of threshold values for background estimation.

        The keys are the ph selections in self.ph_streams and the values
        are 1-D arrays of size nch.
        """
        Th_us = {}
        for ph_sel in self.ph_streams:
            th_us = np.zeros(self.nch)
            for ich, ph in enumerate(self.iter_ph_times(ph_sel=ph_sel)):
                if ph.size > 0:
                    bg_rate, _ = bg.exp_fit(ph, tail_min_us=tail_min_us0)
                    th_us[ich] = 1e6*F_bg/bg_rate
            Th_us[ph_sel] = th_us
        # Save the input used to generate Th_us
        self.add(bg_auto_th_us0=tail_min_us0, bg_auto_F_bg=F_bg)
        return Th_us

    def _get_bg_th_arrays(self, tail_min_us):
        """Return a dict of threshold values for background estimation.

        The keys are the ph selections in self.ph_streams and the values
        are 1-D arrays of size nch.
        """
        n_streams = len(self.ph_streams)

        if np.size(tail_min_us) == 1:
            tail_min_us = np.repeat(tail_min_us, n_streams)
        elif np.size(tail_min_us) == n_streams:
            tail_min_us = np.asarray(tail_min_us)
        elif np.size(tail_min_us) != n_streams:
            raise ValueError, 'Wrong tail_min_us length (%d).' %\
                    len(tail_min_us)
        Th_us = {}
        for i, key in enumerate(self.ph_streams):
            Th_us[key] = np.ones(self.nch)*tail_min_us[i]
        # Save the input used to generate Th_us
        self.add(bg_th_us_user=tail_min_us)
        return Th_us

    def _clean_bg_data(self):
        """Remove background fields specific of only one fit type.

        Computing background with manual or 'auto' threshold resut in
        different sets of attributes being saved. This method removes these
        attributes and should be called before recomputing the background
        to avoid having old stale attributes of a previous background fit.
        """
        # Attributes specific of manual or 'auto' bg fit
        field_list = ['bg_auto_th_us0', 'bg_auto_F_bg', 'bg_th_us_user']
        for field in field_list:
            if field in self:
                self.delete(field)

    def _get_num_periods(self, time_s):
        """Return the number of periods using `time_s` as period duration.
        """
        t_max_mch = np.array([ph[-1] for ph in self.iter_ph_times()])
        # Take the ceil to have at least 1 periods
        # Take the min to avoid having ch with 0 photons in the last period
        nperiods = np.ceil(t_max_mch*self.clk_p/time_s).min().astype('int32')
        # Discard last period if shorter than 0.1s
        if nperiods > 1:
            t_last_period = [t_max_i*self.clk_p - (nperiods - 1)*time_s
                                for t_max_i in t_max_mch]
            if np.min(t_last_period) < 0.1:
                nperiods -= 1
        return int(nperiods)

    def calc_bg(self, fun, time_s=60, tail_min_us=500, F_bg=2,
                error_metrics='KS'):
        """Compute time-dependent background rates for all the channels.

        Compute background rates for donor, acceptor and both detectors.
        The rates are computed every `time_s` seconds, allowing to
        track possible variations during the measurement.

        Arguments:
            fun (function): function for background estimation (example
                `bg.exp_fit`)
            time_s (float, seconds): compute background each time_s seconds
            tail_min_us (float, tuple or string): min threshold in us for
                photon waiting times to use in background estimation.
                If float is the same threshold for 'all', DD, AD and AA photons
                and for all the channels.
                If a 3 or 4 element tuple, each value is used for 'all', DD, AD
                or AA photons, same value for all the channels.
                If 'auto', the threshold is computed for each stream ('all',
                DD, DA, AA) and for each channel as `bg_F * rate_ml0`.
                `rate_ml0` is an initial estimation of the rate performed using
                :func:`bg.exp_fit` and a fixed threshold (default 250us).
            F_bg (float): when `tail_min_us` is 'auto', is the factor by which
                the initial background estimation if multiplied to compute the
                threshold.
            error_metrics (string): Specifies the error metric to use.
                See :func:`fretbursts.background.exp_fit` for more details.

        The background estimation functions are defined in the module
        `background` (conventionally imported as `bg`).

        Example:
            Compute background with `bg.exp_fit`, every 20s, with a
            threshold of 200us for photon waiting times::

               d.calc_bg(bg.exp_fit, time_s=20, tail_min_us=200)

        Returns:
            None, all the results are saved in the object itself.
        """
        pprint(" - Calculating BG rates ... ")
        self._clean_bg_data()

        if tail_min_us == 'auto':
            bg_auto_th = True
            Th_us = self._get_auto_bg_th_arrays(F_bg=F_bg)
        else:
            bg_auto_th = False
            Th_us = self._get_bg_th_arrays(tail_min_us)

        kwargs = dict(clk_p=self.clk_p, error_metrics=error_metrics)
        nperiods = self._get_num_periods(time_s)
        bg_time_clk = time_s/self.clk_p

        BG, BG_dd, BG_ad, BG_da, BG_aa, Lim, Ph_p = [], [], [], [], [], [], []
        rate_m, rate_dd, rate_ad, rate_da, rate_aa = [], [], [], [], []
        BG_err, BG_dd_err, BG_ad_err, BG_da_err, BG_aa_err = [], [], [], [], []
        for ich, ph_ch in enumerate(self.iter_ph_times()):
            th_us_ch_all = Th_us[Ph_sel('all')][ich]
            th_us_ch_dd = Th_us[Ph_sel(Dex='Dem')][ich]
            th_us_ch_ad = Th_us[Ph_sel(Dex='Aem')][ich]
            if self.ALEX:
                th_us_ch_da = Th_us[Ph_sel(Aex='Dem')][ich]
                th_us_ch_aa = Th_us[Ph_sel(Aex='Aem')][ich]

            dd_mask = self.get_ph_mask(ich, ph_sel=Ph_sel(Dex='Dem'))
            ad_mask = self.get_ph_mask(ich, ph_sel=Ph_sel(Dex='Aem'))
            if self.ALEX:
                da_mask = self.get_ph_mask(ich, ph_sel=Ph_sel(Aex='Dem'))
                aa_mask = self.get_ph_mask(ich, ph_sel=Ph_sel(Aex='Aem'))

            lim, ph_p = [], []
            bg, bg_dd, bg_ad, bg_da, bg_aa = [zeros(nperiods) for _ in range(5)]
            zeros_list = [zeros(nperiods) for _ in range(5)]
            bg_err, bg_dd_err, bg_ad_err, bg_da_err, bg_aa_err = zeros_list
            for ip in xrange(nperiods):
                i0 = 0 if ip == 0 else i1           # pylint: disable=E0601
                i1 = (ph_ch < (ip+1)*bg_time_clk).sum()
                lim.append((i0, i1-1))
                ph_p.append((ph_ch[i0], ph_ch[i1-1]))

                ph_i = ph_ch[i0:i1]
                bg[ip], bg_err[ip] = fun(ph_i, tail_min_us=th_us_ch_all,
                                         **kwargs)

                # This supports cases of D-only or A-only timestamps
                # where self.A_em[ich] is a bool and not a bool-array
                # In this case, either `dd_mask` or `ad_mask` is
                # slice(None) (all-elements selection)
                if type(dd_mask) is slice and dd_mask == slice(None):
                    bg_dd[ip], bg_dd_err[ip] = bg[ip], bg_err[ip]
                    continue
                if type(ad_mask) is slice and ad_mask == slice(None):
                    bg_ad[ip], bg_ad_err[ip] = bg[ip], bg_err[ip]
                    continue

                dd_mask_i = dd_mask[i0:i1]
                if dd_mask_i.any():
                    bg_dd[ip], bg_dd_err[ip] = fun(ph_i[dd_mask_i],
                                       tail_min_us=th_us_ch_dd, **kwargs)

                ad_mask_i = ad_mask[i0:i1]
                if ad_mask_i.any():
                    bg_ad[ip], bg_ad_err[ip] = fun(ph_i[ad_mask_i],
                                       tail_min_us=th_us_ch_ad, **kwargs)

                if self.ALEX and aa_mask.any():
                    da_mask_i = da_mask[i0:i1]
                    bg_da[ip], bg_da_err[ip] = fun(ph_i[da_mask_i],
                                       tail_min_us=th_us_ch_da, **kwargs)
                    aa_mask_i = aa_mask[i0:i1]
                    bg_aa[ip], bg_aa_err[ip] = fun(ph_i[aa_mask_i],
                                       tail_min_us=th_us_ch_aa, **kwargs)

            Lim.append(lim);     Ph_p.append(ph_p)
            BG.append(bg);       BG_err.append(bg_err)
            BG_dd.append(bg_dd); BG_dd_err.append(bg_dd_err)
            BG_ad.append(bg_ad); BG_ad_err.append(bg_ad_err)
            BG_da.append(bg_da); BG_da_err.append(bg_da_err)
            BG_aa.append(bg_aa); BG_aa_err.append(bg_aa_err)
            rate_m.append(bg.mean())
            rate_dd.append(bg_dd.mean())
            rate_ad.append(bg_ad.mean())
            if self.ALEX:
                rate_da.append(bg_da.mean())
                rate_aa.append(bg_aa.mean())
        self.add(bg=BG, bg_dd=BG_dd, bg_ad=BG_ad, bg_da=BG_da, bg_aa=BG_aa,
                 bg_err=BG_err, bg_dd_err=BG_dd_err, bg_ad_err=BG_ad_err,
                 bg_da_err=BG_da_err, bg_aa_err=BG_aa_err,
                 Lim=Lim, Ph_p=Ph_p, nperiods=nperiods,
                 bg_fun=fun, bg_fun_name=fun.__name__,
                 bg_time_s=time_s, bg_ph_sel=Ph_sel('all'), rate_m=rate_m,
                 rate_dd=rate_dd, rate_ad=rate_ad,
                 rate_da=rate_da, rate_aa=rate_aa,
                 bg_th_us=Th_us, bg_auto_th=bg_auto_th)
        pprint("[DONE]\n")

    def recompute_bg_lim_ph_p(self, ph_sel=Ph_sel(Dex='Dem'), mute=False):
        """Recompute self.Lim and selp.Ph_p relative to ph selection `ph_sel`
        `ph_sel` is a Ph_sel object selecting the timestamps in which self.Lim
        and self.Ph_p are being computed.
        """
        if self.bg_ph_sel == ph_sel: return

        pprint(" - Recomputing background limits for %s ... " % \
                str(ph_sel), mute)
        bg_time_clk = self.bg_time_s/self.clk_p
        Lim, Ph_p = [], []
        for ph_ch, lim in zip(self.iter_ph_times(ph_sel), self.Lim):
            lim, ph_p = [], []
            for ip in xrange(self.nperiods):
                i0 = 0 if ip == 0 else i1  # pylint: disable=E0601
                i1 = (ph_ch < (ip+1)*bg_time_clk).sum()
                lim.append((i0, i1-1))
                ph_p.append((ph_ch[i0], ph_ch[i1-1]))
            Lim.append(lim)
            Ph_p.append(ph_p)
        self.add(Lim=Lim, Ph_p=Ph_p, bg_ph_sel=ph_sel)
        pprint("[DONE]\n", mute)

    ##
    # Burst analysis methods
    #
    def _calc_burst_period(self):
        """Compute for each burst the "period" `bp`.
        Periods are times intervals on which the BG is computed.
        """
        P = []
        for b, lim in zip(self.mburst, self.Lim):
            p = zeros(b.shape[0], dtype=np.int16)
            if b.size > 0:
                bis = b_istart(b)
                for i, (l0, l1) in enumerate(lim):
                    p[(bis >= l0)*(bis <= l1)] = i
            P.append(p)
        self.add(bp=P)

    def _param_as_mch_array(self, par):
        """Regardless of `par` size, return an arrays with size == nch.
        if `par` is scalar the arrays repeats the calar multiple times
        if `par is a list/array must be of length `nch`.
        """
        assert size(par) == 1 or size(par) == self.nch
        return np.repeat(par, self.nch) if size(par) == 1 else np.asarray(par)

    def bg_from(self, ph_sel):
        """Return the background rates for the specified photon selection.
        """
        ph_sel = self._check_ph_sel(ph_sel)

        if ph_sel == Ph_sel('all'):
            return self.bg

        BG = {Ph_sel(Dex='Dem'): self.bg_dd,
              Ph_sel(Dex='Aem'): self.bg_ad}
        if self.ALEX:
            bg_Dex = [bg_dd + bg_ad for bg_dd, bg_ad in
                            zip(self.bg_dd, self.bg_ad)]
            bg_Aex = [bg_da + bg_aa for bg_da, bg_aa in
                            zip(self.bg_da, self.bg_aa)]
            bg_Dem = [bg_dd + bg_da for bg_dd, bg_da in
                            zip(self.bg_dd, self.bg_da)]
            bg_Aem = [bg_ad + bg_aa for bg_ad, bg_aa in
                            zip(self.bg_ad, self.bg_aa)]
            bg_noDA = [bg_dd + bg_ad + bg_aa for bg_dd, bg_ad, bg_aa in
                            zip(self.bg_dd, self.bg_ad, self.bg_aa)]
            BG.update(**{Ph_sel(Aex='Aem'): self.bg_aa,
                         Ph_sel(Aex='Dem'): self.bg_da,
                         Ph_sel(Dex='DAem'): bg_Dex,
                         Ph_sel(Aex='DAem'): bg_Aex,
                         Ph_sel(Dex='Dem', Aex='Dem'): bg_Dem,
                         Ph_sel(Dex='Aem', Aex='Aem'): bg_Aem,
                         Ph_sel(Dex='DAem', Aex='Aem'): bg_noDA})
        if ph_sel not in BG:
            raise NotImplementedError('Photong selection not implemented.')
        return BG[ph_sel]

    def _calc_T(self, m, P, F=1., ph_sel=Ph_sel('all')):
        """If P is None use F, otherwise uses both P *and* F (F defaults to 1).
        """
        # Regardless of F and P sizes, FF and PP are arrays with size == nch
        FF = self._param_as_mch_array(F)
        PP = self._param_as_mch_array(P)
        if P is None:
            find_T = lambda m, Fi, Pi, bg: 1.*m/(bg*Fi) # NOTE: ignoring P_i
        else:
            if F != 1:
                print "WARNING: BS prob. th. with modified BG rate (F=%.1f)"%F
            find_T = lambda m, Fi, Pi, bg: find_optimal_T_bga(bg*Fi, m, 1-Pi)
        TT, T, rate_th = [], [], []
        bg_bs = self.bg_from(ph_sel)
        for bg_ch, F_ch, P_ch in zip(bg_bs, FF, PP):
            # All "T" are in seconds
            Tch = find_T(m, F_ch, P_ch, bg_ch)
            TT.append(Tch)
            T.append(Tch.mean())
            rate_th.append(np.mean(m/Tch))
        self.add(TT=TT, T=T, bg_bs=bg_bs, FF=FF, PP=PP, F=F, P=P,
                 rate_th=rate_th)

    def _burst_search_rate(self, m, L, min_rate_cps, ph_sel=Ph_sel('all'),
                           verbose=True, pure_python=False):
        """Compute burst search using a fixed minimum photon rate.

        Arguments:
            min_rate_cps (float or array): minimum photon rate for burst start
                if array if one value per channel.
        """
        bsearch = _get_bsearch_func(pure_python=pure_python)

        Min_rate_cps = self._param_as_mch_array(min_rate_cps)
        mburst = []
        T_clk = (1.*m/Min_rate_cps)/self.clk_p
        for ich, (ph, t_clk) in enumerate(zip(self.iter_ph_times(ph_sel),
                                          T_clk)):
            label = '%s CH%d' % (ph_sel, ich+1) if verbose else None
            mb = bsearch(ph, L, m, t_clk, label=label, verbose=verbose)
            mburst.append(mb)
        self.add(mburst=mburst, min_rate_cps=Min_rate_cps, T=T_clk*self.clk_p)

    def _burst_search_TT(self, m, L, ph_sel=Ph_sel('all'), verbose=True,
                         pure_python=False, mute=False):
        """Compute burst search with params `m`, `L` on ph selection `ph_sel`

        Requires the list of arrays `self.TT` with the max time-thresholds in
        the different burst periods for each channel (use `._calc_T()`).
        """
        bsearch = _get_bsearch_func(pure_python=pure_python)

        self.recompute_bg_lim_ph_p(ph_sel=ph_sel, mute=mute)
        MBurst = []
        label = ''
        for ich, (ph, T) in enumerate(zip(self.iter_ph_times(ph_sel),
                                          self.TT)):
            MB = []
            Tck = T/self.clk_p
            for ip, (l0, l1) in enumerate(self.Lim[ich]):
                if verbose:
                    label='%s CH%d-%d' % (ph_sel, ich+1, ip)
                mb = bsearch(ph[l0:l1+1], L, m, Tck[ip], label=label,
                        verbose=verbose)
                if mb.size > 0: # if we found at least one burst
                    mb[:, iistart] += l0
                    mb[:, iiend] += l0
                    MB.append(mb)
            if len(MB) > 0:
                MBurst.append(np.vstack(MB))
            else:
                MBurst.append(np.array([]))
        self.add(mburst=MBurst)
        if ph_sel != Ph_sel('all'):
            # Convert the burst data to be relative to ph_times_m.
            # Convert both Lim/Ph_p and mburst, as they are both needed
            # to compute `.bp`.
            self.recompute_bg_lim_ph_p(ph_sel=Ph_sel('all'), mute=mute)
            self._fix_mburst_from(ph_sel=ph_sel, mute=mute)

    def _fix_mburst_from(self, ph_sel, mute=False):
        """Convert burst data from any ph_sel to 'all' timestamps selection.
        """
        assert type(ph_sel) is Ph_sel and ph_sel != Ph_sel('all')
        pprint(' - Fixing  burst data to refer to ph_times_m ... ', mute)
        old_MBurst = [mb.copy() for mb in self.mburst]

        # Note that mburst is modified in-place
        it_ph_masks = self.iter_ph_masks(ph_sel=ph_sel)
        ph_sizes = [ph.size for ph in self.iter_ph_times()]
        for mburst, ph_size, mask in zip(self.mburst, ph_sizes, it_ph_masks):
            index = np.arange(ph_size, dtype=np.int32)
            mburst[:, iistart] = index[mask][mburst[:, iistart]]
            mburst[:, iiend] = index[mask][mburst[:, iiend]]
            mburst[:, inum_ph] = mburst[:, iiend] - mburst[:, iistart] + 1

        for mb, old_mb in zip(self.mburst, old_MBurst):
            assert (mb[:, iistart] >= old_mb[:, iistart]).all()
            assert (mb[:, iiend] >= old_mb[:, iiend]).all()
            assert (mb[:, inum_ph] >= old_mb[:, inum_ph]).all()
        pprint('[DONE]\n', mute)

    def burst_search_t(self, L=10, m=10, P=None, F=6., min_rate_cps=None,
            nofret=False, max_rate=False, dither=False, ph_sel=Ph_sel('all'),
            verbose=False, mute=False, pure_python=False):
        """Performs a burst search with specified parameters.

        This method performs a sliding-window burst search without
        binning the timestamps. The burst starts when the rate of `m`
        photons is above a minimum rate, and stops when the rate falls below
        the threshold.

        The minimum rate can be explicitly specified (`min_rate`) or computed
        as a function of the background rate (using `F` or `P`).

        Parameters:
            m (int): number of consecutive photons used to compute the
                photon rate.
            L (int): minimum number of photons in burst
            P (float): threshold for burst detection expressed as a
                probability that a detected bursts is not due to a Poisson
                background. If not None, `P` overrides `F`. Note that the
                background process is experimentally super-Poisson so this
                probability is not physically meaningful.
            F (float): if `P` is None: min.rate/bg.rate ratio for burst search
                else: `P` refers to a Poisson rate = bg_rate*F
            min_rate (float or list/array): min. rate in cps for burst start.
                If not None has the precedence over `P` and `F`.
                If non-scalar, contains one rate per each channel.
            ph_sel (Ph_sel object): object defining the photon selection
                used for burst search. Default: all photons.
                See :class:`fretbursts.ph_sel.Ph_sel` for details.
            pure_python (bool): if True, uses the pure python functions even
                when the optimized Cython functions are available.

        Note:
            when using `P` or `F` the background rates are needed, so
            `.calc_bg()` must be called before the burst search.

        Example:
            d.burst_search_t(L=10, m=10, F=6)

        Returns:
            None, all the results are saved in the object.
        """
        pprint(" - Performing burst search (verbose=%s) ..." % verbose, mute)
        # Erase any previous burst data
        self.delete_burst_data()

        if min_rate_cps is not None:
            self._burst_search_rate(m=m, L=L, min_rate_cps=min_rate_cps,
                                    ph_sel=ph_sel, verbose=verbose,
                                    pure_python=pure_python)
        else:
            # Compute TT
            self._calc_T(m=m, P=P, F=F, ph_sel=ph_sel)
            # Use TT and compute mburst
            self._burst_search_TT(L=L, m=m, ph_sel=ph_sel, verbose=verbose,
                                  pure_python=pure_python, mute=mute)
        pprint("[DONE]\n", mute)

        pprint(" - Calculating burst periods ...", mute)
        self._calc_burst_period()                       # writes bp
        pprint("[DONE]\n", mute)

        self.add(m=m, L=L, ph_sel=ph_sel)  # P and F are saved in _calc_T()

        # The correction flags are both set here and in calc_ph_num() so that
        # they are always consistent. Case 1: we perform only burst search
        # (with no call to calc_ph_num). Case 2: we re-call calc_ph_num()
        # without doing a new burst search
        self.add(bg_corrected=False, leakage_corrected=False,
                 dir_ex_corrected=False, dithering=False)

        if not nofret:
            pprint(" - Counting D and A ph and calculating FRET ... \n", mute)
            self.calc_fret(count_ph=True, corrections=True, dither=dither,
                           mute=mute, pure_python=pure_python)
            pprint("   [DONE Counting D/A]\n", mute)
        if max_rate:
            pprint(" - Computing max rates in burst ...", mute)
            self.calc_max_rate(m=m)
            pprint("[DONE]\n", mute)

    def calc_ph_num(self, alex_all=False, pure_python=False):
        """Computes number of D, A (and AA) photons in each burst.

        Arguments:
            alex_all (bool): if True and self.ALEX is True, computes also the
                donor channel photons during acceptor excitation (`nda`)
            pure_python (bool): if True, uses the pure python functions even
                when the optimized Cython functions are available.

        Returns:
            Saves `nd`, `na`, `nt` (and eventually `naa`, `nda`) in self.
            Returns None.
        """
        mch_count_ph_in_bursts = _get_mch_count_ph_in_bursts_func(pure_python)

        if not self.ALEX:
            nt = [b_size(b).astype(float) if b.size > 0 else np.array([])\
                        for b in self.mburst]
            A_em = [self.get_A_em(ich) for ich in xrange(self.nch)]
            if type(A_em[0]) is slice:
                # This to support the case of A-only or D-only data
                n0 = [np.zeros(mb.shape[0]) for mb in self.mburst]
                if A_em[0] == slice(None):
                    nd, na = n0, nt    # A-only case
                elif A_em[0] == slice(0):
                    nd, na = nt, n0    # D-only case
            else:
                # This is the usual case with photons in both D and A channel
                na = mch_count_ph_in_bursts(self.mburst, Mask=A_em)
                nd = [t - a for t, a in zip(nt, na)]
            assert (nt[0] == na[0] + nd[0]).all()
        if self.ALEX:
            Mask = [d_em*d_ex for d_em, d_ex in zip(self.D_em, self.D_ex)]
            nd = mch_count_ph_in_bursts(self.mburst, Mask)

            Mask = [a_em*d_ex for a_em, d_ex in zip(self.A_em, self.D_ex)]
            na = mch_count_ph_in_bursts(self.mburst, Mask)

            Mask = [a_em*a_ex for a_em, a_ex in zip(self.A_em, self.A_ex)]
            naa = mch_count_ph_in_bursts(self.mburst, Mask)
            self.add(naa=naa)

            if alex_all:
                Mask = [d_em*a_ex for d_em, a_ex in zip(self.D_em, self.A_ex)]
                nda = mch_count_ph_in_bursts(self.mburst, Mask)
                self.add(nda=nda)

            nt = [d+a+aa for d, a, aa in zip(nd, na, naa)]
            assert (nt[0] == na[0] + nd[0] + naa[0]).all()
        self.add(nd=nd, na=na, nt=nt,
                 bg_corrected=False, leakage_corrected=False,
                 dir_ex_corrected=False, dithering=False)

    def fuse_bursts(self, ms=0, process=True, mute=False):
        """Return a new :class:`Data` object with nearby bursts fused together.

        Arguments:
            ms (float): fuse all burst separated by less than `ms` millisecs.
                If < 0 no burst is fused. Note that with ms = 0, overlapping
                bursts are fused.
            process (bool): if True (default), reprocess the burst data in
                the new object applying corrections and computing FRET.
            mute (bool): if True suppress any printed output.

        """
        if ms < 0: return self
        mburst = mch_fuse_bursts(self.mburst, ms=ms, clk_p=self.clk_p)
        new_d = Data(**self)
        for k in ['E', 'S', 'nd', 'na', 'naa', 'nt', 'lsb', 'bp']:
            if k in new_d: new_d.delete(k)
        new_d.add(bg_corrected=False, leakage_corrected=False,
                  dir_ex_corrected=False, dithering=False)
        new_d.add(mburst=mburst, fuse=ms)
        if 'bg' in new_d: new_d._calc_burst_period()
        if process:
            pprint(" - Counting D and A ph and calculating FRET ... \n", mute)
            new_d.calc_fret(count_ph=True, corrections=True,
                            dither=self.dithering, mute=mute)
            pprint("   [DONE Counting D/A and FRET]\n", mute)
        return new_d

    ##
    # Corrections methods
    #
    def background_correction_t(self, relax_nt=False, mute=False):
        """Apply background correction to burst sizes (nd, na,...)
        """
        if self.bg_corrected: return -1
        pprint("   - Applying background correction.\n", mute)
        self.add(bg_corrected=True)
        for ich, mb in enumerate(self.mburst):
            if mb.size == 0: continue  # if no bursts skip this ch
            period = self.bp[ich]
            nd, na, bg_d, bg_a, width = self.expand(ich, width=True)
            nd -= bg_d
            na -= bg_a
            if relax_nt:
                # This does not guarantee that nt = nd + na
                self.nt[ich] -= self.bg[ich][period] * width
            else:
                self.nt[ich] = nd + na
            if self.ALEX:
                self.naa[ich] -= self.bg_aa[ich][period] * width
                if 'nda' in self:
                    self.nda[ich] -= self.bg_da[ich][period] * width
                self.nt[ich] += self.naa[ich]

    def leakage_correction(self, mute=False):
        """Apply leakage correction to burst sizes (nd, na,...)
        """
        if self.leakage_corrected: return -1
        pprint("   - Applying leakage correction.\n", mute)
        assert (size(self.leakage) == 1) or (size(self.leakage) == self.nch)
        Lk = self.get_leakage_array()
        for i in range(self.nch):
            if self.na[i].size == 0: continue  # if no bursts skip this ch
            self.na[i] -= self.nd[i]*Lk[i]
            self.nt[i] = self.nd[i] + self.na[i]
            if self.ALEX: self.nt[i] += self.naa[i]
        self.add(leakage_corrected=True)

    def direct_excitation_correction(self, mute=False):
        """Apply direct excitation correction to bursts (ALEX-only).

        The applied correction is: na -= naa*dir_ex
        """
        if self.dir_ex_corrected: return -1
        pprint("   - Applying direct excitation correction.\n", mute)
        for i in range(self.nch):
            if self.na[i].size == 0: continue  # if no bursts skip this ch
            self.na[i] -= self.naa[i]*self.dir_ex
            self.nt[i] = self.nd[i] + self.na[i]
            if self.ALEX: self.nt[i] += self.naa[i]
        self.add(dir_ex_corrected=True)

    def dither(self, lsb=2, mute=False):
        """Add dithering (uniform random noise) to burst sizes (nd, na,...).
        `lsb` is the amplitude of dithering (if lsb=2 then noise is in -1..1).
        """
        if self.dithering: return -1
        pprint("   - Applying burst-size dithering.\n", mute)
        self.add(dithering=True)
        for nd, na in zip(self.nd, self.na):
            nd += lsb*(np.random.rand(nd.size)-0.5)
            na += lsb*(np.random.rand(na.size)-0.5)
        if self.ALEX:
            for naa in self.naa:
                naa += lsb*(np.random.rand(naa.size)-0.5)
            if 'nda' in self:
                for nda in self.nda:
                    nda += lsb*(np.random.rand(nda.size)-0.5)
        self.add(lsb=lsb)

    def calc_chi_ch(self):
        """Calculate the gamma correction prefactor factor `chi_ch` (array).
        `chi_ch` is a ch-dependent prefactor for gamma used to correct
        the dispersion of fitted E peaks (`E_fit`).
        `chi_ch` is returned, to apply the correction use `update_chi_ch()`.
        See also notebook: Gamma corrections.
        """
        if 'E_fit' not in self:
            print "ERROR: E_fit values not found. Call a `.fit_E_*` first."
            return
        EE = self.E_fit.mean()  # Mean E value among the CH
        chi_ch = (1/EE - 1)/(1/self.E_fit - 1)
        return chi_ch

    def corrections(self, mute=False):
        """Apply corrections on burst-counts: nd, na, nda, naa.

        The corrections are: background, leakage (or bleed-through) and
        direct excitation (dir_ex).
        """
        self.background_correction_t(mute=mute)
        self.leakage_correction(mute=mute)
        if 'dir_ex' in self and self.ALEX:
            self.direct_excitation_correction(mute=mute)

    def _update_corrections(self):
        """Recompute corrections whose flag is True.

        Check the flags .bg_corrected, .leakage_corrected, .dir_ex_corrected,
        dithering and recompute the corresponding correction if the flag
        is True (i.e. if the correction was already applied).

        Allows to recompute only the corrections the are already applied.
        """
        if 'mburst' not in self: return  # no burst search performed yet

        old_bg_corrected = self.bg_corrected
        old_leakage_corrected = self.leakage_corrected
        old_dir_ex_corrected = self.dir_ex_corrected
        old_dithering = self.dithering
        self.calc_ph_num()       # recompute uncorrected na, nd, nda, naa
        if old_bg_corrected:
            self.background_correction_t()
        if old_leakage_corrected:
            self.leakage_correction()
        if old_dir_ex_corrected:
            self.direct_excitation_correction()
        if old_dithering:
            self.dither(self.lsb)
        # Recompute E and S with no corrections (because already applied)
        self.calc_fret(count_ph=False, corrections=False)

    def update_bt(self, BT):
        """Deprecated. Use .update_leakage() instead.
        """
        print 'WARNING: The method .update_bt() is deprecated. '
        print 'Use .update_leakage() instead.'
        self.update_leakage(BT)

    def update_leakage(self, leakage):
        """Apply/update leakage (or bleed-through) correction.
        """
        self.add(leakage=leakage, leakage_corrected=True)
        self._update_corrections()

    def update_dir_ex(self, dir_ex):
        """Apply/update direct excitation correction with value `dir_ex`.
        """
        self.add(dir_ex=dir_ex, dir_ex_corrected=True)
        self._update_corrections()

    def update_chi_ch(self, chi_ch):
        """Change the `chi_ch` value and recompute FRET."""
        self.add(chi_ch=chi_ch)
        self.calc_fret(corrections=False)

    def update_gamma(self, gamma):
        """Change the `gamma` value and recompute FRET."""
        self.add(gamma=gamma)
        self.calc_fret(corrections=False)

    def get_gamma_array(self):
        """Get the array of gamma values (one per ch).
        Use this function to obtain an array of gamma values
        regardless of the actual `gamma` (that can be scalar).
        Gamma values are multiplied by `chi_ch`.
        """
        assert (size(self.gamma) == 1) or (size(self.gamma) == self.nch)
        gamma = self.gamma
        G = np.r_[[gamma]*self.nch] if np.size(gamma) == 1 else gamma
        G *= self.chi_ch
        return G

    def get_leakage_array(self):
        """Get the array of leakage coefficients, one per ch.

        This function always returns an array of leakage values regardless
        whether `leakage` is scalar or array.

        Each element of the returned array is multiplied by `chi_ch`.
        """
        lk_size = size(self.leakage)
        assert (lk_size == 1) or (lk_size == self.nch)
        Lk = np.r_[[self.leakage]*self.nch] if lk_size == 1 else self.leakage
        Lk *= self.chi_ch
        return Lk

    def calc_sbr(self, ph_sel=Ph_sel('all'), gamma=1.):
        """Return Signal-to-Background Ratio (SBR) for each burst.

        Arguments:
            ph_sel (Ph_sel object): object defining the photon selection
                for which to compute the sbr. Changes the photons used for
                burst size and the corresponding background rate. Valid values
                here are Ph_sel('all'), Ph_sel(Dex='Dem'), Ph_sel(Dex='Aem').
                See :class:`fretbursts.ph_sel.Ph_sel` for details.
            gamma (float): gamma value used to compute corrected burst size
                in the case `ph_sel` is Ph_sel('all'). Ignored otherwise.
        Returns:
            A list of arrays (one per channel) with one value per burst.
            The list is also saved in `sbr` attribute.
        """
        sbr = []
        for ich, mb in enumerate(self.mburst):
            if mb.size == 0:
                sbr.append([])
                continue  # if no bursts skip this ch
            nd, na, bg_d, bg_a = self.expand(ich)
            nt = self.burst_sizes_ich(ich=ich, gamma=gamma)

            signal = {Ph_sel('all'): nt,
                      Ph_sel(Dex='Dem'): nd, Ph_sel(Dex='Aem'): na}

            background = {Ph_sel('all'): bg_d + bg_a,
                          Ph_sel(Dex='Dem'): bg_d, Ph_sel(Dex='Aem'): bg_a}

            sbr.append(1.*signal[ph_sel]/background[ph_sel])
        self.add(sbr=sbr)
        return sbr

    ##
    # FRET and stoichiometry methods
    #
    def calc_fret(self, count_ph=False, corrections=True, dither=False,
                  mute=False, pure_python=False):
        """Compute FRET (and stoichiometry if ALEX) for each burst.

        This is an high-level functions that can be run after burst search.
        By default, it will count Donor and Acceptor photons, perform
        corrections (background, bleed-through), and compute gamma-corrected
        FRET efficiencies (and stoichiometry if ALEX).

        Arguments:
            count_ph (bool): if True (default), calls :meth:`calc_ph_num` to
                counts Donor and Acceptor photons in each bursts
            corrections (bool):  if True (default), applies background and
                bleed-through correction to burst data
            dither (bool): whether to apply dithering to burst size.
                Default False.
            mute (bool): whether to mute all the printed output. Default False.
            pure_python (bool): if True, uses the pure python functions even
                when the optimized Cython functions are available.

        Returns:
            None, all the results are saved in the object.
        """
        if count_ph:
            self.calc_ph_num(pure_python=pure_python, alex_all=True)
        if dither:
            self.dither(mute=mute)
        if corrections:
            self.corrections(mute=mute)
        self.calculate_fret_eff()
        if self.ALEX:
            self.calculate_stoich()
            self.calc_alex_hist()
        for fitter in ['E_fitter', 'S_fitter']:
            if fitter in self:
                self.delete(fitter)

    def calculate_fret_eff(self):
        """Compute FRET efficiency (`E`) for each burst."""
        G = self.get_gamma_array()
        E = [1.*na/(g*nd+na) for nd, na, g in zip(self.nd, self.na, G)]
        self.add(E=E)

    def calculate_stoich(self):
        """Compute "stochiometry" (the `S` parameter) for each burst."""
        G = self.get_gamma_array()
        S = [1.0*(g*d+a)/(g*d+a+aa) for d, a, aa, g in
                zip(self.nd, self.na, self.naa, G)]
        self.add(S=S)

    def calc_alex_hist(self, bin_step=0.05):
        """Compute the ALEX histogram with given bin width `bin_step`"""
        self.add(bin_step=bin_step)
        ES_hist_tot = [ES_histog(E, S, bin_step) for E, S in
                                                        zip(self.E, self.S)]
        E_bins, S_bins = ES_hist_tot[0][1], ES_hist_tot[0][2]
        ES_hist = [h[0] for h in ES_hist_tot]
        E_ax = E_bins[:-1] + 0.5*bin_step
        S_ax = S_bins[:-1] + 0.5*bin_step
        self.add(ES_hist=ES_hist, E_bins=E_bins, S_bins=S_bins,
                 E_ax=E_ax, S_ax=S_ax)

    ##
    # Information methods
    #
    def status(self, add="", noname=False):
        """Return a string with burst search, corrections and selection info.
        """
        name = "" if noname else self.name()
        s = name
        if 'L' in self: # burst search has been done
            if 'min_rate_cps' in self:
                s += " BS_%s L%d m%d MR%d" % (self.ph_sel, self.L, self.m,
                                              np.mean(self.min_rate_cps*1e-3))
            else:
                s += " BS_%s L%d m%d P%s F%.1f" % \
                        (self.ph_sel, self.L, self.m, self.P, np.mean(self.F))
        if 'gamma' in self: s += " G%.3f" % np.mean(self.gamma)
        if 'bg_fun' in self: s += " BG%s" % self.bg_fun.__name__[:-4]
        if 'bg_time_s' in self: s += "-%ds" % self.bg_time_s
        if 'fuse' in self: s += " Fuse%.1fms" % self.fuse
        if 'bg_corrected' in self and self.bg_corrected:
            s += " bg"
        if 'leakage_corrected' in self and self.leakage_corrected:
            s += " Lk%.3f" % np.mean(self.leakage)
        if 'dir_ex_corrected' in self and self.dir_ex_corrected:
            s += " dir%.1f" % (self.dir_ex*100)
        if 'dithering' in self and self.dithering:
            s += " Dith%d" % self.lsb
        if 's' in self: s += ' '.join(self.s)
        return s + add

    def name(self):
        """Return short filename (last subfolder + fname with no extension)"""
        name = basename = os.path.splitext(os.path.basename(self.fname))[0]
        last_dir = os.path.basename(os.path.dirname(self.fname))
        if last_dir is not '':
            name = '_'.join([last_dir, basename])
        return name


    def Name(self, add=""):
        """Return short filename + status information."""
        n = self.status(add=add)
        return n

    def __repr__(self):
        return self.status()

    def stats(self, string=False):
        """Print common statistics (BG rates, #bursts, mean size, ...)"""
        s = print_burst_stats(self)
        if string:
            return s
        else:
            print s

    ##
    # FRET fitting methods
    #
    def fit_E_m(self, E1=-1, E2=2, weights='size', gamma=1.):
        """Fit E in each channel with the mean using bursts in [E1,E2] range.

        Note:
            This two fitting are equivalent (but the first is much faster)::

                fit_E_m(weights='size')
                fit_E_minimize(kind='E_size', weights='sqrt')

            However `fit_E_minimize()` does not provide a model curve.
        """
        Mask = Sel_mask(self, select_bursts.E, E1=E1, E2=E2)

        fit_res, fit_model_F = zeros((self.nch, 2)), zeros(self.nch)
        for ich, (nd, na, E, mask) in enumerate(
                                        zip(self.nd, self.na, self.E, Mask)):
            w = fret_fit.get_weights(nd[mask], na[mask], weights, gamma)
            # Compute weighted mean
            fit_res[ich, 0] = np.dot(w, E[mask])/w.sum()
            # Compute weighted variance
            fit_res[ich, 1] = np.sqrt(
                    np.dot(w, (E[mask] - fit_res[ich,0])**2)/w.sum())
            fit_model_F[ich] = 1.*mask.sum()/mask.size

        fit_model = lambda x, p: SS.norm.pdf(x, p[0], p[1])
        self.add(fit_E_res=fit_res, fit_E_name='Moments',
                E_fit=fit_res[:,0], fit_E_curve=True, fit_E_E1=E1, fit_E_E2=E2,
                fit_E_model=fit_model, fit_E_model_F=fit_model_F)
        self.fit_E_calc_variance()
        return self.E_fit

    def fit_E_ML_poiss(self, E1=-1, E2=2, method=1, **kwargs):
        """ML fit for E modeling size ~ Poisson, using bursts in [E1,E2] range.
        """
        assert method in [1, 2, 3]
        fit_fun = {1: fret_fit.fit_E_poisson_na, 2: fret_fit.fit_E_poisson_nt,
                   3: fret_fit.fit_E_poisson_nd}
        Mask = Sel_mask(self, select_bursts.E, E1=E1, E2=E2)
        fit_res = zeros(self.nch)
        for ich, mask in zip(xrange(self.nch), Mask):
            nd, na, bg_d, bg_a = self.expand(ich)
            bg_x = bg_d if method == 3 else bg_a
            fit_res[ich] = fit_fun[method](nd[mask], na[mask],
                    bg_x[mask], **kwargs)
        self.add(fit_E_res=fit_res, fit_E_name='MLE: na ~ Poisson',
                E_fit=fit_res, fit_E_curve=False, fit_E_E1=E1, fit_E_E2=E2)
        self.fit_E_calc_variance()
        return self.E_fit

    def fit_E_ML_binom(self, E1=-1, E2=2, **kwargs):
        """ML fit for E modeling na ~ Binomial, using bursts in [E1,E2] range.
        """
        Mask = Sel_mask(self, select_bursts.E, E1=E1, E2=E2)
        fit_res = np.array([fret_fit.fit_E_binom(_d[mask], _a[mask], **kwargs)
                for _d, _a, mask in zip(self.nd, self.na, Mask)])
        self.add(fit_E_res=fit_res, fit_E_name='MLE: na ~ Binomial',
                E_fit=fit_res, fit_E_curve=False, fit_E_E1=E1, fit_E_E2=E2)
        self.fit_E_calc_variance()
        return self.E_fit

    def fit_E_minimize(self, kind='slope', E1=-1, E2=2, **kwargs):
        """Fit E using method `kind` ('slope' or 'E_size') and bursts in [E1,E2]
        If `kind` is 'slope' the fit function is fret_fit.fit_E_slope()
        If `kind` is 'E_size' the fit function is fret_fit.fit_E_E_size()
        Additional arguments in `kwargs` are passed to the fit function.
        """
        assert kind in ['slope', 'E_size']
        # Build a dictionary fun_d so we'll call the function fun_d[kind]
        fun_d = dict(slope=fret_fit.fit_E_slope, E_size=fret_fit.fit_E_E_size)
        Mask = Sel_mask(self, select_bursts.E, E1=E1, E2=E2)
        fit_res = np.array([fun_d[kind](nd[mask], na[mask], **kwargs)
                for nd, na, mask in zip(self.nd,self.na,Mask)])
        fit_name = dict(slope='Linear slope fit', E_size='E_size fit')
        self.add(fit_E_res=fit_res, fit_E_name=fit_name[kind],
                E_fit=fit_res, fit_E_curve=False, fit_E_E1=E1, fit_E_E2=E2)
        self.fit_E_calc_variance()
        return self.E_fit

    def fit_E_two_gauss_EM(self, fit_func=two_gaussian_fit_EM,
                           weights='size', gamma=1., **kwargs):
        """Fit the E population to a Gaussian mixture model using EM method.
        Additional arguments in `kwargs` are passed to the fit_func().
        """
        fit_res = zeros((self.nch, 5))
        for ich, (nd, na, E) in enumerate(zip(self.nd, self.na, self.E)):
            w = fret_fit.get_weights(nd, na, weights, gamma)
            fit_res[ich, :] = fit_func(E, weights=w, **kwargs)
        self.add(fit_E_res=fit_res, fit_E_name=fit_func.__name__,
                E_fit=fit_res[:,2], fit_E_curve=True,
                fit_E_model=two_gauss_mix_pdf,
                fit_E_model_F=np.repeat(1, self.nch))
        return self.E_fit

    def fit_E_generic(self, E1=-1, E2=2, fit_fun=two_gaussian_fit_hist,
            weights=None, gamma=1., **fit_kwargs):
        """Fit E in each channel with `fit_fun` using burst in [E1,E2] range.
        All the fitting functions are defined in
        :mod:`fretbursts.fit.gaussian_fitting`.

        Parameters:
            weights (string or None): specifies the type of weights
                If not None `weights` will be passed to
                `fret_fit.get_weights()`. `weights` can be not-None only when
                using fit functions that accept weights (the ones ending in
                `_hist` or `_EM`)
            gamma (float): passed to `fret_fit.get_weights()` to compute
                weights

        All the additional arguments are passed to `fit_fun`. For example `p0`
        or `mu_fix` can be passed (see `fit.gaussian_fitting` for details).

        Note:
            Use this method for CDF/PDF or hist fitting.
            For EM fitting use :meth:`fit_E_two_gauss_EM()`.
        """
        if fit_fun.__name__.startswith("gaussian_fit"):
            fit_model = lambda x, p: SS.norm.pdf(x, p[0], p[1])
            if 'mu0' not in fit_kwargs: fit_kwargs.update(mu0=0.5)
            if 'sigma0' not in fit_kwargs: fit_kwargs.update(sigma0=0.3)
            iE, nparam = 0, 2
        elif fit_fun.__name__ == "two_gaussian_fit_hist_min_ab":
            fit_model = two_gauss_mix_ab
            if 'p0' not in fit_kwargs:
                fit_kwargs.update(p0=[0, .05, 0.5, 0.6, 0.1, 0.5])
            iE, nparam = 3, 6
        elif fit_fun.__name__.startswith("two_gaussian_fit"):
            fit_model = two_gauss_mix_pdf
            if 'p0' not in fit_kwargs:
                fit_kwargs.update(p0=[0, .05, 0.6, 0.1, 0.5])
            iE, nparam = 2, 5
        else:
            raise ValueError, "Fitting function not recognized."

        Mask = Sel_mask(self, select_bursts.E, E1=E1, E2=E2)

        fit_res, fit_model_F = zeros((self.nch, nparam)), zeros(self.nch)
        for ich, (nd, na, E, mask) in enumerate(
                                        zip(self.nd, self.na, self.E, Mask)):
            if '_hist' in fit_fun.__name__ or '_EM' in fit_fun.__name__:
                if weights is None:
                    w = None
                else:
                    w = fret_fit.get_weights(nd[mask], na[mask], weights, gamma)
                fit_res[ich, :] = fit_fun(E[mask], weights=w, **fit_kwargs)
            else:
                # Non-histogram fits (PDF/CDF) do not support weights
                fit_res[ich, :] = fit_fun(E[mask], **fit_kwargs)
            fit_model_F[ich] = 1.*mask.sum()/mask.size

        # Save enough info to generate a fit plot (see hist_fret in burst_plot)
        self.add(fit_E_res=fit_res, fit_E_name=fit_fun.__name__,
                E_fit=fit_res[:,iE], fit_E_curve=True, fit_E_E1=E1,fit_E_E2=E2,
                fit_E_model=fit_model, fit_E_model_F=fit_model_F,
                fit_E_weights=weights, fit_E_gamma=gamma,
                fit_E_kwargs=fit_kwargs)
        return self.E_fit

    def fit_from(self, D):
        """Copy fit results from another Data() variable.
        Now that the fit methods accept E1,E1 parameter this probabily useless.
        """
        fit_data = ['fit_E_res', 'fit_E_name', 'E_fit', 'fit_E_curve',
                'fit_E_E1', 'fit_E_E2=E2', 'fit_E_model', 'fit_E_model_F',
                'fit_guess', 'fit_fix']  # NOTE Are these last two still used ?
        for name in fit_data:
            if name in D:
                self[name] = D[name]
                setattr(self, name, self[name])
        # Deal with the normalization to the number of bursts
        self.add(fit_model_F=r_[[1.*old_E.size/new_E.size \
                for old_E,new_E in zip(D.E, self.E)]])

    def fit_E_calc_variance(self, weights='sqrt', dist='DeltaE',
            E_fit=None, E1=-1, E2=2):
        """Compute several versions of WEIGHTED std.dev. of the E estimator.
        `weights` are multiplied *BEFORE* squaring the distance/error
        `dist` can be 'DeltaE' or 'SlopeEuclid'

        Note:
            This method is still experimental
        """
        assert dist in ['DeltaE', 'SlopeEuclid']
        if E_fit is None:
            E_fit = self.E_fit
            E1 = self.fit_E_E1 if 'fit_E_E1' in self else -1
            E2 = self.fit_E_E2 if 'fit_E_E2' in self else 2
        else:
            # If E_fit is not None the specified E1,E2 range is used
            if E1 < 0 and E2 > 1:
                pprint('WARN: E1 < 0 and E2 > 1 (wide range of E eff.)\n')
        if size(E_fit) == 1 and self.nch > 0:
            E_fit = np.repeat(E_fit, self.nch)
        assert size(E_fit) == self.nch

        E_sel = [Ei[(Ei>E1)*(Ei<E2)] for Ei in self.E]
        Mask = Sel_mask(self, select_bursts.E, E1=E1, E2=E2)

        E_var, E_var_bu, E_var_ph = \
                zeros(self.nch), zeros(self.nch), zeros(self.nch)
        for i, (Ech, nt, mask) in enumerate(zip(E_sel, self.nt, Mask)):
            nt_s = nt[mask]
            nd_s, na_s = self.nd[i][mask], self.na[i][mask]
            w = fret_fit.get_weights(nd_s, na_s, weights)
            info_ph = nt_s.sum()
            info_bu = nt_s.size

            if dist == 'DeltaE':
                distances = (Ech - E_fit[i])
            elif dist == 'SlopeEuclid':
                distances = fret_fit.get_dist_euclid(nd_s, na_s, E_fit[i])

            residuals = distances * w
            var = np.mean(residuals**2)
            var_bu = np.mean(residuals**2)/info_bu
            var_ph = np.mean(residuals**2)/info_ph
            #lvar = np.mean(log(residuals**2))
            #lvar_bu = np.mean(log(residuals**2)) - log(info_bu)
            #lvar_ph = np.mean(log(residuals**2)) - log(info_ph)
            E_var[i], E_var_bu[i], E_var_ph[i] = var, var_bu, var_ph
            assert (-np.isnan(E_var[i])).all() # check there is NO NaN
        self.add(E_var=E_var, E_var_bu=E_var_bu, E_var_ph=E_var_ph)
        return E_var


