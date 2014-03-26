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
import cPickle as pickle
import numpy as np
from numpy import zeros, size, r_
import scipy.stats as SS

from utils.misc import pprint, clk_to_s, shorten_fname
from poisson_threshold import find_optimal_T_bga
import fret_fit

from burstsearch.burstsearchlib import (
        itstart, iwidth, inum_ph, iistart, iiend, itend,
        ba_pure_o, mch_count_ph_in_bursts, 
        b_start, b_end, b_width,
        b_istart, b_iend, b_size, 
        b_rate, 
        b_separation)

try:
    from burstsearch.c_burstsearch import ba_c#, ba_pure_c
    ba = ba_c
    print " - Optimized burst search (cython) loaded."
except:
    ba = ba_pure_o
    print " - Fallback to pure python burst search."
try:
    from burstsearch.c_burstsearch import c_mch_count_ph_in_bursts
    mch_count_ph_in_bursts = c_mch_count_ph_in_bursts
    print " - Optimized ph_count (cython) loaded."
except:
    print " - Fallback to pure python ph_count."

import background as bg
import select_bursts
from fit.gaussian_fitting import (gaussian_fit_hist,
                                  gaussian_fit_cdf,
                                  two_gaussian_fit_hist,
                                  two_gaussian_fit_hist_min,
                                  two_gaussian_fit_hist_min_ab,
                                  two_gaussian_fit_EM,
                                  two_gauss_mix_pdf,
                                  two_gauss_mix_ab,)


def deprecate(function, old_name, new_name):
    def deprecated_function(*args, **kwargs):
        pprint("Function %s is deprecated, use %s instead.\n" %\
                (old_name, new_name))
        res = function(*args, **kwargs)
        return res
    return deprecated_function

# Redefine some old functions that have been renamed so old scripts will not
# break but will print a warning
bg_calc_exp = deprecate(bg.exp_fit, 'bg_calc_exp', 'bg.exp_fit')
bg_calc_exp_cdf = deprecate(bg.exp_fit, 'bg_calc_exp_cdf', 'bg.exp_cdf_fit')

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

    ##
    ## Alternative implementation (always applies corrections)
    ##
    # NOTE that boolean masking implies numpy array copy
    # On the contrary slicing only makes a new view of the array
    #new_mburst = [mburst[mask] for mburst,mask in zip(d_orig.mburst,Masks)]
    #ds.add(mburst=new_mburst)
    #ds.calc_fret(count_ph=True, corrections=True, dither=d_orig.dithering)
    fields = ['mburst', 'nd', 'na', 'nt', 'bp']
    if d_orig.ALEX: fields += ['naa']
    if hasattr(d_orig, 'max_rate'): fields += ['max_rate']
    for k in fields:
        # Reading hint: think of d[k] as another name for (for ex.) d.nd, etc...

        # Make new lists to contain the filtered data
        # (otherwise changing ds.nd[0] changes also d_orig.nd[0])
        ds[k] = [np.array([])]*d_orig.nch
        # Assign also the attribute to maintain the same object interface
        setattr(ds,k,ds[k])
        # Assign the new data
        for i,mask in enumerate(Masks):
            if d_orig[k][i].size == 0: continue # -> no bursts in current ch
            # Note that boolean masking implies numpy array copy
            # On the contrary slicing only makes a new view of the array
            ds[k][i] = d_orig[k][i][mask]
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
ph_rate = lambda m, ph: 1.*m/(ph[m-1:] - ph[:ph.size-m+1])   # rate
ph_rate_t = lambda m, ph: 0.5*(ph[m-1:] + ph[:ph.size-m+1])  # time for rate

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
def b_rate_max(b, ph, m=3):
    """Returns the max (m-photons) rate reached inside each bursts.
    """
    PHB = [ph[ bu[iistart] : bu[iiend]+1 ] for bu in b]
    rates_max = np.array([ph_rate(m=m, ph=phb).max() for phb in PHB])
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
    """Returns the all the m-ph-rates of burst #bi."""
    return ph_rate(m, b_ph_times(b, bi, ph, pad=pad))


def find_burst(bursts, size, width_ms, clk_p=12.5e-9):
    """Find b_index of burst(s) of given size AND width."""
    width = (width_ms*1e-3)/clk_p
    th = 0.01e-3/clk_p # 800clk or 10us @ clk_p=12.5e-9s
    return np.where((b_size(bursts) == size)*(abs(b_width(bursts)-width) < th))

def b_fuse(mburst, ms=0, clk_p=12.5e-9):
    """Fuse bursts separated by less than `ms` (milli-secs).
    
    Note:
        This function will fuse only 2 consecutive bursts. If there are groups
        of 3 or more consecutive bursts all closer that `ms` the function
        must be run several times to fuse all them together.
    
    Parameters:
        mburst (2D array): Nx6 array of burst data, one row per burst
            See `burstseach.burstseachlib.py` for details.
        ms (float):
            minimum waiting time between bursts (in millisec). Burst closer 
            than that will be fuse in asingle burst.
    
    Returns:
        new_mbursts (2D array): new array of burst data
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
    #pprint(" - Fused %4d of %5d bursts (%.1f%%).\n" %\
    #        (first_burst.sum(), mburst.shape[0],
    #         100.*first_burst.sum()/mburst.shape[0]))
    return new_burst[reorder, :]

def mch_fuse_bursts(MBurst, ms=0, clk_p=12.5e-9):
    """Multi-ch version of `fuse_bursts`. `MBurst` is a list of arrays.
    """
    mburst = [b.copy() for b in MBurst] # safety copy
    new_mburst = []
    ch = 0
    for mb in mburst:
        ch += 1
        print " - - - - - CHANNEL %2d - - - - " % ch
        if mb.size == 0: continue
        z = 0
        init_nburst = mb.shape[0]
        new_nburst, nburst = 0, 1  # starting condition
        while (new_nburst < nburst):
            z += 1
            nburst = mb.shape[0]
            mb = b_fuse(mb, ms=ms, clk_p=clk_p)
            new_nburst = mb.shape[0]
        new_mburst.append(mb)
        delta_b = init_nburst-nburst
        pprint(" --> [CH %d] END Fused %d bursts (%.1f%%, %d iter)\n\n" %\
                (ch, delta_b, 100.*delta_b/init_nburst, z))
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

def ES_histo(E, S, bin_step=0.05, E_bins=None, S_bins=None):
    """Returns 2D (ALEX) histogram and bins of bursts (E,S).
    """
    if E_bins is None:
        E_bins = np.arange(-0.2, 1.2+1e-4, bin_step)
    if S_bins is None:
        S_bins = np.arange(-0.2, 1.2+1e-4, bin_step)
    H, E_bins, S_bins = np.histogram2d(E, S, bins=[E_bins, S_bins])
    return H, E_bins, S_bins

def gamma_correct_E(Er, gamma):
    """Apply gamma correction to the uncorrected FRET `Er`."""
    Er = np.asarray(Er)
    return Er/(gamma-gamma*Er+Er)

def gamma_uncorrect_E(E, gamma):
    """Reverse gamma correction and return uncorrected FRET."""
    E = np.asarray(E)
    return gamma*E/(1 - E + gamma*E)

def delta(x):
    """Return x.max() - x.min()"""
    return x.max() - x.min()


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
    In the follow all the standard data-attributes are listed.
    
    Note: 
        Attributes of type "*list*" contain one element per channel.
        Each element, in turn, can be an array. For example `.ph_times_m[i]`
        is the array of timestamps for channel `i`; or `.nd[i]` is the array
        of donor counts in each burst for channel `i`.

    **Measurement attributes**
    
    Attributes:
        fname (string): measuremets file name
        nch (int): number of channels
        clk_p (float): clock period in seconds for timestamps in `ph_times_m`
        ph_times_m (list): list of timestamp arrays (int64). Each array 
            contains all the timestamps (donor+acceptor) in one channel.
        A_em (list): list of boolean arrays marking acceptor timestamps. Each 
            array is a boolean mask for the corresponding ph_times_m array.
        BT (float or array of floats): bleedthrough or leakage fraction.
            May be scalar or same size as nch.
        gamma (float or array of floats): gamma factor.
            May be scalar or same size as nch.
        A_em (list of boolean arrays):
            boolean mask for `.ph_times_m[i]` for acceptor emission
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
    
    These attributes contain the estimated background rate. Each arribute is
    a list (one element per channel) of arrays. Each array contains several
    rates computed every `time_s` seconds of measurements. `time_s` is an 
    argument passed to `.calc_bg()` method. Each `time_s` measurement slice
    is here called **`period`**.

    Attributes:
        bg (list of arrays):  total background (donor + acceptor) during donor 
            excitation
        bg_dd (list of arrays): background of donor emission during donor 
            excitation
        bg_ad (list of arrays): background of acceptor emission during donor 
            excitation
        bg_aa (list of arrays): background of acceptor emission during 
            acceptor excitation
        nperiods (int): number of periods in which timestamps are splitted for 
            background calculation
        bg_fun (function): function used to compute the background rates
        Lim (list): for each ch, is a list of index pairs for `.ph_times_m[i]` 
            for **first** and **last** photon in each period.
        Ph_p (list): for each ch, is a list of timestamps pairs for 
            **first** and **last** photon of each period.

    Other attributes (per-ch mean of `bg`, `bg_dd`, `bg_ad` and `bg_aa`)::
    
        rate_m: array of bg rates for D+A channel pairs (ex. 4 for 4 spots)
        rate_dd: array of bg rates for D em (and D ex if ALEX)
        rate_da: array of bg rates for A em (and D ex if ALEX)
        rate_aa: array of bg rates for A em and A ex (only for ALEX)

    **Burst search parameters (user input)**

    Attributes:    
        ph_sel (string): valid values 'DA' (default), 'D' or 'A'
            type of ph selection for burst search (donor, acceptor or both)
        m (int): number of consecutive timestamps used to compute the 
            local rate during burst search
        L (int): min. number of photons for a burst to be identified and saved
        P (float, probability): valid values [0..1].
            Probability that a burst-start is due to a Poisson background.
            The emplyed Poisson rate is the one computed by `.calc_bg()`.
        F (float): `(F * background_rate)` is the minimum rate for burst-start

    **Burst search data (available after burst search)**
    
    In the following, when not specified, in parameters marked as (list of 
    arrays) each array contains one element per bursts. `mburst` contain one 
    "row" per burst. `TT` arrays contain one element per `period` (see above:
    background attributes).
    
    Attributes:
        mburst (list of arrays): list of 2-D arrays containing burst data. 
            Each row contains data for 1 burst. Each column contains burst 
            fields like burst start, end, duration, size and indexes.
            The proper way to access these fields is through the
            function b_* (ex: b_start, b_end, b_width, b_size, etc...).
            
        TT (list of arrays): list of arrays of *T* values (in sec.). A *T*
            value is the maximum delay between `m` photons to have a 
            burst-start. Each channels has an array of *T* values, one for 
            each background "period" (see above).
        T (array): per-channel mean of `TT` 
    
        nd, na (list of arrays): number of donor or acceptor photons during 
            donor excitation in each burst
        nt (list of arrays): total number photons (nd+na)
        naa (list of arrays): number of acceptor photons in each bursts
            during acceptor excitation **[ALEX only]**
        bp (list of arrays): time period for each burst. Same shape as `nd`. 
            This is needed to identify the background rate for each burst.
        bg_bs (list): background rates used for threshold computation in burst 
            search (is a reference to `bg`, `bg_dd` or `bg_ad`).
    
        fuse (None or float): if not None, the burst seaparation in ms below
            which bursts have been fused (see `.fuse_bursts()`).
    
        E (list): FRET efficiency value for each burst:
                    E = na/(na + gamma*nd).
        S (list): stochiometry value for each burst:
                    S = (gamma*nd + na) /(gamma*nd + na + naa)
    """
    def __init__(self, **kwargs):
        # Default values
        init_kw = dict(ALEX=False, BT=0., gamma=1, chi_ch=1., s=[])
        # Override with user data
        init_kw.update(**kwargs)
        DataContainer.__init__(self, **init_kw)

    ##
    # Infrastructure methods: they return a new Data object
    #
    def copy(self):
        """Copy data in a new object. All arrays copied except for ph_times_m
        """
        print 'Deep copy executed.'
        new_d = Data(**self) # this make a shallow copy (like a pointer)

        ## Deep copy (not just reference) or array data
        for k in ['mburst', 'nd', 'na', 'nt', 'naa', 'E', 'S']:
            # Making sure k is defined
            if k in self:
                # Make the new list with a copy of the original arrays
                new_d[k] = [self[k][i].copy() for i in range(self.nch)]
                # Set the attribute: new_d.k = new_d[k]
                setattr(new_d, k, new_d[k])
        return new_d

    def iter_ph_times(self, ph_sel='DA'):
        """Iterator that returns the arrays of timestamps in `.ph_times_m`.
        This method always returns in-memory arrays, even when ph_times_m
        is a disk-baked list of arrays.
        `ph_sel` can be 'DA', D' or 'A' and allows to retrive timstamps for
        donor-acceptor, donor or acceptor emission.
        """
        for ich in xrange(self.nch):
            yield self.get_ph_times(ich, ph_sel=ph_sel)

    def get_ph_times(self, ich, ph_sel='DA'):
        """Returns timestamps array for channel `ich`.
        This method always returns in-memory arrays, even when ph_times_m
        is a disk-baked list of arrays.
        `ph_sel` can be 'DA', D' or 'A' and allows to retrive timstamps for
        donor-acceptor, donor or acceptor emission.
        """
        assert ph_sel in ['DA', 'D', 'A', 'AA']
        ph = self.ph_times_m[ich]
        if type(self.ph_times_m) is not list:
            ph = ph[:]

        if ph_sel == 'DA':
            return ph
        elif ph_sel == 'D':
            if self.ALEX:
                return ph[self.get_D_ex(ich)*self.get_D_em(ich)]
            else:
                return ph[self.get_D_em(ich)]
        elif ph_sel == 'A':
            if self.ALEX:
                return ph[self.get_D_ex(ich)*self.get_A_em(ich)]
            else:
                return ph[self.get_A_em(ich)]
        elif ph_sel == 'AA':
            if not self.ALEX:
                raise ValueError('Selection "AA" is valid only for ALEX data')
            return ph[self.A_ex[ich]*self.A_em[ich]]

    def _get_mask(self, ich, mask_name, negate=False):
        """Get the bool array `mask_name` for channel `ich`.
        If the bool array is a scalar return a slice (full or empty)
        """
        mask = np.asarray(getattr(self, mask_name)[ich])
        if negate:
            mask = np.logical_not(mask)
        if len(mask.shape) == 0:
            mask = slice(None) if mask else slice(0)
        return mask

    def get_A_em(self, ich):
        """Returns a mask to select photons detected in the acceptor ch."""
        return self._get_mask(ich, 'A_em')

    def get_D_em(self, ich):
        """Returns a mask to select photons detected in the donor ch."""
        return self._get_mask(ich, 'A_em', negate=True)

    def get_A_ex(self, ich):
        """Returns a mask to select photons in acceptor-excitation periods."""
        return self._get_mask(ich, 'A_ex')

    def get_D_ex(self, ich):
        """Returns a mask to select photons in donor-excitation periods."""
        return self._get_mask(ich, 'D_ex')

    def slice_ph(self, time_s1=5, time_s2=None, s='slice'):
        """Return a new Data object with ph in [`time_s1`,`time_s2`] (seconds)
        """
        if time_s2 is None: time_s2 = self.time_max()
        t1_clk, t2_clk = time_s1/self.clk_p, time_s2/self.clk_p
        assert np.array([t1_clk < ph.max() for ph in self.ph_times_m]).all()

        masks = ((ph > t1_clk)*(ph <= t2_clk) for ph in self.iter_ph_times())

        sliced_d = Data()
        slice_fields = ['ph_times_m', 'A_em', 'D_em', 'A_ex', 'D_ex']
        for name in slice_fields:
            if name in self:
                sliced_d[name] = [data[m] for data, m in zip(self[name], masks)]
                setattr(sliced_d, name, sliced_d[name])
        copy_fields = ['fname', 'nch', 'clk_p', 'BT', 'gamma', 'ALEX',
                'alex_period', 'D_ON', 'A_ON']
        for name in copy_fields:
            if name in self:
                sliced_d[name] = self[name]
                setattr(sliced_d, name, sliced_d[name])
        # Start timestamps from 0 to avoid problems with BG calc
        for ich in range(self.nch):
            ph_i = sliced_d.get_ph_times(ich)
            ph_i -= t1_clk
        sliced_d.s.append(s)
        return sliced_d

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

    def collapse(self):
        """Returns an object with 1-ch data joining the multi-ch data.
        """
        dc = Data(**self)

        mburst = np.vstack(self.mburst)
        burst_start = b_start(mburst)
        sort_index = burst_start.argsort()
        dc.add(mburst=[mburst[sort_index]])

        ich_burst = [i*np.ones(nb) for i, nb in enumerate(self.num_bu())]
        dc.add(ich_burst=np.hstack(ich_burst)[sort_index])

        for k in ['nd', 'na', 'naa', 'nt', 'bp']:
            if k in self:
                dc[k] = [np.hstack(self[k])[sort_index]]
                setattr(dc, k, dc[k])
        dc.add(nch=1)
        dc.add(chi_ch=1.)
        dc.update_gamma(np.mean(self.get_gamma_array()))
        return dc

    ##
    # Utility methods
    #
    def get_params(self):
        """Returns a plain dict containing only parameters and no arrays.
        This can be used as a summary of data analysis parameters.
        An addtional keys `name' and `Names` are added with values
        from `.name()` and `.Name()`.
        """
        p_names = ['fname', 'clk_p', 'nch', 'ph_sel', 'L', 'm', 'F', 'P',
                'BT', 'gamma', 'bg_time_s', 'nperiods',
                'rate_dd', 'rate_ad', 'rate_aa', 'rate_m', 'T', 'rate_th',
                'bg_corrected', 'bt_corrected', 'dithering', #'PP', 'TT',
                'chi_ch', 's', 'ALEX']
        p_dict = dict(self)
        for name in p_dict.keys():
            if name not in p_names:
                p_dict.pop(name)
        p_dict.update(name=self.name(), Name=self.Name())
        return p_dict

    def expand(self, ich, width=False):
        """Return nd, na, bg_d, bg_a (and optionally width) for bursts in `ich`.
        """
        period = self.bp[ich]
        w = b_width(self.mburst[ich])*self.clk_p
        bg_a = self.bg_ad[ich][period]*w
        bg_d = self.bg_dd[ich][period]*w
        if width:
            return self.nd[ich], self.na[ich], bg_d, bg_a, w
        else:
            return self.nd[ich], self.na[ich], bg_d, bg_a

    def time_max(self):
        """Return the measurement time (last photon) in seconds."""
        if 'ph_times_m' in self:
            return max([t[-1]*self.clk_p for t in self.ph_times_m])
        elif 'ph_times_t' in self:
            return self.ph_times_t.max()
        elif 'mburst' in self:
            return max([b_end(mb)[-1]*self.clk_p for mb in self.mburst])
        else:
            raise ValueError("No timestamps or bursts found.")

    def num_bu(self):
        """Return an array with number of bursts for each channel."""
        return np.r_[[mb.shape[0] for mb in self.mburst]]

    def ph_select(self):
        """Return masks of ph inside bursts."""
        self.ph_in_burst = mch_ph_select(self.ph_times_m, self.mburst)

    def cal_max_rate(self, m):
        """Compute the max m-photon rate reached in each burst."""
        Max_Rate = [b_rate_max(mb, ph, m=m)/self.clk_p
                for ph, mb in zip(self.iter_ph_times(), self.mburst)]
        self.add(max_rate=Max_Rate)

    ##
    # Background analysis methods
    #
    def calc_bg_cache(self, fun, time_s=60, **kwargs):
        """Same as `calc_bg` but try to load from cache. Caches results.
        Example:
            d.calc_bg_cache(bg.exp_fit, time_s=20, tail_min_us=200)
        """
        kw_str = ''.join(['_'+k+str(kwargs[k]) for k in kwargs])
        fname = self.fname+'_'+fun.__name__+str(time_s)+kw_str+'.bg'
        try:
            var = pickle.load(open(fname,'rb'))
            if not 'ph_sel' in var:
                var.update(ph_sel='DA')
                pprint(' - Added ph_sel manually (old bg cache file).\n')
            self.add(**var)
            pprint(" - BG rates loaded from cache (%s).\n" % fname)
        except:
            self.calc_bg(fun, time_s, **kwargs)
            keys = ['bg', 'bg_dd', 'bg_ad', 'bg_aa', 'rate_m', 'Lim', 'Ph_p',
                    'nperiods', 'bg_fun', 'bg_time_s', 'ph_sel',
                    'rate_dd', 'rate_ad', 'rate_aa']
            D = dict([(k, self[k]) for k in keys])
            pprint(" - Pickling BG ... ")
            pickle.dump(D, open(fname,'wb'), -1)
            pprint("DONE\n")

    def calc_bg(self, fun, time_s=60, **kwargs):
        """Compute time-dependent background rates for all the channels.
        
        Compute background rates for donor, acceptor and both detectors.
        The rates are computed every `time_s` seconds, allowing to 
        track possible variations during the measurement.
        
        Arguments:
            fun (function): function for background estimation (example 
                `bg.exp_fit`)
            time_s (float, seconds): compute backgound each time_s seconds
            kwargs: additional arguments to be passed to `fun`.
        
        The background estimation functions are defined in the module
        `background` (imported as `bg` in burstlib).
        
        Example:
            Compute background with `bg.exp_fit`, every 20s, with a
            threshold of 200us for photon waiting time::
            
               d.calc_bg(bg.exp_fit, time_s=20, tail_min_us=200)
            
        Returns:
            None, all the results are saved in the object.
        """
        pprint(" - Calculating BG rates ... ")

        if 'tail_min_us' in kwargs:
            tail_min_us = kwargs.pop('tail_min_us')
            if np.size(tail_min_us) == 1:
                th_us = thd_us = tha_us = thaa_us = tail_min_us
            elif np.size(tail_min_us) == 3:
                th_us, thd_us, tha_us = tail_min_us
            elif  np.size(tail_min_us) == 4:
                th_us, thd_us, tha_us, thaa_us = tail_min_us
            else:
                raise ValueError('`tail_min_us` must be of size 1, 3 or 4.')
        else:
            th_us = thd_us = tha_us = thaa_us = 300
            print 'Warning: no threshold specified, using default (300)'

        kwargs.update(clk_p=self.clk_p)
        time_clk = time_s/self.clk_p
        BG, BG_dd, BG_ad, BG_aa, Lim, Ph_p = [], [], [], [], [], []
        rate_m, rate_dd, rate_ad, rate_aa = [], [], [], []
        BG_err, BG_dd_err, BG_ad_err, BG_aa_err = [], [], [], []
        for ich, ph in enumerate(self.iter_ph_times()):
            nperiods = int(np.ceil(ph[-1]/time_clk))
            if not self.ALEX:
                ad_mask = self.A_em[ich]
                if type(ad_mask) is bool:
                    dd_mask =  not ad_mask
                else:
                    dd_mask = -ad_mask
            else:
                dd_mask = self.D_em[ich]*self.D_ex[ich]
                ad_mask = self.A_em[ich]*self.D_ex[ich]
                aa_mask = self.A_em[ich]*self.A_ex[ich]

            lim, ph_p = [], []
            bg, bg_dd, bg_ad, bg_aa = [zeros(nperiods) for _ in xrange(4)]
            zeros_list = [zeros(nperiods) for _ in xrange(4)]
            bg_err, bg_dd_err, bg_ad_err, bg_aa_err = zeros_list
            for ip in xrange(nperiods):
                i0 = 0 if ip == 0 else i1           # pylint: disable=E0601
                i1 = (ph < (ip+1)*time_clk).sum()

                ph_i = ph[i0:i1]
                bg[ip], bg_err[ip] = fun(ph_i, tail_min_us=th_us, **kwargs)

                # This to support cases of D-only or A-only timestamps
                # where self.A_em[ich] is a bool and not a bool-array
                if type(dd_mask) is bool:
                    if dd_mask: bg_dd[ip], bg_dd_err[ip] = bg[ip], bg_err[ip]
                else:
                    dd_mask_i = dd_mask[i0:i1]
                    if dd_mask_i.any():
                        bg_dd[ip], bg_dd_err[ip] = fun(
                                ph_i[dd_mask_i], tail_min_us=thd_us, **kwargs)

                if type(ad_mask) is bool:
                    if ad_mask: bg_ad[ip], bg_ad_err[ip] = bg[ip], bg_err[ip]
                else:
                    ad_mask_i = ad_mask[i0:i1]
                    if ad_mask_i.any():
                        bg_ad[ip], bg_ad_err[ip] = fun(
                                ph_i[ad_mask_i], tail_min_us=tha_us, **kwargs)

                if self.ALEX and aa_mask.any():
                    aa_mask_i = aa_mask[i0:i1]
                    bg_aa[ip], bg_aa_err[ip] = fun(
                                ph_i[aa_mask_i], tail_min_us=thaa_us, **kwargs)

                lim.append((i0, i1-1))
                ph_p.append((ph[i0], ph[i1-1]))
            
            Lim.append(lim); Ph_p.append(ph_p)            
            BG.append(bg); BG_err.append(bg_err)
            BG_dd.append(bg_dd); BG_dd_err.append(bg_dd_err)
            BG_ad.append(bg_ad); BG_ad_err.append(bg_ad_err)
            BG_aa.append(bg_aa); BG_aa_err.append(bg_aa_err)
            rate_m.append(bg.mean())
            rate_dd.append(bg_dd.mean())
            rate_ad.append(bg_ad.mean())
            if self.ALEX: rate_aa.append(bg_aa.mean())
        self.add(bg=BG, bg_dd=BG_dd, bg_ad=BG_ad, bg_aa=BG_aa, bg_err=BG_err, 
                 bg_dd_err=BG_dd_err, bg_ad_err=BG_ad_err, bg_aa_err=BG_aa_err, 
                 Lim=Lim, Ph_p=Ph_p, nperiods=nperiods, bg_fun=fun, 
                 bg_time_s=time_s, ph_sel='DA', rate_m=rate_m, 
                 rate_dd=rate_dd, rate_ad=rate_ad, rate_aa=rate_aa)
        pprint("[DONE]\n")

    def recompute_bg_lim_ph_p(self, ph_sel='D'):
        """Recompute self.Lim and selp.Ph_p relative to ph selection `ph_sel`
        `ph_sel` (can be 'D','A', or 'DA') selects the timestamps (donor,
        acceptor or both) on which self.Lim and selp.Ph_p are computed.
        """
        assert ph_sel in ['DA', 'D', 'A']
        if 'ph_sel' in self and self.ph_sel == ph_sel: return

        pprint(" - Recomputing limits for current ph selection (%s) ... " % \
                ph_sel)
        bg_time_clk = self.bg_time_s/self.clk_p
        Lim, Ph_p = [], []
        for ph_x, lim in zip(self.iter_ph_times(ph_sel), self.Lim):
            lim, ph_p = [], []
            for ip in xrange(self.nperiods):
                i0 = 0 if ip == 0 else i1  # pylint: disable=E0601
                i1 = (ph_x < (ip+1)*bg_time_clk).sum()
                lim.append((i0, i1-1))
                ph_p.append((ph_x[i0], ph_x[i1-1]))
            Lim.append(lim)
            Ph_p.append(ph_p)
        self.add(Lim=Lim, Ph_p=Ph_p, ph_sel=ph_sel)
        pprint("[DONE]\n")

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

    def _calc_T(self, m, P, F=1., ph_sel='DA'):
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
        BG = {'DA': self.bg, 'D': self.bg_dd, 'A': self.bg_ad}
        for bg_ch, F_ch, P_ch in zip(BG[ph_sel], FF, PP):
            # All "T" are in seconds            
            Tch = find_T(m, F_ch, P_ch, bg_ch)
            TT.append(Tch)
            T.append(Tch.mean())
            rate_th.append(np.mean(m/Tch))
        self.add(TT=TT, T=T, bg_bs=BG[ph_sel], FF=FF, PP=PP, F=F, P=P, 
                 rate_th=rate_th)

    def _select_ph_times(self, ph_sel):
        """Return a list of ph_times with photons specified by `ph_sel`.
        `ph_sel` is a string that can be 'D', 'A' or 'DA'.
        """
        assert ph_sel in ['DA', 'D', 'A']
        if ph_sel == 'DA':
            PH = self.ph_times_m
        elif ph_sel == 'D':
            PH = [p[-a] for p, a in zip(self.iter_ph_times(), self.A_em)]
        elif ph_sel == 'A':
            PH = [p[a] for p, a in zip(self.iter_ph_times(), self.A_em)]
        return PH

    def _burst_search_rate(self, m, L, min_rate_cps, ph_sel='DA', 
                           verbose=True):
        """Experimental burst search with one fixed rate.
        May break stuff downstream!!
        """
        Min_rate_cps = self._param_as_mch_array(min_rate_cps)
        mburst = []
        T_clk = (1.*m/Min_rate_cps)/self.clk_p
        for ich, (ph, t_clk) in enumerate(zip(self.iter_ph_times(ph_sel), 
                                          T_clk)):
            label = '%s CH%d' % (ph_sel, ich+1) if verbose else None
            mb = ba(ph, L, m, t_clk, label=label, verbose=verbose) 
            mburst.append(mb)
        self.add(mburst=mburst, min_rate_cps=Min_rate_cps, T=T_clk*self.clk_p)

    def _burst_search_TT(self, m, L, ph_sel='DA', verbose=True):
        """Compute burst search with params `m`, `L` on ph selection `ph_sel`
        
        Requires a list of arrays `self.TT` with the max time-thresholds in
        the different burst periods for each channel (use `._calc_T()`).
        """
        self.recompute_bg_lim_ph_p(ph_sel=ph_sel)
        MBurst = []
        label = ''
        for ich, (ph, T) in enumerate(zip(self.iter_ph_times(ph_sel),
                                          self.TT)):
            MB = []
            Tck = T/self.clk_p
            for ip, (l0, l1) in enumerate(self.Lim[ich]):
                if verbose:
                    label='%s CH%d-%d' % (ph_sel, ich+1, ip)
                mb = ba(ph[l0:l1+1], L, m, Tck[ip], label=label,
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
        if ph_sel != 'DA':
            # Convert the burst data to be relative to ph_times_m.
            # Convert both Lim/Ph_p and mburst, as they are both needed
            # to compute `.bp`.
            self.recompute_bg_lim_ph_p(ph_sel='DA')
            self._fix_mburst_from(ph_sel=ph_sel)

    def _fix_mburst_from(self, ph_sel):
        """Convert burst data from 'D' or 'A' timestamps to 'DA' timestamps
        """
        assert ph_sel in ['A', 'D']
        pprint(' - Fixing  burst data to refer to ph_times_m ... ')
        Mask = self.A_em if ph_sel == 'A' else [-a for a in self.A_em]
        old_MBurst = [mb.copy() for mb in self.mburst]
        # Note that mburst is modified in-place
        for mburst, ph_times, mask in \
                    zip(self.mburst, self.iter_ph_times(), Mask):
            index = np.arange(ph_times.size, dtype=np.int32)
            mburst[:, iistart] = index[mask][mburst[:, iistart]]
            mburst[:, iiend] = index[mask][mburst[:, iiend]]
            mburst[:, inum_ph] = mburst[:, iiend] - mburst[:, iistart] + 1

        for mb, old_mb in zip(self.mburst, old_MBurst):
            assert (mb[:, iistart] >= old_mb[:, iistart]).all()
            assert (mb[:, iiend] >= old_mb[:, iiend]).all()
            assert (mb[:, inum_ph] >= old_mb[:, inum_ph]).all()
        pprint('[DONE]\n')

    def burst_search_t(self, L=10, m=10, P=0.95, F=1., min_rate_cps=None, 
            nofret=False, max_rate=False, dither=False, ph_sel='DA', 
            verbose=False, mute=False):
        """Performs a burst search with specified parameters.
        
        This method performs a sliding-window burst search and does not 
        require binning the timestamps. The burst starts when the rate of `m`
        photons is above a minimum rate, and stops when the rate falls below
        the threshold.
        
        The minimum rate can be explicitly specified (`min_rate`) or computed
        as a function of the background rate (using `F` or `P`).
        
        Parameters:        
            m (int): number of of consecutive photons used to compute the 
                photon rate.
            L (int): minimum number of photons in burst
            P (float): Probability to detect non-Poisson bursts. If not None 
                this overrides `F`.
            F (float): if `P` is None: min.rate/bg.rate ratio for burst search
                else: `P` refers to a Poisson rate = bg_rate*F
            min_rate (float or list/array): min. rate in cps for burst start.
                If not None has the precedence over `P` and `F`.
                If non-scalar, contains one rate per each channel.
            ph_sel (string): can be 'DA', 'D' or 'A'. Specifies the photon 
                selection on which to perform the burst search: 'DA' for 
                donor+acceptor, 'D' or 'A' for donor or acceptor photons.
        
        Note: 
            when using `P` or `F` the background rates are needed, so 
            `.calc_bg()` must be called before the burst search.
              
        Example:
            d.burst_search_t(L=10, m=10, F=6, P=None, ph_sel='DA')
            
        Returns:
            None, all the results are saved in the object.
        """
        pprint(" - Performing burst search (verbose=%s) ..." % verbose, mute)
        if min_rate_cps is not None:
            self._burst_search_rate(m=m, L=L, min_rate_cps=min_rate_cps, 
                                    ph_sel=ph_sel, verbose=verbose)
        else:
            # Compute TT
            self._calc_T(m=m, P=P, F=F, ph_sel=ph_sel)
            # Use TT and compute mburst
            self._burst_search_TT(L=L, m=m, ph_sel=ph_sel, verbose=verbose)
        pprint("[DONE]\n", mute)

        pprint(" - Calculating burst periods ...", mute)
        self._calc_burst_period()                       # writes bp
        pprint("[DONE]\n", mute)

        self.add(m=m, L=L)  # P and F are saved in _calc_T()
        self.add(bg_corrected=False, bt_corrected=False, dithering=False)
        for k in ['E', 'S', 'nd', 'na', 'naa', 'nt', 'fuse', 'lsb']:
            if k in self: self.delete(k)
        if not nofret:
            pprint(" - Counting D and A ph and calculating FRET ... \n", mute)
            self.calc_fret(count_ph=True, corrections=True, dither=dither, 
                           mute=mute)
            pprint("   [DONE Counting D/A]\n", mute)
        if max_rate:
            pprint(" - Computing max rates in burst ...", mute)
            self.cal_max_rate(m=3)
            pprint("[DONE]\n", mute)

    def cal_ph_num(self):
        """After burst search computes number of D and A ph in each burst.
        """
        if not self.ALEX:
            A_em = self.A_em
            if type(A_em[0]) is bool:
                bool_funcs = [np.ones if v else np.zeros for v in A_em]
                A_em = [bfunc(ph.shape, dtype=bool) for ph, bfunc in
                                    zip(self.ph_times_m, bool_funcs)]
            na = mch_count_ph_in_bursts(self.mburst, Mask=A_em)
            nt = [b[:, inum_ph].astype(float) if b.size > 0 else np.array([])\
                    for b in self.mburst]
            nd = [t-a for t, a in zip(nt, na)]
            assert (nt[0] == na[0] + nd[0]).all()
        if self.ALEX:
            Mask = [d_em*d_ex for d_em, d_ex in zip(self.D_em, self.D_ex)]
            nd = mch_count_ph_in_bursts(self.mburst, Mask)

            Mask = [a_em*d_ex for a_em, d_ex in zip(self.A_em, self.D_ex)]
            na = mch_count_ph_in_bursts(self.mburst, Mask)

            Mask = [a_em*a_ex for a_em, a_ex in zip(self.A_em, self.A_ex)]
            naa = mch_count_ph_in_bursts(self.mburst, Mask)
            self.add(naa=naa)

            nt = [d+a+aa for d, a, aa in zip(nd, na, naa)]
            assert (nt[0] == na[0] + nd[0] + naa[0]).all()
        self.add(nd=nd, na=na, nt=nt,
                bt_corrected=False, bg_corrected=False, dithering=False)

    def fuse_bursts(self, ms=0, process=True):
        """Return a new Data() fusing bursts separated by less than ms.
        If process==True, it applies corrections and (re)compute FRET.
        """
        if ms < 0: return self
        mburst = mch_fuse_bursts(self.mburst, ms=ms, clk_p=self.clk_p)
        new_d = Data(**self)
        for k in ['E', 'S', 'nd', 'na', 'naa', 'nt', 'lsb', 'bp']:
            if k in new_d: new_d.delete(k)
        new_d.add(bt_corrected=False, bg_corrected=False, dithering=False)
        new_d.add(mburst=mburst, fuse=ms)
        if 'bg' in new_d: new_d._calc_burst_period()
        if process:
            pprint(" - Counting D and A ph and calculating FRET ... \n")
            new_d.calc_fret(count_ph=True, corrections=True,
                            dither=self.dithering)
            pprint("   [DONE Counting D/A and FRET]\n")
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
            width = mb[:, iwidth]*self.clk_p
            period = self.bp[ich]
            self.nd[ich] -= self.bg_dd[ich][period] * width
            self.na[ich] -= self.bg_ad[ich][period] * width
            if relax_nt:
                # This does not guarantee that nt = nd + na
                self.nt[ich] -= self.bg[ich][period] * width
            else:
                self.nt[ich] = self.nd[ich] + self.na[ich]
            if self.ALEX:
                self.naa[ich] -= self.bg_aa[ich][period] * width
                self.nt[ich] += self.naa[ich]

    def bleed_through_correction(self, mute=False):
        """Apply bleed-through correction to burst sizes (nd, na,...)
        """
        if self.bt_corrected: return -1
        pprint("   - Applying leakage correction.\n", mute)
        assert (size(self.BT) == 1) or (size(self.BT) == self.nch)
        BT = self.get_BT_array()
        for i in range(self.nch):
            if self.na[i].size == 0: continue  # if no bursts skip this ch
            self.na[i] -= self.nd[i]*BT[i]
            self.nt[i] = self.nd[i] + self.na[i]
            if self.ALEX: self.nt[i] += self.naa[i]
        self.add(bt_corrected=True)

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
                naa += lsb*(np.random.rand(na.size)-0.5)
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
        """Apply both background BG and bleed-through (BT) corrections."""
        if 'bg' in self: self.background_correction_t(mute=mute)
        else: self.background_correction()
        self.bleed_through_correction(mute=mute)

    def update_bt(self, BT):
        """Change the bleed-through value `BT` and recompute FRET."""
        if not self.bt_corrected:
            # if not previously BT-corrected just apply BT correction
            # NOTE: previous BG correction or dithering is maintained
            self.add(BT=BT)
            self.bleed_through_correction()
        else:
            # if already BT-corrected need to recompute na,nd,nt to avoid
            # overcorrection
            self.add(BT=BT, bt_corrected=False)
            old_bg_corrected = self.bg_corrected
            old_dithering = self.dithering
            self.cal_ph_num()       # recomupte uncorrected na,nd...
            if old_bg_corrected: self.background_correction_t()
            self.bleed_through_correction()
            if old_dithering: self.dither(self.lsb)
        # Recompute FRET with no corrections (because already applied)
        self.calc_fret(count_ph=False, corrections=False)

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

    def get_BT_array(self):
        """Get the array of bleed-through coefficients (one per ch).
        Use this function to obtain an array of BT values
        regardless of the actual `BT` (that can be scalar).
        BT values are multiplied by `chi_ch`.
        """
        assert (size(self.BT) == 1) or (size(self.BT) == self.nch)
        BT = np.r_[[self.BT]*self.nch] if size(self.BT) == 1 else self.BT
        BT *= self.chi_ch
        return BT

    def calc_sbr(self):
        """Return Signal-to-background ratio for each burst.
        """
        sbr = []
        for ich, mb in enumerate(self.mburst):
            if mb.size == 0: 
                sbr.append([])                
                continue  # if no bursts skip this ch
            width = b_width(mb)*self.clk_p
            raw_sizes = b_size(mb)
            period = self.bp[ich]
            bg_bursts = self.bg_dd[ich][period] * width
            bg_bursts += self.bg_ad[ich][period] * width
            if self.ALEX:
                bg_bursts += self.bg_aa[ich][period] * width
            sbr.append(1.*raw_sizes/bg_bursts - 1)
        self.add(sbr=sbr)
    ##
    # FRET and stochiometry methods
    #
    def calc_fret(self, count_ph=False, corrections=True, dither=False, 
                  mute=False):
        """Compute FRET (and Stoichiometry if ALEX) for each burst."""
        if count_ph: self.cal_ph_num()
        if dither: self.dither(mute=mute)
        if corrections: self.corrections(mute=mute)
        self.calculate_fret_eff()
        if self.ALEX:
            self.calculate_stoich()
            self.calc_alex_hist()

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
        AH = [ES_histo(E, S, bin_step)[0] for E, S in zip(self.E, self.S)]
        _, E_bins, S_bins = ES_histo(E, S, bin_step)
        E_ax = E_bins[:-1] + 0.5*bin_step
        S_ax = S_bins[:-1] + 0.5*bin_step
        self.add(AH=AH, E_bins=E_bins, S_bins=S_bins, E_ax=E_ax, S_ax=S_ax)

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
        if 'bg_fun' in self: s += " BG%s" % self.bg_fun.__name__[8:]
        if 'bg_time_s' in self: s += "-%d" % self.bg_time_s
        if 'fuse' in self: s += " Fuse%.1fms" % self.fuse
        if 'bt_corrected' in self and self.bt_corrected:
            s += " BT%.3f" % np.mean(self.BT)
        if 'bg_corrected' in self and self.bg_corrected:
            s += " bg"
        if 'dithering' in self and self.dithering:
            s += " Dith%d" % self.lsb
        if 's' in self: s += ' '.join(self.s)
        return s +add

    def name(self):
        """Return short filename"""
        return shorten_fname(self.fname)[:-len('.dat')].replace("/","_")

    def Name(self, add=""):
        """Return short filename + status information."""
        n = self.status(add=add).replace(os.path.sep, '_').replace('/', '_')
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
        All the fitting functions are defined in `fit.gaussian_fitting`.
        
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
            For EM fitting use `fit_E_two_gauss_EM()`.
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
                fit_E_model=fit_model, fit_E_model_F=fit_model_F)
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


