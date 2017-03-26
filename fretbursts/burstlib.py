#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2013-2016 The Regents of the University of California,
#               Antonino Ingargiola <tritemio@gmail.com>
#
"""
This module contains all the main FRETBursts analysis functions.

`burstslib.py` defines the fundamental object `Data()` that contains both the
experimental data (attributes) and the high-level analysis routines (methods).

Furthermore it loads all the remaining **FRETBursts** modules (except for
`loaders.py`).

For usage example see the IPython Notebooks in sub-folder "notebooks".
"""

from __future__ import print_function, absolute_import, division
from future.utils import raise_from
from builtins import range, zip

import os
import hashlib
import numpy as np
import copy
from numpy import zeros, size, r_
import scipy.stats as SS

from .utils.misc import pprint, clk_to_s, deprecate
from .poisson_threshold import find_optimal_T_bga
from . import fret_fit
from . import bg_cache
from .ph_sel import Ph_sel
from .fretmath import gamma_correct_E, gamma_uncorrect_E

from .phtools import burstsearch as bslib
from .phtools.burstsearch import (
    # Burst search function
    bsearch,
    # Photon counting function,
    mch_count_ph_in_bursts
)
from .phtools import phrates
from . import background as bg
from . import select_bursts
from . import fit
from .fit.gaussian_fitting import (gaussian_fit_hist,
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
bg_calc_exp_cdf = deprecate(bg.exp_cdf_fit, 'bg_calc_exp_cdf', 'bg.exp_cdf_fit')


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

def isarray(obj):
    """Test if the object support the array interface.

    Returns True for numpy arrays and pandas sequences.
    """
    return hasattr(obj, '__array__')


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  BURST SELECTION FUNCTIONS
#

def Sel(d_orig, filter_fun, negate=False, nofret=False, **kwargs):
    """Uses `filter_fun` to select a sub-set of bursts from `d_orig`.

    This function is deprecated. Use :meth:`Data.select_bursts` instead.
    """
    d_sel = d_orig.select_bursts(filter_fun, negate=negate,
                                 computefret=not nofret,
                                 **kwargs)
    return d_sel


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Bursts and Timestamps utilities
#
def get_alex_fraction(on_range, alex_period):
    """Get the fraction of period beween two numbers indicating a range.
    """
    assert len(on_range) == 2
    if on_range[0] < on_range[1]:
        fraction = (on_range[1] - on_range[0]) / alex_period
    else:
        fraction = (alex_period + on_range[1] - on_range[0]) / alex_period
    return fraction

def top_tail(nx, a=0.1):
    """Return for each ch the mean size of the top `a` fraction.
    nx is one of nd, na, nt from Data() (list of burst size in each ch).
    """
    assert a > 0 and a < 1
    return np.r_[[n[n > n.max() * (1 - a)].mean() for n in nx]]


##
# Per-burst quatitites from ph-data arrays (timestamps, lifetime, etc..)
#

def _excitation_width(excitation_range, alex_period):
    """Returns duration of alternation period outside selected excitation.
    """
    if excitation_range[1] > excitation_range[0]:
        return alex_period - excitation_range[1] + excitation_range[0]
    elif excitation_range[1] < excitation_range[0]:
        return excitation_range[0] - excitation_range[1]

def _ph_times_compact(ph_times_sel, alex_period, excitation_width):
    """Compact ph_times inplace by removing gaps between alternation periods.

    Arguments:
        ph_times_sel (array): array of timestamps from one alternation period.
        alex_period (scalar): period of alternation in timestamp units.
        excitation_width (float): fraction of `alex_period` covered by
            current photon selection.

    Returns nothing, ph_times is modified in-place.
    """
    # The formula is
    #
    #   gaps = (ph_times_sel // alex_period)*excitation_width
    #   ph_times_sel = ph_times_sel - gaps
    #
    # As a memory optimization the `-gaps` array is reused inplace
    times_minusgaps = (ph_times_sel // alex_period) * (-1 * excitation_width)
    # The formula is ph_times_sel = ph_times_sel - "gaps"
    times_minusgaps += ph_times_sel
    return times_minusgaps

def iter_bursts_start_stop(bursts):
    """Iterate over (start, stop) indexes to slice photons for each burst.
    """
    arr_istart = bursts.istart
    arr_istop = bursts.istop + 1
    for istart, istop in zip(arr_istart, arr_istop):
        yield istart, istop

def iter_bursts_ph(ph_data, bursts, mask=None, compact=False,
                   alex_period=None, excitation_width=None):
    """Iterator over arrays of photon-data for each burst.

    Arguments:
        ph_data (1D array): array of photon-data (timestamps, nanotimes).
        bursts (Bursts object): bursts computed from `ph`.
        mask (boolean mask or None): if not None, is a boolean mask
            to select photons in `ph_data` (for example Donor-ch photons).
        compact (bool): if True, a photon selection of only one excitation
            period is required and the timestamps are "compacted" by
            removing the "gaps" between each excitation period.
        alex_period (scalar): period of alternation in timestamp units.
        excitation_width (float): fraction of `alex_period` covered by
            current photon selection.

    Yields an array with a selection of "photons" for each burst.
    """
    if isinstance(mask, slice) and mask == slice(None):
        mask = None
    if compact:
        assert alex_period is not None
        assert excitation_width is not None
        assert mask is not None
    for start, stop in iter_bursts_start_stop(bursts):
        ph = ph_data[start:stop]
        if mask is not None:
            ph = ph[mask[start:stop]]
        if compact:
            ph = _ph_times_compact(ph, alex_period, excitation_width)
        yield ph

def bursts_ph_list(ph_data, bursts, mask=None):
    """Returna list of ph-data for each burst.

    ph_data can be either the timestamp array on which the burst search
    has been performed or any other array with same size (boolean array,
    nanotimes, etc...)
    """
    return [ph for ph in iter_bursts_ph(ph_data, bursts, mask=mask)]

def burst_ph_stats(ph_data, bursts, func=np.mean, func_kw=None, **kwargs):
    """Reduce burst photons (timestamps, nanotimes) to a scalar using `func`.

    Arguments
        ph_data (1D array): array of photon-data (timestamps, nanotimes).
        bursts (Bursts object): bursts computed from `ph`.
        func (callable): function that takes the burst photon timestamps
            as first argument and returns a scalar.
        func_kw (callable): additional arguments in `func` beyond photon-data.
        **kwargs: additional arguments passed to :func:`iter_bursts_ph`.

    Return
        Array one element per burst.
    """
    if func_kw is None:
        func_kw = {}
    burst_stats = []
    for burst_ph in iter_bursts_ph(ph_data, bursts, **kwargs):
        burst_stats.append(func(burst_ph, **func_kw))
    return np.asfarray(burst_stats)  # NOTE: asfarray converts None to nan


def ph_in_bursts_mask(ph_data_size, bursts):
    """Return bool mask to select all "ph-data" inside any burst."""
    mask = zeros(ph_data_size, dtype=bool)
    for start, stop in iter_bursts_start_stop(bursts):
        mask[start:stop] = True
    return mask


def fuse_bursts_direct(bursts, ms=0, clk_p=12.5e-9, verbose=True):
    """Fuse bursts separated by less than `ms` (milli-seconds).

    This function is a direct implementation using a single loop.
    For a faster implementation see :func:`fuse_bursts_iter`.

    Parameters:
        bursts (BurstsGap object): bursts to be fused.
            See `phtools.burstsearch` for details.
        ms (float): minimum waiting time between bursts (in millisec).
            Bursts closer than that will be fused in a single burst.
        clk_p (float): clock period or timestamp units in seconds.
        verbose (bool): if True print a summary of fused bursts.

    Returns:
        A BurstsGap object containing the new fused bursts.
    """
    max_delay_clk = (ms * 1e-3) / clk_p

    fused_bursts_list = []
    fused_burst = None
    for burst1, burst2 in zip(bursts[:-1], bursts[1:]):
        if fused_burst is not None:
            burst1c = fused_burst
        else:
            burst1c = bslib.BurstGap.from_burst(burst1)

        separation = burst2.start - burst1c.stop
        if separation <= max_delay_clk:
            gap = burst2.start - burst1c.stop
            gap_counts = burst2.istart - burst1c.istop - 1
            if burst1c.istop >= burst2.istart:
                gap = 0
                gap_counts = 0

            fused_burst = bslib.BurstGap(
                start = burst1c.start,
                istart = burst1c.istart,
                stop = burst2.stop,
                istop = burst2.istop,
                gap = burst1c.gap + gap,
                gap_counts = burst1c.gap_counts + gap_counts)
        else:
            if fused_burst is not None:
                fused_bursts_list.append(fused_burst)
                fused_burst = None
            else:
                fused_bursts_list.append(bslib.BurstGap.from_burst(burst1c))

    # Append the last bursts (either a fused or an isolated one)
    if fused_burst is not None:
        fused_bursts_list.append(fused_burst)
    else:
        fused_bursts_list.append(bslib.BurstGap.from_burst(burst2))

    fused_bursts = bslib.BurstsGap.from_list(fused_bursts_list)

    init_num_bursts = bursts.num_bursts
    delta_b = init_num_bursts - fused_bursts.num_bursts
    pprint(" --> END Fused %d bursts (%.1f%%)\n\n" %
           (delta_b, 100 * delta_b / init_num_bursts), mute=not verbose)
    return fused_bursts


def fuse_bursts_iter(bursts, ms=0, clk_p=12.5e-9, verbose=True):
    """Fuse bursts separated by less than `ms` (milli-secs).

    This function calls iteratively :func:`b_fuse` until there are no more
    bursts to fuse. For a slower but more readable version see
    :func:`fuse_bursts_direct`.

    Parameters:
        bursts (BurstsGap object): bursts to be fused.
            See `phtools.burstsearch` for details.
        ms (float): minimum waiting time between bursts (in millisec).
            Bursts closer than that will be fused in a single burst.
        clk_p (float): clock period or timestamp units in seconds.
        verbose (bool): if True print a summary of fused bursts.

    Returns:
        A BurstsGap object containing the new fused bursts.
    """
    init_nburst = bursts.num_bursts
    bursts = bslib.BurstsGap(bursts.data)
    z = 0
    new_nburst, nburst = 0, 1  # starting condition
    while new_nburst < nburst:
        z += 1
        nburst = bursts.num_bursts
        bursts = b_fuse(bursts, ms=ms, clk_p=clk_p)
        new_nburst = bursts.num_bursts
    delta_b = init_nburst - nburst
    pprint(" --> END Fused %d bursts (%.1f%%, %d iter)\n\n" %
           (delta_b, 100 * delta_b / init_nburst, z), mute=not verbose)
    return bursts


def b_fuse(bursts, ms=0, clk_p=12.5e-9):
    """Fuse bursts separated by less than `ms` (milli-secs).

    This is a low-level function which fuses pairs of consecutive
    bursts separated by less than `ms` millisec.
    If there are 3 or more consecutive bursts separated by less than `ms`
    only the first 2 are fused.
    See :func:`fuse_bursts_iter` or :func:`fuse_bursts_direct` for
    higher level functions.

    Parameters:
        bursts (BurstsGap object): bursts to be fused.
            See `phtools.burstsearch` for details.
        ms (float): minimum waiting time between bursts (in millisec).
            Bursts closer than that will be fused in a single burst.
        clk_p (float): clock period or timestamp units in seconds.

    Returns:
        A BurstsGap object containing the new fused bursts.
    """
    max_delay_clk = (ms * 1e-3) / clk_p
    # Nearby bursts masks
    delays_below_th = (bursts.separation <= max_delay_clk)
    if not np.any(delays_below_th):
        return bursts

    buffer_mask = np.hstack([(False,), delays_below_th, (False,)])
    first_bursts = buffer_mask[1:]
    second_bursts = buffer_mask[:-1]

    # Keep only the first pair in case of more than 2 consecutive bursts
    first_bursts ^= (second_bursts * first_bursts)
    # note that previous in-place operation also modifies `second_bursts`

    both_bursts = first_bursts + second_bursts

    # istart is from the first burst, istop is from the second burst
    fused_bursts1 = bursts[first_bursts]
    fused_bursts2 = bursts[second_bursts]

    # Compute gap and gap_counts
    gap = fused_bursts2.start - fused_bursts1.stop
    gap_counts = fused_bursts2.istart - fused_bursts1.istop - 1  # yes it's -1
    overlaping = fused_bursts1.istop >= fused_bursts2.istart
    gap[overlaping] = 0
    gap_counts[overlaping] = 0

    # Assign the new burst data
    # fused_bursts1 has alredy the right start and istart
    fused_bursts1.istop = fused_bursts2.istop
    fused_bursts1.stop = fused_bursts2.stop
    fused_bursts1.gap += gap
    fused_bursts1.gap_counts += gap_counts

    # Join fused bursts with the remaining bursts
    new_burst = fused_bursts1.join(bursts[~both_bursts], sort=True)
    return new_burst


def mch_fuse_bursts(MBurst, ms=0, clk_p=12.5e-9, verbose=True):
    """Multi-ch version of `fuse_bursts`. `MBurst` is a list of Bursts objects.
    """
    mburst = [b.copy() for b in MBurst]  # safety copy
    new_mburst = []
    ch = 0
    for mb in mburst:
        ch += 1
        pprint(" - - - - - CHANNEL %2d - - - - \n" % ch, not verbose)
        if mb.num_bursts == 0:
            continue
        new_bursts = fuse_bursts_iter(mb, ms=ms, clk_p=clk_p, verbose=verbose)
        new_mburst.append(new_bursts)
    return new_mburst


def burst_stats(mburst, clk_p):
    """Compute average duration, size and burst-delay for bursts in mburst.
    """
    width_stats = np.array([[b.width.mean(), b.width.std()] for b in mburst
                            if b.num_bursts > 0]).T
    height_stats = np.array([[b.counts.mean(), b.counts.std()] for b in mburst
                             if b.num_bursts > 0]).T
    mean_burst_delay = np.array([b.separation.mean() for b in mburst
                                 if b.num_bursts > 0])
    return (clk_to_s(width_stats, clk_p) * 1e3, height_stats,
            clk_to_s(mean_burst_delay, clk_p))


def print_burst_stats(d):
    """Print some bursts statistics."""
    nch = len(d.mburst)
    width_ms, height, delays = burst_stats(d.mburst, d.clk_p)
    s = "\nNUMBER OF BURSTS: m = %d, L = %d" % (d.m, d.L)
    s += "\nPixel:          "+"%7d "*nch % tuple(range(1, nch+1))
    s += "\n#:              "+"%7d "*nch % tuple([b.num_bursts for b in d.mburst])
    s += "\nT (us) [BS par] "+"%7d "*nch % tuple(np.array(d.T)*1e6)
    s += "\nBG Rat T (cps): "+"%7d "*nch % tuple(d.bg_mean[Ph_sel('all')])
    s += "\nBG Rat D (cps): "+"%7d "*nch % tuple(d.bg_mean[Ph_sel(Dex='Dem')])
    s += "\nBG Rat A (cps): "+"%7d "*nch % tuple(d.bg_mean[Ph_sel(Dex='Aem')])
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
    if isinstance(mask, slice):
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

    def delete(self, *args, **kwargs):
        """Delete an element (attribute and/or dict entry). """
        warning = kwargs.get('warning', True)
        for name in args:
            try:
                self.pop(name)
            except KeyError:
                if warning:
                    print(' WARNING: Name %s not found (dict).' % name)
            try:
                delattr(self, name)
            except AttributeError:
                if warning:
                    print(' WARNING: Name %s not found (attr).' % name)


class Data(DataContainer):
    """
    Container for all the information (timestamps, bursts) of a dataset.

    Data() contains all the information of a dataset (name, timestamps, bursts,
    correction factors) and provides several methods to perform analysis
    (background estimation, burst search, FRET fitting, etc...).

    When loading a measurement file a Data() object is created by one
    of the loader functions in `loaders.py`. Data() objects can be also
    created with :meth:`Data.copy`, :meth:`Data.fuse_bursts()` or
    :meth:`Data.select_bursts`.

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

    The background is computed with :meth:`Data.calc_bg`
    and is estimated in chunks of equal duration called *background periods*.
    Estimations are performed in each spot and photon stream.
    The following attributes contain the estimated background rate.

    Attributes:
        bg (dict): background rates for the different photon streams,
            channels and background periods. Keys are `Ph_sel` objects
            and values are lists (one element per channel) of arrays (one
            element per background period) of background rates.
        bg_mean (dict): mean background rates across the entire measurement
            for the different photon streams and channels. Keys are `Ph_sel`
            objects and values are lists (one element per channel) of
            background rates.
        nperiods (int): number of periods in which timestamps are split for
            background calculation
        bg_fun (function): function used to compute the background rates
        Lim (list): each element of this list is a list of index pairs for
            `.ph_times_m[i]` for **first** and **last** photon in each period.
        Ph_p (list): each element in this list is a list of timestamps pairs
            for **first** and **last** photon of each period.
        bg_ph_sel (Ph_sel object): photon selection used by Lim and Ph_p.
            See :mod:`fretbursts.ph_sel` for details.
        Th_us (dict): thresholds in us used to select the tail of the
            interphoton delay distribution. Keys are `Ph_sel` objects
            and values are lists (one element per channel) of arrays (one
            element per background period).

    Additionlly, there are a few deprecated attributes (`bg_dd`, `bg_ad`,
    `bg_da`, `bg_aa`, `rate_dd`, `rate_ad`, `rate_da`, `rate_aa` and `rate_m`)
    which will be removed in a future version.
    Please use :attr:`Data.bg` and :attr:`Data.bg_mean` instead.

    **Burst search parameters (user input)**

    These are the parameters used to perform the burst search
    (see :meth:`burst_search`).

    Attributes:
        ph_sel (Ph_sel object): photon selection used for burst search.
            See :mod:`fretbursts.ph_sel` for details.
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
        mburst (list of Bursts objects): list Bursts() one element per channel.
            See :class:`fretbursts.phtools.burstsearch.Bursts`.

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
    ph_fields = ['ph_times_m', 'nanotimes', 'particles',
                 'A_em', 'D_em', 'A_ex', 'D_ex']

    # Attribute names containing background data.
    # The attribute `bg` is a dict with photon-selections as keys and
    # list of arrays as values. Each list contains one element per channel and
    # each array one element per background period.
    # The attributes `.Lim`  and `.Ph_p` are lists with one element per channel.
    # Each element is a lists-of-tuples (one tuple per background period).
    # These attributes do not exist before computing the background.
    bg_fields = ['bg', 'Lim', 'Ph_p']

    # Attribute names containing per-burst data.
    # Each attribute is a list (1 element per ch) of arrays (1 element
    # per burst).
    # They do not necessarly exist. For example 'naa' exists only for ALEX
    # data. Also none of them exist before performing a burst search.
    burst_fields = ['E', 'S', 'mburst', 'nd', 'na', 'nt', 'bp', 'nda', 'naa',
                    'max_rate', 'sbr']

    # Quantities (scalars or arrays) defining the current set of bursts
    burst_metadata = ['m', 'L', 'T', 'TT', 'F', 'FF', 'P', 'PP', 'rate_th',
                      'bg_bs', 'ph_sel', 'bg_corrected', 'leakage_corrected',
                      'dir_ex_corrected', 'dithering', 'fuse', 'lsb']

    # List of photon selections on which the background is computed
    _ph_streams = [Ph_sel('all'), Ph_sel(Dex='Dem'), Ph_sel(Dex='Aem'),
                   Ph_sel(Aex='Dem'), Ph_sel(Aex='Aem')]

    @property
    def ph_streams(self):
        if self.ALEX:
            return self._ph_streams
        else:
            return [Ph_sel('all'), Ph_sel(Dex='Dem'), Ph_sel(Dex='Aem')]

    def __init__(self, leakage=0., gamma=1., dir_ex=0., **kwargs):
        # Default values
        init_kw = dict(ALEX=False, _leakage=float(leakage), _gamma=float(gamma),
                       _dir_ex=float(dir_ex), _beta=1., _chi_ch=1., s=[])
        # Override with user data
        init_kw.update(**kwargs)
        DataContainer.__init__(self, **init_kw)

    def __getattr__(self, name):
        """Single-channel shortcuts for per-channel fields.

        Appending a '_' to a per-channel field avoids specifying the channel.
        For example use d.nd_ instead if d.nd[0].
        """
        msg_missing_attr = "'%s' object has no attribute '%s'" %\
                           (self.__class__.__name__, name)
        if name.startswith('_') or not name.endswith('_'):
            raise AttributeError(msg_missing_attr)

        field = name[:-1]
        try:
            value = self.__getitem__(field)
        except KeyError:
            raise AttributeError(msg_missing_attr)
        else:
            # Support lists, tuples and object with array interface
            if isinstance(value, (list, tuple)) or isarray(value):
                if len(value) == self.nch:
                    return value[0]
            raise ValueError('Name "%s" is not a per-channel field.' % field)

    def copy(self, mute=False):
        """Copy data in a new object. All arrays copied except for ph_times_m
        """
        pprint('Deep copy executed.\n', mute)
        new_d = Data(**self)  # this make a shallow copy (like a pointer)

        # Deep copy (not just reference) or array data
        for field in self.burst_fields + self.bg_fields:
            # Making sure k is defined
            if field in self:

                # Make a deepcopy of the per-channel lists
                new_d[field] = copy.deepcopy(self[field])

                # Set the attribute: new_d.k = new_d[k]
                setattr(new_d, field, new_d[field])
        return new_d

    ##
    # Methods for photon timestamps (ph_times_m) access
    #
    def ph_times_hash(self, hash_name='md5', hexdigest=True):
        """Return an hash for the timestamps arrays.
        """
        m = hashlib.new(hash_name)
        for ph in self.iter_ph_times():
            if isinstance(ph, np.ndarray):
                m.update(ph.data)
            else:
                # TODO Handle ph_times in PyTables files
                raise NotImplementedError
        if hexdigest:
            return m.hexdigest()
        else:
            return m

    @property
    def ph_data_sizes(self):
        """Array of total number of photons (ph-data) for each channel.
        """
        if not hasattr(self, '_ph_data_sizes'):
            # This works both for numpy arrays and pytables arrays
            self._ph_data_sizes = np.array([ph.shape[0] for ph in
                                            self.ph_times_m])
        return self._ph_data_sizes

    def _fix_ph_sel(self, ph_sel):
        """For non-ALEX data fix Aex to allow stable comparison."""
        msg = 'Photon selection must be of type `Ph_sel` (it was `%s` instead).'
        assert isinstance(ph_sel, Ph_sel), (msg % type(ph_sel))
        if self.ALEX or ph_sel.Dex != 'DAem':
            return ph_sel
        else:
            return Ph_sel(Dex=ph_sel.Dex, Aex='DAem')

    def _is_allph(self, ph_sel):
        """Return whether a photon selection `ph_sel` covers all photon."""
        if self.ALEX:
            return ph_sel == Ph_sel(Dex='DAem', Aex='DAem')
        else:
            return ph_sel.Dex == 'DAem'

    def get_ph_mask(self, ich=0, ph_sel=Ph_sel('all')):
        """Returns a mask for `ph_sel` photons in channel `ich`.

        The masks are either boolean arrays or slices (full or empty). In
        both cases they can be used to index the timestamps of the
        corresponding channel.

        Arguments:
            ph_sel (Ph_sel object): object defining the photon selection.
                See :mod:`fretbursts.ph_sel` for details.
        """
        assert isinstance(ich, int)

        if self._is_allph(ph_sel):
            # Note that slice(None) is equivalent to [:].
            # Also, numpy arrays are not copied when sliced.
            # So getting all photons with this mask is efficient
            # Note: the drawback is that the slice cannot be indexed
            #       (where a normal boolean array would)
            return slice(None)

        # Handle the case when A_em contains slice objects
        if isinstance(self.A_em[ich], slice):
            if self.A_em[ich] == slice(None):
                if ph_sel.Dex == 'Dem':
                    return slice(0)
                if ph_sel.Dex == 'Aem':
                    return slice(None)
            elif self.A_em[ich] == slice(0):
                if ph_sel.Dex == 'Dem':
                    return slice(None)
                if ph_sel.Dex == 'Aem':
                    return slice(0)
            else:
                msg = 'When a slice, A_em can only be slice(None) or slice(0).'
                raise NotImplementedError(msg)

        # Base selections
        elif ph_sel == Ph_sel(Dex='Dem'):
            return self.get_D_em_D_ex(ich)
        elif ph_sel == Ph_sel(Dex='Aem'):
            return self.get_A_em_D_ex(ich)
        elif ph_sel == Ph_sel(Aex='Dem'):
            return self.get_D_em(ich) * self.get_A_ex(ich)
        elif ph_sel == Ph_sel(Aex='Aem'):
            return self.get_A_em(ich) * self.get_A_ex(ich)

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
            return self.get_D_ex(ich) + self.get_A_em(ich) * self.get_A_ex(ich)

        else:
            raise ValueError('Photon selection not implemented.')

    def iter_ph_masks(self, ph_sel=Ph_sel('all')):
        """Iterator returning masks for `ph_sel` photons.

        Arguments:
            ph_sel (Ph_sel object): object defining the photon selection.
                See :mod:`fretbursts.ph_sel` for details.
        """
        for ich in range(self.nch):
            yield self.get_ph_mask(ich, ph_sel=ph_sel)

    def get_ph_times(self, ich=0, ph_sel=Ph_sel('all'), compact=False):
        """Returns the timestamps array for channel `ich`.

        This method always returns in-memory arrays, even when ph_times_m
        is a disk-backed list of arrays.

        Arguments:
            ph_sel (Ph_sel object): object defining the photon selection.
                See :mod:`fretbursts.ph_sel` for details.
            compact (bool): if True, a photon selection of only one excitation
                period is required and the timestamps are "compacted" by
                removing the "gaps" between each excitation period.
        """
        ph = self.ph_times_m[ich]

        # If not a list is an on-disk array, we need to load it
        if not isinstance(ph, np.ndarray):
            if hasattr(self, '_ph_cache') and self._ph_cache_ich == ich:
                ph = self._ph_cache
            else:
                ph = ph.read()
                self._ph_cache = ph
                self._ph_cache_ich = ich

        ph = ph[self.get_ph_mask(ich, ph_sel=ph_sel)]
        if compact:
            ph = self._ph_times_compact(ph, ph_sel)
        return ph

    def iter_ph_times(self, ph_sel=Ph_sel('all'), compact=False):
        """Iterator that returns the arrays of timestamps in `.ph_times_m`.

        Arguments:
            Same arguments as :meth:`get_ph_mask` except for `ich`.
        """
        for ich in range(self.nch):
            yield self.get_ph_times(ich, ph_sel=ph_sel, compact=compact)

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
            return self.get_D_em(ich) * self.get_D_ex(ich)
        else:
            return self.get_D_em(ich)

    def get_A_em_D_ex(self, ich=0):
        """Returns a mask of acceptor photons during donor-excitation."""
        if self.ALEX:
            return self.get_A_em(ich) * self.get_D_ex(ich)
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

        if isinstance(mask, slice) and mask == slice(None):
            ph_times_period = ph_times[period_slice]
        else:
            ph_times_period = ph_times[period_slice][mask[period_slice]]
        return ph_times_period

    def _assert_compact(self, ph_sel):
        msg = ('Option compact=True requires a photon selection \n'
               'from a single excitation period (either Dex or Aex).')
        if not self.ALEX:
            raise ValueError('Option compact=True requires ALEX data.')
        if ph_sel.Dex is not None and ph_sel.Aex is not None:
            raise ValueError(msg)

    def _excitation_width(self, ph_sel, ich=0):
        """Returns duration of alternation period outside selected excitation.
        """
        self._assert_compact(ph_sel)
        if ph_sel.Aex is None:
            excitation_range = self._D_ON_multich[ich]
        elif ph_sel.Dex is None:
            excitation_range = self._A_ON_multich[ich]
        return _excitation_width(excitation_range, self.alex_period)

    def _ph_times_compact(self, ph, ph_sel):
        """Return timestamps in one excitation period with "gaps" removed.

        It takes timestamps in the specified alternation period and removes
        gaps due to time intervals outside the alternation period selection.
        This allows to correct the photon rates distorsion due to alternation.

        Arguments:
            ph (array): timestamps array from which gaps have to be removed.
                This array **is modified in-place**.
            ph_sel (Ph_sel object): photon selection to be compacted.
                Note that only one excitation must be specified, but the
                emission can be 'Dem', 'Aem' or 'DAem'.
                See :mod:`fretbursts.ph_sel` for details.

        Returns:
            Array of timestamps in one excitation periods with "gaps" removed.
        """
        excitation_width = self._excitation_width(ph_sel)
        return _ph_times_compact(ph, self.alex_period, excitation_width)

    def _get_tuple_multich(self, name):
        """Get a n-element tuple field in multi-ch format (1 row per ch)."""
        field = np.array(self[name])
        if field.ndim == 1:
            field = np.repeat([field], self.nch, axis=0)
        return field

    @property
    def _D_ON_multich(self):
        return self._get_tuple_multich('D_ON')

    @property
    def _A_ON_multich(self):
        return self._get_tuple_multich('A_ON')

    @property
    def _det_donor_accept_multich(self):
        return self._get_tuple_multich('det_donor_accept')

    ##
    # Methods and properties for burst-data access
    #
    @property
    def num_bursts(self):
        """Array of number of bursts in each channel."""
        return np.array([bursts.num_bursts for bursts in self.mburst])

    @property
    def burst_widths(self):
        """List of arrays of burst duration in seconds. One array per channel.
        """
        return [bursts.width * self.clk_p for bursts in self.mburst]

    def burst_sizes_ich(self, ich=0, gamma=1., add_naa=False,
                        beta=1., donor_ref=True):
        """Return gamma corrected burst sizes for channel `ich`.

        If `donor_ref == True` (default) the gamma corrected burst size is
        computed according to::

            1)    nd + na / gamma

        Otherwise, if `donor_ref == False`, the gamma corrected burst size is::

            2)    nd * gamma  + na

        With the definition (1) the corrected burst size is equal to the raw
        burst size for zero-FRET or D-only bursts (that's why is `donor_ref`).
        With the definition (2) the corrected burst size is equal to the raw
        burst size for 100%-FRET bursts.

        In an ALEX measurement, use `add_naa = True` to add counts from
        AexAem stream to the returned burst size. The argument `gamma` and
        `beta` are used to correctly scale `naa` so that it become commensurate
        with the Dex corrected burst size. In particular, when using
        definition (1) (i.e. `donor_ref = True`), the total burst size is::

            (nd + na/gamma) + naa / (beta * gamma)

        Conversely, when using definition (2) (`donor_ref = False`), the
        total burst size is::

            (nd * gamma + na) + naa / beta

        Arguments:
            ich (int): the spot number, only relevant for multi-spot.
                In single-spot data there is only one channel (`ich=0`)
                so this argument may be omitted. Default 0.
            add_naa (boolean): when True, add a term for AexAem photons when
                computing burst size. Default False.
            gamma (float): coefficient for gamma correction of burst
                sizes. Default: 1. For more info see explanation above.
            beta (float): beta correction factor used for the AexAem term
                of the burst size. Default 1. If `add_naa = False` or
                measurement is not ALEX this argument is ignored.
                For more info see explanation above.

        Returns
            Array of burst sizes for channel `ich`.

        See also :meth:`fretbursts.burstlib.Data.get_naa_corrected`.
        """
        if donor_ref:
            burst_size = self.nd[ich] + self.na[ich] / gamma
        else:
            burst_size = self.nd[ich] * gamma + self.na[ich]

        if self.ALEX and add_naa:
            kws = dict(ich=ich, gamma=gamma, beta=beta, donor_ref=donor_ref)
            burst_size += self.get_naa_corrected(**kws)
        return burst_size

    def get_naa_corrected(self, ich=0, gamma=1., beta=1., donor_ref=True):
        """Return corrected naa array for channel `ich`.

        Arguments:
            ich (int): the spot number, only relevant for multi-spot.
            gamma (floats): gamma-factor to use in computing the corrected naa.
            beta (float): beta-factor to use in computing the corrected naa.
            donor_ref (bool): Select the convention for `naa` correction.
                If True (default), uses `naa / (beta * gamma)`. Otherwise,
                uses `naa / beta`. A consistent convention should be used
                for the corrected Dex burst size in order to make it
                commensurable with naa.

        See also :meth:`fretbursts.burstlib.Data.burst_sizes_ich`.
        """
        if donor_ref:
            naa_term = self.naa[ich] / (gamma * beta)
        else:
            naa_term = self.naa[ich] / beta
        return naa_term

    def burst_sizes(self, gamma=1., add_naa=False, beta=1., donor_ref=True):
        """Return gamma corrected burst sizes for all the channel.

        Compute burst sizes by calling :meth:`burst_sizes_ich` for each
        channel. See :meth:`burst_sizes_ich` for a description of the arguments.

        Returns
            List of arrays of burst sizes, one array per channel.
        """
        kwargs = dict(gamma=gamma, add_naa=add_naa, beta=beta,
                      donor_ref=donor_ref)
        bsize_list = [self.burst_sizes_ich(ich, **kwargs) for ich in
                      range(self.nch)]
        return np.array(bsize_list)

    def iter_bursts_ph(self, ich=0):
        """Iterate over (start, stop) indexes to slice photons for each burst.
        """
        for istart, istop in iter_bursts_start_stop(self.mburst[ich]):
            yield istart, istop

    def bursts_slice(self, N1=0, N2=-1):
        """Return new Data object with bursts between `N1` and `N2`
        `N1` and `N2` can be scalars or lists (one per ch).
        """
        if np.isscalar(N1): N1 = [N1] * self.nch
        if np.isscalar(N2): N2 = [N2] * self.nch
        assert len(N1) == len(N2) == self.nch
        d = Data(**self)
        d.add(mburst=[b[n1:n2].copy() for b, n1, n2 in zip(d.mburst, N1, N2)])
        d.add(nt=[nt[n1:n2] for nt, n1, n2 in zip(d.nt, N1, N2)])
        d.add(nd=[nd[n1:n2] for nd, n1, n2 in zip(d.nd, N1, N2)])
        d.add(na=[na[n1:n2] for na, n1, n2 in zip(d.na, N1, N2)])
        if self.ALEX:
            d.add(naa=[aa[n1:n2] for aa, n1, n2 in zip(d.naa, N1, N2)])
        d.calc_fret()  # recalculate fret efficiency
        return d

    def delete_burst_data(self):
        """Erase all the burst data"""
        for name in self.burst_fields + self.burst_metadata:
            if name in self:
                self.delete(name)
        for name in ('E_fitter', 'S_fitter'):
            if hasattr(self, name):
                delattr(self, name)

    ##
    # Methods for high-level data transformation
    #
    def slice_ph(self, time_s1=0, time_s2=None, s='slice'):
        """Return a new Data object with ph in [`time_s1`,`time_s2`] (seconds)

        If ALEX, this method must be called right after
        :func:`fretbursts.loader.alex_apply_periods` (with `delete_ph_t=True`)
        and before any background estimation or burst search.
        """
        if time_s2 is None:
            time_s2 = self.time_max
        if time_s2 >= self.time_max and time_s1 <= 0:
            return self.copy()
        assert time_s1 < self.time_max

        t1_clk, t2_clk = int(time_s1 / self.clk_p), int(time_s2 / self.clk_p)
        masks = [(ph >= t1_clk) * (ph < t2_clk) for ph in self.iter_ph_times()]

        new_d = Data(**self)
        for name in self.ph_fields:
            if name in self:
                new_d[name] = [a[mask] for a, mask in zip(self[name], masks)]
                setattr(new_d, name, new_d[name])
        new_d.delete_burst_data()

        # Shift timestamps to start from 0 to avoid problems with BG calc
        for ich in range(self.nch):
            ph_i = new_d.get_ph_times(ich)
            ph_i -= t1_clk
        new_d.s.append(s)

        # Delete eventual cached properties
        for attr in ['_time_min', '_time_max']:
            if hasattr(new_d, attr):
                delattr(new_d, attr)
        return new_d

    def collapse(self, update_gamma=True):
        """Returns an object with 1-spot data joining the multi-spot data.

        The argument `update_gamma` (bool, default True) allows to avoid
        recomputing gamma as the average of the original gamma. This flag
        should be always True. Set False only for testing/debugging.
        """
        dc = Data(**self)

        bursts = bslib.Bursts.merge(self.mburst, sort=False)
        # Sort by start times, and when equal by stop times
        indexsort = np.lexsort((bursts.stop, bursts.start))
        dc.add(mburst=[bursts[indexsort]])

        ich_burst = [i * np.ones(nb) for i, nb in enumerate(self.num_bursts)]
        dc.add(ich_burst=np.hstack(ich_burst)[indexsort])

        for name in self.burst_fields:
            if name in self and name is not 'mburst':
                # Concatenate arrays along axis = 0
                value = [np.concatenate(self[name])[indexsort]]
                dc.add(**{name: value})
        dc.add(nch=1)
        dc.add(_chi_ch=1.)
        # NOTE: Updating gamma has the side effect of recomputing E
        #       (and S if ALEX). We need to update gamma because, in general,
        #       gamma can be an array with a value for each ch.
        if update_gamma:
            dc._update_gamma(np.mean(self.get_gamma_array()))
        return dc

    ##
    # Utility methods
    #
    def get_params(self):
        """Returns a plain dict containing only parameters and no arrays.
        This can be used as a summary of data analysis parameters.
        Additional keys `name' and `Names` are added with values
        from `.name` and `.Name()`.
        """
        p_names = ['fname', 'clk_p', 'nch', 'ph_sel', 'L', 'm', 'F', 'P',
                   '_leakage', '_dir_ex', '_gamma', 'bg_time_s',
                   'T', 'rate_th',
                   'bg_corrected', 'leakage_corrected', 'dir_ex_corrected',
                   'dithering', '_chi_ch', 's', 'ALEX']
        p_dict = dict(self)
        for name in p_dict.keys():
            if name not in p_names:
                p_dict.pop(name)
        p_dict.update(name=self.name, Name=self.Name(), bg_mean=self.bg_mean,
                      nperiods=self.nperiods)
        return p_dict

    def expand(self, ich=0, alex_naa=False, width=False):
        """Return per-burst D and A sizes (nd, na) and their background counts.

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
        w = self.mburst[ich].width * self.clk_p
        bg_a = self.bg[Ph_sel(Dex='Aem')][ich][period] * w
        bg_d = self.bg[Ph_sel(Dex='Dem')][ich][period] * w
        res = [self.nd[ich], self.na[ich]]
        if self.ALEX and alex_naa:
            bg_aa = self.bg[Ph_sel(Aex='Aem')][ich][period] * w
            res.extend([self.naa[ich], bg_d, bg_a, bg_aa])
        else:
            res.extend([bg_d, bg_a])
        if width:
            res.append(w)
        return res

    @property
    def time_max(self):
        """The last recorded time in seconds."""
        if not hasattr(self, '_time_max'):
            self._time_max = self._time_reduce(last=True, func=max)
        return self._time_max

    @property
    def time_min(self):
        """The first recorded time in seconds."""
        if not hasattr(self, '_time_min'):
            self._time_min = self._time_reduce(last=False, func=min)
        return self._time_min

    def _time_reduce(self, last=True, func=max):
        """Return first or last timestamp per-ch, reduced with `func`.
        """
        idx = -1 if last else 0

        # Get either ph_times_m or ph_times_t
        ph_times = None
        for ph_times_name in ['ph_times_m', 'ph_times_t']:
            try:
                ph_times = self[ph_times_name]
            except KeyError:
                pass

        if ph_times is not None:
            # This works with both numpy arrays and pytables arrays
            time = func(t[idx] for t in ph_times if t.shape[0] > 0)
        elif 'mburst' in self:
            if last:
                time = func(bursts[idx].stop for bursts in self.mburst)
            else:
                time = func(bursts[idx].start for bursts in self.mburst)
        else:
            raise ValueError("No timestamps or bursts found.")

        return time * self.clk_p

    def ph_in_bursts_mask_ich(self, ich=0, ph_sel=Ph_sel('all')):
        """Return mask of all photons inside bursts for channel `ich`.

        Returns
            Boolean array for photons in channel `ich` and photon
            selection `ph_sel` that are inside any burst.
        """
        bursts_mask = ph_in_bursts_mask(self.ph_data_sizes[ich],
                                        self.mburst[ich])
        if self._is_allph(ph_sel):
            return bursts_mask
        else:
            ph_sel_mask = self.get_ph_mask(ich=ich, ph_sel=ph_sel)
            return ph_sel_mask * bursts_mask

    def ph_in_bursts_ich(self, ich=0, ph_sel=Ph_sel('all')):
        """Return timestamps of photons inside bursts for channel `ich`.

        Returns
            Array of photon timestamps in channel `ich` and photon
            selection `ph_sel` that are inside any burst.
        """
        ph_all = self.get_ph_times(ich=ich)
        bursts_mask = self.ph_in_bursts_mask_ich(ich, ph_sel)
        return ph_all[bursts_mask]

    ##
    # Background analysis methods
    #
    def _obsolete_bg_attr(self, attrname, ph_sel):
        print('The Data.%s attribute is deprecated. Please use '
              'Data.bg(%s) instead.' % (attrname, repr(ph_sel)))
        bg_attrs = ('bg_dd', 'bg_ad', 'bg_da', 'bg_aa')
        bg_mean_attrs = ('rate_m', 'rate_dd', 'rate_ad', 'rate_da', 'rate_aa')
        assert attrname in bg_attrs or attrname in bg_mean_attrs
        if attrname in bg_attrs:
            bg_field = 'bg'
        elif attrname in bg_mean_attrs:
            bg_field = 'bg_mean'
        try:
            value = getattr(self, bg_field)[ph_sel]
        except AttributeError as e:
            # This only happens when trying to access 'bg' because
            # 'bg_mean' raises RuntimeError when missing.
            msg = 'No attribute `%s` found. Please compute background first.'
            raise_from(RuntimeError(msg % bg_field), e)
        return value

    @property
    def rate_m(self):
        return self._obsolete_bg_attr('rate_m', Ph_sel('all'))

    @property
    def rate_dd(self):
        return self._obsolete_bg_attr('rate_dd', Ph_sel(Dex='Dem'))

    @property
    def rate_ad(self):
        return self._obsolete_bg_attr('rate_ad', Ph_sel(Dex='Aem'))

    @property
    def rate_da(self):
        return self._obsolete_bg_attr('rate_da', Ph_sel(Aex='Dem'))

    @property
    def rate_aa(self):
        return self._obsolete_bg_attr('rate_aa', Ph_sel(Aex='Aem'))

    @property
    def bg_dd(self):
        return self._obsolete_bg_attr('bg_dd', Ph_sel(Dex='Dem'))

    @property
    def bg_ad(self):
        return self._obsolete_bg_attr('bg_ad', Ph_sel(Dex='Aem'))

    @property
    def bg_da(self):
        return self._obsolete_bg_attr('bg_da', Ph_sel(Aex='Dem'))

    @property
    def bg_aa(self):
        return self._obsolete_bg_attr('bg_aa', Ph_sel(Aex='Aem'))

    def calc_bg_cache(self, fun, time_s=60, tail_min_us=500, F_bg=2,
                      error_metrics=None, fit_allph=True,
                      recompute=False):
        """Compute time-dependent background rates for all the channels.

        This version is the cached version of :meth:`calc_bg`.
        This method tries to load the background data from a cache file.
        If a saved background data is not found, it computes
        the background and stores it to disk.

        The arguments are the same as :meth:`calc_bg` with the only addition
        of `recompute` (bool) to force a background recomputation even if
        a cached version is found.

        Form more details on the other arguments see :meth:`calc_bg`.
        """
        bg_cache.calc_bg_cache(self, fun, time_s=time_s,
                               tail_min_us=tail_min_us, F_bg=F_bg,
                               error_metrics=error_metrics, fit_allph=fit_allph,
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
                    th_us[ich] = 1e6 * F_bg / bg_rate
            Th_us[ph_sel] = th_us
        # Save the input used to generate Th_us
        self.add(bg_auto_th_us0=tail_min_us0, bg_auto_F_bg=F_bg)
        return Th_us

    def _get_bg_th_arrays(self, tail_min_us, nperiods):
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
            raise ValueError('Wrong tail_min_us length (%d).' %
                             len(tail_min_us))
        th_us = {}
        for i, key in enumerate(self.ph_streams):
            th_us[key] = np.ones(nperiods) * tail_min_us[i]
        # Save the input used to generate Th_us
        self.add(bg_th_us_user=tail_min_us)
        return th_us

    def _clean_bg_data(self):
        """Remove background fields specific of only one fit type.

        Computing background with manual or 'auto' threshold results in
        different sets of attributes being saved. This method removes these
        attributes and should be called before recomputing the background
        to avoid having old stale attributes of a previous background fit.
        """
        # Attributes specific of manual or 'auto' bg fit
        field_list = ['bg_auto_th_us0', 'bg_auto_F_bg', 'bg_th_us_user']
        for field in field_list:
            if field in self:
                self.delete(field)
        if hasattr(self, '_bg_mean'):
            delattr(self, '_bg_mean')

    def _get_num_periods(self, time_s):
        """Return the number of periods using `time_s` as period duration.
        """
        duration = self.time_max - self.time_min
        # Take the ceil to have at least 1 periods
        nperiods = np.ceil(duration / time_s)
        # Discard last period if negligibly small to avoid problems with
        # background fit with very few photons.
        if nperiods > 1:
            last_period = self.time_max - time_s * (nperiods - 1)
            # Discard last period if smaller than 3% of the bg period
            if last_period < time_s * 0.03:
                nperiods -= 1
        return int(nperiods)

    def calc_bg(self, fun, time_s=60, tail_min_us=500, F_bg=2,
                error_metrics=None, fit_allph=True):
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
            fit_allph (bool): if True (default) the background for the
                all-photon is fitted. If False it is computed as the sum of
                backgrounds in all the other streams.

        The background estimation functions are defined in the module
        `background` (conventionally imported as `bg`).

        Example:
            Compute background with `bg.exp_fit` (inter-photon delays MLE
            tail fitting), every 30s, with automatic tail-threshold::

               d.calc_bg(bg.exp_fit, time_s=20, tail_min_us='auto')

        Returns:
            None, all the results are saved in the object itself.
        """
        pprint(" - Calculating BG rates ... ")
        self._clean_bg_data()
        kwargs = dict(clk_p=self.clk_p, error_metrics=error_metrics)
        nperiods = self._get_num_periods(time_s)
        streams_noall = [s for s in self.ph_streams if s != Ph_sel('all')]

        bg_auto_th = tail_min_us == 'auto'
        if bg_auto_th:
            tail_min_us0 = 250
            self.add(bg_auto_th_us0=tail_min_us0, bg_auto_F_bg=F_bg)
            auto_th_kwargs = dict(clk_p=self.clk_p, tail_min_us=tail_min_us0)
            th_us = {}
            for key in self.ph_streams:
                th_us[key] = np.zeros(nperiods)
        else:
            th_us = self._get_bg_th_arrays(tail_min_us, nperiods)

        Lim, Ph_p = [], []
        BG, BG_err = [], []
        Th_us = []
        for ich, ph_ch in enumerate(self.iter_ph_times()):
            masks = {sel: self.get_ph_mask(ich, ph_sel=sel)
                     for sel in self.ph_streams}

            bins = ((np.arange(nperiods + 1) * time_s + self.time_min) /
                    self.clk_p)
            # Note: histogram bins are half-open, e.g. [a, b)
            counts, _ = np.histogram(ph_ch, bins=bins)
            lim, ph_p = [], []
            bg = {sel: np.zeros(nperiods) for sel in self.ph_streams}
            bg_err = {sel: np.zeros(nperiods) for sel in self.ph_streams}
            i1 = 0
            for ip in range(nperiods):
                i0 = i1
                i1 += counts[ip]
                lim.append((i0, i1 - 1))
                ph_p.append((ph_ch[i0], ph_ch[i1 - 1]))
                ph_i = ph_ch[i0:i1]

                if fit_allph:
                    sel = Ph_sel('all')
                    if bg_auto_th:
                        _bg, _ = fun(ph_i, **auto_th_kwargs)
                        th_us[sel][ip] = 1e6 * F_bg / _bg
                    bg[sel][ip], bg_err[sel][ip] = \
                        fun(ph_i, tail_min_us=th_us[sel][ip], **kwargs)

                for sel in streams_noall:
                    # This supports cases of D-only or A-only timestamps
                    # where self.A_em[ich] is a bool and not a bool-array
                    # In this case, the mask of either DexDem or DexAem is
                    # slice(None) (all-elements selection).
                    if isinstance(masks[sel], slice):
                        if masks[sel] == slice(None):
                            bg[sel][ip] = bg[Ph_sel('all')][ip]
                            bg_err[sel][ip] = bg_err[Ph_sel('all')][ip]
                        continue
                    else:
                        ph_i_sel = ph_i[masks[sel][i0:i1]]

                    if ph_i_sel.size > 0:
                        if bg_auto_th:
                            _bg, _ = fun(ph_i_sel, **auto_th_kwargs)
                            th_us[sel][ip] = 1e6 * F_bg / _bg
                        bg[sel][ip], bg_err[sel][ip] = \
                            fun(ph_i_sel, tail_min_us=th_us[sel][ip], **kwargs)

            if not fit_allph:
                bg[Ph_sel('all')] += sum(bg[s] for s in streams_noall)
                bg_err[Ph_sel('all')] += sum(bg_err[s] for s in streams_noall)
            Lim.append(lim)
            Ph_p.append(ph_p)
            BG.append(bg)
            BG_err.append(bg_err)
            Th_us.append(th_us)

        # Make Dict Of Lists (DOL) from Lists of Dicts
        BG_dol, BG_err_dol, Th_us_dol = {}, {}, {}
        for sel in self.ph_streams:
            BG_dol[sel] = [bg_ch[sel] for bg_ch in BG]
            BG_err_dol[sel] = [err_ch[sel] for err_ch in BG_err]
            Th_us_dol[sel] = [th_ch[sel] for th_ch in Th_us]

        self.add(bg=BG_dol, bg_err=BG_err_dol, bg_th_us=Th_us_dol,
                 Lim=Lim, Ph_p=Ph_p,
                 bg_fun=fun, bg_fun_name=fun.__name__,
                 bg_time_s=time_s, bg_ph_sel=Ph_sel('all'),
                 bg_auto_th=bg_auto_th,  # bool, True if the using auto-threshold
                 )
        pprint("[DONE]\n")

    @property
    def nperiods(self):
        return len(self.bg[Ph_sel('all')][0])

    @property
    def bg_mean(self):
        if 'bg' not in self:
            raise RuntimeError('No background found, compute it first.')
        if not hasattr(self, '_bg_mean'):
            self._bg_mean = {k: [bg_ch.mean() for bg_ch in bg_ph_sel]
                             for k, bg_ph_sel in self.bg.items()}
        return self._bg_mean

    def recompute_bg_lim_ph_p(self, ph_sel, mute=False):
        """Recompute self.Lim and selp.Ph_p relative to ph selection `ph_sel`
        `ph_sel` is a Ph_sel object selecting the timestamps in which self.Lim
        and self.Ph_p are being computed.
        """
        ph_sel = self._fix_ph_sel(ph_sel)
        if self.bg_ph_sel == ph_sel:
            return

        pprint(" - Recomputing background limits for %s ... " %
               str(ph_sel), mute)
        bg_time_clk = self.bg_time_s / self.clk_p
        Lim, Ph_p = [], []
        for ph_ch, lim in zip(self.iter_ph_times(ph_sel), self.Lim):
            bins = np.arange(self.nperiods + 1) * bg_time_clk
            # Note: histogram bins are half-open, e.g. [a, b)
            counts, _ = np.histogram(ph_ch, bins=bins)
            lim, ph_p = [], []
            i1 = 0
            for ip in range(self.nperiods):
                i0 = i1
                i1 += counts[ip]
                lim.append((i0, i1 - 1))
                ph_p.append((ph_ch[i0], ph_ch[i1-1]))
            Lim.append(lim)
            Ph_p.append(ph_p)
        self.add(Lim=Lim, Ph_p=Ph_p, bg_ph_sel=ph_sel)
        pprint("[DONE]\n", mute)

    ##
    # Burst analysis methods
    #
    def _calc_burst_period(self):
        """Compute for each burst the "background period" `bp`.
        Background periods are the time intervals on which the BG is computed.
        """
        P = []
        for b, lim in zip(self.mburst, self.Lim):
            p = zeros(b.num_bursts, dtype=np.int16)
            if b.num_bursts > 0:
                istart = b.istart
                for i, (l0, l1) in enumerate(lim):
                    p[(istart >= l0) * (istart <= l1)] = i
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
        ph_sel = self._fix_ph_sel(ph_sel)
        if ph_sel in self.ph_streams:
            return self.bg[ph_sel]
        elif ph_sel == Ph_sel(Dex='DAem'):
            sel = Ph_sel(Dex='Dem'), Ph_sel(Dex='Aem')
            bg = [b1 + b2 for b1, b2 in zip(self.bg[sel[0]], self.bg[sel[1]])]
        elif ph_sel == Ph_sel(Aex='DAem'):
            sel = Ph_sel(Aex='Dem'), Ph_sel(Aex='Aem')
            bg = [b1 + b2 for b1, b2 in zip(self.bg[sel[0]], self.bg[sel[1]])]
        elif ph_sel == Ph_sel(Dex='Dem', Aex='Dem'):
            sel = Ph_sel(Dex='Dem'), Ph_sel(Aex='Dem')
            bg = [b1 + b2 for b1, b2 in zip(self.bg[sel[0]], self.bg[sel[1]])]
        elif ph_sel == Ph_sel(Dex='Aem', Aex='Aem'):
            sel = Ph_sel(Dex='Aem'), Ph_sel(Aex='Aem')
            bg = [b1 + b2 for b1, b2 in zip(self.bg[sel[0]], self.bg[sel[1]])]
        elif ph_sel == Ph_sel(Dex='DAem', Aex='Aem'):
            sel = (Ph_sel(Dex='Dem'), Ph_sel(Dex='Aem'), Ph_sel(Aex='Aem'))
            bg = [b1 + b2 + b3 for b1, b2, b3 in
                  zip(self.bg[sel[0]], self.bg[sel[1]], self.bg[sel[2]])]
        else:
            raise NotImplementedError('Photon selection %s not implemented.' %
                                      ph_sel)
        return bg


    def _calc_T(self, m, P, F=1., ph_sel=Ph_sel('all'), c=-1):
        """If P is None use F, otherwise uses both P *and* F (F defaults to 1).

        When P is None, compute the time lag T for burst search according to::

            T = (m - 1 - c) / (F * bg_rate)

        """
        # Regardless of F and P sizes, FF and PP are arrays with size == nch
        FF = self._param_as_mch_array(F)
        PP = self._param_as_mch_array(P)
        if P is None:
            # NOTE: the following lambda ignores Pi
            find_T = lambda m, Fi, Pi, bg: (m - 1 - c) / (bg * Fi)
        else:
            if F != 1:
                print("WARNING: BS prob. th. with modified BG rate (F=%.1f)"
                      % F)
            find_T = lambda m, Fi, Pi, bg: find_optimal_T_bga(bg*Fi, m, 1-Pi)
        TT, T, rate_th = [], [], []
        bg_bs = self.bg_from(ph_sel)
        for bg_ch, F_ch, P_ch in zip(bg_bs, FF, PP):
            # All "T" are in seconds
            Tch = find_T(m, F_ch, P_ch, bg_ch)
            TT.append(Tch)
            T.append(Tch.mean())
            rate_th.append(np.mean(m / Tch))
        self.add(TT=TT, T=T, bg_bs=bg_bs, FF=FF, PP=PP, F=F, P=P,
                 rate_th=rate_th)

    def _burst_search_rate(self, m, L, min_rate_cps, c=-1, ph_sel=Ph_sel('all'),
                           compact=False, index_allph=True, verbose=True,
                           pure_python=False):
        """Compute burst search using a fixed minimum photon rate.

        The burst starts when, for `m` consecutive photons::

            (m - 1 - c) / (t[last] - t[first]) >= min_rate_cps

        Arguments:
            min_rate_cps (float or array): minimum photon rate for burst start
                if array if one value per channel.
        """
        bsearch = _get_bsearch_func(pure_python=pure_python)

        Min_rate_cps = self._param_as_mch_array(min_rate_cps)
        mburst = []
        T_clk = (m - 1 - c) / Min_rate_cps / self.clk_p
        for ich, t_clk in enumerate(T_clk):
            ph_bs = ph = self.get_ph_times(ich=ich, ph_sel=ph_sel)
            if compact:
                ph_bs = self._ph_times_compact(ph, ph_sel)
            label = '%s CH%d' % (ph_sel, ich + 1) if verbose else None
            burstarray = bsearch(ph_bs, L, m, t_clk, label=label, verbose=verbose)
            if burstarray.size > 1:
                bursts = bslib.Bursts(burstarray)
                if compact:
                    bursts.recompute_times(ph, out=bursts)
            else:
                bursts = bslib.Bursts.empty()
            mburst.append(bursts)
        self.add(mburst=mburst, rate_th=Min_rate_cps, T=T_clk * self.clk_p)
        if ph_sel != Ph_sel('all') and index_allph:
            self._fix_mburst_from(ph_sel=ph_sel)

    def _burst_search_TT(self, m, L, ph_sel=Ph_sel('all'), verbose=True,
                         compact=False, index_allph=True, pure_python=False,
                         mute=False):
        """Compute burst search with params `m`, `L` on ph selection `ph_sel`

        Requires the list of arrays `self.TT` with the max time-thresholds in
        the different burst periods for each channel (use `._calc_T()`).
        """
        bsearch = _get_bsearch_func(pure_python=pure_python)

        self.recompute_bg_lim_ph_p(ph_sel=ph_sel, mute=mute)
        MBurst = []
        label = ''
        for ich, T in enumerate(self.TT):
            ph_bs = ph = self.get_ph_times(ich=ich, ph_sel=ph_sel)
            if compact:
                ph_bs = self._ph_times_compact(ph, ph_sel)
            burstarray_ch_list = []
            Tck = T / self.clk_p

            for ip, (l0, l1) in enumerate(self.Lim[ich]):
                if verbose:
                    label = '%s CH%d-%d' % (ph_sel, ich + 1, ip)
                burstarray = bsearch(ph_bs, L, m, Tck[ip], slice_=(l0, l1 + 1),
                                     label=label, verbose=verbose)
                if burstarray.size > 1:
                    burstarray_ch_list.append(burstarray)

            if len(burstarray_ch_list) > 0:
                data = np.vstack(burstarray_ch_list)
                bursts = bslib.Bursts(data)
                if compact:
                    bursts.recompute_times(ph, out=bursts)
            else:
                bursts = bslib.Bursts.empty()
            MBurst.append(bursts)

        self.add(mburst=MBurst)
        if ph_sel != Ph_sel('all') and index_allph:
            # Convert the burst data to be relative to ph_times_m.
            # Convert both Lim/Ph_p and mburst, as they are both needed
            # to compute `.bp`.
            self.recompute_bg_lim_ph_p(ph_sel=Ph_sel('all'), mute=mute)
            self._fix_mburst_from(ph_sel=ph_sel, mute=mute)

    def _fix_mburst_from(self, ph_sel, mute=False):
        """Convert burst data from any ph_sel to 'all' timestamps selection.
        """
        assert isinstance(ph_sel, Ph_sel) and not self._is_allph(ph_sel)
        pprint(' - Fixing  burst data to refer to ph_times_m ... ', mute)

        for bursts, mask in zip(self.mburst,
                                self.iter_ph_masks(ph_sel=ph_sel)):
            bursts.recompute_index_expand(mask, out=bursts)

        pprint('[DONE]\n', mute)

    def burst_search(self, L=None, m=10, F=6., P=None, min_rate_cps=None,
                     ph_sel=Ph_sel('all'), compact=False, index_allph=True,
                     c=-1, computefret=True, max_rate=False, dither=False,
                     pure_python=False, verbose=False, mute=False):
        """Performs a burst search with specified parameters.

        This method performs a sliding-window burst search without
        binning the timestamps. The burst starts when the rate of `m`
        photons is above a minimum rate, and stops when the rate falls below
        the threshold. The result of the burst search is stored in the
        `mburst` attribute (a list of Bursts objects, one per channel)
        containing start/stop times and indexes. By default, after burst
        search, this method computes donor and acceptor counts, it applies
        burst corrections (background, leakage, etc...) and computes
        E (and S in case of ALEX). You can skip these steps by passing
        `computefret=False`.

        The minimum rate can be explicitly specified with the `min_rate_cps`
        argument, or computed as a function of the background rate with the
        `F` argument.

        Parameters:
            m (int): number of consecutive photons used to compute the
                photon rate. Typical values 5-20. Default 10.
            L (int or None): minimum number of photons in burst. If None
                (default) L = m is used.
            F (float): defines how many times higher than the background rate
                is the minimum rate used for burst search
                (`min rate = F * bg. rate`), assuming that `P = None` (default).
                Typical values are 3-9. Default 6.
            P (float): threshold for burst detection expressed as a
                probability that a detected bursts is not due to a Poisson
                background. If not None, `P` overrides `F`. Note that the
                background process is experimentally super-Poisson so this
                probability is not physically very meaningful. Using this
                argument is discouraged.
            min_rate_cps (float or list/array): minimum rate in cps for burst
                start. If not None, it has the precedence over `P` and `F`.
                If non-scalar, contains one rate per each multispot channel.
                Typical values range from 20e3 to 100e3.
            ph_sel (Ph_sel object): defines the "photon selection" (or stream)
                to be used for burst search. Default: all photons.
                See :mod:`fretbursts.ph_sel` for details.
            compact (bool): if True, a photon selection of only one excitation
                period is required and the timestamps are "compacted" by
                removing the "gaps" between each excitation period.
            index_allph (bool): if True (default), the indexes of burst start
                and stop (`istart`, `istop`) are relative to the full
                timestamp array. If False, the indexes are relative to
                timestamps selected by the `ph_sel` argument.
            c (float): correction factor used in the rate vs time-lags relation.
                `c` affects the computation of the burst-search parameter `T`.
                When `F` is not None, `T = (m - 1 - c) / (F * bg_rate)`.
                When using `min_rate_cps`, `T = (m - 1 - c) / min_rate_cps`.
            computefret (bool): if True (default) compute donor and acceptor
                counts, apply corrections (background, leakage, direct
                excitation) and compute E (and S). If False, skip all these
                steps and stop just after the initial burst search.
            max_rate (bool): if True compute the max photon rate inside each
                burst using the same `m` used for burst search. If False
                (default) skip this step.
            dither (bool): if True applies dithering corrections to burst
                counts. Default False. See :meth:`Data.dither`.
            pure_python (bool): if True, uses the pure python functions even
                when optimized Cython functions are available.

        Note:
            when using `P` or `F` the background rates are needed, so
            `.calc_bg()` must be called before the burst search.

        Example:
            d.burst_search(m=10, F=6)

        Returns:
            None, all the results are saved in the `Data` object.
        """
        ph_sel = self._fix_ph_sel(ph_sel)
        if compact:
            self._assert_compact(ph_sel)
        pprint(" - Performing burst search (verbose=%s) ..." % verbose, mute)
        # Erase any previous burst data
        self.delete_burst_data()
        if L is None:
            L = m
        if min_rate_cps is not None:
            # Saves rate_th in self
            self._burst_search_rate(m=m, L=L, min_rate_cps=min_rate_cps, c=c,
                                    ph_sel=ph_sel, compact=compact,
                                    index_allph=index_allph,
                                    verbose=verbose, pure_python=pure_python)
        else:
            # Compute TT, saves P and F in self
            self._calc_T(m=m, P=P, F=F, ph_sel=ph_sel, c=c)
            # Use TT and compute mburst
            self._burst_search_TT(L=L, m=m, ph_sel=ph_sel, compact=compact,
                                  index_allph=index_allph, verbose=verbose,
                                  pure_python=pure_python, mute=mute)
        pprint("[DONE]\n", mute)

        pprint(" - Calculating burst periods ...", mute)
        self._calc_burst_period()                       # writes bp
        pprint("[DONE]\n", mute)

        # (P, F) or rate_th are saved in _calc_T() or _burst_search_rate()
        self.add(m=m, L=L, ph_sel=ph_sel)

        # The correction flags are both set here and in calc_ph_num() so that
        # they are always consistent. Case 1: we perform only burst search
        # (with no call to calc_ph_num). Case 2: we re-call calc_ph_num()
        # without doing a new burst search
        self.add(bg_corrected=False, leakage_corrected=False,
                 dir_ex_corrected=False, dithering=False)
        self._burst_search_postprocess(
            computefret=computefret, max_rate=max_rate, dither=dither,
            pure_python=pure_python, mute=mute)

    def _burst_search_postprocess(self, computefret, max_rate, dither,
                                  pure_python, mute):
        if computefret:
            pprint(" - Counting D and A ph and calculating FRET ... \n", mute)
            self.calc_fret(count_ph=True, corrections=True, dither=dither,
                           mute=mute, pure_python=pure_python)
            pprint("   [DONE Counting D/A]\n", mute)
        if max_rate:
            pprint(" - Computing max rates in burst ...", mute)
            self.calc_max_rate(m=self.m)
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
            nt = [b.counts.astype(float) if b.num_bursts > 0 else np.array([])
                  for b in self.mburst]
            A_em = [self.get_A_em(ich) for ich in range(self.nch)]
            if isinstance(A_em[0], slice):
                # This is to support the case of A-only or D-only data
                n0 = [np.zeros(mb.num_bursts) for mb in self.mburst]
                if A_em[0] == slice(None):
                    nd, na = n0, nt    # A-only case
                elif A_em[0] == slice(0):
                    nd, na = nt, n0    # D-only case
            else:
                # This is the usual case with photons in both D and A channels
                na = mch_count_ph_in_bursts(self.mburst, A_em)
                nd = [t - a for t, a in zip(nt, na)]
            assert (nt[0] == na[0] + nd[0]).all()
        if self.ALEX:
            # The "new style" would be:
            #Mask = [m for m in self.iter_ph_masks(Ph_sel(Dex='Dem'))]
            Mask = [d_em * d_ex for d_em, d_ex in zip(self.D_em, self.D_ex)]
            nd = mch_count_ph_in_bursts(self.mburst, Mask)

            Mask = [a_em * d_ex for a_em, d_ex in zip(self.A_em, self.D_ex)]
            na = mch_count_ph_in_bursts(self.mburst, Mask)

            Mask = [a_em * a_ex for a_em, a_ex in zip(self.A_em, self.A_ex)]
            naa = mch_count_ph_in_bursts(self.mburst, Mask)
            self.add(naa=naa)

            if alex_all:
                Mask = [d_em * a_ex for d_em, a_ex in zip(self.D_em, self.A_ex)]
                nda = mch_count_ph_in_bursts(self.mburst, Mask)
                self.add(nda=nda)

            nt = [d + a + aa for d, a, aa in zip(nd, na, naa)]
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
        if ms < 0:
            return self
        mburst = mch_fuse_bursts(self.mburst, ms=ms, clk_p=self.clk_p)
        new_d = Data(**self)
        for k in ['E', 'S', 'nd', 'na', 'naa', 'nt', 'lsb', 'bp']:
            if k in new_d:
                new_d.delete(k)
        new_d.add(bg_corrected=False, leakage_corrected=False,
                  dir_ex_corrected=False, dithering=False)
        new_d.add(mburst=mburst, fuse=ms)
        if 'bg' in new_d:
            new_d._calc_burst_period()
        if process:
            pprint(" - Counting D and A ph and calculating FRET ... \n", mute)
            new_d.calc_fret(count_ph=True, corrections=True,
                            dither=self.dithering, mute=mute)
            pprint("   [DONE Counting D/A and FRET]\n", mute)
        return new_d


    ##
    # Burst selection and filtering
    #
    def select_bursts(self, filter_fun, negate=False, computefret=True,
                      args=None, **kwargs):
        """Return an object with bursts filtered according to `filter_fun`.

        This is the main method to select bursts according to different
        criteria. The selection rule is defined by the selection function
        `filter_fun`. FRETBursts provides a several predefined selection
        functions see :ref:`burst_selection`. New selection
        functions can be defined and passed to this method to implement
        arbitrary selection rules.

        Arguments:
            filter_fun (fuction): function used for burst selection
            negate (boolean): If True, negates (i.e. take the complementary)
                of the selection returned by `filter_fun`. Default `False`.
            computefret (boolean): If True (default) recompute donor and
                acceptor counts, corrections and FRET quantities (i.e. E, S)
                in the new returned object.
            args (tuple or None): positional arguments for `filter_fun()`

        kwargs:
            Additional keyword arguments passed to `filter_fun()`.

        Returns:
            A new :class:`Data` object containing only the selected bursts.

        Note:
            In order to save RAM, the timestamp arrays (`ph_times_m`)
            of the new Data() points to the same arrays of the original
            Data(). Conversely, all the bursts data (`mburst`, `nd`, `na`,
            etc...) are new distinct objects.
        """
        Masks, str_sel = self.select_bursts_mask(filter_fun, negate=negate,
                                                 return_str=True, args=args,
                                                 **kwargs)
        d_sel = self.select_bursts_mask_apply(Masks, computefret=computefret,
                                              str_sel=str_sel)
        return d_sel

    def select_bursts_mask(self, filter_fun, negate=False, return_str=False,
                           args=None, **kwargs):
        """Returns mask arrays to select bursts according to `filter_fun`.

        The function `filter_fun` is called to compute the mask arrays for
        each channel.

        This method is useful when you want to apply a selection from one
        object to a second object. Otherwise use :meth:`Data.select_bursts`.

        Arguments:
            filter_fun (fuction): function used for burst selection
            negate (boolean): If True, negates (i.e. take the complementary)
                of the selection returned by `filter_fun`. Default `False`.
            return_str: if True return, for each channel, a tuple with
                a bool array and a string that can be added to the measurement
                name to indicate the selection. If False returns only
                the bool array. Default False.
            args (tuple or None): positional arguments for `filter_fun()`

        kwargs:
            Additional keyword arguments passed to `filter_fun()`.

        Returns:
            A list of boolean arrays (one per channel) that define the burst
            selection. If `return_str` is True returns a list of tuples, where
            each tuple is a bool array and a string.

        See also:
            :meth:`Data.select_bursts`, :meth:`Data.select_bursts_mask_apply`
        """
        # Create the list of bool masks for the bursts selection
        if args is None:
            args = tuple()
        M = [filter_fun(self, i, *args, **kwargs) for i in range(self.nch)]
        # Make sure the selection function has the right return signature
        msg = 'The second argument returned by `%s` must be a string.'
        assert np.all([isinstance(m[1], str) for m in M]), msg % filter_fun
        # Make sure all boolean masks have the right size
        msg = ("The size of boolean masks returned by `%s` needs to match "
               "the number of bursts.")
        assert np.all([m[0].size == n for m, n in zip(M, self.num_bursts)]), (
            msg % filter_fun)
        Masks = [-m[0] if negate else m[0] for m in M]
        str_sel = M[0][1]
        if return_str:
            return Masks, str_sel
        else:
            return Masks

    def select_bursts_mask_apply(self, masks, computefret=True, str_sel=''):
        """Returns a new Data object with bursts selected according to `masks`.

        This method select bursts using a list of boolean arrays as input.
        Since the user needs to create the boolean arrays first, this method
        is useful when experimenting with new selection criteria that don't
        have a dedicated selection function. Usually, however, it is easier
        to select bursts through :meth:`Data.select_bursts` (using a
        selection function).

        Arguments:
            masks (list of arrays): each element in this list is a boolean
                array that selects bursts in a channel.
            computefret (boolean): If True (default) recompute donor and
                acceptor counts, corrections and FRET quantities (i.e. E, S)
                in the new returned object.

        Returns:
            A new :class:`Data` object containing only the selected bursts.

        Note:
            In order to save RAM, the timestamp arrays (`ph_times_m`)
            of the new Data() points to the same arrays of the original
            Data(). Conversely, all the bursts data (`mburst`, `nd`, `na`,
            etc...) are new distinct objects.

        See also:
            :meth:`Data.select_bursts`, :meth:`Data.select_mask`
        """
        # Attributes of ds point to the same objects of self
        ds = Data(**self)

        ##Copy the per-burst fields that must be filtered
        used_fields = [field for field in Data.burst_fields if field in self]
        for name in used_fields:

            # Recreate the current attribute as a new list to avoid modifying
            # the old list that is also in the original object.
            # The list is initialized with empty arrays because this is the
            # valid value when a ch has no bursts.
            ds.add(**{name: [np.array([])] * self.nch})

            # Assign the new data
            for ich, mask in enumerate(masks):
                if self[name][ich].size == 0:
                    continue  # -> no bursts in ch
                # Note that boolean masking implies numpy array copy
                # On the contrary slicing only makes a new view of the array
                ds[name][ich] = self[name][ich][mask]

        # Recompute E and S
        if computefret:
            ds.calc_fret(count_ph=False)
        # Add the annotation about the filter function
        ds.s = list(self.s + [str_sel])  # using append would modify also self
        return ds

    ##
    # Burst corrections
    #
    def background_correction(self, relax_nt=False, mute=False):
        """Apply background correction to burst sizes (nd, na,...)
        """
        if self.bg_corrected:
            return -1
        pprint("   - Applying background correction.\n", mute)
        self.add(bg_corrected=True)
        for ich, bursts in enumerate(self.mburst):
            if bursts.num_bursts == 0:
                continue  # if no bursts skip this ch
            period = self.bp[ich]
            nd, na, bg_d, bg_a, width = self.expand(ich, width=True)
            nd -= bg_d
            na -= bg_a
            if relax_nt:
                # This does not guarantee that nt = nd + na
                self.nt[ich] -= self.bg_from(Ph_sel('all'))[ich][period] * width
            else:
                self.nt[ich] = nd + na
            if self.ALEX:
                bg_aa = self.bg_from(Ph_sel(Aex='Aem'))
                self.naa[ich] -= bg_aa[ich][period] * width
                if 'nda' in self:
                    bg_da = self.bg_from(Ph_sel(Aex='Dem'))
                    self.nda[ich] -= bg_da[ich][period] * width
                self.nt[ich] += self.naa[ich]

    def leakage_correction(self, mute=False):
        """Apply leakage correction to burst sizes (nd, na,...)
        """
        if self.leakage_corrected:
            return -1
        pprint("   - Applying leakage correction.\n", mute)
        Lk = self.get_leakage_array()
        for i, num_bursts in enumerate(self.num_bursts):
            if num_bursts == 0:
                continue  # if no bursts skip this ch
            self.na[i] -= self.nd[i] * Lk[i]
            self.nt[i] = self.nd[i] + self.na[i]
            if self.ALEX:
                self.nt[i] += self.naa[i]
        self.add(leakage_corrected=True)

    def direct_excitation_correction(self, mute=False):
        """Apply direct excitation correction to bursts (ALEX-only).

        The applied correction is: na -= naa*dir_ex
        """
        if self.dir_ex_corrected:
            return -1
        pprint("   - Applying direct excitation correction.\n", mute)
        for i, num_bursts in enumerate(self.num_bursts):
            if num_bursts == 0:
                continue  # if no bursts skip this ch
            self.na[i] -= self.naa[i] * self.dir_ex
            self.nt[i] = self.nd[i] + self.na[i]
            if self.ALEX:
                self.nt[i] += self.naa[i]
        self.add(dir_ex_corrected=True)

    def dither(self, lsb=2, mute=False):
        """Add dithering (uniform random noise) to burst counts (nd, na,...).

        The dithering amplitude is the range -0.5*lsb .. 0.5*lsb.
        """
        if self.dithering:
            return -1
        pprint("   - Applying burst-size dithering.\n", mute)
        self.add(dithering=True)
        for nd, na in zip(self.nd, self.na):
            nd += lsb * (np.random.rand(nd.size) - 0.5)
            na += lsb * (np.random.rand(na.size) - 0.5)
        if self.ALEX:
            for naa in self.naa:
                naa += lsb * (np.random.rand(naa.size) - 0.5)
            if 'nda' in self:
                for nda in self.nda:
                    nda += lsb * (np.random.rand(nda.size) - 0.5)
        self.add(lsb=lsb)

    def calc_chi_ch(self):
        """Calculate the gamma correction prefactor factor `chi_ch` (array).

        `chi_ch` is a ch-dependent prefactor for gamma used to correct
        the dispersion of fitted E peaks (`E_fit`).
        This method returns `chi_ch`. To apply the correction assign the array
        to the `Data.chi_ch` attribute.
        """
        if 'E_fit' not in self:
            print("ERROR: E_fit values not found. Call a `.fit_E_*` first.")
            return

        EE = self.E_fit.mean()  # Mean E value among the CH
        chi_ch = (1/EE - 1) / (1/self.E_fit - 1)
        return chi_ch

    def corrections(self, mute=False):
        """Apply corrections on burst-counts: nd, na, nda, naa.

        The corrections are: background, leakage (or bleed-through) and
        direct excitation (dir_ex).
        """
        self.background_correction(mute=mute)
        self.leakage_correction(mute=mute)
        if self.ALEX:
            self.direct_excitation_correction(mute=mute)

    def _update_corrections(self):
        """Recompute corrections whose flag is True.

        Checks the flags .bg_corrected, .leakage_corrected, .dir_ex_corrected,
        .dithering and recomputes the correction if the corresponding flag
        is True (i.e. if the correction was already applied).
        Note that this method is not used for gamma and beta corrections
        because these do not affect the `nd`, `na` and `naa` quantities but
        are only applied when computing E, S and corrected size.

        Differently from :meth:`corrections`, this allows to recompute
        corrections that have already been applied.
        """
        if 'mburst' not in self:
            return  # no burst search performed yet
        old_bg_corrected = self.bg_corrected
        old_leakage_corrected = self.leakage_corrected
        old_dir_ex_corrected = self.dir_ex_corrected
        old_dithering = self.dithering
        self.calc_ph_num()       # recompute uncorrected na, nd, nda, naa
        if old_bg_corrected:
            self.background_correction()
        if old_leakage_corrected:
            self.leakage_correction()
        if old_dir_ex_corrected:
            self.direct_excitation_correction()
        if old_dithering:
            self.dither(self.lsb)
        # Recompute E and S with no corrections (because already applied)
        self.calc_fret(count_ph=False, corrections=False)

    @property
    def leakage(self):
        """Spectral leakage (bleed-through) of D emission in the A channel.
        """
        return self._leakage

    @leakage.setter
    def leakage(self, leakage):
        self._update_leakage(leakage)

    def _update_leakage(self, leakage):
        """Apply/update leakage (or bleed-through) correction.
        """
        assert (np.size(leakage) == 1) or (np.size(leakage) == self.nch)
        self.add(_leakage=np.asfarray(leakage), leakage_corrected=True)
        self._update_corrections()

    @property
    def dir_ex(self):
        """Direct excitation correction factor."""
        return self._dir_ex

    @dir_ex.setter
    def dir_ex(self, value):
        self._update_dir_ex(value)

    def _update_dir_ex(self, dir_ex):
        """Apply/update direct excitation correction with value `dir_ex`.
        """
        assert np.size(dir_ex) == 1
        self.add(_dir_ex=float(dir_ex), dir_ex_corrected=True)
        self._update_corrections()

    @property
    def beta(self):
        """Beta factor used to correct S (compensates Dex and Aex unbalance).
        """
        return self._beta

    @beta.setter
    def beta(self, value):
        self._update_beta(value)

    def _update_beta(self, beta):
        """Change the `beta` value and recompute E and S."""
        assert np.size(beta) == 1
        self.add(_beta=float(beta))
        if 'mburst' in self:
            # Recompute E and S and delete fitter objects
            self.calc_fret(corrections=False)

    @property
    def chi_ch(self):
        """Per-channel relative gamma factor."""
        return self._chi_ch

    @chi_ch.setter
    def chi_ch(self, value):
        self._update_chi_ch(value)

    def _update_chi_ch(self, chi_ch):
        """Change the `chi_ch` value and recompute E and S."""
        msg = 'chi_ch is a per-channel correction and must have size == nch.'
        assert np.size(chi_ch) == self.nch, ValueError(msg)
        self.add(_chi_ch=np.asfarray(chi_ch))
        if 'mburst' in self:
            # Recompute E and S and delete fitter objects
            self.calc_fret(corrections=False)

    @property
    def gamma(self):
        """Gamma correction factor (compensates DexDem and DexAem unbalance).
        """
        return self._gamma

    @gamma.setter
    def gamma(self, value):
        self._update_gamma(value)

    def _update_gamma(self, gamma):
        """Change the `gamma` value and recompute E and S."""
        assert (np.size(gamma) == 1) or (np.size(gamma) == self.nch)
        self.add(_gamma=np.asfarray(gamma))
        if 'mburst' in self:
            # Recompute E and S and delete fitter objects
            self.calc_fret(corrections=False)

    def get_gamma_array(self):
        """Get the array of gamma factors, one per ch.

        It always returns an array of gamma factors regardless of
        whether `self.gamma` is scalar or array.

        Each element of the returned array is multiplied by `chi_ch`.
        """
        gamma = self.gamma
        G = np.repeat(gamma, self.nch) if np.size(gamma) == 1 else gamma
        G *= self.chi_ch
        return G

    def get_leakage_array(self):
        """Get the array of leakage coefficients, one per ch.

        It always returns an array of leakage coefficients regardless of
        whether `self.leakage` is scalar or array.

        Each element of the returned array is multiplied by `chi_ch`.
        """
        leakage = self.leakage
        Lk = np.r_[[leakage] * self.nch] if np.size(leakage) == 1 else leakage
        Lk *= self.chi_ch
        return Lk

    ##
    # Methods to compute burst quantities: FRET, S, SBR, max_rate, etc ...
    #
    def calc_sbr(self, ph_sel=Ph_sel('all'), gamma=1.):
        """Return Signal-to-Background Ratio (SBR) for each burst.

        Arguments:
            ph_sel (Ph_sel object): object defining the photon selection
                for which to compute the sbr. Changes the photons used for
                burst size and the corresponding background rate. Valid values
                here are Ph_sel('all'), Ph_sel(Dex='Dem'), Ph_sel(Dex='Aem').
                See :mod:`fretbursts.ph_sel` for details.
            gamma (float): gamma value used to compute corrected burst size
                in the case `ph_sel` is Ph_sel('all'). Ignored otherwise.
        Returns:
            A list of arrays (one per channel) with one value per burst.
            The list is also saved in `sbr` attribute.
        """
        ph_sel = self._fix_ph_sel(ph_sel)
        sbr = []
        for ich, mb in enumerate(self.mburst):
            if mb.num_bursts == 0:
                sbr.append([])
                continue  # if no bursts skip this ch
            nd, na, bg_d, bg_a = self.expand(ich)
            nt = self.burst_sizes_ich(ich=ich, gamma=gamma)

            signal = {Ph_sel('all'): nt,
                      Ph_sel(Dex='Dem'): nd, Ph_sel(Dex='Aem'): na}

            background = {Ph_sel('all'): bg_d + bg_a,
                          Ph_sel(Dex='Dem'): bg_d, Ph_sel(Dex='Aem'): bg_a}

            sbr.append(signal[ph_sel] / background[ph_sel])
        self.add(sbr=sbr)
        return sbr

    def calc_burst_ph_func(self, func, func_kw, ph_sel=Ph_sel('all'),
                           compact=False, ich=0):
        """Evaluate a scalar function from photons in each burst.

        This method allow calling an arbitrary function on the photon
        timestamps of each burst. For example if `func` is `np.mean` it
        computes the mean time in each bursts.

        Arguments:
            func (callable): function that takes as first argument an array of
                timestamps for one burst.
            func_kw (callable): additional arguments to be passed  `func`.
            ph_sel (Ph_sel object): object defining the photon selection.
                See :mod:`fretbursts.ph_sel` for details.
            compact (bool): if True, a photon selection of only one excitation
                period is required and the timestamps are "compacted" by
                removing the "gaps" between each excitation period.

        Returns:
            A list (on element per channel) array. The array size is equal to
            the number of bursts in the corresponding channel.
        """
        if compact:
            self._assert_compact(ph_sel)

        kwargs = dict(func=func, func_kw=func_kw, compact=compact)
        if self.ALEX:
            kwargs.update(alex_period=self.alex_period)
        if compact:
            kwargs.update(excitation_width=self._excitation_width(ph_sel))

        results_mch = [burst_ph_stats(ph, bursts, mask=mask, **kwargs)
                       for ph, mask, bursts in
                       zip(self.iter_ph_times(),
                           self.iter_ph_masks(ph_sel=ph_sel),
                           self.mburst)]
        return results_mch

    def calc_max_rate(self, m, ph_sel=Ph_sel('all'), compact=False,
                      c=phrates.default_c):
        """Compute the max m-photon rate reached in each burst.

        Arguments:
            m (int): number of timestamps to use to compute the rate.
                As for burst search, typical values are 5-20.
            ph_sel (Ph_sel object): object defining the photon selection.
                See :mod:`fretbursts.ph_sel` for details.
            c (float): this parameter is used in the definition of the
                rate estimator which is `(m - 1 - c) / t[last] - t[first]`.
                For more details see :func:`.phtools.phrates.mtuple_rates`.
        """
        ph_sel = self._fix_ph_sel(ph_sel)
        Max_Rate = self.calc_burst_ph_func(func=phrates.mtuple_rates_max,
                                           func_kw=dict(m=m, c=c),
                                           ph_sel=ph_sel, compact=compact)
        Max_Rate = [mr / self.clk_p - bg[bp] for bp, bg, mr in
                    zip(self.bp, self.bg_from(ph_sel), Max_Rate)]
        params = dict(m=m, ph_sel=ph_sel, compact=compact)
        self.add(max_rate=Max_Rate, max_rate_params=params)

    def calc_fret(self, count_ph=False, corrections=True, dither=False,
                  mute=False, pure_python=False):
        """Compute FRET (and stoichiometry if ALEX) for each burst.

        This is an high-level functions that can be run after burst search.
        By default, it will count Donor and Acceptor photons, perform
        corrections (background, leakage), and compute gamma-corrected
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
        self._calculate_fret_eff()
        if self.ALEX:
            self._calculate_stoich()
            #self._calc_alex_hist()

        for attr in ('ES_binwidth', 'ES_hist', 'E_fitter', 'S_fitter'):
            # E_fitter and S_fitter are only attributes
            # so we cannot use the membership syntax (attr in self)
            if hasattr(self, attr):
                self.delete(attr, warning=False)

    def _calculate_fret_eff(self):
        """Compute FRET efficiency (`E`) for each burst."""
        G = self.get_gamma_array()
        E = [na / (g*nd + na) for nd, na, g in zip(self.nd, self.na, G)]
        self.add(E=E)

    def _calculate_stoich(self):
        """Compute "stoichiometry" (the `S` parameter) for each burst."""
        G = self.get_gamma_array()
        S = [(g*d + a) / (g*d + a + aa/self.beta) for d, a, aa, g in
             zip(self.nd, self.na, self.naa, G)]
        self.add(S=S)

    def _calc_alex_hist(self, binwidth=0.05):
        """Compute the ALEX histogram with given bin width `bin_step`"""
        if 'ES_binwidth' in self and self.ES_binwidth == binwidth:
            return

        ES_hist_tot = [ES_histog(E, S, binwidth) for E, S in
                       zip(self.E, self.S)]
        E_bins, S_bins = ES_hist_tot[0][1], ES_hist_tot[0][2]
        ES_hist = [h[0] for h in ES_hist_tot]
        E_ax = E_bins[:-1] + 0.5 * binwidth
        S_ax = S_bins[:-1] + 0.5 * binwidth
        self.add(ES_hist=ES_hist, E_bins=E_bins, S_bins=S_bins,
                 E_ax=E_ax, S_ax=S_ax, ES_binwidth=binwidth)

    ##
    # Methods for measurement info
    #
    def status(self, add="", noname=False):
        """Return a string with burst search, corrections and selection info.
        """
        name = "" if noname else self.name
        s = name
        if 'L' in self:  # burst search has been done
            if 'rate_th' in self:
                s += " BS_%s L%d m%d MR%d" % (self.ph_sel, self.L, self.m,
                                              np.mean(self.rate_th) * 1e-3)
            else:
                P_str = '' if self.P is None else ' P%s' % self.P
                s += " BS_%s L%d m%d F%.1f%s" % \
                     (self.ph_sel, self.L, self.m, np.mean(self.F), P_str)
        s += " G%.3f" % np.mean(self.gamma)
        if 'bg_fun' in self: s += " BG%s" % self.bg_fun.__name__[:-4]
        if 'bg_time_s' in self: s += "-%ds" % self.bg_time_s
        if 'fuse' in self: s += " Fuse%.1fms" % self.fuse
        if 'bg_corrected' in self and self.bg_corrected:
            s += " bg"
        if 'leakage_corrected' in self and self.leakage_corrected:
            s += " Lk%.3f" % np.mean(self.leakage*100)
        if 'dir_ex_corrected' in self and self.dir_ex_corrected:
            s += " dir%.1f" % (self.dir_ex*100)
        if 'dithering' in self and self.dithering:
            s += " Dith%d" % self.lsb
        if 's' in self: s += ' '.join(self.s)
        return s + add

    @property
    def name(self):
        """Measurement name: last subfolder + file name with no extension."""
        if not hasattr(self, '_name'):
            basename = str(os.path.splitext(os.path.basename(self.fname))[0])
            name = basename
            last_dir = str(os.path.basename(os.path.dirname(self.fname)))
            if len(last_dir) > 0:
                name = '_'.join([last_dir, basename])
            self.add(_name=name)
        return self._name

    @name.setter
    def name(self, value):
        self.add(_name=value)

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
            print(s)

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
        Mask = self.select_bursts_mask(select_bursts.E, E1=E1, E2=E2)

        fit_res, fit_model_F = zeros((self.nch, 2)), zeros(self.nch)
        for ich, (nd, na, E, mask) in enumerate(zip(
                self.nd, self.na, self.E, Mask)):
            w = fret_fit.get_weights(nd[mask], na[mask],
                                     weights=weights, gamma=gamma)
            # Compute weighted mean
            fit_res[ich, 0] = np.dot(w, E[mask])/w.sum()
            # Compute weighted variance
            fit_res[ich, 1] = np.sqrt(
                np.dot(w, (E[mask] - fit_res[ich, 0])**2)/w.sum())
            fit_model_F[ich] = mask.sum()/mask.size

        fit_model = lambda x, p: SS.norm.pdf(x, p[0], p[1])
        self.add(fit_E_res=fit_res, fit_E_name='Moments',
                 E_fit=fit_res[:, 0], fit_E_curve=True, fit_E_E1=E1,
                 fit_E_E2=E2, fit_E_model=fit_model,
                 fit_E_model_F=fit_model_F)
        self.fit_E_calc_variance()
        return self.E_fit

    def fit_E_ML_poiss(self, E1=-1, E2=2, method=1, **kwargs):
        """ML fit for E modeling size ~ Poisson, using bursts in [E1,E2] range.
        """
        assert method in [1, 2, 3]
        fit_fun = {1: fret_fit.fit_E_poisson_na, 2: fret_fit.fit_E_poisson_nt,
                   3: fret_fit.fit_E_poisson_nd}
        Mask = self.select_bursts_mask(select_bursts.E, E1=E1, E2=E2)
        fit_res = zeros(self.nch)
        for ich, mask in zip(range(self.nch), Mask):
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
        Mask = self.select_bursts_mask(select_bursts.E, E1=E1, E2=E2)
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
        fun_d = dict(slope=fret_fit.fit_E_slope,
                     E_size=fret_fit.fit_E_E_size)
        Mask = self.select_bursts_mask(select_bursts.E, E1=E1, E2=E2)
        fit_res = np.array([fun_d[kind](nd[mask], na[mask], **kwargs)
                            for nd, na, mask in
                            zip(self.nd, self.na, Mask)])
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
            w = fret_fit.get_weights(nd, na, weights=weights, gamma=gamma)
            fit_res[ich, :] = fit_func(E, weights=w, **kwargs)
        self.add(fit_E_res=fit_res, fit_E_name=fit_func.__name__,
                 E_fit=fit_res[:, 2], fit_E_curve=True,
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
            raise ValueError("Fitting function not recognized.")

        Mask = self.select_bursts_mask(select_bursts.E, E1=E1, E2=E2)

        fit_res, fit_model_F = zeros((self.nch, nparam)), zeros(self.nch)
        for ich, (nd, na, E, mask) in enumerate(zip(
                self.nd, self.na, self.E, Mask)):
            if '_hist' in fit_fun.__name__ or '_EM' in fit_fun.__name__:
                if weights is None:
                    w = None
                else:
                    w = fret_fit.get_weights(nd[mask], na[mask],
                                             weights=weights, gamma=gamma)
                fit_res[ich, :] = fit_fun(E[mask], weights=w, **fit_kwargs)
            else:
                # Non-histogram fits (PDF/CDF) do not support weights
                fit_res[ich, :] = fit_fun(E[mask], **fit_kwargs)
            fit_model_F[ich] = mask.sum()/mask.size

        # Save enough info to generate a fit plot (see hist_fret in burst_plot)
        self.add(fit_E_res=fit_res, fit_E_name=fit_fun.__name__,
                 E_fit=fit_res[:, iE], fit_E_curve=True, fit_E_E1=E1,
                 fit_E_E2=E2, fit_E_model=fit_model,
                 fit_E_model_F=fit_model_F, fit_E_weights=weights,
                 fit_E_gamma=gamma, fit_E_kwargs=fit_kwargs)
        return self.E_fit

    def fit_from(self, D):
        """Copy fit results from another Data() variable.
        Now that the fit methods accept E1,E1 parameter this probabily useless.
        """
        # NOTE Are 'fit_guess' and 'fit_fix' still used ?
        fit_data = ['fit_E_res', 'fit_E_name', 'E_fit', 'fit_E_curve',
                    'fit_E_E1', 'fit_E_E2=E2', 'fit_E_model',
                    'fit_E_model_F', 'fit_guess', 'fit_fix']
        for name in fit_data:
            if name in D:
                self[name] = D[name]
                setattr(self, name, self[name])
        # Deal with the normalization to the number of bursts
        self.add(fit_model_F=r_[[old_E.size/new_E.size \
                                 for old_E, new_E in zip(D.E, self.E)]])

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

        E_sel = [Ei[(Ei > E1)*(Ei < E2)] for Ei in self.E]
        Mask = self.select_bursts_mask(select_bursts.E, E1=E1, E2=E2)

        E_var, E_var_bu, E_var_ph = \
            zeros(self.nch), zeros(self.nch), zeros(self.nch)
        for i, (Ech, nt, mask) in enumerate(zip(E_sel, self.nt, Mask)):
            nt_s = nt[mask]
            nd_s, na_s = self.nd[i][mask], self.na[i][mask]
            w = fret_fit.get_weights(nd_s, na_s, weights=weights)
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
