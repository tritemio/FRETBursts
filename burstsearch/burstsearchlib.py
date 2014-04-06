#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
Core burst search and photon counting functions.

All the burst search functions return a 2-D array (burst array) of shape Nx6,
where N is the number of bursts. The 6 columns contain the burst data.

The burst array can be indexed of sliced along the first dimension (row wise
or axis=0) to select one or more bursts. However it is discouraged to access
specific burst fields in the second dimension (columns or axis=1). The proper
way to access burst data is using the burst-data functions (starting with
`b_`). These functions are both clearer and less bug-prone than using column
index. In other terms the second dimension of burst-array should be considered
opaque.

The list of burst-data functions is: :func:`b_start`, :func:`b_end`,
:func:`b_size`, :func:`b_width`, :func:`b_istart`, :func:`b_istart`,
:func:`b_iend`, :func:`b_rate`, :func:`b_separation`.

For example, assume a burst array `mburst`. To take a slice of only the first
10 bursts you can do::

    mburst2 = mburst[:10]   # new array with burst data of the first 10 bursts

Obtain the burst start of all the bursts::

    b_start(mbursts)

Obtain the burst size (number of photons) for 10-th to 20-th burst::

    b_size(mbursts[10:20])

"""

import numpy as np
from utils.misc import pprint

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  GLOBAL VARIABLES
#

# Define name for the 6 column indexes in the Nx6 array containing the bursts
itstart, iwidth, inum_ph, iistart, iiend, itend = 0, 1, 2, 3, 4, 5

# Quick functions for burst data
def b_start(b):
    """Time of 1st ph in burst"""
    return b[:, itstart]

def b_end(b):
    """Time of last ph in burst"""
    return b[:, itend]

def b_width(b):
    """Burst width in clk cycles"""
    return b[:, iwidth]

def b_istart(b):
    """Index of 1st ph in burst"""
    return b[:, iistart]

def b_iend(b):
    """Index of last ph in burst"""
    return b[:, iiend]

def b_size(b):
    """Number of ph in the burst"""
    return b[:, inum_ph]

def b_rate(b):
    """Mean rate of ph in burst"""
    return 1.*b_size(b)/b_width(b)

def b_separation(b):
    """Separation between nearby bursts"""
    return b_start(b)[1:] - b_end(b)[:-1]

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  LOW-LEVEL BURST SEARCH FUNCTIONS
#

def bsearch_py(t, L, m, T, label='Burst search', verbose=True):
    """Sliding window burst search. Pure python implementation.

    Finds bursts in the array `t` (int64). A burst starts when the photon rate
    is above a minimum threshold, and ends when the rate falls below the same
    threshold. The rate-threshold is defined by the ratio `m`/`T` (`m` photons
    in a time interval `T`). A burst is discarded if it has less than `L`
    photons.

    Arguments:
        t (array, int64): array of timestamps on which to perform the search
        L (int): minimum number of photons in a bursts. Bursts with size
            (or counts) < L are discarded.
        m (int): number of consecutive photons used to compute the rate.
        T (float): max time separation of `m` photons to be inside a burst
        label (string): a label printed when the function is called
        verbose (bool): if False, the function does not print anything.

    Returns:
        2D array of burst data, one row per burst, shape (N, 6), type int64.
        Each row contains (in the order): time of burst start, burst duration,
        number of photons, index of start time, time of burst end.
        To extract burst information it's safer to use the utility functions
        `b_*` (i.e. :func:`b_start`, :func:`b_size`, :func:`b_width`, etc...).
    """
    if verbose: pprint('Python search (v): %s\n' % label)
    bursts = []
    in_burst = False
    above_min_rate = (t[m-1:]-t[:t.size-m+1])<=T
    for i in xrange(t.size-m+1):
        if above_min_rate[i]:
            if not in_burst:
                i_start = i
                in_burst = True
        elif in_burst:
            in_burst = False
            i_end = i+m-1
            if i_end - i_start >= L:
                burst_start, burst_end = t[i_start], t[i_end-1]
                bursts.append([burst_start, burst_end-burst_start,
                    i_end-i_start, i_start, i_end-1, burst_end])
    return np.array(bursts, dtype=np.int64)


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Functions to count D and A ph in bursts
#

def count_ph_in_bursts_py(bursts, mask):
    """Counts number of photons in each burst counting only photons in `mask`.

    This function takes a burst-array and a boolean mask (photon selection)
    and computes the number of photons selected by the mask.
    It is used, for example, to count donor and acceptor photons in each burst.

    This is a reference implementation. In practice the multi-channel
    is always used instead (see :func:`mch_count_ph_in_bursts_py`).

    Arguments:
        Mburst (list of 2D arrays, int64): a list of burst-arrays, one per ch.
        Mask (list of 1D boolean arrays): a list of photon masks (one per ch),
            For each channel, the boolean mask must be of the same size of the
            timestamp array used for burst search.

    Returns:
        A list of 1D arrays, each containing the number of photons in the
        photon selection mask.
    """
    assert (itstart, iwidth, inum_ph, iistart, iiend, itend) == (0,1,2,3,4,5)
    num_ph = np.zeros(bursts.shape[0], dtype=np.int16)
    for i,b in enumerate(bursts):
        # Counts photons between start and end of current burst b
        num_ph[i] = mask[b[iistart]:b[iiend]+1].sum()
    return num_ph

def mch_count_ph_in_bursts_py(Mburst, Mask):
    """Counts number of photons in each burst counting only photons in `Mask`.

    This multi-channel function takes a list of burst-arrays and a list of
    photon masks and compute the number of photons selected by the mask.
    It is used, for example, to count donor and acceptor photons in each burst.

    Arguments:
        Mburst (list of 2D arrays, int64): a list of burst-arrays, one per ch.
        Mask (list of 1D boolean arrays): a list of photon masks (one per ch),
            For each channel, the boolean mask must be of the same size of the
            timestamp array used for burst search.

    Returns:
        A list of 1D arrays, each containing the number of photons in the
        photon selection mask.
    """
    assert (itstart, iwidth, inum_ph, iistart, iiend, itend) == (0,1,2,3,4,5)
    Num_ph = []
    for bursts, mask in zip(Mburst, Mask):
        num_ph = np.zeros(bursts.shape[0], dtype=np.int16)
        for i,b in enumerate(bursts):
            # Counts photons between start and end of current burst b
            num_ph[i] = mask[b[iistart]:b[iiend]+1].sum()
            #count_ph_in_bursts(bursts, mask).astype(float)
        Num_ph.append(num_ph.astype(float))
    return Num_ph


##
#  Try to import the optimized Cython functions
#

try:
    from burstsearchlib_c import bsearch_c
    bsearch = bsearch_c
    print " - Optimized (cython) burst search loaded."
except ImportError:
    bsearch = bsearch_py
    print " - Fallback to pure python burst search."

try:
    from burstsearchlib_c import mch_count_ph_in_bursts_c
    mch_count_ph_in_bursts = mch_count_ph_in_bursts_c
    print " - Optimized (cython) photon counting loaded."
except ImportError:
    mch_count_ph_in_bursts = mch_count_ph_in_bursts_py
    print " - Fallback to pure python photon counting."

