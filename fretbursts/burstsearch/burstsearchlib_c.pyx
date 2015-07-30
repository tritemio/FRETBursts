#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
Optimized version of burst search functions to be compiled in C with Cython.
To compile run::

    python setup.py build_ext --inplace

"""

import sys
import numpy as np
cimport numpy as np


cdef pprint(s):
    """Print immediately, even if inside a busy loop."""
    sys.stdout.write(s)
    sys.stdout.flush()

def bsearch_c(np.int64_t[:] t, np.int16_t L, np.int16_t m, np.float64_t T,
              label='burst search', verbose=True):
    """Sliding window burst search. Cython implementation (fastest version).

    Finds bursts in the array `t` (int64). A burst starts when the photon rate
    is above a minimum threshold, and ends when the rate falls below the same
    threshold. The rate-threshold is defined by the ratio `m`/`T` (`m` photons
    in a time interval `T`). A burst is discarded if it has less than `L`
    photons.

    Arguments:
        t (array, int64): 1D array of timestamps on which to perform the search
        L (int16): minimum number of photons in a bursts. Bursts with size
            (or counts) < L are discarded.
        m (int16): number of consecutive photons used to compute the rate.
        T (float64): max time separation of `m` photons to be inside a burst
        label (string): a label printed when the function is called
        verbose (bool): if False, the function does not print anything.

    Returns:
        2D array of burst data, one row per burst, shape (N, 6), type int64.
        Each row contains (in the order): time of burst start, burst duration,
        number of photons, index of start time, time of burst end.
        To extract burst information it's safer to use the utility functions
        `b_*` (i.e. :func:`b_start`, :func:`b_size`, :func:`b_width`, etc...).
    """
    cdef int i, i_start, i_end
    cdef np.int64_t burst_start, burst_end, t1, t2
    cdef np.int8_t in_burst = False

    if verbose: pprint('C Burst search: %s\n' % label)
    bursts = []

    for i in xrange(t.size-m+1):
        if (t[i+m-1] - t[i]) <= T:
            if not in_burst:
                i_start = i
                in_burst = 1
        elif in_burst:
            in_burst = 0
            i_end = i+m-1
            if i_end - i_start >= L:
                burst_start, burst_end = t[i_start], t[i_end-1]
                bursts.append([burst_start, burst_end-burst_start,
                    i_end-i_start, i_start, i_end-1, burst_end])
    return np.array(bursts, dtype=np.int64)


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Functions to count D and A ph in bursts
#
def mch_count_ph_in_bursts_c(mburst, Mask):
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
    cdef np.int8_t itstart, iwidth, inum_ph, iistart, iiend, itend
    cdef np.ndarray[np.int64_t, ndim=2] bursts
    cdef np.ndarray[np.int8_t, ndim=1] mask
    cdef np.ndarray[np.int16_t, ndim=1] num_ph
    cdef np.ndarray[np.int64_t, ndim=1] b
    cdef int i,ii
    cdef np.ndarray m
    Mask_c = [m.astype(np.int8) for m in Mask]
    itstart, iwidth, inum_ph, iistart, iiend, itend = 0,1,2,3,4,5
    Num_ph = []
    for bursts, mask in zip(mburst, Mask_c):
        num_ph = np.zeros(bursts.shape[0], dtype=np.int16)
        for i,b in enumerate(bursts):
            # Counts donors between start and end of current burst b
            #num_ph[i] = mask[b[iistart]:b[iiend]+1].sum()
            for ii in xrange(b[iistart], b[iiend]+1):
                num_ph[i] += mask[ii]
            #count_ph_in_bursts(bursts, mask).astype(float)
        Num_ph.append(num_ph.astype(np.float))
    return Num_ph
