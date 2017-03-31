#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014-2016 The Regents of the University of California,
#               Antonino Ingargiola <tritemio@gmail.com>
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

def bsearch_c(np.int64_t[:] times, np.int16_t L, np.int16_t m,
              np.float64_t T, slice_=None,
              label='Burst search', verbose=True):
    """Sliding window burst search. Pure python implementation.

    Finds bursts in the array `times` (int64). A burst starts when the photon rate
    is above a minimum threshold, and ends when the rate falls below the same
    threshold. The rate-threshold is defined by the ratio `m`/`T` (`m` photons
    in a time interval `T`). A burst is discarded if it has less than `L`
    photons.

    Arguments:
        times (array, int64): array of timestamps on which to perform the search
        L (int): minimum number of photons in a bursts. Bursts with size
            (or counts) < L are discarded.
        m (int): number of consecutive photons used to compute the rate.
        T (float): max time separation of `m` photons to be inside a burst
        slice_ (tuple): 2-element tuple used to slice times
        label (string): a label printed when the function is called
        verbose (bool): if False, the function does not print anything.

    Returns:
        Array of burst data Nx4, type int64.
    """
    cdef np.int64_t islice1, islice2, i, i_start, i_stop
    cdef np.uint8_t in_burst

    if verbose:
        pprint('Python search (v): %s\n' % label)
    if slice_ is not None:
        islice1 = slice_[0]
        islice2 = slice_[1]
    else:
        islice1 = 0
        islice2 = times.size

    bursts = []
    in_burst = 0
    for i in range(islice1, islice2 - m + 1):
        if (times[i+m-1] - times[i]) <= T:
            if not in_burst:
                in_burst = 1
                i_start = i
        elif in_burst:
            in_burst = 0
            i_stop = i + m - 2
            if i_stop - i_start + 1 >= L:
                bursts.append((i_start, i_stop,
                               times[i_start], times[i_stop]))

    if in_burst:
        i_stop = i + m - 1
        if i_stop - i_start + 1 >= L:
            bursts.append((i_start, i_stop, times[i_start], times[i_stop]))

    return np.array(bursts, dtype='int64')


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Functions to count D and A ph in bursts
#
def count_ph_in_bursts_c(np.int64_t[:,:] burstdata, np.uint8_t[:] mask):
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
    cdef np.int64_t i, ii
    cdef np.int32_t[:] num_ph

    num_ph = np.zeros(burstdata.shape[0], dtype=np.int32)
    for i in range(burstdata.shape[0]):
        for ii in range(burstdata[i, 0], burstdata[i, 1]+1):
            num_ph[i] += mask[ii]
    return num_ph

def mch_count_ph_in_bursts_c(mburst, masks):
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
    num_ph_list = []
    for bursts, mask in zip(mburst, masks):
        num_ph_list.append(
            np.asfarray(count_ph_in_bursts_c(bursts.data, mask.view(np.uint8))))
    return num_ph_list
