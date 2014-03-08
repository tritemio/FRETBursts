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

def pprint(s):
    """Print immediately, even if inside a busy loop."""
    sys.stdout.write(s)
    sys.stdout.flush()

def ba_c(np.ndarray[np.int64_t, ndim=1] t, L, m, T, label='burst search',
	 verbose=True):
    """FIFO burst search (T in clk periods). Optimized implementation."""
    cdef int i, i_start, i_end
    cdef np.int64_t burst_start, burst_end
    assert t.dtype == np.int64
    #cdef char in_burst
    #assert t.dtype == DTYPE
    if verbose: pprint('C Burst search: %s\n' % label)
    bursts = []
    in_burst = False
    cdef np.ndarray[np.int8_t, ndim=1] above_min_rate
    above_min_rate = np.int8(((t[m-1:]-t[:t.size-m+1])<=T))
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

def ba_pure_c(np.ndarray[np.int64_t, ndim=1] t,
		int L, int m, np.float64_t T, label='Burst search'):
    """FIFO burst search (T in clk periods). Reference implementation."""
    cdef int i, i_start, i_end, iburst
    cdef np.int64_t burst_start, burst_end
    cdef char in_burst
    assert t.dtype == np.int64
    pprint('C Burst search (pure): %s\n'%label)
    #bursts = []
    cdef np.ndarray[np.int64_t, ndim=2] bursts = np.zeros((3e3, 6), 
		    dtype=np.int64)
    iburst = 0
    in_burst = 0
    for i in xrange(t.size-m+1):
        i_end = i+m-1
        if t[i_end]-t[i] <= T:
            if not in_burst:
                i_start = i
                in_burst = 1
        elif in_burst:
            # Note that i_end is the index of the last ph in the current time 
            # window, while the last ph in a burst is (i_end-1). 
            # Note however that the number of ph in a burst is (i_end-i_start),
            # not (i_end-1)-i_start as may erroneously appears at first glace.
            in_burst = 0
            if i_end - i_start >= L:
                burst_start, burst_end = t[i_start], t[i_end-1]
                #bursts.append([burst_start,burst_end-burst_start,
                #    i_end-i_start, i_start, i_end-1, burst_end])
                bursts[iburst][0] = burst_start
                bursts[iburst][0] = burst_end-burst_start
                bursts[iburst][0] = i_end-i_start
                bursts[iburst][0] = i_end-1
                bursts[iburst][0] = burst_end
                bursts[iburst][0] = burst_start
                iburst += 1
    return np.array(bursts, dtype=np.int64)[:iburst]

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#  Functions to count D and A ph in bursts
#

def c_mch_count_ph_in_bursts(mburst, Mask):
    """Counts num ph in each burst counting only ph in Mask (multi-ch version).
    
    mburst: is a list of burst data (Nx4 array), one per ch.
    Mask: is a list of ph masks (one per ch), same size as ph_time_m used for
          burst search.
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
            for ii in range(b[iistart],b[iiend]+1):
                num_ph[i] += mask[ii]
            #count_ph_in_bursts(bursts, mask).astype(float)
        Num_ph.append(num_ph.astype(np.float))
    return Num_ph

