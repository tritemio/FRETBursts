#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
Core burst search and photon counting functions.

All the burst search functions return a 2-D array of shape Nx6, where N is 
the number of bursts (burst array). The 6 columns contain burst data 
(time start, width, number of photons, index of time start, index of time end, 
time end).

The burst array can be indexed of sliced along the first dimension (row wise
or axis=0) to select one or more bursts. However to access specific burst 
fields in the second dimension (columns or axis=1) the b_* utility functions 
should be used. These is both clearer and less bug-prone than using column
index to access burst data.

For example, assume a burst array `mburst`. To take a slice of only the first 
10 bursts you can do::

    mburst2 = mburst[:10]   # new array with burst data of the first 10 bursts
    
Obtain the burst start fo all the bursts::
    
    b_start(mbursts)

Obtain the burst size (number of photons) for burst 10 to 20::

    b_size(mbursts[10:20])
    
"""

import numpy as np
from utils.misc import pprint

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#  GLOBAL VARIABLES
#

# Define name for the 6 column indexes in the Nx6 array containing the bursts
itstart, iwidth, inum_ph, iistart, iiend, itend = 0, 1, 2, 3, 4, 5

# Quick functions for bursts start and end
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

def wba(t,L,m,T):
    """Sliding window burst analysis. (Classical)"""
    bursts = []
    in_burst = False
    for i in range(t.size):
        ph_in_T = ((t>=t[i])*(t<t[i]+T)).sum()
        if ph_in_T > m:
            if not in_burst:
                burst_start = t[i]
                i_start = i
                in_burst = True
        elif in_burst:
            in_burst = False
            burst_end = t[i-1+m]
            if i-1+m - i_start >= L:
                bursts.append([burst_start,burst_end-burst_start, 
                    i-1+m-i_start])
    return np.array(bursts)

def ba_pure(t,L,m,T,label='Burst search'):
    """FIFO burst search (T in clk periods). Reference implementation."""
    pprint('Python search: %s\n'%label)
    bursts = []
    in_burst = False
    for i in xrange(t.size-m+1):
        i_end = i+m-1
        if t[i_end]-t[i] <= T:
            if not in_burst:
                i_start = i
                in_burst = True
        elif in_burst:
            # Note that i_end is the index of the last ph in the current time 
            # window, while the last ph in a burst is (i_end-1). 
            # Note however that the number of ph in a burst is (i_end-i_start),
            # not (i_end-1)-i_start as may erroneously appears at first glace.
            in_burst = False
            if i_end - i_start >= L:
                burst_start, burst_end = t[i_start], t[i_end-1]
                bursts.append([burst_start,burst_end-burst_start,
                    i_end-i_start, i_start, i_end-1, burst_end])
    return np.array(bursts, dtype=np.int64)

def ba_pure_o(t, L, m, T, label='Burst search', verbose=True):
    """FIFO burst search (T in clk periods). Optimized implementation."""
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

def ba_pure_a(t,L,m,T,label='Burst search'):
    """FIFO burst search (T in clk periods). Optimized implementation."""
    pprint('Python search (a): %s\n'%label)
    bursts = np.zeros((3e3, 6))
    iburst = 0
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
                bursts[iburst] = (burst_start, burst_end-burst_start,
                    i_end-i_start, i_start, i_end-1, burst_end)
                iburst += 1
    return np.array(bursts, dtype=np.int64)[:iburst]


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#  Functions to count D and A ph in bursts
#

def count_ph_in_bursts(bursts, mask):
    """Return num. ph in each burst considering only ph in mask.

    bursts: Nx4 array containing burst data
    mask: 1-D bool array (same size as the ph array), true if a ph is "counted"
    """
    assert (itstart, iwidth, inum_ph, iistart, iiend, itend) == (0,1,2,3,4,5)
    num_ph = np.zeros(bursts.shape[0], dtype=np.int16)
    for i,b in enumerate(bursts):
        # Counts donors between start and end of current burst b
        num_ph[i] = mask[b[iistart]:b[iiend]+1].sum()
    return num_ph

def mch_count_ph_in_bursts(mburst, Mask):
    """Counts num ph in each burst counting only ph in Mask (multi-ch version).
    
    mburst: is a list of burst data (Nx4 array), one per ch.
    Mask: is a list of ph masks (one per ch), same size as ph_time_m used for
          burst search.
    """
    assert (itstart, iwidth, inum_ph, iistart, iiend, itend) == (0,1,2,3,4,5)
    Num_ph = []
    for bursts, mask in zip(mburst, Mask):
        num_ph = np.zeros(bursts.shape[0], dtype=np.int16)
        for i,b in enumerate(bursts):
            # Counts donors between start and end of current burst b
            num_ph[i] = mask[b[iistart]:b[iiend]+1].sum()
            #count_ph_in_bursts(bursts, mask).astype(float)
        Num_ph.append(num_ph.astype(float))
    return Num_ph

