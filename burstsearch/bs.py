
from numpy import *
from poisson_threshold import find_optimal_T
from utils import pprint

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#  GLOBAL VARIABLES
#

# Define name for the 6 column indexes in the Nx6 array containing the bursts
itstart, iwidth, inum_ph, iistart, iiend, itend = 0,1,2,3,4,5

# Quick functions for bursts start and end
b_start = lambda b: b[:,itstart]            # time of 1st ph in burst
b_end = lambda b: b[:,itend]                # time of last ph in burst
b_width = lambda b: b[:,iwidth]             # burst width in clk cycles
b_istart = lambda b: b[:,iistart]           # index of 1st ph in burst
b_iend = lambda b: b[:,iiend]               # index of last ph in burst
b_size = lambda b: b[:,inum_ph]             # number of ph in the burst
b_rate = lambda b: 1.*b_size(b)/b_width(b)  # mean rate of ph in burst

# Separation between nearby bursts
b_separation = lambda b: b[1:,itstart] - b_end(b)[:-1]

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
    return array(bursts)

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
    return array(bursts, dtype=int64)

def ba_pure_o(t,L,m,T,label='Burst search'):
    """FIFO burst search (T in clk periods). Optimized implementation."""
    pprint('Python search (v): %s\n'%label)
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
    return array(bursts, dtype=int64)

def ba_pure_a(t,L,m,T,label='Burst search'):
    """FIFO burst search (T in clk periods). Optimized implementation."""
    pprint('Python search (a): %s\n'%label)
    bursts = zeros((3e3, 6))
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
    return array(bursts, dtype=int64)[:iburst]

def ba_cpp(t,L,m,T, label='Burst search', fast=False):
    """FIFO burst search (T in clk periods). Optimized. Prints what it does.

    If fast==True the C++ burst search is performed.

    [OBSOLETE] C++ code is difficult to maintain. Now the same high
    performances are achieved through cython (see .pyx file).
    """
    if fast: 
        bursts = c_findbursts(t, L, m, T)
    else:
        pprint('Python search: %s\n'%label)
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
    return array(bursts, dtype=int64)
    # NOTE: Returns int64 and not uint64 so that aritmethics does not
    # overflows (ex burst separation < 0)


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#  Functions to count D and A ph in bursts
#

def count_ph_in_bursts(bursts, mask):
    """Return num. ph in each burst considering only ph in mask.

    bursts: Nx4 array containing burst data
    mask: 1-D bool array (same size as the ph array), true if a ph is "counted"
    """
    assert (itstart, iwidth, inum_ph, iistart, iiend, itend) == (0,1,2,3,4,5)
    num_ph = zeros(bursts.shape[0], dtype=int16)
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
        num_ph = zeros(bursts.shape[0], dtype=int16)
        for i,b in enumerate(bursts):
            # Counts donors between start and end of current burst b
            num_ph[i] = mask[b[iistart]:b[iiend]+1].sum()
            #count_ph_in_bursts(bursts, mask).astype(float)
        Num_ph.append(num_ph.astype(float))
    return Num_ph

