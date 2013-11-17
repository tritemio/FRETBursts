#
# WARNING: No speed increase with these functions. Use pure python version.
#
# Optimized version of burst fuse functions to be compiled in C with cython.
# To compile run: python setup.py build_ext --inplace
#

import sys
import numpy as NP
cimport numpy as NP

cdef NP.int8_t itstart, iwidth, inum_ph, iistart, iiend, itend
itstart, iwidth, inum_ph, iistart, iiend, itend = 0,1,2,3,4,5

def pprint(s):
    """Print immediately, even if inside a busy loop."""
    sys.stdout.write(s)
    sys.stdout.flush()

# Quick functions for bursts start and end
cdef NP.ndarray[NP.int64_t, ndim=1] b_start(NP.ndarray[NP.int64_t, ndim=2] b): 
    """time of 1st ph in burst"""
    return b[:,itstart]            
cdef NP.ndarray[NP.int64_t, ndim=1] b_end(NP.ndarray[NP.int64_t, ndim=2] b): 
    """time of last ph in burst"""
    return b[:,itend]             
cdef NP.ndarray[NP.int64_t, ndim=1] b_width(NP.ndarray[NP.int64_t, ndim=2] b): 
    """burst width in clk cycles"""
    return b[:,iwidth]          
cdef NP.ndarray[NP.int64_t, ndim=1] b_istart(NP.ndarray[NP.int64_t, ndim=2] b): 
    """index of 1st ph in burst"""
    return b[:,iistart]        
cdef NP.ndarray[NP.int64_t, ndim=1] b_iend(NP.ndarray[NP.int64_t, ndim=2] b): 
    """index of last ph in burst"""
    return b[:,iiend]            
cdef NP.ndarray[NP.int64_t, ndim=1] b_size(NP.ndarray[NP.int64_t, ndim=2] b): 
    """number of ph in the burst"""
    return b[:,inum_ph]          

# Separation between nearby bursts
b_separation = lambda b: b[1:,itstart] - b_end(b)[:-1]

cdef NP.ndarray[NP.int64_t, ndim=2] c_b_fuse(
		NP.ndarray[NP.int64_t, ndim=1] ph, 
		NP.ndarray[NP.int64_t, ndim=2] b, 
		NP.float_t ms, 
		NP.float_t clk_p):
    """Fuse touching or nearby bursts in the Nx6 burst array 'b'."""
    max_delay_clk = (ms*1e-3)/clk_p
    # Nearby bursts masks
    delays = (b_separation(b) <= max_delay_clk)
    first_burst = NP.hstack([delays, (False,)])
    second_burst = NP.hstack([(False,), delays])
    # Maintain just the 1st in case there were more than 2 consecutive bursts
    first_burst -= (second_burst*first_burst)
    second_burst[1:] = first_burst[:-1]
    both_burst = first_burst + second_burst

    # istart is from the first bursts, iend is from the second burst
    fused_burst1 = b[first_burst,:] # slicing makes a copy
    fused_burst2 = b[second_burst,:]
    
    #fused_burst1 = b[first_burst] # pure bool mask, no copy, b will be changed
    #fused_burst2 = b[second_burst]
    
    num_ph = b_size(fused_burst1) + b_size(fused_burst2)
    overlap = b_iend(fused_burst1) - b_istart(fused_burst2) + 1
    overlap[overlap < 0] = 0
    num_ph -= overlap
    #assert (num_ph <= (b_size(fused_burst1) + b_size(fused_burst2))).all()
    
    width = b_width(fused_burst1) + b_width(fused_burst2)
    t_overlap = b_end(fused_burst1) - b_start(fused_burst2)
    #assert (t_overlap[overlap > 0] >= 0).all()
    width[overlap > 0] -= t_overlap[overlap > 0]
    #assert (width <= (b_width(fused_burst1) + b_width(fused_burst2))).all()
    
    # Assign the new burst data
    # fused_burst1 has alredy the right tstart and istart
    fused_burst1[:,inum_ph] = num_ph
    fused_burst1[:,iwidth] = width
    fused_burst1[:,iiend] = b_iend(fused_burst2)
    fused_burst1[:,itend] = b_end(fused_burst2)
    
    new_burst = NP.vstack([fused_burst1, b[-both_burst,:]])
    reorder = new_burst[:,itstart].argsort()
    pprint(" - (C) Fused %4d of %5d bursts (%.1f%%).\n" %\
            (first_burst.sum(), b.shape[0], 100.*first_burst.sum()/b.shape[0]))
    return new_burst[reorder,:]


def c_mch_fuse_bursts(ph_times_m, MBurst, ms=0, clk_p=12.5e-9):
    mburst = [b.copy() for b in MBurst] # safety copy
    new_mburst = []
    ch = 0
    for ph, mb in zip(ph_times_m, mburst):
        ch += 1
        print " (C) - - - - - CHANNEL %2d - - - - " % ch
        if mb.size == 0: continue
        z = 0
        new_nburst, nburst = 0, 1  # starting condition
        while (new_nburst < nburst):
            z += 1
            nburst = mb.shape[0]
            mb = c_b_fuse(ph, mb, ms, clk_p)
            new_nburst = mb.shape[0]
        new_mburst.append(mb)
        pprint(" --> [CH %d] (C) END Fusing (%d iterations)\n\n" % (ch,z))
    return new_mburst

