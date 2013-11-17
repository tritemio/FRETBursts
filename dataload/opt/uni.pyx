import sys
import numpy as NP
cimport numpy as NP

def unwind_uni_o(NP.ndarray[NP.int32_t, ndim=1] times, 
        NP.ndarray[NP.uint8_t, ndim=1] det, 
        NP.int32_t nch=8, times_nbit=28, debug=True):
    """64bit conversion and merging of corresponding D/A channels.
    """
    cdef NP.int32_t ii, i_d, i_a, sd, sa
    cdef NP.int32_t ich, det_d, det_a
    cdef NP.ndarray[NP.int64_t, ndim=1] T, t_d, t_a
    cdef NP.ndarray[NP.int8_t, ndim=1] A
    cdef NP.ndarray[NP.int32_t, ndim=1] t_d1, t_a1
    cdef NP.int64_t v

    diff = lambda a: a[1:]-a[:-1]

    num_spad = nch*2
    expo = (2**times_nbit)
    ph_times_m, A_det = [[]]*nch, [[]]*nch
    for ich in xrange(nch):
        det_d, det_a = ich+1, ich+1+nch
        t_d1, t_a1 = times[det == det_d], times[det == det_a]
        t_d = t_d1.astype(NP.int64) + \
                NP.hstack([0,NP.cumsum((diff(t_d1)<0),dtype=NP.int16)])*expo
        t_a = t_a1.astype(NP.int64) + \
                NP.hstack([0,NP.cumsum((diff(t_a1)<0),dtype=NP.int16)])*expo
        del t_d1, t_a1

        ## Merge sort
        T = NP.zeros(t_d.size+t_a.size, dtype=NP.int64)
        A = NP.zeros(T.size, dtype=NP.int8)
        i_d, i_a = 0,0
        sd, sa = t_d.size, t_a.size
        for ii in range(T.size):
            if i_a == sa or (i_d < sd and (t_d[i_d] <= t_a[i_a])):
                v = t_d[i_d]
                i_d += 1
            #if i_a == t_a.size or (i_d < t_d.size and (t_d[i_d] <= t_a[i_a])):
            else:
                v = t_a[i_a]
                i_a +=1
                A[ii] = 1
            T[ii] = v
        ph_times_m[ich] = T
        A_det[ich] = A	
    return ph_times_m, A_det

