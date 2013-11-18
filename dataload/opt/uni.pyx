import sys
import numpy as np
cimport numpy as np

def unwind_uni_o(np.ndarray[np.int32_t, ndim=1] times, 
        np.ndarray[np.uint8_t, ndim=1] det, 
        np.int32_t nch=8, times_nbit=28, debug=True):
    """64bit conversion and merging of corresponding D/A channels.
    """
    cdef np.int32_t ii, i_d, i_a, sd, sa
    cdef np.int32_t ich, det_d, det_a
    cdef np.ndarray[np.int64_t, ndim=1] T, t_d, t_a
    cdef np.ndarray[np.int8_t, ndim=1] A
    cdef np.ndarray[np.int32_t, ndim=1] t_d1, t_a1
    cdef np.int64_t v

    cumsum, hstack, int16, int64 = np.cumsum, np.hstack, np.int16, np.int64
    diff = lambda a: a[1:]-a[:-1]

    num_spad = nch*2
    ts_max = (2**times_nbit)
    ph_times_m, A_det = [[]]*nch, [[]]*nch
    for ich in xrange(nch):
        det_d, det_a = ich+1, ich+1+nch
        t_d1, t_a1 = times[det == det_d], times[det == det_a]
        t_d = t_d1.astype(int64)
        t_d += hstack([0, cumsum((diff(t_d1)<0), dtype=int16)])*ts_max
        t_a = t_a1.astype(int64) + \
        t_a += hstack([0, cumsum((diff(t_a1)<0), dtype=int16)])*ts_max
        del t_d1, t_a1

        ## Merge sort
        T = np.zeros(t_d.size+t_a.size, dtype=int64)
        A = np.zeros(T.size, dtype=np.int8)
        i_d, i_a = 0,0
        sd, sa = t_d.size, t_a.size
        for ii in xrange(T.size):
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

