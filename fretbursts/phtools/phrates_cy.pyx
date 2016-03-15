"""
Cython version of functions to compute photon rates.
"""

from __future__ import division
import numpy as np
cimport numpy as np
from libc.math cimport exp, fabs


ctypedef np.int64_t DTYPE_t


cdef _kde_gaussian_cy(DTYPE_t[:] timestamps, DTYPE_t tau, DTYPE_t[:] time_axis):
    """Computes Gaussian KDE for `timestamps` evaluated at `time_axis`.
    """
    cdef np.int64_t timestamps_size, ipos, ineg, it, itx
    cdef np.float64_t[:] rates
    cdef np.float64_t tau2
    cdef DTYPE_t tau_lim, t

    timestamps_size = timestamps.size
    rates = np.zeros((time_axis.size,), dtype=np.float64)
    tau_lim = 3 * tau   # 3 tau = 99.7 % of the Gaussian
    tau2 = 2 * (tau**2)

    ipos, ineg = 0, 0  # indexes for timestamps
    for it in range(time_axis.size):
        t = time_axis[it]
        while ipos < timestamps_size and timestamps[ipos] - t < tau_lim:
            ipos += 1
        while ineg < timestamps_size and t - timestamps[ineg] > tau_lim:
            ineg += 1

        for itx in range(ineg, ipos):
            rates[it] += exp(-((timestamps[itx] - t)**2)/tau2)
    return rates

cdef _kde_laplace_cy(DTYPE_t[:] timestamps, DTYPE_t tau, DTYPE_t[:] time_axis):
    """Computes exponential KDE for `timestamps` evaluated at `time_axis`.
    """
    cdef np.int64_t timestamps_size, ipos, ineg, it, itx
    cdef np.float64_t[:] rates
    cdef DTYPE_t tau_lim, t

    timestamps_size = timestamps.size
    rates = np.zeros((time_axis.size,), dtype=np.float64)
    tau_lim = 5 * tau

    ipos, ineg = 0, 0  # indexes for timestamps
    for it in range(time_axis.size):
        t = time_axis[it]
        while ipos < timestamps_size and timestamps[ipos] - t < tau_lim:
            ipos += 1
        while ineg < timestamps_size and t - timestamps[ineg] > tau_lim:
            ineg += 1

        for itx in range(ineg, ipos):
            rates[it] += exp(-fabs(timestamps[itx] - t)/tau)
    return rates


cdef _kde_rect_cy(DTYPE_t[:] timestamps, DTYPE_t tau, DTYPE_t[:] time_axis):
    """Computes rectangular KDE for `timestamps` evaluated at `time_axis`.
    """
    cdef np.int64_t timestamps_size, ipos, ineg, it
    cdef np.float64_t[:] rates
    cdef DTYPE_t tau_lim, t

    timestamps_size = timestamps.size
    rates = np.zeros((time_axis.size,), dtype=np.float64)
    tau_lim = tau // 2

    ipos, ineg = 0, 0  # indexes for timestamps
    for it in range(time_axis.size):
        t = time_axis[it]
        while ipos < timestamps_size and timestamps[ipos] - t < tau_lim:
            ipos += 1
        while ineg < timestamps_size and t - timestamps[ineg] > tau_lim:
            ineg += 1

        rates[it] = ipos - ineg
    return rates

def kde_gaussian_cy(timestamps, tau, time_axis=None):
    """Computes Gaussian KDE for `timestamps` evaluated at `time_axis`.
    """
    if time_axis is None:
        time_axis = timestamps
    return np.asarray(_kde_gaussian_cy(timestamps, tau, time_axis))

def kde_laplace_cy(timestamps, tau, time_axis=None):
    """Computes exponential KDE for `timestamps` evaluated at `time_axis`.
    """
    if time_axis is None:
        time_axis = timestamps
    return np.asarray(_kde_laplace_cy(timestamps, tau, time_axis))

def kde_rect_cy(timestamps, tau, time_axis=None):
    """Computes rectangular KDE for `timestamps` evaluated at `time_axis`.
    """
    if time_axis is None:
        time_axis = timestamps
    return np.asarray(_kde_rect_cy(timestamps, tau, time_axis))
