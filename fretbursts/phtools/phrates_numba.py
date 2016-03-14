"""
Numba-optimized version of functions to compute photon rates.
"""

from __future__ import division
import numpy as np
import numba
from math import exp, fabs


@numba.jit
def _kde_laplace_numba(timestamps, tau, time_axis=None):
    """Computes exponential KDE using numba optimization.
    """
    if time_axis is None:
        time_axis = timestamps
    t_size = time_axis.size
    timestamps_size = timestamps.size
    rates = np.zeros((t_size,), dtype=np.float64)
    tau_lim = 5 * tau

    ipos, ineg = 0, 0  # indexes for timestamps
    for it, t in enumerate(time_axis):
        while ipos < timestamps_size and timestamps[ipos] - t < tau_lim:
            ipos += 1
        while ineg < timestamps_size and t - timestamps[ineg] > tau_lim:
            ineg += 1

        for itx in range(ineg, ipos):
            rates[it] += exp(-fabs(timestamps[itx] - t)/tau)

    return rates


@numba.jit
def _kde_gaussian_numba(timestamps, tau, time_axis=None):
    """Computes Gaussian KDE using numba optimization.
    """
    if time_axis is None:
        time_axis = timestamps
    timestamps_size = timestamps.size
    rates = np.zeros((time_axis.size,), dtype=np.float64)
    tau_lim = 3 * tau   # 3 tau = 99.7 % of the Gaussian
    tau2 = 2 * (tau**2)

    ipos, ineg = 0, 0  # indexes for timestamps
    for it, t in enumerate(time_axis):
        while ipos < timestamps_size and timestamps[ipos] - t < tau_lim:
            ipos += 1
        while ineg < timestamps_size and t - timestamps[ineg] > tau_lim:
            ineg += 1

        for itx in range(ineg, ipos):
            rates[it] += exp(-((timestamps[itx] - t)**2)/tau2)

    return rates


@numba.jit
def _kde_rect_numba(timestamps, tau, time_axis=None):
    """Computes rectangular KDE using numba optimization.
    """
    if time_axis is None:
        time_axis = timestamps
    timestamps_size = timestamps.size
    rates = np.zeros((time_axis.size,), dtype=np.float64)
    tau_lim = tau / 2

    ipos, ineg = 0, 0  # indexes for timestamps
    for it, t in enumerate(time_axis):
        while ipos < timestamps_size and timestamps[ipos] - t < tau_lim:
            ipos += 1
        while ineg < timestamps_size and t - timestamps[ineg] > tau_lim:
            ineg += 1

        rates[it] = ipos - ineg

    return rates

@numba.jit
def _kde_laplace_self_numba(ph, tau):
    """Computes exponential KDE for `timestamps` evaluated at `time_axis`.
    """
    ph_size = ph.size
    ipos, ineg = 0, 0
    rates = np.zeros((ph_size,), dtype=np.float64)
    nph = np.zeros((ph_size,), dtype=np.int16)
    tau_lim = 5*tau
    for i, t in enumerate(ph):
        while ipos < ph_size and ph[ipos] - t < tau_lim:
            ipos += 1

        while t - ph[ineg] > tau_lim:
            ineg += 1

        for itx in range(ineg, ipos):
            rates[i] += exp(-fabs(ph[itx]-t)/tau)
            nph[i] += 1

    return rates, nph
