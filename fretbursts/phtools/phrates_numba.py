"""
Numba-optimized version of functions to compute KDE-based photon rates.
"""

from __future__ import division
import numpy as np
import numba
from math import exp, fabs


@numba.jit
def kde_laplace_numba(timestamps, tau, time_axis=None):
    """Computes exponential KDE for `timestamps` evaluated at `time_axis`.
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
def kde_gaussian_numba(timestamps, tau, time_axis=None):
    """Computes Gaussian KDE for `timestamps` evaluated at `time_axis`.
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
def kde_rect_numba(timestamps, tau, time_axis=None):
    """Computes rectangular KDE for `timestamps` evaluated at `time_axis`.
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

##
# "self" functions: evaluating KDE on the same position as the timestamps
#
@numba.jit
def kde_laplace_self_numba(ph, tau):
    """Computes exponential KDE for `timestamps` evaluated at `timestamps`.
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

##
# Special functions
#
@numba.jit
def kde_laplace_nph(timestamps, tau, time_axis=None):
    """Computes exponential KDE for `timestamps` evaluated at `time_axis`.

    Computes KDE rates of `timestamps` and number of photon used to compute
    each rate. Number of photons are the one in the 10*tau range around the
    current time.

    The kernel used is a symmetric-exponential (i.e. laplace distribution)::

        kernel = exp( -|t - t0| / tau)

    The rate is computed for each time in `time_axis`.
    When ``time_axis`` is None them ``timestamps`` is used also as time axis.

    Arguments:
        timestamps (array): arrays of photon timestamps
        tau (float): time constant of the exponential kernel
        time_axis (array or None): array of time points where the rate is
            computed. If None, uses `timestamps` as time axis.

    Returns:
        2-element tuple containing

        - **rates** (*array*): the unnormalized rates (just the sum of the
          exponential kernels). To obtain rates in Hz divide the
          array by `2*tau` (or other conventional `x*tau` duration).
        - **nph** (*array*): number of photons in -5*tau..5*tau window
          for each timestamp. Proportional to the rate computed
          with KDE and rectangular kernel.
    """
    if time_axis is None:
        time_axis = timestamps
    t_size = time_axis.size
    timestamps_size = timestamps.size
    rates = np.zeros((t_size,), dtype=np.float64)
    nph = np.zeros((t_size,), dtype=np.int16)
    tau_lim = 5 * tau

    ipos, ineg = 0, 0  # indexes for timestamps
    for it, t in enumerate(time_axis):
        while ipos < timestamps_size and timestamps[ipos] - t < tau_lim:
            ipos += 1
        while ineg < timestamps_size and t - timestamps[ineg] > tau_lim:
            ineg += 1

        for itx in range(ineg, ipos):
            rates[it] += exp(-fabs(timestamps[itx] - t)/tau)
            nph[it] += 1

    return rates, nph
