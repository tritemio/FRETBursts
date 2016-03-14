
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014-2016 Antonino Ingargiola <tritemio@gmail.com>
#
"""
This module provides functions to compute photon rates from timestamps
arrays. Different methods to compute rates are implemented:

1. Consecutive set of *m* timestamps ("sliding m-tuple") [TODO, move from burstlib]
2. KDE-based methods with Gaussian or Laplace distribution or rectangular
   kernels.

It also provide a few specific functions useful for the implementation
of the 2CDE method [1].

Note:
    When using of "sliding m-tuple" method (1), rates can be only
    computed for each consecutive set of *m* timestamps. The time-axis can be
    computed from the mean timestamp in each m-tuple.

    When using the KDE method, rates can be computed at any time point.
    Practically, the time points at which rates are computed are timestamps
    (in a photon stream). In other words, we don't normally use a uniformly
    sampled time axis but we use a timestamps array as time axis for the rate.

    Note that computing rates with a fixed sliding time window and sampling
    the function by centering the window on each timestamp is equivalent to
    a KDE-based rate computation using a rectangular kernel.

References:
    [1] Tomov et al. "Disentangling Subpopulations in Single-Molecule FRET ..."
        Biophys. J. (2012). doi:10.1016/j.bpj.2011.11.4025.

"""

import numpy as np
import numba
from math import exp, fabs


##
# Functions to compute rates using m-tuple of photon timestamps
#
def mtuple_delays(ph, m):
    """Compute array of m-photons delays of size ph.size - m + 1."""
    return ph[m-1:] - ph[:ph.size-m+1]

def mtuple_delays_min(ph, m):
    """Compute the min m-photons delay in `ph`."""
    if ph.size < m:
        return None
    else:
        return mtuple_delays(ph=ph, m=m).min()

def mtuple_rates(ph, m):
    """Compute array of m-photons rates of size ph.size - m + 1."""
    return m/(ph[m-1:] - ph[:ph.size-m+1])

def mtuple_rates_t(ph, m):
    """Compute mean time for each rate computed by `mtuple_rates`."""
    return 0.5*(ph[m-1:] + ph[:ph.size-m+1])  # time for rate

def mtuple_rates_max(ph, m):
    """Compute max m-photon rate in `ph`."""
    if ph.size < m:
        return None
    else:
        return mtuple_rates(ph=ph, m=m).max()


##
# General purpose functions to compute rates
#
@numba.jit
def _kde_laplace_numba(timestamps, tau, time_axis=None):
    """Computes exponential KDE for `timestamps` evaluated at `time_axis`.

    Computes KDE rates of `timestamps` using a symmetric-exponential kernel
    (i.e. laplace distribution)::

        kernel = exp( -|t - t0| / tau)

    The rate is computed for each time in `time_axis`.
    When ``time_axis`` is None them ``timestamps`` is used also as time axis.
    For a similar function returning also the number of photons used to
    compute each rate see :func:`kde_laplace_nph`.

    Arguments:
        timestamps (array): arrays of photon timestamps
        tau (float): time constant of the exponential kernel
        time_axis (array or None): array of time points where the rate is
            computed. If None, uses `timestamps` as time axis.

    Returns:
        rates (array): non-normalized rates (just the sum of the
        exponential kernels). To obtain rates in Hz divide the
        array by `2*tau` (or other conventional x*tau duration).
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
    """Computes Gaussian KDE for `timestamps` evaluated at `time_axis`.

    Computes KDE rates of `timestamps` using a Gaussian kernel.

    The rate is computed for each time in `time_axis`.
    When ``time_axis`` is None them ``timestamps`` is used also as time axis.

    Arguments:
        timestamps (array): arrays of photon timestamps
        tau (float): sigma of the Gaussian kernel
        time_axis (array or None): array of time points where the rate is
            computed. If None, uses `timestamps` as time axis.

    Returns:
        rates (array): non-normalized rates (just the sum of the
        Gaussian kernels). To obtain rates in Hz divide the
        array by `2.5*tau`.
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
    """Computes KDE with rect kernel for `timestamps` evaluated at `time_axis`.

    Computes KDE rates of `timestamps` using a rectangular kernel.

    The rate is computed for each time in `time_axis`.
    When ``time_axis`` is None them ``timestamps`` is used also as time axis.

    Arguments:
        timestamps (array): arrays of photon timestamps
        tau (float): duration of the rectangular kernel
        time_axis (array or None): array of time points where the rate is
            computed. If None, uses `timestamps` as time axis.

    Returns:
        rates (array): non-normalized rates (just the sum of the
        rectangular kernels). To obtain rates in Hz divide the
        array by `tau`.
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

kde_laplace = _kde_laplace_numba
kde_gaussian = _kde_gaussian_numba
kde_rect = _kde_rect_numba

##
# Functions evaluating rates at the same location as timestamps
#
def _kde_laplace_self(ph, tau):
    """Computes exponential KDE for each photon in `ph`.

    This function computes the rate of timestamps in `ph`
    using a KDE and symmetric-exponential kernel (i.e. laplace distribution)::

        kernel = exp( -|t - t0| / tau)

    The rate is evaluated for each element in `ph` (that's why name ends
    with ``_self``).

    Arguments:
        ph (array): arrays of photon timestamps
        tau (float): time constant of the exponential kernel

    Returns:
        rates, nph (arrays): the unnormalized rates (just the sum of the
        exponential kernels). To obtain rates in Hz divide the
        array by `2*tau` (or other conventional x*tau duration).
        nph is the number of photons in -5*tau..5*tau window
        for each timestamp. Proportional to the rate computed
        with KDE and rectangular kernel.
        """
    ph_size = ph.size
    ipos, ineg = 0, 0
    rates = np.zeros((ph_size,), dtype=np.float64)
    nph = np.zeros((ph_size,), dtype=np.int16)
    tau_lim = 5*tau
    for i, t in enumerate(ph):
        # Increment ipos until falling out of N*tau (tau_lim)
        # ipos is the first value *outside* the limit
        while ipos < ph_size and ph[ipos] - t < tau_lim:
            ipos += 1

        # Increment ineg until falling inside N*tau (tau_lim)
        # ineg is the first value *inside* the limit
        while t - ph[ineg] > tau_lim:
            ineg += 1

        delta_t = np.abs(ph[ineg:ipos] - t)
        rates[i] = np.exp(-delta_t / tau).sum()
        nph[i] = ipos - ineg
    return rates, nph


@numba.jit
def _kde_laplace_self_numba(ph, tau):
    """Numba version of `kde_laplace_self()`
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
# Custom functions needed for 2CDE
#
@numba.jit
def kde_laplace_nph(timestamps, tau, time_axis=None):
    """Computes exponential KDE for `timestamps` evaluated at `time_axis`.

    Computes KDE rates of `timestamps` and number of photon used to compute
    each rate. Number of photons are the one in the 10*tau range around the
    current time. This function is used when computing the 2CDE estimator
    nbKDE which requires both the photon rates the number of photons.
    For a similar function computing only the rates see :func:`kde_laplace`.

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
        rates, nph (arrays): the unnormalized rates (just the sum of the
        exponential kernels). To obtain rates in Hz divide the
        array by `2*tau` (or other conventional x*tau duration).
        nph is the number of photons in -5*tau..5*tau window
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

@numba.jit
def kde_nbKDE(timestamps, tau):
    """Computes nbKDE for `timestamps`. See Tomov et al. BJ 2012.
    """
    timestamps_size = timestamps.size
    rates = np.zeros((timestamps_size,), dtype=np.float64)
    nph = np.zeros((timestamps_size,), dtype=np.int16)
    tau_lim = 5 * tau

    ipos, ineg = 0, 0  # indexes for timestamps
    for it, t in enumerate(timestamps):
        while ipos < timestamps_size and timestamps[ipos] - t < tau_lim:
            ipos += 1
        while t - timestamps[ineg] > tau_lim:
            ineg += 1

        for itx in range(ineg, ipos):
            rates[it] += exp(-fabs(timestamps[itx] - t)/tau)
            nph[it] += 1

    nbkde = (1 + 2/nph) * (rates - 1)
    return nbkde, rates, nph
