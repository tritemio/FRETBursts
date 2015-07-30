#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
The module :mod:`burstsearch.burstsearchlib` provides the low-level (or core)
burst search and photon counting functions.

All the burst search functions return a 2-D array (burst array) of shape Nx6,
where N is the number of bursts. The 6 columns contain the burst data.

In order to select one or more bursts, the burst array can be indexed or
sliced along the first dimension (row wise or axis=0). The second dimension
(columns) of the burst-array should be considered opaque, and the burst-data
functions (starting with `b_`) should be used to access it (instead of the
indexing the column). Using the `b_*` functions is both clearer and less
bug-prone than using a column index.

The list of burst-data functions is: :func:`b_start`, :func:`b_end`,
:func:`b_size`, :func:`b_width`, :func:`b_istart`, :func:`b_istart`,
:func:`b_iend`, :func:`b_rate`, :func:`b_separation`. You can note that they
are more than 6 because some of them return additional quantities derived
from the burst array columns.

As an example, assume having a burst array `mburst`. To take a slice of only
the first 10 bursts you can do::

    mburst10 = mburst[:10]   # new array with burst data of the first 10 bursts

To obtain the burst start of all the bursts::

    b_start(mbursts)

To obtain the burst size (number of photons) for the 10-th to 20-th burst::

    b_size(mbursts[10:20])

"""

from __future__ import division, print_function
from builtins import range, zip

import numpy as np
from fretbursts.utils.misc import pprint


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  GLOBAL VARIABLES
#

# Define name for the 6 column indexes in the Nx6 array containing the bursts
itstart, iwidth, inum_ph, iistart, iiend, itend = 0, 1, 2, 3, 4, 5

# Quick functions for burst data (burst array)
def b_start(bursts):
    """Time of 1st ph in burst"""
    return bursts[:, itstart]

def b_end(bursts):
    """Time of last ph in burst"""
    return bursts[:, itend]

def b_width(bursts):
    """Burst width in clk cycles"""
    return bursts[:, iwidth]

def b_istart(bursts):
    """Index of 1st ph in burst"""
    return bursts[:, iistart]

def b_iend(bursts):
    """Index of last ph in burst"""
    return bursts[:, iiend]

def b_size(bursts):
    """Number of ph in the burst"""
    return bursts[:, inum_ph]

def b_ph_rate(bursts):
    """Photon rate in burst (tot size/duration)"""
    return b_size(bursts) / b_width(bursts)

def b_separation(bursts):
    """Separation between nearby bursts"""
    return b_start(bursts)[1:] - b_end(bursts)[:-1]

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  LOW-LEVEL BURST SEARCH FUNCTIONS
#

def bsearch_py(t, L, m, T, label='Burst search', verbose=True):
    """Sliding window burst search. Pure python implementation.

    Finds bursts in the array `t` (int64). A burst starts when the photon rate
    is above a minimum threshold, and ends when the rate falls below the same
    threshold. The rate-threshold is defined by the ratio `m`/`T` (`m` photons
    in a time interval `T`). A burst is discarded if it has less than `L`
    photons.

    Arguments:
        t (array, int64): array of timestamps on which to perform the search
        L (int): minimum number of photons in a bursts. Bursts with size
            (or counts) < L are discarded.
        m (int): number of consecutive photons used to compute the rate.
        T (float): max time separation of `m` photons to be inside a burst
        label (string): a label printed when the function is called
        verbose (bool): if False, the function does not print anything.

    Returns:
        2D array of burst data, one row per burst, shape (N, 6), type int64.
        Each row contains (in the order): time of burst start, burst duration,
        number of photons, index of start time, time of burst end.
        To extract burst information it's safer to use the utility functions
        `b_*` (i.e. :func:`b_start`, :func:`b_size`, :func:`b_width`, etc...).
    """
    if verbose: pprint('Python search (v): %s\n' % label)
    bursts = []
    in_burst = False
    above_min_rate = (t[m-1:] - t[:t.size-m+1]) <= T
    for i in range(t.size-m+1):
        if above_min_rate[i]:
            if not in_burst:
                i_start = i
                in_burst = True
        elif in_burst:
            # Note that i_end is the index of the last ph in the current time
            # window, while the last ph in a burst is (i_end-1).
            # Note however that the number of ph in a burst is (i_end-i_start),
            # not (i_end-1)-i_start as may erroneously appears at first glace.
            in_burst = False
            i_end = i + m - 1
            if i_end - i_start >= L:
                burst_start, burst_end = t[i_start], t[i_end-1]
                bursts.append([burst_start, burst_end-burst_start,
                               i_end-i_start, i_start, i_end-1, burst_end])
    return np.array(bursts, dtype=np.int64)


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Functions to count D and A photons in bursts
#

def count_ph_in_bursts_py(bursts, mask):
    """Counts number of photons in each burst counting only photons in `mask`.

    This function takes a burst-array and a boolean mask (photon selection)
    and computes the number of photons selected by the mask.
    It is used, for example, to count donor and acceptor photons
    in each burst.

    This is a reference implementation. In practice the multi-channel
    is always used instead (see :func:`mch_count_ph_in_bursts_py`).

    Arguments:
        bursts (2D array, int64): the burst-array.
        mask (1D boolean array): the photon mask. The boolean mask must be
            of the same size of the timestamp array used for burst search.

    Returns:
        A 1D array containing the number of photons in each burst
        counting only photons in the selection mask.
    """
    num_ph = np.zeros(bursts.shape[0], dtype=np.int16)
    for i, burst in enumerate(bursts):
        # Counts photons between start and end of current `burst`
        num_ph[i] = mask[burst[iistart] : burst[iiend]+1].sum()
    return num_ph

def mch_count_ph_in_bursts_py(Mburst, Mask):
    """Counts number of photons in each burst counting only photons in `Mask`.

    This multi-channel function takes a list of burst-arrays and a list of
    photon masks and compute the number of photons selected by the mask.
    It is used, for example, to count donor and acceptor photons in
    each burst.

    Arguments:
        Mburst (list of 2D arrays, int64): a list of burst-arrays, one per ch.
        Mask (list of 1D boolean arrays): a list of photon masks (one per ch),
            For each channel, the boolean mask must be of the same size of the
            timestamp array used for burst search.

    Returns:
        A list of 1D array, each containing the number of photons
        in each burst counting only photons in the selection mask.
    """
    Num_ph = []
    for bursts, mask in zip(Mburst, Mask):
        num_ph = np.zeros(bursts.shape[0], dtype=np.int16)
        for i, burst in enumerate(bursts):
            # Counts photons between start and end of current `burst`
            num_ph[i] = mask[burst[iistart] : burst[iiend]+1].sum()
            #count_ph_in_bursts(bursts, mask).astype(float)
        Num_ph.append(num_ph.astype(float))
    return Num_ph


##
#  Try to import the optimized Cython functions
#

try:
    from burstsearchlib_c import bsearch_c
    bsearch = bsearch_c
    print(" - Optimized (cython) burst search loaded.")
except ImportError:
    bsearch = bsearch_py
    print(" - Fallback to pure python burst search.")

try:
    from burstsearchlib_c import mch_count_ph_in_bursts_c
    mch_count_ph_in_bursts = mch_count_ph_in_bursts_c
    print(" - Optimized (cython) photon counting loaded.")
except ImportError:
    mch_count_ph_in_bursts = mch_count_ph_in_bursts_py
    print(" - Fallback to pure python photon counting.")

##
#  Additional functions processing burst data

def burst_and(bursts_d, bursts_a):
    """From 2 burst arrays return bursts defined as intersection (AND rule).

    The two input burst-arrays come from 2 different burst searches.
    Returns new bursts representing the overlapping bursts in the 2 inputs
    with start and stop defined as intersection (or AND) operator.

    The format of both input and output arrays is "burst-array" as returned
    by :func:`bsearch_py`.

    Arguments:
        bursts_d (array): burst array-1
        bursts_a (array): burst array 2. The number of burst in each of the
            input array can be different.

    Returns:
        Burst-array representing the intersection (AND) of overlapping bursts.
    """
    bstart_d, bend_d = b_start(bursts_d), b_end(bursts_d)
    bstart_a, bend_a = b_start(bursts_a), b_end(bursts_a)

    bursts = []
    burst0 = np.zeros(6, dtype=np.int64)
    i_d, i_a = 0, 0
    while i_d < bursts_d.shape[0] and i_a < bursts_a.shape[0]:
        # Skip any disjoint burst
        if bend_a[i_a] < bstart_d[i_d]:
            i_a += 1
            continue
        if bend_d[i_d] < bstart_a[i_a]:
            i_d += 1
            continue

        # Assign start and stop according the AND rule
        if bstart_a[i_a] < bstart_d[i_d] < bend_a[i_a]:
            start_burst = bursts_d[i_d]
        elif bstart_d[i_d] < bstart_a[i_a] < bend_d[i_d]:
            start_burst = bursts_a[i_a]

        if bend_d[i_d] < bend_a[i_a]:
            end_burst = bursts_d[i_d]
            i_d += 1
        else:
            end_burst = bursts_a[i_a]
            i_a += 1

        burst = burst0.copy()
        burst[itstart] = start_burst[itstart]
        burst[iistart] = start_burst[iistart]
        burst[itend] = end_burst[itend]
        burst[iiend] = end_burst[iiend]

        # Compute new width and size
        burst[iwidth] = burst[itend] - burst[itstart]
        burst[inum_ph] = burst[iiend] - burst[iistart] + 1

        bursts.append(burst)

    return np.vstack(bursts)

def recompute_burst_times(bursts, times):
    """Recomputes start, stop and width applying index data to `times`.

    This function computes burst start, stop and width using the index of
    timestamps in `bursts` and and using `times` as timestamps array.

    Arguments:
        bursts (array): input burst array
        times (array): array of photon timestamps

    Returns:
        A new burst array with recomputed "time" data.
    """
    newbursts = bursts.copy()
    for i, burst in enumerate(bursts):
        newbursts[i] = burst
        newbursts[i, itstart] = times[burst[iistart]]
        newbursts[i, itend] = times[burst[iiend]]
        newbursts[i, iwidth] = newbursts[i, itend] - newbursts[i, itstart]
    return newbursts

class Bursts():
    """A container for burst data.

    This class provides a container for burst data. It provides a large
    set of attributes (`start`, `stop`, `size`, etc...) that can be
    accessed to obtain burst data. Only a few fundamental attributes are
    stored, the others are comuputed on-fly using python properties.

    Some basic methods for burst manipulation are provided.
    `recompute_times` recompute start and stop times using the current
    start and stop index and a new timestamps array passed as argument.
    `recompute_index_*` recompute start and stop index to refer to an expanded
    or reduced timestamp selection.

    Other methods are:

    - `and_gate` computing burst intersection
    - `or_gate` (TODO, should it be added?): computing burst union
    - `fuse_bursts` (TODO, should it be added?)

    """
    _i_istart, _i_istop, _i_tstart, _i_tstop = 0, 1, 2, 3

    def __init__(self, data):
        assert data.shape[1] == 4
        self.data = np.atleast_2d(data)

    def copy(self):
        return Bursts(self.data.copy())

    @property
    def num_bursts(self):
        return self.data.shape[0]

    def __getitem__(self, i):
        return Bursts(self.data[i])

    def __iter__(self):
        for i in range(self.num_bursts):
            yield self[i]

    @property
    def start(self):
        """Time of 1st ph in each burst"""
        return self.data[:, Bursts._i_tstart]

    @start.setter
    def start(self, value):
        self.data[:, Bursts._i_tstart] = value

    @property
    def stop(self):
        """Time of last ph in each burst"""
        return self.data[:, Bursts._i_tstop]

    @stop.setter
    def stop(self, value):
        self.data[:, Bursts._i_tstop] = value

    @property
    def istart(self):
        """Index of 1st ph in each burst"""
        return self.data[:, Bursts._i_istart]

    @istart.setter
    def istart(self, value):
        self.data[:, Bursts._i_istart] = value

    @property
    def istop(self):
        """Index of last ph in each burst"""
        return self.data[:, Bursts._i_istop]

    @istop.setter
    def istop(self, value):
        self.data[:, Bursts._i_istop] = value

    @property
    def width(self):
        return self.stop - self.start

    @property
    def size(self):
        return self.istop - self.istart + 1

    @property
    def ph_rate(self):
        """Photon rate in burst (tot size/duration)"""
        return self.size / self.width

    @property
    def separation(self):
        """Separation between nearby bursts"""
        return self.start[1:] - self.stop[:-1]

    def recompute_times(self, times):
        """Recomputes start, stop times applying index data to `times`.

        This method computes burst start, stop using the index of
        timestamps in `bursts` and and using `times` as timestamps array.

        Arguments:
            times (array): array of photon timestamps

        Returns:
            A new Bursts object with recomputed start/stop times.
        """
        newbursts = self.copy()
        newbursts.start = times[self.istart]
        newbursts.stop = times[self.istart]
        return newbursts

    def recompute_index_expand(self, mask):
        """Recompute istart and istop from selection `mask` to full timestamps.

        This method returns a new Bursts object with same start and stop times
        and recomputed istart and istop. Old istart, istop are assumed to
        be index of a reduced array `timestamps[mask]`. New istart, istop
        are computed to be index of a "full" timestamps array of size
        `mask.size`.

        Arguments:
            mask (bool array): boolean mask defining the timestamp selection
                on which the old istart and istop were computed.

        Returns:
            A new Bursts object with recomputed istart/istop.
        """
        newbursts = self.copy()
        index = np.arange(mask.size, dtype=np.int32)
        newbursts.istart = index[mask][self.istart]
        newbursts.istop = index[mask][self.istop]
        return newbursts

    def recompute_index_reduce(self, mask):
        """Recompute istart and istop on reduced timestamps selection `mask`.

        This method returns a new Bursts object with same start and stop times
        and recomputed istart and istop. Old istart, istop are assumed to
        be index of a "full" timestamps array of size `mask.size`. New istart,
        istop are computed to be index of a reduced array `timestamps[mask]`.

        Note that it is required that all the start and stop times are
        also contained in the reduced timestamps selection.

        Arguments:
            mask (bool array): boolean mask defining the timestamp selection
                on which the new istart and istop are computed.

        Returns:
            A new Bursts object with recomputed istart/istop times.
        """
        newbursts = self.copy()
        ## Untested, to be checked
        newbursts.istart = np.nonzero(mask[self.istart])[0]
        newbursts.istop = np.nonzero(mask[self.istop])[0]

        # Check that we are not missing any start or stop index
        assert newbursts.istart.size == self.istart.size
        assert newbursts.istop.size == self.istop.size

        return newbursts

    def and_gate(self, bursts2):
        """From 2 burst arrays return bursts defined as intersection (AND rule).

        The two input burst-arrays come from 2 different burst searches.
        Returns new bursts representing the overlapping bursts in the 2 inputs
        with start and stop defined as intersection (or AND) operator.

        Both input and output are Bursts objects.

        Arguments:
            bursts_a (Bursts object): second set of bursts to be intersected
                with bursts in self. The number of bursts in self and
                `bursts_a` can be different.

        Returns:
            Bursts object containing intersections (AND) of overlapping bursts.
        """
        bstart_d, bend_d = self.start, self.end
        bstart_a, bend_a = bursts2.start, bursts2.end

        bursts = []
        i_d, i_a = 0, 0
        while i_d < self.num_bursts and i_a < bursts2.num_bursts:
            # Skip any disjoint burst
            if bend_a[i_a] < bstart_d[i_d]:
                i_a += 1
                continue
            if bend_d[i_d] < bstart_a[i_a]:
                i_d += 1
                continue

            # Assign start and stop according the AND rule
            if bstart_a[i_a] < bstart_d[i_d] < bend_a[i_a]:
                start_burst = self[i_d]
            elif bstart_d[i_d] < bstart_a[i_a] < bend_d[i_d]:
                start_burst = bursts2[i_a]

            if bend_d[i_d] < bend_a[i_a]:
                end_burst = self[i_d]
                i_d += 1
            else:
                end_burst = bursts2[i_a]
                i_a += 1

            burst = start_burst.copy()
            burst.stop = end_burst.stop
            burst.istop = end_burst.istop

            bursts.append(burst)

        return Bursts(np.vstack(bursts.data))

