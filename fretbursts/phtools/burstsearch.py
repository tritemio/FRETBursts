#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014-2016 The Regents of the University of California,
#               Antonino Ingargiola <tritemio@gmail.com>
#
"""
The module :mod:`phtools.burstsearch` provides the low-level (or core)
burst search and photon counting functions.
This module also provides :class:`Bursts`, a container for a set of bursts.
:class:`Bursts` provides attributes for the main burst quatitites (`istart`,
`istop`, `start`, `stop`, `counts`,  `width`, etc...). It implements the
iterator interface (iterate burst by burst). Moreover :class:`Bursts` can
be indexed (`[]`, i.e. `getitem` interface) supporting the same indexing as a
numpy 1-D array.

The burst search functions return a 2-D array (burst array) of shape Nx4,
where N is the number of bursts. This array can used to build a `Bursts` object
using::

    Bursts(bursts_array)

As an example, let assume having a burst array `bursts`. To take a slice of
only the first 10 bursts you can do::

    bursts10 = bursts[:10]   # new Bursts object with the first 10 bursts

To obtain the burst start of all the bursts::

    bursts.start

To obtain the burst counts (number of photons) for the 10-th to 20-th burst::

    bursts[10:20].counts

For efficiency, when iterating over `Bursts` the returned burst is a
named tuple :class:`Burst`, which implements the same attributes as `Bursts`
(istart, istop, start, stop, counts and width).
This results in faster iteration and attribute access than using `Bursts`
objects with only one burst.

Three methods allow to transform Bursts to refer to a new timestamps array:

- :meth:`Bursts.recompute_times`
- :meth:`Bursts.recompute_index_expand`
- :meth:`Bursts.recompute_index_reduce`

Finally, in order to support fusion of consecutive bursts, we provide the class
:class:`BurstsGap` (and single-burst version :class:`BurstGap`) which add the
attributes `gap` and `gap_counts` that contains the duration and the number
of photons in gaps inside a burst. The attribute `width` is the total burst
duration minus `gap`, while `counts` is the total number of photons minus
photons falling inside gaps (gaps are open intervals, do not include edges).
"""

from __future__ import division, print_function
from builtins import range, zip

from collections import namedtuple
import numpy as np
import pandas as pd

from fretbursts.utils.misc import pprint

pd.set_option('display.max_rows', 10)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  LOW-LEVEL BURST SEARCH FUNCTIONS
#

def bsearch_py(times, L, m, T, slice_=None,
               label='Burst search', verbose=True):
    """Sliding window burst search. Pure python implementation.

    Finds bursts in the array `time` (int64). A burst starts when the photon rate
    is above a minimum threshold, and ends when the rate falls below the same
    threshold. The rate-threshold is defined by the ratio `m`/`T` (`m` photons
    in a time interval `T`). A burst is discarded if it has less than `L`
    photons.

    Arguments:
        times (array, int64): array of timestamps on which to perform the search
        L (int): minimum number of photons in a bursts. Bursts with size
            (or counts) < L are discarded.
        m (int): number of consecutive photons used to compute the rate.
        T (float): max time separation of `m` photons to be inside a burst
        slice_ (tuple): 2-element tuple used to slice times
        label (string): a label printed when the function is called
        verbose (bool): if False, the function does not print anything.

    Returns:
        Array of burst data Nx4, type int64.
        Column order is: istart, istop, start, stop.
    """
    if verbose:
        pprint('Python search (v): %s\n' % label)

    i_time0 = 0
    if slice_ is not None:
        times = times[slice_[0]:slice_[1]]
        i_time0 = slice_[0]

    bursts = []
    above_min_rate = (times[m-1:] - times[:times.size-m+1]) <= T

    it = enumerate(above_min_rate)
    for i, above_min_rate_ in it:
        if not above_min_rate_:
            continue

        i_start = i
        for i, above_min_rate_ in it:
            if not above_min_rate_:
                break
        i_stop = i + m - 2  # index of last ph in burst
        if i_stop - i_start + 1 >= L:
            bursts.append((i_start, i_stop, times[i_start], times[i_stop]))

    if above_min_rate_:
        # Correct burst-stop off by 1 when last burst does not finish
        i_stop += 1
        if i_stop - i_start + 1 >= L:
            bursts.pop()
            bursts.append((i_start, i_stop, times[i_start], times[i_stop]))

    bursts = np.array(bursts, dtype='int64')
    if bursts.size > 0:
        # Add slice offset to istart and istop
        bursts[:, :2] += i_time0
    return bursts


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Functions to count D and A photons in bursts
#

def count_ph_in_bursts(bursts, mask):
    """Counts number of photons in each burst counting only photons in `mask`.

    This function takes a :class:`Bursts` object and a boolean mask (photon
    selection) and computes the number of photons selected by the mask.
    It is used, for example, to count donor and acceptor photons
    in each burst.

    For a multi-channel version see :func:`mch_count_ph_in_bursts_py`.

    Arguments:
        bursts (Bursts object): the bursts used as input
        mask (1D boolean array): the photon mask. The boolean mask must be
            of the same size of the timestamp array used for burst search.

    Returns:
        A 1D array containing the number of photons in each burst
        counting only photons in the selection mask.
    """
    num_ph = np.zeros(bursts.num_bursts, dtype=np.int32)
    for i, burst in enumerate(bursts):
        # Counts photons between start and end of current `burst`
        num_ph[i] = mask[burst.istart : burst.istop + 1].sum()
    return num_ph


def mch_count_ph_in_bursts_py(Mburst, Mask):
    """Counts number of photons in each burst counting only photons in `Mask`.

    This multi-channel function takes a list of a :class:`Bursts` objects and
    photon masks and computes the number of photons selected by the mask
    in each channel.

    It is used, for example, to count donor and acceptor photons in
    each burst.

    For a single-channel version see :func:`count_ph_in_bursts_py`.

    Arguments:
        Mburst (list Bursts objects): a list of bursts collections, one per ch.
        Mask (list of 1D boolean arrays): a list of photon masks (one per ch),
            For each channel, the boolean mask must be of the same size of the
            timestamp array used for burst search.

    Returns:
        A list of 1D array, each containing the number of photons
        in each burst counting only photons in the selection mask.
    """
    return [np.asfarray(count_ph_in_bursts(bursts, mask))
            for bursts, mask in zip(Mburst, Mask)]


##
#  Try to import the optimized Cython functions
#

try:
    from burstsearch_c import bsearch_c
    bsearch = bsearch_c
    print(" - Optimized (cython) burst search loaded.")
except ImportError:
    bsearch = bsearch_py
    print(" - Fallback to pure python burst search.")

try:
    from burstsearch_c import mch_count_ph_in_bursts_c
    mch_count_ph_in_bursts = mch_count_ph_in_bursts_c
    print(" - Optimized (cython) photon counting loaded.")
except ImportError:
    mch_count_ph_in_bursts = mch_count_ph_in_bursts_py
    print(" - Fallback to pure python photon counting.")


class Burst(namedtuple('Burst', ['istart', 'istop', 'start', 'stop'])):
    """Container for a single burst."""
    __slots__ = ()

    @property
    def width(self):
        """Burst duration in timestamps unit."""
        return self.stop - self.start

    @property
    def counts(self):
        """Number of photons in the burst."""
        return self.istop - self.istart + 1

    @property
    def ph_rate(self):
        """Photon rate in the burst (total photon counts/duration)."""
        return self.counts / self.width


class BurstGap(namedtuple('BurstGap',
                          'istart istop start stop gap gap_counts')):
    __slots__ = ()

    @staticmethod
    def from_burst(burst):
        """Build a `BurstGap` from a :class:`Burst` object."""
        return BurstGap(start=burst.start, stop=burst.stop,
                        istart=burst.istart, istop=burst.istop,
                        gap=0, gap_counts=0)

    @property
    def width(self):
        """Burst duration in timestamps unit, minus `gap` time."""
        return self.stop - self.start - self.gap

    @property
    def counts(self):
        """Number of photons in the burst, minus `gap_counts`."""
        return self.istop - self.istart + 1 - self.gap_counts


class Bursts(object):
    """A container for burst data.

    This class provides a container for burst data. It has a
    set of attributes (`start`, `stop`, `istart`, `istop`, `counts`, `width`,
    `ph_rate`, `separation`) that can be accessed to obtain burst data.
    Only a few fundamental attributes are stored, the others are comuputed
    on-fly using python properties.

    Other attributes are `dataframe` (a `pandas.DataFrame` with the complete
    burst data), `num_bursts` (the number of bursts).

    `Bursts` objects can be built from a list of single :class:`Burst` objects
    by using the method :meth:`Bursts.from_list`, or from 2D arrays
    containing bursts data (one row per burst; columns: istart, istop, start,
    stop) such as the ones returned by burst search functions (e.g.
    :func:`bsearch_py`).

    `Bursts` object are iterable, yielding one burst a time (:class:`Burst`
    objects). `Bursts` can be compared for equality (with `==`) and copied
    (:meth:`Bursts.copy`).

    Additionally basic methods for burst manipulation are provided:

    - `recompute_times` recompute start and stop times using the current
      start and stop index and a new timestamps array passed as argument.
    - `recompute_index_*` recompute start and stop indexes to refer to an
      expanded or reduced timestamp selection.

    Other methods are:

    - `and_gate` computing burst intersection with a second set of bursts.
      Used to implement the dual-channel burst search (DCBS).

    Methods that may be implemented in the future:

    - `or_gate`: computing union with a second set of bursts.
    - `fuse_bursts`: fuse nearby bursts.

    """
    _i_istart, _i_istop, _i_start, _i_stop = 0, 1, 2, 3
    _ncols = 4

    def __init__(self, burstarray):
        if burstarray.ndim == 1:
            assert burstarray.size == self._ncols
            burstarray = burstarray[np.newaxis, :]
        self.data = burstarray

    @classmethod
    def empty(cls, num_bursts=0):
        """Return an empty `Bursts()` object."""
        return cls(np.zeros((num_bursts, cls._ncols),
                            dtype=np.int64))

    @classmethod
    def from_list(cls, bursts_list):
        """Build a new `Bursts()` object from a list of :class:`Burst`.
        """
        bursts = cls.empty(len(bursts_list))
        for i, burst in enumerate(bursts_list):
            bursts.istart[i], bursts.istop[i] = burst.istart, burst.istop
            bursts.start[i], bursts.stop[i] = burst.start, burst.stop
        return bursts

    @classmethod
    def merge(cls, list_of_bursts, sort=False):
        """Merge `Bursts` in `list_of_bursts`, returning a new `Bursts` object.
        """
        mergedata = np.vstack([b.data for b in list_of_bursts])
        if sort:
            # Sort by start times, and when equal by stop times
            indexsort = np.lexsort((mergedata[:, 3], mergedata[:, 2]))
            mergedata = mergedata[indexsort]
        return cls(mergedata)

    #
    # Basic interface
    #
    def copy(self):
        """Return a new copy of current `Bursts` object."""
        return self.__class__(self.data.copy())

    def __getitem__(self, i):
        return self.__class__(self.data[i])

    def __iter__(self):
        for bdata in self.data:
            # Bursts has a fast __init__() but relatively slow
            # attribute access:
            #yield Bursts(bdata)

            # Burst (namedtuple) has fast attribute access (>10x faster)
            # and __init__ as fast as Bursts.__init__ when avoiding the *bdata
            # syntax (which creates an additional tuple).
            #yield Burst(*bdata)
            yield Burst(bdata[0], bdata[1], bdata[2], bdata[3])

    def __eq__(self, other_bursts):
        return np.all(self.data == other_bursts.data)

    @property
    def num_bursts(self):
        """Number of bursts."""
        return self.data.shape[0]

    @property
    def size(self):
        """Number of bursts. Used for compatibility with ndarray.size.
        Use `Bursts.num_bursts` preferentially.
        """
        return self.num_bursts

    @property
    def dataframe(self):
        """A `pandas.DataFrame` containing burst data, one row per burst.
        """
        return pd.DataFrame(self.data,
                            columns=['istart', 'istop', 'start', 'stop'])

    def __repr__(self):
        return self.dataframe.__repr__()

    def _repr_html_(self):
        return self.dataframe._repr_html_()

    def join(self, bursts, sort=False):
        """Join the current `Bursts` object with another one. Returns a copy.
        """
        return self.merge([self, bursts], sort=sort)

    #
    # Burst data attributes/properties
    #
    def _set_data(self, column, values):
        """This method allows subclasses to easily extend the setters.
        """
        self.data[:, column] = values

    @property
    def start(self):
        """Time of 1st ph in each burst"""
        return self.data[:, Bursts._i_start]

    @start.setter
    def start(self, value):
        self._set_data(Bursts._i_start, value)

    @property
    def stop(self):
        """Time of last ph in each burst"""
        return self.data[:, Bursts._i_stop]

    @stop.setter
    def stop(self, value):
        self._set_data(Bursts._i_stop, value)

    @property
    def istart(self):
        """Index of 1st ph in each burst"""
        return self.data[:, Bursts._i_istart]

    @istart.setter
    def istart(self, value):
        self._set_data(Bursts._i_istart, value)

    @property
    def istop(self):
        """Index of last ph in each burst"""
        return self.data[:, Bursts._i_istop]

    @istop.setter
    def istop(self, value):
        self._set_data(Bursts._i_istop, value)

    @property
    def width(self):
        """Burst duration in timestamps units."""
        return self.stop - self.start

    @property
    def counts(self):
        """Number of photons in each burst."""
        return self.istop - self.istart + 1

    @property
    def ph_rate(self):
        """Photon rate in burst (tot size/duration)"""
        return self.counts / self.width

    @property
    def separation(self):
        """Separation between nearby bursts"""
        return self.start[1:] - self.stop[:-1]

    #
    # Burst manipulation methods
    #
    def recompute_times(self, times, out=None):
        """Recomputes start, stop times using timestamps from a new array.

        This method computes burst start, stop using the index of timestamps
        from the current object and timestamps from the passed array `times`.

        This is useful, for example, when burst search is computed on a
        "compacted" timestamps array (i.e. removing the gaps outside the
        alternation period in usALEX experiments), and afterwards the "real"
        start and stop times needs to be recomputed.

        Arguments:
            times (array): array of photon timestamps
            out (None or Bursts): if None (default), do computations on a copy
                of the current object. Otherwise, modify the `Bursts` object
                passed (can be used for in-place operations).

        Returns:
            `Bursts` object with recomputed start/stop times.
        """
        if out is None:
            out = self.copy()
        out.start = times[self.istart]
        out.stop = times[self.istop]
        return out

    def recompute_index_expand(self, mask, out=None):
        """Recompute istart and istop from selection `mask` to full timestamps.

        This method returns a new Bursts object with recomputed istart and
        istop. Old istart, istop are assumed to be index of a reduced array
        `timestamps[mask]`. New istart, istop are computed to be index of
        a "full" timestamps array of size `mask.size`.

        This is useful when performing burst search on a timestamps selection
        and we want to transform the burst data to use the index of the "full"
        timestamps array.

        Arguments:
            mask (bool array): boolean mask defining the timestamps selection
                on which the old istart and istop were computed.
            out (None or Bursts): if None (default), do computations on a copy
                of the current object. Otherwise, modify the `Bursts` object
                passed (can be used for in-place operations).

        Returns:
            `Bursts` object with recomputed istart/istop.
        """
        if out is None:
            out = self.copy()
        index = np.arange(mask.size, dtype=np.int32)
        out.istart = index[mask][self.istart]
        out.istop = index[mask][self.istop]
        return out

#    def recompute_index_reduce2(self, times_reduced, out=None):
#        """Recompute istart and istop on reduced timestamps `times_reduced`.
#
#        Extremely inefficient (but very simple!) version of
#        `Bursts.recompute_index` used for testing.
#        """
#        if out is None:
#            out = self.copy()
#
#        for i, burst in enumerate(self):
#            # The first index ([0]) accesses the tuple returned by nonzero.
#            # The second index ([0] or [-1]) accesses the array inside the
#            # tuple. This array can have size > 1 when burst start or stop
#            # happens on a repeated timestamp.
#            out[i].istart = np.nonzero(times_reduced == burst.start)[0][0]
#            out[i].istop = np.nonzero(times_reduced == burst.stop)[0][-1]
#        return out

    def recompute_index_reduce(self, times_reduced, out=None):
        """Recompute istart and istop on reduced timestamps `times_reduced`.

        This method returns a new Bursts object with same start and stop times
        and recomputed istart and istop. Old istart, istop are assumed to
        be index of a "full" timestamps array of size `mask.size`. New istart,
        istop are computed to be index of the reduced timestamps array
        `timestamps_reduced`.

        Note: it is required that all the start and stop times are
        also contained in the reduced timestamps selection.

        This method is the inverse of :meth:`recompute_index_expand`.

        Arguments:
            times_reduced (array): array of selected timestamps used to
                compute the new istart and istop. This array needs to be
                a sub-set of the original timestamps array.
            out (None or Bursts): if None (default), do computations on a copy
                of the current object. Otherwise, modify the `Bursts` object
                passed (can be used for in-place operations).

        Returns:
            `Bursts` object with recomputed istart/istop times.
        """
        if out is None:
            out = self.copy()

        # Go through the timestamps searching for start
        # and stop of each burst in order
        it = 0
        for ib, burst in enumerate(self):
            while times_reduced[it] != burst.start:
                it += 1
            out[ib].istart = it
            it_saved = it + 1

            while times_reduced[it] != burst.stop:
                it += 1
            # in case of repeated timestamps, istop needs to point
            # to the last of the repeats
            while it < times_reduced.size and times_reduced[it] == burst.stop:
                it += 1
            out[ib].istop = it - 1

            # Done with current burst, before starting a new burst,
            # reset `it` to `istart+1`
            it = it_saved
        return out

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
        start_d, stop_d = self.start, self.stop
        start_a, stop_a = bursts2.start, bursts2.stop

        bursts = []
        i_d, i_a = 0, 0
        while i_d < self.num_bursts and i_a < bursts2.num_bursts:
            # Skip any disjoint burst
            if stop_a[i_a] < start_d[i_d]:
                i_a += 1
                continue
            if stop_d[i_d] < start_a[i_a]:
                i_d += 1
                continue

            # Assign start and stop according the AND rule
            if start_a[i_a] < start_d[i_d] < stop_a[i_a]:
                start_burst = self[i_d]
            elif start_d[i_d] < start_a[i_a] < stop_d[i_d]:
                start_burst = bursts2[i_a]

            if stop_d[i_d] < stop_a[i_a]:
                end_burst = self[i_d]
                i_d += 1
            else:
                end_burst = bursts2[i_a]
                i_a += 1

            burst = start_burst.copy()
            burst.stop = end_burst.stop
            burst.istop = end_burst.istop

            bursts.append(burst)

        return Bursts.from_list(bursts)


class BurstsGap(Bursts):
    """A container for bursts with optional gaps.

    This class extend Bursts adding the attributes/properties `gap`
    (a duration) and `gap_counts` (counts in gap) that allow accounting
    for gaps inside bursts.

    """
    _i_gap, _i_gap_counts = 4, 5
    _ncols = 6

    def __init__(self, data):
        if data.shape[1] == Bursts._ncols:
            datag = np.zeros((data.shape[0], self._ncols), dtype=np.int64)
            datag[:, :Bursts._ncols] = data
            data = datag
        super(BurstsGap, self).__init__(data)

    @classmethod
    def from_list(cls, bursts_list):
        """Build a new `BurstsGap()` from a list of :class:`BurstGap`.
        """
        bursts = super(cls, cls).from_list(bursts_list)
        for i, burst in enumerate(bursts_list):
            bursts.gap[i], bursts.gap_counts[i] = burst.gap, burst.gap_counts
        return bursts

    def __iter__(self):
        for bdata in self.data:
           yield BurstGap(bdata[0], bdata[1], bdata[2], bdata[3],
                          bdata[4], bdata[5])

    @property
    def gap(self):
        """Time gap inside a burst"""
        return self.data[:, BurstsGap._i_gap]

    @gap.setter
    def gap(self, value):
        self._set_data(BurstsGap._i_gap, value)

    @property
    def gap_counts(self):
        """Number of photons falling inside gaps of each burst."""
        return self.data[:, BurstsGap._i_gap_counts]

    @gap_counts.setter
    def gap_counts(self, value):
        self._set_data(BurstsGap._i_gap_counts, value)

    @property
    def width(self):
        """Burst duration in timestamps units, minus the `gap` time.
        """
        return self.stop - self.start - self.gap

    @property
    def counts(self):
        """Number of photons in each burst, minus the `gap_counts`.
        """
        return self.istop - self.istart + 1 - self.gap_counts
