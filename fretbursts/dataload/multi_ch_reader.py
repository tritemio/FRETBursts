#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014-2015 Antonino Ingargiola <tritemio@gmail.com>
#
"""
This module contains routines to load and preprocess timestamp data
saved by the NI board (multi-ch data) using two 32bit words (timestamp
and detector) for each photon.

This data is produced by the 4 and 8-spot smFRET setup, first generation.
"""

from __future__ import absolute_import, division
from builtins import range, zip

import os
import numpy as np

from ..utils.misc import pprint


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  DATA LOADING
#

def read_int32_int32_file(fname, n_bytes_to_read=-1):
    """Read the data file with 32+32 bit format (int32 version)."""
    try:
        f = open(fname, 'rb')
    except IOError:
        fname += '.dat'
        f = open(fname, 'rb')

    # Reading the header
    lines = [f.readline() for _ in range(3)]
    words_per_photon = lines[1].split()[-1]
    assert words_per_photon == b'2'

    #  Reading data in int32
    bytes_in_file = os.path.getsize(fname) - f.tell()
    if n_bytes_to_read < 4:
        n_bytes_to_read = bytes_in_file
    N_bytes = (int(min(n_bytes_to_read, bytes_in_file)) // 4) * 4
    data = np.ndarray(shape=(N_bytes // 4,), dtype='>i4',
                      buffer=f.read(N_bytes))
    detector = data[::2] + 1
    ph_times = (data[1::2] - data[1])
    assert ((detector < 17) * (detector >= 0)).all()
    return ph_times, detector.astype('uint8')


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  DATA CONVERSION
#

def swap_donor_acceptor(detectors, nch=4):
    """Swap the donor and the acceptor channels."""
    donors = detectors > nch
    acceptors = -donors
    det_d = detectors[donors]
    det_a = detectors[acceptors]
    detectors[donors] = det_d - nch
    detectors[acceptors] = det_a + nch
    return detectors

def load_data_ordered16(fname, n_bytes_to_read=-1, nch=8, swap_D_A=False,
                        remap_D=False, remap_A=False, mute=False):
    """Load data, unroll the 32bit overflow and order in increasing order."""
    pprint(' - Loading data "%s" ... ' % fname, mute)
    ph_times, detector = read_int32_int32_file(fname, n_bytes_to_read)
    pprint(" [DONE]\n", mute)
    pprint(" - Processing data ... ", mute)
    if remap_D:
        pprint("\n   - Inverting DONOR pixels order ... ", mute)
        detector[detector <= nch] = nch+1-detector[detector <= nch]
        pprint(" [DONE]\n", mute)
    if remap_A:
        pprint("\n   - Inverting ACCEPTOR pixels order ... ", mute)
        detector[detector > nch] = 2*nch+1-(detector[detector > nch] - nch)
        pprint(" [DONE]\n", mute)
    if swap_D_A:
        pprint("\n   - Swapping D and A channels ... ", mute)
        detector = swap_donor_acceptor(detector, nch=8)
        pprint(" [DONE]\n", mute)
    ph_times_m, red, ph_times_ma = unwind_uni(ph_times, detector)
    pprint("   [DONE Processing]\n", mute)

    return ph_times_m, red, ph_times_ma

def unwind_uni(times, det, nch=8, times_nbit=28, debug=True):
    """64bit conversion and merging of corresponding D/A channels."""
    diff = lambda a: a[1:]-a[:-1]
    ts_max = 2**times_nbit
    num_spad = nch*2

    times_ma = [times[det==d].astype('int64') for d in range(1, num_spad+1)]
    for i, t in enumerate(times_ma):
        ## This debug check is a tautology (cumsum of a bool mask >= 0)
        #if debug: assert (cumsum((diff(t) < 0)) >= 0).all()
        t[1:] += np.cumsum((diff(t) < 0), dtype='int64')*ts_max
        if debug: assert (np.diff(t) > 0).all()
        #print i, (t < 0).sum(), (diff(t) < 0).sum(), "\n"

    ph_times_m, red = nch*[0], nch*[0]
    bones = lambda n: np.ones(n, dtype=bool)
    bzeros = lambda n: np.zeros(n, dtype=bool)
    for i in range(nch):
        # Merge channels and sort
        ph_times_m[i] = np.hstack([times_ma[i], times_ma[i+nch]])
        index_sort = ph_times_m[i].argsort()
        red[i] = np.hstack([bzeros(times_ma[i].size),
                            bones(times_ma[i+nch].size)])
        red[i] = red[i][index_sort]
        ph_times_m[i] = ph_times_m[i][index_sort]
    return ph_times_m, red, times_ma


def unwind_uni_c(times, det, nch=8, times_nbit=28, debug=True):
    """64bit conversion and merging of corresponding D/A channels.
    This version is identical to the cython version but uses numpy sort.
    """
    cumsum, hstack, int16, int64 = np.cumsum, np.hstack, np.int16, np.int64
    diff = lambda a: a[1:]-a[:-1]
    bones = lambda n: np.ones(n, dtype=bool)
    bzeros = lambda n: np.zeros(n, dtype=bool)

    num_spad = nch*2
    ts_max = 2**times_nbit
    ph_times_m, A_det = [[]]*nch, [[]]*nch
    for ich in range(nch):
        det_d, det_a = ich+1, ich+1+nch
        t_d1, t_a1 = times[det == det_d], times[det == det_a]
        t_d = t_d1.astype(int64)
        t_d += hstack([0, cumsum((diff(t_d1)<0), dtype=int16)])*ts_max
        t_a = t_a1.astype(int64)
        t_a += hstack([0, cumsum((diff(t_a1)<0), dtype=int16)])*ts_max
        del t_d1, t_a1

        T = hstack([t_d, t_a])
        index_sort = T.argsort(kind='mergesort')

        ph_times_m[ich] = T[index_sort]
        A_det[ich] = hstack([bzeros(t_d.size), bones(t_a.size)])[index_sort]
        del index_sort, T
    return ph_times_m, A_det

def unwind_uni_o(times, det, nch=8, times_nbit=28, debug=True):
    """64bit conversion and merging of corresponding D/A channels.
    This version is about 10% faster than unwind_uni().
    """
    diff = lambda a: a[1:]-a[:-1]
    bones = lambda n: np.ones(n, dtype=bool)
    bzeros = lambda n: np.zeros(n, dtype=bool)

    num_spad = nch*2
    expo = (2**times_nbit)
    ph_times_m, A_det = [[]]*nch, [[]]*nch
    for ich in range(nch):
        det_d, det_a = ich+1, ich+1+nch
        #a d_mask, a_mask = (det == det_d), (det == det_a)

        #t_d, t_a = times[d_mask].astype(int64), times[a_mask].astype(int64)
        #t_d[1:] += cumsum((diff(t_d) < 0), dtype=int64)*expo
        #t_a[1:] += cumsum((diff(t_a) < 0), dtype=int64)*expo

        #a t_d, t_a = times[d_mask], times[a_mask]
        t_d, t_a = times[det == det_d], times[det == det_a]
        # overflows must be < 2^16
        step_d = np.cumsum((diff(t_d) < 0), dtype='int16')
        step_a = np.cumsum((diff(t_a) < 0), dtype='int16')

        t1 = np.hstack([t_d, t_a]).astype('int64')
        t2 = np.hstack([0, step_d, 0, step_a])
        #ph_times_m[ich] = NE.evaluate("t1+t2*expo") # this is slower!
        ph_times_m[ich] = t1 + t2*expo
        del t1, t2
        #ph_times_m[ich] = concatenate([t_d,t_a]).astype(int64) + \
        #        hstack([0,step_d,0,step_a])*expo
        #ph_times_m[ich] = concatenate([t_d, t_a])
        index_sort = ph_times_m[ich].argsort(kind='mergesort')
        ph_times_m[ich] = ph_times_m[ich][index_sort]
        A_det[ich] = np.hstack([bzeros(t_d.size), bones(t_a.size)])[index_sort]
        del index_sort
    return ph_times_m, A_det


if __name__ == '__main__':
    ## Some test

    #fname = gui_fname()
    #ph_times, det = read_int32_int32_file(fname)

    ph_m, a_em, t_ma = unwind_uni(ph_times, det)    # plain
    ph_mc, a_emc = unwind_uni_c(ph_times, det)      # 15% faster than plain
    ph_mo, a_emo = unwind_uni_o(ph_times, det)      # 15% faster than plain

    print([(ph == phc).all() for ph, phc in zip(ph_m, ph_mc)])
    print([(ae == aec).all() for ae, aec in zip(a_em, a_emc)])

    ## NOTE: write a compare function that takes into account same timestamps
    ##       in donor and acceptor ch.
