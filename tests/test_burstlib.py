# -*- coding: utf-8 -*-
"""
FRETBursts - a confocal single-molecule FRET burst analysis software.

Copyright (C) 2011-2014 Antonino Ingargiola tritemio@gmail.com

This module contains automated test for FRETBursts.
"""

import pytest

import numpy as np
from loaders import load_multispot8
from fretbursts_path_def import data_dir
import background as bg
import burstlib as bl


@pytest.fixture(scope="module")
def data():
    fn = "12d_New_30p_320mW_steer_3.dat"
    dir_ = "2013-05-15/"
    fname = data_dir+dir_+fn
    
    BT = 0.038
    gamma = 0.43
    
    d = load_multispot8(fname=fname, BT=BT, gamma=gamma)
    d.calc_bg(bg.exp_fit, time_s=30, tail_min_us=300)
    d.burst_search_t(L=10, m=10, P=None, F=7, ph_sel='DA')
    return d


##
# Test functions
#

def test_b_end_b_iend(data):
    """Test coherence between b_end() and b_iend()"""
    d = data    
    for i in xrange(d.nch):
        assert (d.ph_times_m[i][bl.b_iend(d.mburst[i])] == \
                bl.b_end(d.mburst[i])
                ).all()

def test_monotonic_burst_start(data):
    """Test for monotonic burst_start."""
    d = data
    for i in xrange(d.nch):
        assert (np.diff(bl.b_start(d.mburst[i])) > 0).all()

def test_monotonic_burst_end(data):
    """Test for monotonic burst_end."""
    d = data
    for i in xrange(d.nch):
        assert (np.diff(bl.b_end(d.mburst[i])) > 0).all()

def test_burst_start_end_size(data):
    """Test consistency between burst start, end and size"""
    d = data
    for mb in d.mburst:
        assert (mb[:, bl.iiend] == mb[:, bl.iistart]+mb[:, bl.inum_ph]-1).all()
        assert (bl.b_iend(mb) == bl.b_istart(mb) + bl.b_size(mb) - 1).all()

def test_burst_fuse_0ms(data):
    """Test that after fusing with ms=0 the sum of bursts sizes is that same
    as the number of ph in bursts (via burst selection).
    """
    d = data
    if not hasattr(d, 'fuse'):
        df = d.fuse_bursts(ms=0)
        for ph, mb in zip(df.ph_times_m, df.mburst):
            m = bl.ph_select(ph, mb)
            assert m.sum() == bl.b_size(mb).sum()
#            
if __name__ == '__main__':
    pytest.main("-x -v tests/test_burstlib.py")
