#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2017 Antonino Ingargiola <tritemio@gmail.com>
#
"""
Unit tests for Burst, Bursts, BurstGap, BurstsGap.

Running the tests requires `py.test`.
"""

from __future__ import division
from builtins import range, zip

import numpy as np
import pytest

from fretbursts.phtools import burstsearch as bs


def test_Burst():
    burst = bs.Burst(0, 1, 100, 200)
    assert burst.counts == 2
    assert burst.width == 100
    assert burst.ph_rate == 2 / 100


def test_BurstGap():
    burst = bs.BurstGap(0, 2, 100, 300, 100, 1)
    assert burst.counts == 2
    assert burst.width == 100
    assert burst.gap == 100
    assert burst.gap_counts == 1


def test_Bursts():
    barray = np.zeros((100, 4), dtype=np.int64)
    barray[:, 0] = np.arange(100) * 10
    barray[:, 1] = barray[:, 0] + 5
    barray[:, 2] = np.arange(100) * 100
    barray[:, 3] = barray[:, 2] + 50

    bursts = bs.Bursts(barray)
    assert bursts.num_bursts == 100
    assert (bursts.counts == 6).all()
    assert (bursts.width == 50).all()

    # Test indexing and slicing
    assert bursts[0] == bs.Bursts(barray[0])
    assert bursts[:1] == bursts[0]
    assert (bursts[:10:2].data == barray[:10:2]).all()


def test_BurstsGap():
    barray = np.zeros((100, 6), dtype=np.int64)
    barray[:, 0] = np.arange(100) * 10
    barray[:, 1] = barray[:, 0] + 6
    barray[:, 2] = np.arange(100) * 100
    barray[:, 3] = barray[:, 2] + 80
    barray[:, 4] = 30
    barray[:, 5] = 1

    bursts = bs.BurstsGap(barray)
    assert bursts.num_bursts == 100
    assert (bursts.counts == 6).all()
    assert (bursts.width == 50).all()

    # Test indexing and slicing
    assert bursts[0] == bs.BurstsGap(barray[0])
    assert bursts[:1] == bursts[0]
    assert (bursts[:10:2].data == barray[:10:2]).all()


if __name__ == '__main__':
    pytest.main("-x -v fretbursts/tests/test_Bursts.py")
