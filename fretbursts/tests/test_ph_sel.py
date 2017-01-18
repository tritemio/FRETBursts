#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
Module containing automated unit tests for FRETBursts.

Running the tests requires `py.test`.
"""

from __future__ import division

from itertools import product
import pytest
from fretbursts.ph_sel import Ph_sel


def test_ph_sel():
    # Assign one excitation
    for i, k in enumerate(('Dex', 'Aex')):
        for v in ('Dem', 'Aem', 'DAem'):
            p = Ph_sel(**{k: v})
            assert p[i] == v
            if v == 'DAem':
                assert str(p) == k      # [Dex|Aex] representations
            else:
                assert str(p) == k + v  # [Dex|Aex][Dem|Aem] representations

    # Assign both excitations
    for dv, av in product(('Dem', 'Aem', 'DAem', None), repeat=2):
        if dv is None and av is None:
            with pytest.raises(ValueError):
                Ph_sel(Dex=dv, Aex=av)
        else:
            p = Ph_sel(Dex=dv, Aex=av)
            assert p.Dex == dv
            assert p.Aex == av
            if dv == av != 'DAem':
                p_str = dv if dv != 'DAem' else 'all'
                assert str(p) == p_str  # [Dem|Aem|all] representations

    # Test some corner cases
    assert Ph_sel('all') == Ph_sel(Dex='DAem', Aex='DAem')
    assert str(Ph_sel(Dex='DAem', Aex='Aem')) == 'DexDAem_AexAem'

    # Test str <-> Ph_sel mapping
    m = Ph_sel._get_str_mapping()
    assert len(set(m.values())) == len(set(m.keys()))
    m_inv = {v: k for k, v in m.items()}
    assert m_inv == Ph_sel._get_str_mapping(invert=True)

    # Test Ph_sel.from_str()
    str_reprs = Ph_sel._get_str_mapping().values()
    for p_str in str_reprs:
        assert Ph_sel.from_str(p_str) == m_inv[p_str]
