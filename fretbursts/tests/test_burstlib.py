#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
Module containing automated unit tests for FRETBursts.

Running the tests requires `py.test`.
"""

import pytest
import numpy as np

from fretbursts import loader
import fretbursts.background as bg
import fretbursts.burstlib as bl
import fretbursts.burstlib_ext as bext
from fretbursts.ph_sel import Ph_sel

# data subdir in the notebook folder
DATASETS_DIR = u'notebooks/data/'


def _alex_process(d):
        loader.alex_apply_period(d)
        d.calc_bg(bg.exp_fit, time_s=30, tail_min_us=300)
        d.burst_search(L=10, m=10, F=7)

def load_dataset_1ch(process=True):
    fn = "0023uLRpitc_NTP_20dT_0.5GndCl.hdf5"
    fname = DATASETS_DIR + fn
    d = loader.hdf5(fname=fname)
    d.add(det_donor_accept=(0, 1), alex_period=4000,
          D_ON=(2850, 580), A_ON=(900, 2580))
    if process:
        _alex_process(d)
    return d

def load_dataset_8ch():
    fn = "12d_New_30p_320mW_steer_3.hdf5"
    fname = DATASETS_DIR + fn
    d = loader.hdf5(fname=fname)
    d.calc_bg(bg.exp_fit, time_s=30, tail_min_us=300)
    d.burst_search(L=10, m=10, F=7)
    return d

@pytest.fixture(scope="module", params=[
                                    load_dataset_1ch,
                                    load_dataset_8ch,
                                    ])
def data(request):
    load_func = request.param
    d = load_func()
    return d


@pytest.fixture(scope="module")
def data_8ch(request):
    d = load_dataset_8ch()
    return d

@pytest.fixture(scope="module")
def data_1ch(request):
    d = load_dataset_1ch()
    return d


##
# List comparison functions
#

def list_equal(list1, list2):
    """Test numerical equality of all the elements in the two lists.
    """
    return np.all([val1 == val2 for val1, val2 in zip(list1, list2)])

def list_array_equal(list1, list2):
    """Test numerical equality between two lists of arrays.
    """
    return np.all([np.all(arr1 == arr2) for arr1, arr2 in zip(list1, list2)])

def list_array_allclose(list1, list2):
    """Test float closeness (np.allclose) between two lists of arrays.
    """
    return np.all([np.allclose(arr1, arr2) for arr1, arr2 in zip(list1, list2)])

##
# Test functions
#

def test_time_min_max():
    """Test time_min and time_max for ALEX data."""
    d = load_dataset_1ch(process=False)
    assert d.time_max == d.ph_times_t.max()*d.clk_p
    assert d.time_min == d.ph_times_t.min()*d.clk_p
    _alex_process(d)
    assert d.time_max == d.ph_times_m[0][-1]*d.clk_p
    assert d.time_min == d.ph_times_m[0][0]*d.clk_p
    d.delete('ph_times_m')
    assert d.time_max == bl.b_end(d.mburst[0])[-1]*d.clk_p
    assert d.time_min == bl.b_start(d.mburst[0])[0]*d.clk_p

def test_time_min_max_multispot(data_8ch):
    """Test time_min and time_max for multi-spot data."""
    d = data_8ch
    assert d.time_max == max(t[-1] for t in d.ph_times_m)*d.clk_p
    assert d.time_min == min(t[0] for t in d.ph_times_m)*d.clk_p

def test_bg_calc(data):
    """Smoke test bg_calc() and test deletion of bg fields.
    """
    data.calc_bg(bg.exp_fit, time_s=30, tail_min_us=300)
    data.calc_bg(bg.exp_fit, time_s=30, tail_min_us='auto', F_bg=1.7)
    assert 'bg_auto_th_us0' in data
    assert 'bg_auto_F_bg' in data
    assert 'bg_th_us_user' not in data
    data.calc_bg(bg.exp_fit, time_s=30, tail_min_us=300)
    assert 'bg_auto_th_us0' not in data
    assert 'bg_auto_F_bg' not in data
    assert 'bg_th_us_user' in data
    data.calc_bg(bg.exp_fit, time_s=30, tail_min_us='auto', F_bg=1.7)

def test_bg_from(data):
    """Test the method .bg_from() for all the ph_sel combinations.
    """
    d = data

    bg = d.bg_from(ph_sel=Ph_sel('all'))
    assert list_array_equal(bg, d.bg)
    bg = d.bg_from(ph_sel=Ph_sel(Dex='Dem'))
    assert list_array_equal(bg, d.bg_dd)
    bg = d.bg_from(ph_sel=Ph_sel(Dex='Aem'))
    assert list_array_equal(bg, d.bg_ad)

    if not d.ALEX:
        bg = d.bg_from(ph_sel=Ph_sel(Dex='DAem'))
        assert list_array_equal(bg, d.bg)
    else:
        bg = d.bg_from(ph_sel=Ph_sel(Aex='Dem'))
        assert list_array_equal(bg, d.bg_da)
        bg = d.bg_from(ph_sel=Ph_sel(Aex='Aem'))
        assert list_array_equal(bg, d.bg_aa)

        bg = d.bg_from(ph_sel=Ph_sel(Dex='DAem'))
        bg_c = [bg1 + bg2 for bg1, bg2 in zip(d.bg_dd, d.bg_ad)]
        assert list_array_equal(bg, bg_c)

        bg = d.bg_from(ph_sel=Ph_sel(Aex='DAem'))
        bg_c = [bg1 + bg2 for bg1, bg2 in zip(d.bg_da, d.bg_aa)]
        assert list_array_equal(bg, bg_c)

        bg = d.bg_from(ph_sel=Ph_sel(Dex='Dem', Aex='Dem'))
        bg_c = [bg1 + bg2 for bg1, bg2 in zip(d.bg_dd, d.bg_da)]
        assert list_array_equal(bg, bg_c)

        bg = d.bg_from(ph_sel=Ph_sel(Dex='Aem', Aex='Aem'))
        bg_c = [bg1 + bg2 for bg1, bg2 in zip(d.bg_ad, d.bg_aa)]
        assert list_array_equal(bg, bg_c)

        bg = d.bg_from(ph_sel=Ph_sel(Dex='DAem', Aex='Aem'))
        bg_c = [bg1 + bg2 + bg3 for bg1, bg2, bg3 in
                        zip(d.bg_dd, d.bg_ad, d.bg_aa)]
        assert list_array_equal(bg, bg_c)


def test_iter_ph_times(data):
    """Test method .iter_ph_times() for all the ph_sel combinations.
    """
    # TODO add all the ph_sel combinations like in test_bg_from()
    d = data

    assert list_array_equal(d.ph_times_m, d.iter_ph_times())

    for ich, ph in enumerate(d.iter_ph_times(Ph_sel(Dex='Dem'))):
        if d.ALEX:
            assert (ph == d.ph_times_m[ich][d.D_em[ich]*d.D_ex[ich]]).all()
        else:
            assert (ph == d.ph_times_m[ich][-d.A_em[ich]]).all()

    for ich, ph in enumerate(d.iter_ph_times(Ph_sel(Dex='Aem'))):
        if d.ALEX:
            assert (ph == d.ph_times_m[ich][d.A_em[ich]*d.D_ex[ich]]).all()
        else:
            assert (ph == d.ph_times_m[ich][d.A_em[ich]]).all()

    if d.ALEX:
        for ich, ph in enumerate(d.iter_ph_times(Ph_sel(Aex='Dem'))):
            assert (ph == d.ph_times_m[ich][d.D_em[ich]*d.A_ex[ich]]).all()
        for ich, ph in enumerate(d.iter_ph_times(Ph_sel(Aex='Aem'))):
            assert (ph == d.ph_times_m[ich][d.A_em[ich]*d.A_ex[ich]]).all()

        for ich, ph in enumerate(d.iter_ph_times(Ph_sel(Dex='DAem'))):
            assert (ph == d.ph_times_m[ich][d.D_ex[ich]]).all()
        for ich, ph in enumerate(d.iter_ph_times(Ph_sel(Aex='DAem'))):
            assert (ph == d.ph_times_m[ich][d.A_ex[ich]]).all()

        for ich, ph in enumerate(d.iter_ph_times(Ph_sel(Dex='Dem', Aex='Dem'))):
            assert (ph == d.ph_times_m[ich][d.D_em[ich]]).all()
        for ich, ph in enumerate(d.iter_ph_times(Ph_sel(Dex='Aem', Aex='Aem'))):
            assert (ph == d.ph_times_m[ich][d.A_em[ich]]).all()

        for ich, ph in enumerate(d.iter_ph_times(
                                    Ph_sel(Dex='DAem', Aex='Aem'))):
            mask = d.D_ex[ich] + d.A_em[ich]*d.A_ex[ich]
            assert (ph == d.ph_times_m[ich][mask]).all()
    else:
        assert list_array_equal(d.iter_ph_times(),
                                d.iter_ph_times(Ph_sel(Dex='DAem')))

def test_get_ph_times_period(data):
    for ich in range(data.nch):
        data.get_ph_times_period(0, ich=ich)
        data.get_ph_times_period(0, ich=ich, ph_sel=Ph_sel(Dex='Dem'))

def test_iter_ph_times_period(data):
    d = data
    for ich in range(data.nch):
        for period, ph_period in enumerate(d.iter_ph_times_period(ich=ich)):
            istart, iend = d.Lim[ich][period]
            assert (ph_period == d.ph_times_m[ich][istart : iend + 1]).all()

        ph_sel = Ph_sel(Dex='Dem')
        mask = d.get_ph_mask(ich=ich, ph_sel=ph_sel)
        for period, ph_period in enumerate(d.iter_ph_times_period(ich=ich,
                                                    ph_sel=ph_sel)):
            istart, iend = d.Lim[ich][period]
            ph_period_test = d.ph_times_m[ich][istart : iend + 1]
            ph_period_test = ph_period_test[mask[istart : iend + 1]]
            assert (ph_period == ph_period_test).all()

def test_burst_search(data):
    """Smoke test and bg_bs check."""
    data.burst_search(L=10, m=10, F=7, ph_sel=Ph_sel(Dex='Dem'))
    assert list_equal(data.bg_bs, data.bg_dd)
    data.burst_search(L=10, m=10, F=7, ph_sel=Ph_sel(Dex='Aem'))
    assert list_equal(data.bg_bs, data.bg_ad)

    if data.ALEX:
        data.burst_search(L=10, m=10, F=7,
                            ph_sel=Ph_sel(Dex='Aem', Aex='Aem'))
        bg_Aem = [b1 + b2 for b1, b2 in zip(data.bg_ad, data.bg_aa)]
        assert list_equal(data.bg_bs, bg_Aem)

    data.burst_search(L=10, m=10, F=7)

def test_burst_sizes(data):
    """Test for .burst_sizes_ich() and burst_sizes()"""
    # Smoke test
    plain_sizes = data.burst_sizes()
    assert len(plain_sizes) == data.nch
    # Test gamma and gamma1 arguments
    bs1 = data.burst_sizes_ich(gamma=0.5)
    bs2 = data.burst_sizes_ich(gamma1=0.5)
    assert np.allclose(bs1, bs2/0.5)
    # Test add_naa
    if data.ALEX:
        bs_no_naa = data.burst_sizes_ich(add_naa=False)
        bs_naa = data.burst_sizes_ich(add_naa=True)
        assert np.allclose(bs_no_naa + data.naa_, bs_naa)

def test_leakage(data):
    """
    Test setting leakage before and after burst search
    """
    # burst search, then set leakage
    data.burst_search()
    data.leakage = 0.04
    na1 = list(data.na)
    # set leakage, then burst search
    data.burst_search()
    na2 = list(data.na)
    assert list_array_equal(na1, na2)

def test_gamma(data):
    """
    Test setting gamma before and after burst search
    """
    # burst search, then set gamma
    data.burst_search()
    data.gamma = 0.5
    E1 = list(data.E)
    # set leakage, then burst search
    data.burst_search()
    E2 = list(data.E)
    assert list_array_equal(E1, E2)

def test_dir_ex(data_1ch):
    """
    Test setting dir_ex before and after burst search
    """
    data = data_1ch
    # burst search, then set dir_ex
    data.burst_search()
    data.dir_ex = 0.05
    na1 = list(data.na)
    # set leakage, then burst search
    data.burst_search()
    na2 = list(data.na)
    assert list_array_equal(na1, na2)

def test_b_functions(data):
    itstart, iwidth, inum_ph, iistart, iiend, itend = 0, 1, 2, 3, 4, 5
    d = data
    for mb in d.mburst:
        assert (bl.b_start(mb) == mb[:, itstart]).all()
        assert (bl.b_end(mb) == mb[:, itend]).all()
        assert (bl.b_width(mb) == mb[:, iwidth]).all()
        assert (bl.b_istart(mb) == mb[:, iistart]).all()
        assert (bl.b_iend(mb) == mb[:, iiend]).all()
        assert (bl.b_size(mb) == mb[:, inum_ph]).all()

        rate = 1.*mb[:, inum_ph]/mb[:, iwidth]
        assert (bl.b_ph_rate(mb) == rate).all()

        separation = mb[1:, itstart] - mb[:-1, itend]
        assert (bl.b_separation(mb) == separation).all()

        assert (bl.b_end(mb) > bl.b_start(mb)).all()

def test_b_end_b_iend(data):
    """Test coherence between b_end() and b_iend()"""
    d = data
    for ph, mb in zip(d.ph_times_m, d.mburst):
        assert (ph[bl.b_iend(mb)] == bl.b_end(mb)).all()

def test_monotonic_burst_start(data):
    """Test for monotonic burst_start."""
    d = data
    for i in xrange(d.nch):
        assert (np.diff(bl.b_start(d.mburst[i])) > 0).all()

def test_monotonic_burst_end(data):
    """Test for monotonic burst_end."""
    d = data
    for mb in d.mburst:
        assert (np.diff(bl.b_end(mb)) > 0).all()

def test_burst_start_end_size(data):
    """Test consistency between burst istart, iend and size"""
    d = data
    for mb in d.mburst:
        size = mb[:, bl.iiend] - mb[:, bl.iistart] + 1
        assert (size == mb[:, bl.inum_ph]).all()
        size2 = bl.b_iend(mb) - bl.b_istart(mb) + 1
        assert (size2 == bl.b_size(mb)).all()

def test_burst_ph_data_functions(data):
    """Tests the functions that operate on per-burst "ph-data".
    """
    d = data
    for bursts, ph, mask in zip(d.mburst, d.iter_ph_times(),
                                d.iter_ph_masks(Ph_sel(Dex='Dem'))):
        bstart = bl.b_start(bursts)
        bend = bl.b_end(bursts)

        for i, (start, stop) in enumerate(bl.iter_bursts_start_stop(bursts)):
            assert ph[start] == bstart[i]
            assert ph[stop-1] == bend[i]

        for i, burst_ph in enumerate(bl.iter_bursts_ph(ph, bursts)):
            assert burst_ph[0] == bstart[i]
            assert burst_ph[-1] == bend[i]

        for i, burst_ph in enumerate(bl.iter_bursts_ph(ph, bursts, mask=mask)):
            if burst_ph.size > 0:
                assert burst_ph[0] >= bstart[i]
                assert burst_ph[-1] <= bend[i]

        stats = bl.burst_ph_stats(ph, bursts, mask=mask)
        assert (stats[~np.isnan(stats)] >= bstart[~np.isnan(stats)]).all()
        assert (stats[~np.isnan(stats)] <= bend[~np.isnan(stats)]).all()

        bistart = bl.b_istart(bursts)
        biend = bl.b_iend(bursts)
        bursts_mask = bl.ph_in_bursts_mask(ph.size, bursts)
        for i, (start, stop) in enumerate(bl.iter_bursts_start_stop(bursts)):
            assert bursts_mask[start:stop].all()
            if start > 0:
                if i > 0 and biend[i-1] < bistart[i] - 1:
                    assert not bursts_mask[start - 1]
            if stop < ph.size:
                if i < bistart.size-1 and bistart[i+1] > biend[i] + 1:
                    assert not bursts_mask[stop]

def test_ph_in_bursts_ich(data):
    """Tests the ph_in_bursts_ich method.
    """
    d = data
    for ich in range(d.nch):
        ph_in_bursts = d.ph_in_bursts_ich(ich)
        ph_in_bursts_dd = d.ph_in_bursts_ich(ich, ph_sel=Ph_sel(Dex='Dem'))
        assert ph_in_bursts_dd.size < ph_in_bursts.size

def test_burst_fuse(data):
    """Test 2 independent implementations of fuse_bursts for consistency.
    """
    d = data
    for mb in d.mburst:
        new_mbursti = bl.fuse_bursts_iter(mb, ms=1)
        new_mburstd = bl.fuse_bursts_direct(mb, ms=1)
        assert (new_mbursti == new_mburstd).all()

def test_burst_fuse_0ms(data):
    """Test that after fusing with ms=0 the sum of bursts sizes is that same
    as the number of ph in bursts (via burst selection).
    """
    d = data

    df = d.fuse_bursts(ms=0)
    for ich, bursts in enumerate(df.mburst):
        mask = bl.ph_in_bursts_mask(df.ph_data_sizes[ich], bursts)
        assert mask.sum() == bl.b_size(bursts).sum()

def test_burst_fuse_separation(data):
    """Test that after fusing bursts the minimum separation is equal
    to the threshold usied during fusing.
    """
    d = data
    fuse_ms = 2
    df = d.fuse_bursts(ms=fuse_ms)
    for mb in df.mburst:
        separation = bl.b_separation(mb)*df.clk_p
        assert separation.min() >= fuse_ms*1e-3

def test_burst_sizes(data):
    """Test that Data.burst_sizes_ich() returns nd + na when gamma = 1.
    """
    d = data
    for ich, (nd, na) in enumerate(zip(d.nd, d.na)):
        burst_size = d.burst_sizes_ich(ich)
        assert (burst_size == nd + na).all()

def test_calc_sbr(data):
    """Smoke test Data.calc_sbr()"""
    data.calc_sbr()

def test_calc_max_rate(data):
    """Smoke test Data.calc_max-rate()"""
    data.calc_max_rate(m=10)

def test_burst_data(data):
    """Smoke test Data.calc_max-rate()"""
    bext.burst_data(data, include_bg=True, include_ph_index=True)

def test_expand(data):
    """Test method `expand()` for `Data()`."""
    d = data
    for ich, mb in enumerate(d.mburst):
        if mb.size == 0: continue  # if no bursts skip this ch
        nd, na, bg_d, bg_a, width = d.expand(ich, width=True)
        width2 = bl.b_width(mb)*d.clk_p
        period = d.bp[ich]
        bg_d2 = d.bg_dd[ich][period] * width2
        bg_a2 = d.bg_ad[ich][period] * width2
        assert (width == width2).all()
        assert (nd == d.nd[ich]).all() and (na == d.na[ich]).all()
        assert (bg_d == bg_d2).all() and (bg_a == bg_a2).all()


def test_burst_corrections(data):
    """Test background and bleed-through corrections."""
    d = data
    d.calc_ph_num(alex_all=True)
    d.corrections()
    leakage = d.get_leakage_array()

    for ich, mb in enumerate(d.mburst):
        if mb.size == 0: continue  # if no bursts skip this ch
        nd, na, bg_d, bg_a, width = d.expand(ich, width=True)
        burst_size_raw = bl.b_size(mb)

        lk = leakage[ich]
        if d.ALEX:
            nda, naa = d.nda[ich], d.naa[ich]
            period = d.bp[ich]
            bg_da = d.bg_da[ich][period]*width
            bg_aa = d.bg_aa[ich][period]*width
            burst_size_raw2 = nd + na + bg_d + bg_a + lk*nd + nda + naa + \
                              bg_da + bg_aa
            assert np.allclose(burst_size_raw, burst_size_raw2)
        else:
            burst_size_raw2 = nd + na + bg_d + bg_a + lk*nd
            assert np.allclose(burst_size_raw, burst_size_raw2)

def test_burst_search_consistency(data):
    """Test consistency of burst data array
    """
    d = data
    for mb, ph in zip(d.mburst, d.iter_ph_times()):
        tot_size = bl.b_size(mb)
        istart, istop = bl.b_istart(mb), bl.b_iend(mb)
        assert np.all(tot_size == istop - istart + 1)
        start, stop, width = bl.b_start(mb), bl.b_end(mb), bl.b_width(mb)
        assert np.all(width == stop - start)
    df = d.fuse_bursts(ms=0)
    for mb, ph in zip(df.mburst, df.iter_ph_times()):
        tot_size = bl.b_size(mb)
        istart, istop = bl.b_istart(mb), bl.b_iend(mb)
        assert np.all(tot_size == istop - istart + 1)
        start, stop, width = bl.b_start(mb), bl.b_end(mb), bl.b_width(mb)
        assert np.all(width == stop - start)
    df = d.fuse_bursts(ms=1)
    for mb, ph in zip(df.mburst, df.iter_ph_times()):
        tot_size = bl.b_size(mb)
        istart, istop = bl.b_istart(mb), bl.b_iend(mb)
        assert np.all(tot_size <= istop - istart + 1)
        start, stop, width = bl.b_start(mb), bl.b_end(mb), bl.b_width(mb)
        assert np.all(width <= stop - start)

def test_burst_size_da(data):
    """Test that nd + na with no corrections is equal to b_size(mburst).
    """
    d = data
    d.calc_ph_num(alex_all=True)
    if d.ALEX:
        for mb, nd, na, naa, nda in zip(d.mburst, d.nd, d.na, d.naa, d.nda):
            tot_size = bl.b_size(mb)
            tot_size2 = nd + na + naa + nda
            assert np.allclose(tot_size, tot_size2)
    else:
        for mb, nd, na in zip(d.mburst, d.nd, d.na):
            tot_size = bl.b_size(mb)
            assert (tot_size == nd + na).all()

def test_collapse(data_8ch):
    """Test the .collapse() method that joins the ch.
    """
    d = data_8ch
    dc1 = d.collapse()
    dc2 = d.collapse(update_gamma=False)

    for name in d.burst_fields:
        if name in d:
            assert np.allclose(dc1[name][0], dc2[name][0])

if __name__ == '__main__':
    pytest.main("-x -v fretbursts/tests/test_burstlib.py")
