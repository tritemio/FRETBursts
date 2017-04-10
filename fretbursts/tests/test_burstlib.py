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
from builtins import range, zip

from collections import namedtuple
import pytest
import numpy as np

try:
    import matplotlib
except ImportError:
    has_matplotlib = False  # OK to run tests without matplotlib
else:
    has_matplotlib = True
    matplotlib.use('Agg')  # but if matplotlib is installed, use Agg

try:
    import numba
except ImportError:
    has_numba = False
else:
    has_numba = True


import fretbursts.background as bg
import fretbursts.burstlib as bl
import fretbursts.burstlib_ext as bext
from fretbursts import loader
from fretbursts import select_bursts
from fretbursts.ph_sel import Ph_sel
from fretbursts.phtools import phrates
if has_matplotlib:
    import fretbursts.burst_plot as bplt


# data subdir in the notebook folder
DATASETS_DIR = u'notebooks/data/'


def _alex_process(d):
    loader.alex_apply_period(d)
    d.calc_bg(bg.exp_fit, time_s=30, tail_min_us=300)
    d.burst_search(L=10, m=10, F=7)

def load_dataset_1ch(process=True):
    fn = "0023uLRpitc_NTP_20dT_0.5GndCl.hdf5"
    fname = DATASETS_DIR + fn
    d = loader.photon_hdf5(fname)
    if process:
        _alex_process(d)
    return d

def load_dataset_8ch():
    fn = "12d_New_30p_320mW_steer_3.hdf5"
    fname = DATASETS_DIR + fn
    d = loader.photon_hdf5(fname)
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


def test_bg_compatlayer_for_obsolete_attrs():
    d = load_dataset_1ch(process=False)
    attrs = ('bg_dd', 'bg_ad', 'bg_da', 'bg_aa',
             'rate_m', 'rate_dd', 'rate_ad', 'rate_da', 'rate_aa')
    for attr in attrs:
        with pytest.raises(RuntimeError):
            getattr(d, attr)
    _alex_process(d)
    for attr in attrs:
        assert isinstance(getattr(d, attr), list)


def test_ph_times_compact(data_1ch):
    """Test calculation of ph_times_compact."""
    def isinteger(x):
        return np.equal(np.mod(x, 1), 0)
    ich = 0
    d = data_1ch

    ph_d = d.get_ph_times(ph_sel=Ph_sel(Dex='DAem'))
    ph_a = d.get_ph_times(ph_sel=Ph_sel(Aex='DAem'))
    ph_dc = d.get_ph_times(ph_sel=Ph_sel(Dex='DAem'), compact=True)
    ph_ac = d.get_ph_times(ph_sel=Ph_sel(Aex='DAem'), compact=True)
    # Test that the difference of ph and ph_compact is multiple of
    # the complementary excitation period duration
    Dex_void = bl._excitation_width(d._D_ON_multich[ich], d.alex_period)
    Aex_void = bl._excitation_width(d._A_ON_multich[ich], d.alex_period)
    assert isinteger((ph_d - ph_dc) / Dex_void).all()
    assert isinteger((ph_a - ph_ac) / Aex_void).all()
    # Test that alternation histogram does not have "gaps" for ph_compact
    bins = np.linspace(0, d.alex_period, num=101)
    hist_dc, _ = np.histogram(ph_dc % d.alex_period, bins=bins)
    hist_ac, _ = np.histogram(ph_ac % d.alex_period, bins=bins)
    assert (hist_dc > 0).all()
    assert (hist_ac > 0).all()


def test_time_min_max():
    """Test time_min and time_max for ALEX data."""
    d = load_dataset_1ch(process=False)
    ich = 0
    assert d.time_max == d.ph_times_t[ich].max() * d.clk_p
    assert d.time_min == d.ph_times_t[ich].min() * d.clk_p
    del d._time_max, d._time_min
    _alex_process(d)
    assert d.time_max == d.ph_times_m[ich][-1] * d.clk_p
    assert d.time_min == d.ph_times_m[ich][0] * d.clk_p
    d.delete('ph_times_m')
    del d._time_max, d._time_min
    assert d.time_max == d.mburst[0].stop[-1] * d.clk_p
    assert d.time_min == d.mburst[0].start[0] * d.clk_p

def test_time_min_max_multispot(data_8ch):
    """Test time_min and time_max for multi-spot data."""
    d = data_8ch
    assert d.time_max == max(t[-1] for t in d.ph_times_m) * d.clk_p
    assert d.time_min == min(t[0] for t in d.ph_times_m) * d.clk_p


def test_bg_calc(data):
    """Smoke test bg_calc() and test deletion of bg fields.
    """
    data.calc_bg(bg.exp_fit, time_s=30, tail_min_us=300)
    assert 'bg_auto_th_us0' not in data
    assert 'bg_auto_F_bg' not in data
    assert 'bg_th_us_user' in data

    data.calc_bg(bg.exp_fit, time_s=30, tail_min_us='auto', F_bg=1.7)
    assert 'bg_auto_th_us0' in data
    assert 'bg_auto_F_bg' in data
    assert 'bg_th_us_user' not in data

    data.calc_bg(bg.exp_fit, time_s=30, tail_min_us='auto', F_bg=1.7,
                 fit_allph=False)
    streams = [s for s in data.ph_streams if s != Ph_sel('all')]
    bg_t = [np.sum(data.bg[s][ich] for s in streams) for ich in range(data.nch)]
    assert list_array_equal(data.bg[Ph_sel('all')], bg_t)


def test_ph_streams(data):
    sel = [Ph_sel('all'), Ph_sel(Dex='Dem'), Ph_sel(Dex='Aem')]
    if data.ALEX:
        sel.extend([Ph_sel(Aex='Aem'), Ph_sel(Aex='Dem')])
    for s in sel:
        assert s in data.ph_streams


def test_bg_from(data):
    """Test the method .bg_from() for all the ph_sel combinations.
    """
    d = data
    for sel in d.ph_streams:
        bg = d.bg_from(ph_sel=sel)
        assert list_array_equal(bg, d.bg[sel])

    if not data.ALEX:
        assert list_array_equal(d.bg_from(Ph_sel('all')),
                                d.bg_from(Ph_sel(Dex='DAem')))
        return

    bg_dd = d.bg_from(ph_sel=Ph_sel(Dex='Dem'))
    bg_ad = d.bg_from(ph_sel=Ph_sel(Dex='Aem'))

    bg = d.bg_from(ph_sel=Ph_sel(Dex='DAem'))
    assert list_array_equal(bg, [b1 + b2 for b1, b2 in zip(bg_dd, bg_ad)])

    bg_aa = d.bg_from(ph_sel=Ph_sel(Aex='Aem'))
    bg_da = d.bg_from(ph_sel=Ph_sel(Aex='Dem'))

    bg = d.bg_from(ph_sel=Ph_sel(Aex='DAem'))
    assert list_array_equal(bg, [b1 + b2 for b1, b2 in zip(bg_aa, bg_da)])

    bg = d.bg_from(ph_sel=Ph_sel(Dex='Dem', Aex='Dem'))
    assert list_array_equal(bg, [b1 + b2 for b1, b2 in zip(bg_dd, bg_da)])

    bg = d.bg_from(ph_sel=Ph_sel(Dex='Aem', Aex='Aem'))
    assert list_array_equal(bg, [b1 + b2 for b1, b2 in zip(bg_ad, bg_aa)])

    bg = d.bg_from(ph_sel=Ph_sel(Dex='DAem'))
    assert list_array_equal(bg, [b1 + b2 for b1, b2 in zip(bg_dd, bg_ad)])

    bg = d.bg_from(ph_sel=Ph_sel(Dex='DAem', Aex='Aem'))
    bg2 = [b1 + b2 + b3 for b1, b2, b3 in zip(bg_dd, bg_ad, bg_aa)]
    assert list_array_equal(bg, bg2)


def test_iter_ph_times(data):
    """Test method .iter_ph_times() for all the ph_sel combinations.
    """
    # TODO add all the ph_sel combinations like in test_bg_from()
    d = data

    assert list_array_equal(d.ph_times_m, d.iter_ph_times())

    for ich, ph in enumerate(d.iter_ph_times(Ph_sel(Dex='Dem'))):
        if d.ALEX:
            assert (ph == d.ph_times_m[ich][d.D_em[ich] * d.D_ex[ich]]).all()
        else:
            assert (ph == d.ph_times_m[ich][-d.A_em[ich]]).all()

    for ich, ph in enumerate(d.iter_ph_times(Ph_sel(Dex='Aem'))):
        if d.ALEX:
            assert (ph == d.ph_times_m[ich][d.A_em[ich] * d.D_ex[ich]]).all()
        else:
            assert (ph == d.ph_times_m[ich][d.A_em[ich]]).all()

    if d.ALEX:
        for ich, ph in enumerate(d.iter_ph_times(Ph_sel(Aex='Dem'))):
            assert (ph == d.ph_times_m[ich][d.D_em[ich] * d.A_ex[ich]]).all()
        for ich, ph in enumerate(d.iter_ph_times(Ph_sel(Aex='Aem'))):
            assert (ph == d.ph_times_m[ich][d.A_em[ich] * d.A_ex[ich]]).all()

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
            mask = d.D_ex[ich] + d.A_em[ich] * d.A_ex[ich]
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
        for period, ph_period in enumerate(
                d.iter_ph_times_period(ich=ich, ph_sel=ph_sel)):
            istart, iend = d.Lim[ich][period]
            ph_period_test = d.ph_times_m[ich][istart : iend + 1]
            ph_period_test = ph_period_test[mask[istart : iend + 1]]
            assert (ph_period == ph_period_test).all()

def test_burst_search_py_cy(data):
    """Test python and cython burst search with background-dependent threshold.
    """
    data.burst_search(pure_python=True)
    mburst1 = [b.copy() for b in data.mburst]
    num_bursts1 = data.num_bursts
    data.burst_search(pure_python=False)
    assert np.all(num_bursts1 == data.num_bursts)
    assert mburst1 == data.mburst
    data.burst_search(L=30, pure_python=True)
    mburst1 = [b.copy() for b in data.mburst]
    num_bursts1 = data.num_bursts
    data.burst_search(L=30, pure_python=False)
    assert np.all(num_bursts1 == data.num_bursts)
    assert mburst1 == data.mburst

def test_burst_search_constant_rates(data):
    """Test python and cython burst search with constant threshold."""
    data.burst_search(min_rate_cps=50e3, pure_python=True)
    assert (data.num_bursts > 0).all()
    mburst1 = [b.copy() for b in data.mburst]
    num_bursts1 = data.num_bursts
    data.burst_search(min_rate_cps=50e3, pure_python=False)
    assert (data.num_bursts > 0).all()
    assert np.all(num_bursts1 == data.num_bursts)
    assert mburst1 == data.mburst

def test_burst_search_L(data):
    """Test burst search with different L arguments."""
    data.burst_search(L=10)
    for bursts in data.mburst:
        assert (bursts.counts >= 10).all()
    num_bursts1 = data.num_bursts
    data.burst_search(L=30)
    for bursts in data.mburst:
        assert (bursts.counts >= 30).all()
    assert np.all(num_bursts1 > data.num_bursts)

def test_burst_search_with_no_bursts(data):
    """Smoke test burst search when some periods have no bursts."""
    # F=600 results in periods with no bursts for the us-ALEX measurement
    # and in no bursts at all for the multi-spot measurements
    data.burst_search(m=10, F=600)

if has_matplotlib:
    def test_stale_fitter_after_burst_search(data):
        """Test that E/S_fitter attributes are deleted on burst search."""
        data.burst_search(L=10, m=10, F=7, ph_sel=Ph_sel(Dex='Dem'))
        bplt.dplot(data, bplt.hist_fret)  # create E_fitter attribute
        if data.ALEX:
            bplt.dplot(data, bplt.hist_S)  # create S_fitter attribute

        data.burst_search(L=10, m=10, F=7, ph_sel=Ph_sel(Dex='Aem'))
        assert not hasattr(data, 'E_fitter')
        if data.ALEX:
            assert not hasattr(data, 'S_fitter')

        bplt.dplot(data, bplt.hist_fret)  # create E_fitter attribute
        if data.ALEX:
            bplt.dplot(data, bplt.hist_S)  # create S_fitter attribute

        data.calc_fret()
        assert not hasattr(data, 'E_fitter')
        if data.ALEX:
            assert not hasattr(data, 'S_fitter')

def test_burst_search(data):
    """Smoke test and bg_bs check."""
    streams = [Ph_sel(Dex='Dem'), Ph_sel(Dex='Aem')]
    if data.ALEX:
        streams.extend([Ph_sel(Dex='Aem', Aex='Aem'), Ph_sel(Dex='DAem')])
    for sel in streams:
        data.burst_search(L=10, m=10, F=7, ph_sel=sel)
        assert list_equal(data.bg_bs, data.bg_from(sel))

    if data.ALEX:
        data.burst_search(m=10, F=7, ph_sel=Ph_sel(Dex='DAem'), compact=True)
    data.burst_search(L=10, m=10, F=7)

def test_burst_search_and_gate(data_1ch):
    """Test consistency of burst search and gate."""
    d = data_1ch
    assert d.ALEX
    d_dex = d.copy()
    d_dex.burst_search(ph_sel=Ph_sel(Dex='DAem'))
    d_aex = d.copy()
    d_aex.burst_search(ph_sel=Ph_sel(Aex='Aem'))
    d_and = bext.burst_search_and_gate(d)
    for bursts_dex, bursts_aex, bursts_and, ph in zip(
            d_dex.mburst, d_aex.mburst, d_and.mburst, d.iter_ph_times()):
        ph_b_mask_dex = bl.ph_in_bursts_mask(ph.size, bursts_dex)
        ph_b_mask_aex = bl.ph_in_bursts_mask(ph.size, bursts_aex)
        ph_b_mask_and = bl.ph_in_bursts_mask(ph.size, bursts_and)
        assert (ph_b_mask_and == ph_b_mask_dex * ph_b_mask_aex).all()

def test_mch_count_ph_num_py_c(data):
    na_py = bl.bslib.mch_count_ph_in_bursts_py(data.mburst, data.A_em)
    na_c = bl.bslib.mch_count_ph_in_bursts_c(data.mburst, data.A_em)
    assert list_array_equal(na_py, na_c)
    assert na_py[0].dtype == np.float64

def test_burst_sizes(data):
    """Test for .burst_sizes_ich() and burst_sizes()"""
    # Smoke test
    plain_sizes = data.burst_sizes()
    assert len(plain_sizes) == data.nch
    # Test gamma and donor_ref arguments
    bs1 = data.burst_sizes_ich(gamma=0.5, donor_ref=True)
    bs2 = data.burst_sizes_ich(gamma=0.5, donor_ref=False)
    assert np.allclose(bs1, bs2 / 0.5)
    # Test add_naa
    if data.ALEX:
        bs_no_naa = data.burst_sizes_ich(add_naa=False)
        bs_naa = data.burst_sizes_ich(add_naa=True)
        assert np.allclose(bs_no_naa + data.naa_, bs_naa)

        # Test beta and donor_ref arguments with gamma=1
        naa1 = data.get_naa_corrected(beta=0.8, donor_ref=True)
        naa2 = data.get_naa_corrected(beta=0.8, donor_ref=False)
        assert np.allclose(naa1, naa2)

        # Test beta and donor_ref arguments with gamma=0.5
        naa1 = data.get_naa_corrected(gamma=0.5, beta=0.8, donor_ref=True)
        naa2 = data.get_naa_corrected(gamma=0.5, beta=0.8, donor_ref=False)
        assert np.allclose(naa1 * 0.5, naa2)

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
    E0 = list(data.E)
    data.gamma = 0.5
    E1 = list(data.E)
    assert not list_array_equal(E0, E1)
    # burst search after setting gamma
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
    na0 = list(data.na)
    data.dir_ex = 0.05
    na1 = list(data.na)
    assert not list_array_equal(na0, na1)
    # burst search after setting dir_ex
    data.burst_search()
    na2 = list(data.na)
    assert list_array_equal(na1, na2)

def test_beta(data_1ch):
    """
    Test setting beta before and after burst search
    """
    data = data_1ch
    # burst search, then set beta
    data.burst_search()
    S0 = list(data.S)
    data.beta = 0.7
    S1 = list(data.S)
    assert not list_array_equal(S0, S1)
    # burst search after setting beta
    data.burst_search()
    S2 = list(data.S)
    assert list_array_equal(S1, S2)

def test_bursts_interface(data):
    d = data
    for b in d.mburst:
        assert (b.start == b.data[:, b._i_start]).all()
        assert (b.stop == b.data[:, b._i_stop]).all()
        assert (b.istart == b.data[:, b._i_istart]).all()
        assert (b.istop == b.data[:, b._i_istop]).all()

        rate = 1.*b.counts/b.width
        assert (b.ph_rate == rate).all()

        separation = b.start[1:] - b.stop[:-1]
        assert (b.separation == separation).all()

        assert (b.stop > b.start).all()

def test_burst_stop_istop(data):
    """Test coherence between b_end() and b_iend()"""
    d = data
    for ph, bursts in zip(d.ph_times_m, d.mburst):
        assert (ph[bursts.istop] == bursts.stop).all()

def test_monotonic_burst_start(data):
    """Test for monotonic burst start times."""
    d = data
    for i in range(d.nch):
        assert (np.diff(d.mburst[i].start) > 0).all()

def test_monotonic_burst_stop(data):
    """Test for monotonic burst stop times."""
    d = data
    for bursts in d.mburst:
        assert (np.diff(bursts.stop) > 0).all()

def test_burst_istart_iend_size(data):
    """Test consistency between burst istart, istop and counts (i.e. size)"""
    d = data
    for bursts in d.mburst:
        counts = bursts.istop - bursts.istart + 1
        assert (counts == bursts.counts).all()

def test_burst_recompute_times(data):
    """Test Bursts.recompute_times method."""
    d = data
    for times, bursts in zip(d.ph_times_m, d.mburst):
        newbursts = bursts.recompute_times(times)
        assert newbursts == bursts

def test_burst_recompute_index(data):
    """Test Bursts.recompute_index_* methods."""
    d = data
    ph_sel = Ph_sel(Dex='Dem')
    d.burst_search(ph_sel=ph_sel, index_allph=True)
    d_sel = d.copy()
    d_sel.burst_search(ph_sel=ph_sel, index_allph=False)
    for times_sel, mask_sel, bursts_sel, times_allph, bursts_allph in zip(
            d.iter_ph_times(ph_sel=ph_sel),
            d.iter_ph_masks(ph_sel=ph_sel),
            d_sel.mburst,
            d.iter_ph_times(),
            d.mburst):
        assert (times_sel[bursts_sel.istart] == bursts_sel.start).all()
        assert (times_sel[bursts_sel.istop] == bursts_sel.stop).all()

        assert (times_allph[bursts_allph.istart] == bursts_allph.start).all()
        assert (times_allph[bursts_allph.istop] == bursts_allph.stop).all()

        # Test individual methods
        bursts_allph2 = bursts_sel.recompute_index_expand(mask_sel)
        assert  bursts_allph2 == bursts_allph
        assert (times_allph[bursts_allph2.istart] == bursts_allph2.start).all()
        assert (times_allph[bursts_allph2.istop] == bursts_allph2.stop).all()

        bursts_sel2 = bursts_allph.recompute_index_reduce(times_sel)
        assert (times_sel[bursts_sel2.istart] == bursts_sel2.start).all()
        assert (times_sel[bursts_sel2.istop] == bursts_sel2.stop).all()
        assert  bursts_sel2 == bursts_sel

        # Test round-trip
        bursts_allph3 = bursts_sel2.recompute_index_expand(mask_sel)
        assert  bursts_allph3 == bursts_allph2
        assert (times_allph[bursts_allph3.istart] == bursts_allph3.start).all()
        assert (times_allph[bursts_allph3.istop] == bursts_allph3.stop).all()

## This test is only used to develop alternative implementations of
## Bursts.recompute_index_reduce() and is normally disabled as it is very slow.
#def test_burst_recompute_index_reduce(data):
#    """Test different versions of Bursts.recompute_index_reduce methods.
#
#    This test is very slow so it's normally disabled.
#    """
#    d = data
#    ph_sel = Ph_sel(Dex='Aem')
#    d.burst_search(ph_sel=ph_sel)
#    d_sel = d.copy()
#    d_sel.burst_search(ph_sel=ph_sel, index_allph=False)
#    for times_sel, bursts_sel, times_allph, bursts_allph in zip(
#            d.iter_ph_times(ph_sel=ph_sel),
#            d_sel.mburst,
#            d.iter_ph_times(),
#            d.mburst):
#        assert (times_allph[bursts_allph.istart] == bursts_allph.start).all()
#        assert (times_allph[bursts_allph.istop] == bursts_allph.stop).all()
#
#        bursts_sel1 = bursts_allph.recompute_index_reduce(times_sel)
#        bursts_sel2 = bursts_allph.recompute_index_reduce2(times_sel)
#        assert  bursts_sel1 == bursts_sel2
#        assert  bursts_sel == bursts_sel1

def test_phrates_mtuple(data):
    d = data
    m = 10
    max_num_ph = 20001
    for ph in d.iter_ph_times():
        phc = ph[:max_num_ph]
        rates = phrates.mtuple_rates(phc, m)
        delays = phrates.mtuple_delays(phc, m)
        t_rates = 0.5 * (phc[m-1:] + phc[:-m+1])
        assert phrates.mtuple_rates_max(phc, m) == rates.max()
        assert phrates.mtuple_delays_min(phc, m) == delays.min()
        assert phrates.default_c == 1
        assert (rates == (m - 1 - phrates.default_c) / delays).all()
        assert (phrates.mtuple_rates_t(phc, m) == t_rates).all()

if has_numba:
    def test_phrates_kde(data):
        d = data
        tau = 5000  # 5000 * 12.5ns = 6.25 us
        for ph in d.iter_ph_times():
            # Test consistency of kde_laplace_nph and (kde_laplace, kde_rect)
            rates = phrates.kde_laplace(ph, tau)
            nrect = phrates.kde_rect(ph, tau*10)
            ratesl, nph = phrates.nb.kde_laplace_nph(ph, tau)
            assert (rates == ratesl).all()
            assert (nph == nrect).all()

            # Test consistency of kde_laplace and _kde_laplace_self_numba
            ratesl2, nph2 = phrates.nb.kde_laplace_self_numba(ph, tau)
            assert (nph2 == nrect).all()
            assert (ratesl2 == rates).all()

            # Smoke test laplace, gaussian, rect with time_axis
            ratesl = phrates.kde_laplace(ph, tau, time_axis=ph+1)
            assert ((ratesl >= 0) * (ratesl < 5e6)).all()
            ratesg = phrates.kde_gaussian(ph, tau, time_axis=ph+1)
            assert ((ratesg >= 0) * (ratesg < 5e6)).all()
            ratesr = phrates.kde_rect(ph, tau, time_axis=ph+1)
            assert ((ratesr >= 0) * (ratesr < 5e6)).all()

    def test_phrates_kde_cy(data):
        d = data
        tau = 5000  # 5000 * 12.5ns = 6.25 us
        for ph in d.iter_ph_times():
            # Test consistency of kde_laplace_nph and (kde_laplace, kde_rect)
            ratesg = phrates.nb.kde_gaussian_numba(ph, tau)
            ratesl = phrates.nb.kde_laplace_numba(ph, tau)
            ratesr = phrates.nb.kde_rect_numba(ph, tau)
            ratesgc = phrates.cy.kde_gaussian_cy(ph, tau)
            rateslc = phrates.cy.kde_laplace_cy(ph, tau)
            ratesrc = phrates.cy.kde_rect_cy(ph, tau)
            assert (ratesg == ratesgc).all()
            assert (ratesl == rateslc).all()
            assert (ratesr == ratesrc).all()

def test_burst_ph_data_functions(data):
    """Tests the functions that iterate or operate on per-burst "ph-data".
    """
    d = data
    for bursts, ph, mask in zip(d.mburst, d.iter_ph_times(),
                                d.iter_ph_masks(Ph_sel(Dex='Dem'))):
        bstart = bursts.start
        bend = bursts.stop

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

        bistart = bursts.istart
        biend = bursts.istop
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
    for bursts in d.mburst:
        new_mbursti = bl.fuse_bursts_iter(bursts, ms=1)
        new_mburstd = bl.fuse_bursts_direct(bursts, ms=1)
        assert new_mbursti == new_mburstd

def test_burst_fuse_0ms(data):
    """Test that after fusing with ms=0 the sum of bursts sizes is that same
    as the number of ph in bursts (via burst selection).
    """
    d = data

    df = d.fuse_bursts(ms=0)
    for ich, bursts in enumerate(df.mburst):
        mask = bl.ph_in_bursts_mask(df.ph_data_sizes[ich], bursts)
        assert mask.sum() == bursts.counts.sum()

def test_burst_fuse_separation(data):
    """Test that after fusing bursts the minimum separation is equal
    to the threshold used during fusing.
    """
    d = data
    fuse_ms = 2
    df = d.fuse_bursts(ms=fuse_ms)
    for bursts in df.mburst:
        separation = bursts.separation*df.clk_p
        assert separation.min() >= fuse_ms*1e-3

def test_calc_sbr(data):
    """Smoke test Data.calc_sbr()"""
    data.calc_sbr()

def test_calc_max_rate(data):
    """Smoke test for Data.calc_max_rate()"""
    data.calc_max_rate(m=10)
    if data.ALEX:
        data.calc_max_rate(m=10, ph_sel=Ph_sel(Dex='DAem'), compact=True)

def test_burst_data(data):
    """Smoke test for bext.burst_data()"""
    bext.burst_data(data, include_bg=True, include_ph_index=True)

def test_print_burst_stats(data):
    """Smoke test for burstlib.print_burst_stats()"""
    bl.print_burst_stats(data)

def test_expand(data):
    """Test method `expand()` for `Data()`."""
    d = data
    for ich, bursts in enumerate(d.mburst):
        if bursts.num_bursts == 0:
            continue  # if no bursts skip this ch
        nd, na, bg_d, bg_a, width = d.expand(ich, width=True)
        width2 = bursts.width * d.clk_p
        period = d.bp[ich]
        bg_d2 = d.bg_from(Ph_sel(Dex='Dem'))[ich][period] * width2
        bg_a2 = d.bg_from(Ph_sel(Dex='Aem'))[ich][period] * width2
        assert (width == width2).all()
        assert (nd == d.nd[ich]).all() and (na == d.na[ich]).all()
        assert (bg_d == bg_d2).all() and (bg_a == bg_a2).all()


def test_burst_corrections(data):
    """Test background and bleed-through corrections."""
    d = data
    d.calc_ph_num(alex_all=True)
    d.corrections()
    leakage = d.get_leakage_array()

    for ich, bursts in enumerate(d.mburst):
        if bursts.num_bursts == 0: continue  # if no bursts skip this ch
        nd, na, bg_d, bg_a, width = d.expand(ich, width=True)
        burst_size_raw = bursts.counts

        lk = leakage[ich]
        if d.ALEX:
            nda, naa = d.nda[ich], d.naa[ich]
            period = d.bp[ich]
            bg_da = d.bg_from(Ph_sel(Aex='Dem'))[ich][period]*width
            bg_aa = d.bg_from(Ph_sel(Aex='Aem'))[ich][period]*width
            burst_size_raw2 = (nd + na + bg_d + bg_a + lk*nd + nda + naa +
                               bg_da + bg_aa)
            assert np.allclose(burst_size_raw, burst_size_raw2)
        else:
            burst_size_raw2 = nd + na + bg_d + bg_a + lk*nd
            assert np.allclose(burst_size_raw, burst_size_raw2)

def test_burst_search_consistency(data):
    """Test consistency of burst data array
    """
    d = data
    for mb, ph in zip(d.mburst, d.iter_ph_times()):
        tot_size = mb.counts
        istart, istop = mb.istart, mb.istop
        assert np.all(tot_size == istop - istart + 1)
        start, stop, width = mb.start, mb.stop, mb.width
        assert np.all(width == stop - start)
    df = d.fuse_bursts(ms=0)
    for mb, ph in zip(df.mburst, df.iter_ph_times()):
        tot_size = mb.counts
        istart, istop = mb.istart, mb.istop
        assert np.all(tot_size == istop - istart + 1)
        start, stop, width = mb.start, mb.stop, mb.width
        assert np.all(width == stop - start)
    df = d.fuse_bursts(ms=1)
    for mb, ph in zip(df.mburst, df.iter_ph_times()):
        tot_size = mb.counts
        istart, istop = mb.istart, mb.istop
        assert np.all(tot_size <= istop - istart + 1)
        start, stop, width = mb.start, mb.stop, mb.width
        assert np.all(width <= stop - start)

def test_E_and_S_with_corrections(data):
    d = data
    gamma = 0.5
    beta = 0.7
    d.gamma = gamma
    d.beta = beta
    for i, (E, nd, na) in enumerate(zip(d.E, d.nd, d.na)):
        assert (E == na / (nd * gamma + na)).all()
        if d.ALEX:
            naa = d.naa[i]
            assert (d.S[i] == (gamma * nd + na) /
                              (gamma * nd + na + naa / beta)).all()


def test_burst_size_da(data):
    """Test that nd + na with no corrections is equal to b_size(mburst).
    """
    d = data
    d.calc_ph_num(alex_all=True)
    if d.ALEX:
        for mb, nd, na, naa, nda in zip(d.mburst, d.nd, d.na, d.naa, d.nda):
            tot_size = mb.counts
            tot_size2 = nd + na + naa + nda
            assert np.allclose(tot_size, tot_size2)
    else:
        for mb, nd, na in zip(d.mburst, d.nd, d.na):
            tot_size = mb.counts
            assert (tot_size == nd + na).all()

def test_burst_selection(data):
    """Smoke test for burst selection methods.
    """
    d = data
    d.select_bursts(select_bursts.size, th1=20, th2=100, add_naa=True)
    d.select_bursts(select_bursts.size, th1=20, th2=100, gamma=0.5)

    M1 = d.select_bursts_mask(select_bursts.consecutive, th1=1e-3, th2=1e4,
                              kind='first')
    M2 = d.select_bursts_mask(select_bursts.consecutive, th1=1e-3, th2=1e4,
                              kind='second')
    Mb = d.select_bursts_mask(select_bursts.consecutive, th1=1e-3, th2=1e4,
                              kind='both')
    Mb2 = [m1 + m2 for m1, m2 in zip(M1, M2)]
    assert list_array_equal(Mb, Mb2)

def test_burst_selection_nocorrections(data):
    """Test burst selection with uncorrected bursts.
    """
    d = data
    d.burst_search(computefret=False)
    d.calc_fret(count_ph=True, corrections=False)
    ds1 = d.select_bursts(select_bursts.size, th1=20, th2=100,
                          computefret=False)
    ds2 = d.select_bursts(select_bursts.size, th1=20, th2=100)
    ds2.calc_ph_num()
    ds2.calc_fret(corrections=False)

    assert list_array_equal(ds1.nd, ds2.nd)
    assert list_array_equal(ds1.na, ds2.na)
    assert list_array_equal(ds1.E, ds2.E)
    if d.ALEX:
        assert list_array_equal(ds1.naa, ds2.naa)
        assert list_array_equal(ds1.E, ds2.E)

def test_burst_selection_ranges(data):
    """Test selection functions having a min-max range.
    """
    d = data
    d.burst_search()
    d.calc_max_rate(m=10, ph_sel=Ph_sel(Dex='DAem'))

    Range = namedtuple('Range', ['min', 'max', 'getter'])

    sel_functions = dict(
        E=Range(0.5, 1, None), nd=Range(30, 40, None), na=Range(30, 40, None),
        time=Range(1, 61, lambda d, ich: d.mburst[ich].start * d.clk_p),
        width=Range(0.5, 1.5, lambda d, ich: d.mburst[ich].width * d.clk_p*1e3),
        peak_phrate=Range(50e3, 150e3, lambda d, ich: d.max_rate[ich]))
    if d.ALEX:
        sel_functions.update(naa=Range(30, 40, None), S=Range(0.3, 0.7, None))

    for func_name, range_ in sel_functions.items():
        func = getattr(select_bursts, func_name)
        getter = range_.getter
        if getter is None:
            getter = lambda d, ich: d[func_name][ich]

        ds = d.select_bursts(func, args=(range_.min, range_.max))
        for ich in range(d.nch):
            selected = getter(ds, ich)
            assert ((selected >= range_.min) * (selected <= range_.max)).all()

def test_join_data(data):
    """Smoke test for bext.join_data() function.
    """
    d = data
    dj = bext.join_data([d, d.copy()])
    assert (dj.num_bursts == 2 * d.num_bursts).all()
    for bursts in dj.mburst:
        assert (np.diff(bursts.start) > 0).all()

def test_collapse(data_8ch):
    """Test the .collapse() method that joins the ch.
    """
    d = data_8ch
    dc1 = d.collapse()
    bursts1 = dc1.mburst[0]
    bursts2 = bl.bslib.Bursts.merge(d.mburst, sort=True)
    assert bursts1 == bursts2
    bursts2 = bl.bslib.Bursts.merge(d.mburst, sort=False)
    indexsort_stop = bursts2.stop.argsort()
    bursts3 = bursts2[indexsort_stop]
    indexsort_start = bursts3.start.argsort()
    bursts4 = bursts3[indexsort_start]
    assert bursts1 == bursts4
    indexsort = np.lexsort((bursts2.stop, bursts2.start))
    for name in d.burst_fields:
        if name not in d or name == 'mburst':
            continue
        newfield = np.hstack(d[name])[indexsort]
        assert np.allclose(dc1[name][0], newfield)

    dc2 = d.collapse(update_gamma=False)
    for name in d.burst_fields:
        if name not in d: continue

        if name == 'mburst':
            assert dc1.mburst[0] == dc2.mburst[0]
        else:
            assert np.allclose(dc1[name][0], dc2[name][0])

if __name__ == '__main__':
    pytest.main("-x -v fretbursts/tests/test_burstlib.py")
