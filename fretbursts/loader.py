#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014-2016 The Regents of the University of California,
#               Antonino Ingargiola <tritemio@gmail.com>
#
"""
The `loader` module contains functions to load each supported data format.
The loader functions load data from a specific format and
return a new :class:`fretbursts.burstlib.Data()` object containing the data.

This module contains the high-level function to load a data-file and
to return a `Data()` object. The low-level functions that perform the binary
loading and preprocessing can be found in the `dataload` folder.
"""

from __future__ import print_function, absolute_import
from builtins import range, zip

import os
import numpy as np
try:
    import cPickle as pickle
except ImportError:
    import pickle
import tables

from phconvert.smreader import load_sm
from .dataload.multi_ch_reader import load_data_ordered16
from .dataload.spcreader import load_spc
from .dataload.manta_reader import (load_manta_timestamps,
                                    load_xavier_manta_data,
                                    get_timestamps_detectors,
                                    # process_timestamps,
                                    process_store,
                                    load_manta_timestamps_pytables)
from .utils.misc import pprint, deprecate
from .burstlib import Data
from . import loader_legacy
import phconvert as phc


def _is_multich(h5data):
    if 'photon_data' in h5data:
        return False
    elif 'photon_data0' in h5data:
        return True
    else:
        msg = 'Cannot find a photon_data group.'
        raise phc.hdf5.Invalid_PhotonHDF5(msg)

def _append_data_ch(d, name, value):
    if name not in d:
        d.add(**{name: [value]})
    else:
        d[name].append(value)

def _load_from_group(d, group, name, dest_name, multich_field=False,
                     ondisk=False, allow_missing=True):
    if allow_missing and name not in group:
        return

    node_value = group._f_get_child(name)
    if not ondisk:
        node_value = node_value.read()
    if multich_field:
        _append_data_ch(d, dest_name, node_value)
    else:
        d.add(**{dest_name: node_value})

def _append_empy_ch(data):
    # Empty channel, fill it with empty arrays
    ph_times = np.array([], dtype='int64')
    _append_data_ch(data, 'ph_times_m', ph_times)

    a_em = np.array([], dtype=bool)
    _append_data_ch(data, 'A_em', a_em)

def _get_measurement_specs(ph_data):
    if 'measurement_specs' not in ph_data:
        # No measurement specs, we will load timestamps and set them all in a
        # conventional photon stream (acceptor emission)
        meas_type = 'smFRET-1color'
        meas_specs = None
    else:
        assert 'measurement_type' in ph_data.measurement_specs
        meas_specs = ph_data.measurement_specs
        meas_type = meas_specs.measurement_type.read().decode()

    if meas_type not in ['smFRET-1color', 'smFRET',
                         'smFRET-usALEX', 'smFRET-nsALEX']:
        raise NotImplementedError('Measurement type "%s" not supported'
                                  ' by FRETBursts.' % meas_type)
    return meas_type, meas_specs

def _load_photon_data_arrays(data, ph_data, meas_type, ondisk=False):
    assert 'timestamps' in ph_data

    # Build mapping to convert Photon-HDF5 to FRETBursts names
    # fields not mapped use the same name on both Photon-HDF5 and FRETBursts
    mapping = {'timestamps': 'ph_times_m',
               'nanotimes': 'nanotimes', 'particles': 'particles'}
    if 'ALEX' in meas_type:
        mapping = {'timestamps': 'ph_times_t', 'detectors': 'det_t',
                   'nanotimes': 'nanotimes_t', 'particles': 'particles_t'}

    # Load all photon-data arrays
    for name in ph_data._v_leaves:
        dest_name = mapping.get(name, name)
        _load_from_group(data, ph_data, name, dest_name=dest_name,
                         multich_field=True, ondisk=ondisk)

    # Timestamps are always present, and their units are always present too
    data.add(clk_p=ph_data.timestamps_specs.timestamps_unit.read())

def _load_nanotimes_specs(data, ph_data):
    nanot_specs = ph_data.nanotimes_specs
    nanotimes_params = {}
    for name in ['tcspc_unit', 'tcspc_num_bins', 'tcspc_range']:
        value = nanot_specs._f_get_child(name).read()
        nanotimes_params.update(**{name: value})
    if 'user' in nanot_specs:
        for name in ['tau_accept_only', 'tau_donor_only',
                     'tau_fret_donor', 'inverse_fret_rate']:
            if name in nanot_specs.user:
                value = nanot_specs.user._f_get_child(name).read()
                nanotimes_params.update(**{name: value})
    _append_data_ch(data, 'nanotimes_params', nanotimes_params)

def _load_alex_periods_donor_acceptor(data, meas_specs):
    # Both us- and ns-ALEX
    try:
        # Try to load alex period definitions
        D_ON = meas_specs.alex_excitation_period1.read()
        A_ON = meas_specs.alex_excitation_period2.read()
    except tables.NoSuchNodeError:
        # But if it fails it's OK, those fields are optional
        print('WARNING: No alternation defintion found.')
    else:
        _append_data_ch(data, 'D_ON', D_ON)
        _append_data_ch(data, 'A_ON', A_ON)

def _compute_acceptor_emission_mask(data, ich, ondisk):
    """For non-ALEX measurements."""
    if data.detectors[ich].dtype.itemsize != 1:
        raise NotImplementedError('Detectors dtype must be 1-byte.')
    donor, accept = data._det_donor_accept_multich[ich]

    # Remove counts not associated with D or A channels
    if not ondisk and len(np.unique(data.detectors[ich])) > 2:
        mask = (data.detectors[ich] == donor) + (data.detectors[ich] == accept)
        data.detectors[ich] = data.detectors[ich][mask]
        data.ph_times_m[ich] = data.ph_times_m[ich][mask]
        if 'nanotimes' in data:
            data.nanotimes[ich] = data.nanotimes[ich][mask]

    # From `detectors` compute boolean mask `A_em`
    if not ondisk and 0 in (accept, donor):
        # In this case we create the boolean mask in-place
        # using the detectors array
        _append_data_ch(data, 'A_em', data.detectors[ich].view(dtype=bool))
        if accept == 0:
            np.logical_not(data.A_em[ich], out=data.A_em[ich])
    else:
        # Create the boolean mask
        _append_data_ch(data, 'A_em', data.detectors[ich][:] == accept)

def _add_usALEX_specs(data, meas_specs):
    try:
        offset = meas_specs.alex_offset.read()
    except tables.NoSuchNodeError:
        print('WARNING: No offset found, assuming offset = 0.')
        offset = 0
    data.add(offset=offset)
    data.add(alex_period=meas_specs.alex_period.read())

def _photon_hdf5_1ch(h5data, data, ondisk=False, nch=1, ich=0):
    data.add(nch=nch)
    ph_data_name = '/photon_data' if nch == 1 else '/photon_data%d' % ich

    # Handle the case of missing channel (e.g. dead pixel)
    if ph_data_name not in h5data:
        _append_empy_ch(data)
        return

    # Load photon_data group and measurement_specs (if present)
    ph_data = h5data._f_get_child(ph_data_name)
    meas_type, meas_specs = _get_measurement_specs(ph_data)

    # Load photon_data arrays
    _load_photon_data_arrays(data, ph_data, meas_type=meas_type, ondisk=ondisk)

    # If nanotimes are present load their specs
    if 'nanotimes' in ph_data:
        _load_nanotimes_specs(data, ph_data)

    # Unless 1-color, load donor and acceptor info
    if meas_type != 'smFRET-1color':
        donor = np.asscalar(meas_specs.detectors_specs.spectral_ch1.read())
        accept = np.asscalar(meas_specs.detectors_specs.spectral_ch2.read())
        _append_data_ch(data, 'det_donor_accept', (donor, accept))

    # Load alternation definition both for ns-ALEX and us-ALEX
    if 'ALEX' in meas_type:
        _load_alex_periods_donor_acceptor(data, meas_specs)

    # Here there are all the special-case for each measurement type
    if meas_type == 'smFRET-1color':
        # Non-FRET or unspecified data, assume all photons are "acceptor"
        _append_data_ch(data, 'A_em', slice(None))

    elif meas_type == 'smFRET':
        _compute_acceptor_emission_mask(data, ich, ondisk=ondisk)

    elif meas_type == 'smFRET-usALEX':
        _add_usALEX_specs(data, meas_specs)

    elif meas_type == 'smFRET-nsALEX':
        data.add(laser_repetition_rate=meas_specs.laser_repetition_rate.read())

    # Set some `data` flags
    data.add(ALEX='ALEX' in meas_type)
    data.add(lifetime='nanotimes' in ph_data)


def _photon_hdf5_multich(h5data, data, ondisk=True):

    ph_times_dict = phc.hdf5.photon_data_mapping(h5data._v_file)
    nch = np.max(list(ph_times_dict.keys())) + 1

    for ich in range(nch):
        _photon_hdf5_1ch(h5data, data, ondisk=ondisk, nch=nch, ich=ich)


def photon_hdf5(filename, ondisk=False, strict=False):
    """Load a data file saved in Photon-HDF5 format version 0.3 or higher.

    Photon-HDF5 is a format for a wide range of timestamp-based
    single molecule data. For more info please see:

    http://photon-hdf5.org/

    Arguments:
        filename (str or pathlib.Path): path of the data file to be loaded.
        ondisk (bool): if True, do not load the timestamps in memory
            using instead references to the HDF5 arrays. Default False.

    Returns:
        :class:`fretbursts.burstlib.Data` object containing the data.
    """
    filename = str(filename)
    assert os.path.isfile(filename), 'File not found.'
    version = phc.hdf5._check_version(filename)
    if version == u'0.2':
        return loader_legacy.hdf5(filename)

    h5file = tables.open_file(filename)
    h5data = h5file.root
    d = Data(fname=filename, data_file=h5data._v_file)

    for grp_name in ['setup', 'sample', 'provenance', 'identity']:
        if grp_name in h5data:
            d.add(**{grp_name:
                     phc.hdf5.dict_from_group(h5data._f_get_child(grp_name))})

    for field_name in ['description', 'acquisition_duration']:
        if field_name in h5data:
            d.add(**{field_name: h5data._f_get_child(field_name).read()})

    if _is_multich(h5data):
        _photon_hdf5_multich(h5data, d, ondisk=ondisk)
    else:
        _photon_hdf5_1ch(h5data, d, ondisk=ondisk)

    return d


##
# Multi-spot loader functions
#
def multispot8(fname, bytes_to_read=-1, swap_D_A=True, leakage=0, gamma=1.):
    """Load a 8-ch multispot file and return a Data() object. Cached version.
    """
    fname_c = fname + '_cache.pickle'
    try:
        var = pickle.load(open(fname_c, 'rb'))
        dx = Data(fname=fname, clk_p=12.5e-9, nch=8, leakage=leakage,
                  gamma=gamma)
        dx.add(ph_times_m=var['ph_times_m'], A_em=var['A_em'], ALEX=False)
        pprint(" - File loaded from cache: %s\n" % fname)
    except IOError:
        dx = multispot8_core(fname, bytes_to_read=bytes_to_read,
                             swap_D_A=swap_D_A, leakage=leakage, gamma=gamma)
        D = {'ph_times_m': dx.ph_times_m, 'A_em': dx.A_em}
        pprint(" - Pickling data ... ")
        pickle.dump(D, open(fname_c, 'wb'), -1)
        pprint("DONE\n")
    return dx

load_multispot8 = deprecate(multispot8, "load_multispot8", "loader.multispot8")

def multispot8_core(fname, bytes_to_read=-1, swap_D_A=True, leakage=0,
                    gamma=1.):
    """Load a 8-ch multispot file and return a Data() object.
    """
    dx = Data(fname=fname, clk_p=12.5e-9, nch=8, leakage=leakage,
              gamma=gamma)
    ph_times_m, A_em, ph_times_det = load_data_ordered16(
        fname=fname, n_bytes_to_read=bytes_to_read, swap_D_A=swap_D_A)
    dx.add(ph_times_m=ph_times_m, A_em=A_em, ALEX=False)
    return dx

def multispot48_simple(fname, leakage=0, gamma=1.,
                       i_start=0, i_stop=None, debug=False):
    """Load a 48-ch multispot file and return a Data() object.
    """
    dx = Data(fname=fname, clk_p=10e-9, nch=48, leakage=leakage, gamma=gamma)
    ph_times_m, big_fifo, ch_fifo = load_manta_timestamps(
        fname, i_start=i_start, i_stop=i_stop, debug=debug)
    A_em = [True] * len(ph_times_m)
    dx.add(ph_times_m=ph_times_m, A_em=A_em, ALEX=False)
    big_fifo_full = np.array([b.any() for b in big_fifo]).any()
    ch_fifo_full = np.array([b.any() for b in ch_fifo]).any()
    if big_fifo_full:
        print('WARNING: Big-FIFO full, flags saved in Data()')
        dx.add(big_fifo=big_fifo)
    if ch_fifo_full:
        print('WARNING: CH-FIFO full, flags saved in Data()')
        dx.add(ch_fifo=ch_fifo)
    return dx

def multispot48(fname, leakage=0, gamma=1., reprocess=False,
                i_start=0, i_stop=None, debug=False):
    """Load a 48-ch multispot file and return a Data() object.
    """
    basename, ext = os.path.splitext(fname)
    fname_h5 = basename + '.hdf5'
    fname_dat = basename + '.dat'

    def load_dat_file():
        pprint(' - Loading DAT file: %s ... ' % fname_dat)
        # Load data from raw file and store it in a HDF5 file
        data = load_xavier_manta_data(fname_dat, i_start=i_start,
                                      i_stop=i_stop, debug=debug)
        pprint('DONE.\n - Extracting timestamps and detectors ... ')
        timestamps, det = get_timestamps_detectors(data, nbits=24)
        pprint('DONE.\n - Processing and storing ... ')
        ph_times_m, big_fifo, ch_fifo = process_store(
            timestamps, det, out_fname=fname_h5, fifo_flag=True, debug=False)
        pprint('DONE.\n')
        return ph_times_m, big_fifo, ch_fifo

    if not (os.path.isfile(fname_dat) or os.path.isfile(fname_h5)):
        raise IOError('Data file "%s" not found' % basename)

    if os.path.exists(fname_h5) and not reprocess:
        # There is a HDF5 file
        try:
            pprint(' - Loading HDF5 file: %s ... ' % fname_h5)
            ph_times_m, big_fifo, ch_fifo = \
                load_manta_timestamps_pytables(fname_h5)
            pprint('DONE.\n')
        except tables.HDF5ExtError:
            pprint('\n  Ops! File may be truncated.\n')
            ph_times_m, big_fifo, ch_fifo = load_dat_file()
    else:
        ph_times_m, big_fifo, ch_fifo = load_dat_file()

    # Current data has only acceptor ch
    A_em = [True] * len(ph_times_m)

    dx = Data(fname=fname, clk_p=10e-9, nch=48, leakage=leakage, gamma=gamma)
    dx.add(ph_times_m=ph_times_m, A_em=A_em, ALEX=False,
           data_file=ph_times_m.data_file, bg_data_file=ph_times_m.data_file)
    big_fifo_full = np.array([b[:].any() for b in big_fifo]).any()
    ch_fifo_full = np.array([b[:].any() for b in ch_fifo]).any()
    if big_fifo_full:
        print('WARNING: Big-FIFO full, flags saved in Data()')
        dx.add(big_fifo=big_fifo)
    if ch_fifo_full:
        print('WARNING: CH-FIFO full, flags saved in Data()')
        dx.add(ch_fifo=ch_fifo)
    return dx


##
# usALEX loader functions
#

# Build masks for the alternating periods
def _select_outer_range(times, period, edges):
    return ((times % period) >= edges[0]) + ((times % period) < edges[1])

def _select_inner_range(times, period, edges):
    return ((times % period) >= edges[0]) * ((times % period) < edges[1])

def _select_range(times, period, edges):
    return _select_inner_range(times, period, edges) if edges[0] < edges[1] \
        else _select_outer_range(times, period, edges)

def usalex(fname, leakage=0, gamma=1., header=None, BT=None):
    """Load usALEX data from a SM file and return a Data() object.

    This function returns a Data() object to which you need to apply
    an alternation selection before performing further analysis (background
    estimation, burst search, etc.).

    The pattern to load usALEX data is the following::

        d = loader.usalex(fname=fname)
        d.add(D_ON=(2850, 580), A_ON=(900, 2580), alex_period=4000)
        plot_alternation_hist(d)

    If the plot looks good, apply the alternation with::

        loader.alex_apply_period(d)

    Now `d` is ready for futher processing such as background estimation,
    burst search, etc...
    """
    if BT is not None:
        print('WARNING: `BT` argument is deprecated, use `leakage` instead.')
        leakage = BT
    if header is not None:
        print('WARNING: `header` argument ignored. '
              '         The header length is now computed automatically.')
    print(" - Loading '%s' ... " % fname)
    ph_times_t, det_t, labels = load_sm(fname, return_labels=True)
    print(" [DONE]\n")

    DONOR_ON = (2850, 580)
    ACCEPT_ON = (930, 2580)
    alex_period = 4000

    dx = Data(fname=fname, clk_p=12.5e-9, nch=1, leakage=leakage, gamma=gamma,
              ALEX=True, lifetime=False,
              D_ON=DONOR_ON, A_ON=ACCEPT_ON, alex_period=alex_period,
              ph_times_t=[ph_times_t], det_t=[det_t], det_donor_accept=(0, 1),
              ch_labels=labels)
    return dx

def usalex_apply_period(d, delete_ph_t=True, remove_d_em_a_ex=False):
    """Applies to the Data object `d` the alternation period previously set.

    Note that you first need to load the data in a variable `d` and then
    set the alternation parameters using `d.add(D_ON=..., A_ON=...)`.

    The typical pattern for loading ALEX data is the following::

        d = loader.photon_hdf5(fname=fname)
        d.add(D_ON=(2850, 580), A_ON=(900, 2580))
        alex_plot_alternation(d)

    If the plot looks good, apply the alternation with::

        loader.alex_apply_period(d)

    Now `d` is ready for futher processing such as background estimation,
    burst search, etc...

    *See also:* :func:`alex_apply_period`.
    """
    ich = 0
    donor_ch, accept_ch = d._det_donor_accept_multich[ich]
    D_ON, A_ON = d._D_ON_multich[ich], d._A_ON_multich[ich]
    # Remove eventual ch different from donor or acceptor
    d_ch_mask_t = (d.det_t[ich] == donor_ch)
    a_ch_mask_t = (d.det_t[ich] == accept_ch)
    valid_mask = d_ch_mask_t + a_ch_mask_t
    ph_times_val = d.ph_times_t[ich][valid_mask]
    d_ch_mask_val = d_ch_mask_t[valid_mask]
    a_ch_mask_val = a_ch_mask_t[valid_mask]
    assert (d_ch_mask_val + a_ch_mask_val).all()
    assert not (d_ch_mask_val * a_ch_mask_val).any()
    if 'offset' in d:
        ph_times_val -= d.offset

    # Build masks for excitation windows
    d_ex_mask_val = _select_range(ph_times_val, d.alex_period, D_ON)
    a_ex_mask_val = _select_range(ph_times_val, d.alex_period, A_ON)
    # Safety check: each ph is either D or A ex (not both)
    assert not (d_ex_mask_val * a_ex_mask_val).any()

    # Select alternation periods, removing transients
    DexAex_mask = d_ex_mask_val + a_ex_mask_val

    # Reduce timestamps and masks to the DexAex_mask selection
    ph_times = ph_times_val[DexAex_mask]
    d_em = d_ch_mask_val[DexAex_mask]
    a_em = a_ch_mask_val[DexAex_mask]
    d_ex = d_ex_mask_val[DexAex_mask]
    a_ex = a_ex_mask_val[DexAex_mask]
    assert d_ex.sum() == d_ex_mask_val.sum()
    assert a_ex.sum() == a_ex_mask_val.sum()

    if remove_d_em_a_ex:
        # Removes donor-ch photons during acceptor excitation
        mask = a_em + d_em * d_ex
        assert (mask == -(a_ex * d_em)).all()

        ph_times = ph_times[mask]
        d_em = d_em[mask]
        a_em = a_em[mask]
        d_ex = d_ex[mask]
        a_ex = a_ex[mask]

    assert d_em.sum() + a_em.sum() == ph_times.size
    assert (d_em + a_em).all()       # masks fill the total array
    assert not (d_em * a_em).any()   # no photon is both D and A
    assert a_ex.size == a_em.size == d_ex.size == d_em.size == ph_times.size
    d.add(ph_times_m=[ph_times],
          D_em=[d_em], A_em=[a_em], D_ex=[d_ex], A_ex=[a_ex])

    if 'particles_t' in d:
        particles_val = d.particles_t[ich][valid_mask]
        particles = particles_val[DexAex_mask]
        d.add(particles=[particles])

    assert d.ph_times_m[ich].size == d.A_em[ich].size

    if delete_ph_t:
        d.delete('ph_times_t')
        d.delete('det_t')
    return d


##
# nsALEX loader functions
#

def nsalex(fname):
    """Load nsALEX data from a SPC file and return a Data() object.

    This function returns a Data() object to which you need to apply
    an alternation selection before performing further analysis (background
    estimation, burst search, etc.).

    The pattern to load nsALEX data is the following::

        d = loader.nsalex(fname=fname)
        d.add(D_ON=(2850, 580), A_ON=(900, 2580))
        alex_plot_alternation(d)

    If the plot looks good apply the alternation with::

        loader.alex_apply_period(d)

    Now `d` is ready for futher processing such as background estimation,
    burst search, etc...
    """
    ph_times_t, det_t, nanotimes = load_spc(fname)

    DONOR_ON = (10, 1500)
    ACCEPT_ON = (2000, 3500)
    nanotimes_nbins = 4095

    dx = Data(fname=fname, clk_p=50e-9, nch=1, ALEX=True, lifetime=True,
              D_ON=DONOR_ON, A_ON=ACCEPT_ON,
              nanotimes_nbins=nanotimes_nbins,
              nanotimes_params=[{'tcspc_num_bins': nanotimes_nbins}],
              ph_times_t=[ph_times_t], det_t=[det_t], nanotimes_t=[nanotimes],
              det_donor_accept=(4, 6))
    return dx

def nsalex_apply_period(d, delete_ph_t=True):
    """Applies to the Data object `d` the alternation period previously set.

    Note that you first need to load the data in a variable `d` and then
    set the alternation parameters using `d.add(D_ON=..., A_ON=...)`.

    The typical pattern for loading ALEX data is the following::

        d = loader.photon_hdf5(fname=fname)
        d.add(D_ON=(2850, 580), A_ON=(900, 2580))
        alex_plot_alternation(d)

    If the plot looks good, apply the alternation with::

        loader.alex_apply_period(d)

    Now `d` is ready for futher processing such as background estimation,
    burst search, etc...

    *See also:* :func:`alex_apply_period`.
    """
    ich = 0
    donor_ch, accept_ch = d._det_donor_accept_multich[ich]
    D_ON_multi, A_ON_multi = d._D_ON_multich[ich], d._A_ON_multich[ich]
    D_ON = [(D_ON_multi[i], D_ON_multi[i+1]) for i in range(0, len(D_ON_multi), 2)]
    A_ON = [(A_ON_multi[i], A_ON_multi[i+1]) for i in range(0, len(A_ON_multi), 2)]

    # Mask for donor + acceptor detectors (discard other detectors)
    d_ch_mask_t = (d.det_t[ich] == donor_ch)
    a_ch_mask_t = (d.det_t[ich] == accept_ch)
    da_ch_mask_t = d_ch_mask_t + a_ch_mask_t

    # Masks for excitation periods
    d_ex_mask_t = np.zeros(d.nanotimes_t[ich].size, dtype='bool')
    for d_on in D_ON:
        d_ex_mask_t += (d.nanotimes_t[ich] > d_on[0]) * (d.nanotimes_t[ich] < d_on[1])

    a_ex_mask_t = np.zeros(d.nanotimes_t[ich].size, dtype='bool')
    for a_on in A_ON:
        a_ex_mask_t += (d.nanotimes_t[ich] > a_on[0]) * (d.nanotimes_t[ich] < a_on[1])

    ex_mask_t = d_ex_mask_t + a_ex_mask_t  # Select only ph during Dex or Aex

    # Total mask: D+A photons, and only during the excitation periods
    mask = da_ch_mask_t * ex_mask_t  # logical AND

    # Apply selection to timestamps and nanotimes
    ph_times = d.ph_times_t[ich][mask]
    nanotimes = d.nanotimes_t[ich][mask]

    # Apply selection to the emission masks
    d_em = d_ch_mask_t[mask]
    a_em = a_ch_mask_t[mask]
    assert (d_em + a_em).all()       # masks fill the total array
    assert not (d_em * a_em).any()   # no photon is both D and A

    # Apply selection to the excitation masks
    d_ex = d_ex_mask_t[mask]
    a_ex = a_ex_mask_t[mask]
    assert (d_ex + a_ex).all()
    assert not (d_ex * a_ex).any()

    d.add(ph_times_m=[ph_times], nanotimes=[nanotimes],
          D_em=[d_em], A_em=[a_em], D_ex=[d_ex], A_ex=[a_ex],)

    if delete_ph_t:
        d.delete('ph_times_t')
        d.delete('det_t')
        d.delete('nanotimes_t')

def alex_apply_period(d, delete_ph_t=True):
    """Apply the ALEX period definition set in D_ON and A_ON attributes.

    This function works both for us-ALEX and ns-ALEX data.

    Note that you first need to load the data in a variable `d` and then
    set the alternation parameters using `d.add(D_ON=..., A_ON=...)`.

    The typical pattern for loading ALEX data is the following::

        d = loader.photon_hdf5(fname=fname)
        d.add(D_ON=(2850, 580), A_ON=(900, 2580))
        alex_plot_alternation(d)

    If the plot looks good, apply the alternation with::

        loader.alex_apply_period(d)

    Now `d` is ready for futher processing such as background estimation,
    burst search, etc...
    """
    assert d.ALEX
    if d.lifetime:
        apply_period_func = nsalex_apply_period
    else:
        apply_period_func = usalex_apply_period
    apply_period_func(d, delete_ph_t=delete_ph_t)
    msg = ('# Total photons (after ALEX selection):  {:10,}\n'
           '#  D  photons in D+A excitation periods: {:10,}\n'
           '#  A  photons in D+A excitation periods: {:10,}\n'
           '# D+A photons in  D  excitation period:  {:10,}\n'
           '# D+A photons in  A  excitation period:  {:10,}\n')
    print(msg.format(d.ph_times_m[0].size, d.D_em[0].sum(), d.A_em[0].sum(),
                     d.D_ex[0].sum(), d.A_ex[0].sum()))
