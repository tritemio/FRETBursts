#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
This module contains functions to load each supported data format.
The loader functions load data from a specific format and
return a new :class:`fretbursts.burstlib.Data()` object containing the data.

This file contains the high-level function to load a data-file and
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

from .dataload.multi_ch_reader import load_data_ordered16
from .dataload.smreader import load_sm
from .dataload.spcreader import load_spc
from .dataload.manta_reader import (load_manta_timestamps,
                                    load_xavier_manta_data,
                                    get_timestamps_detectors,
                                    #process_timestamps,
                                    process_store,
                                    load_manta_timestamps_pytables)
from .utils.misc import pprint, deprecate
from .burstlib import Data
from .hdf5 import hdf5_data_map

from phconvert.hdf5 import dict_from_group
import phconvert as phc


def _append_data_ch(d, name, value):
    if name not in d:
        d.add(**{name: [value]})
    else:
        d[name].append(value)

def _load_from_group(d, group, name, dest_name, ich=None,
                     ondisk=False, allow_missing=True):
    if name not in group:
        return

    node_value = group._f_get_child(name)
    if not ondisk:
        node_value = node_value.read()
    if ich is None:
        d.add(**{dest_name: node_value})
    else:
        _append_data_ch(d, dest_name, node_value)

def _is_multich(h5data):
    if 'photon_data' in h5data:
        return False
    elif 'photon_data0' in h5data:
        return True
    else:
        msg = 'Cannot find a photon_data group.'
        raise phc.hdf5.Invalid_PhotonHDF5(msg)

def _photon_hdf5_1ch(h5data, data, ondisk=False):
    data.add(nch=1)

    ph_data = h5data.photon_data
    if 'measurement_specs' not in ph_data:
        raise NotImplementedError('FRETBursts requires a "measurement_specs" '
                                  ' to interpret a Photon-HDF5 file.')

    assert 'measurement_type' in ph_data.measurement_specs
    meas_specs = ph_data.measurement_specs
    meas_type = meas_specs.measurement_type.read().decode()
    if meas_type not in ['smFRET', 'smFRET-usALEX', 'smFRET-nsALEX']:
        raise NotImplementedError('Meaurement type "%s" not supported'
                                  ' by FRETBursts.' % meas_type)

    ## Load the photon-data arrays
    if  meas_type == 'smFRET':
        mapping = {'timestamps': 'ph_times_m', 'detectors': 'A_em',
                   'nanotimes': 'nanotimes', 'particles': 'particles'}
        ich = 0  # Created a 1-element list for each field
    else:
        mapping = {'timestamps': 'ph_times_t', 'detectors': 'det_t',
                   'nanotimes': 'nanotimes_t', 'particles': 'particles_t'}
        ich = None  # don't warp the arrays in a list

    for name, dest_name in mapping.items():
        _load_from_group(data, ph_data, name, dest_name=dest_name, ich=ich,
                         ondisk=ondisk)

    ## Load the other parameters
    donor = meas_specs.detectors_specs.spectral_ch1.read(),
    accept = meas_specs.detectors_specs.spectral_ch2.read(),

    data.add(
        clk_p = ph_data.timestamps_specs.timestamps_unit.read(),
        det_donor_accept = (donor, accept))

    if meas_type == 'smFRET':
        # We still have to convert A_em to a boolean mask
        # because until now it is just the detector list

        # Make sure there are at most 2 detectors
        assert len(np.unique(data.A_em[0])) <= 2
        if accept and not donor:
            # In this case we can convert without a copy
            data.add(A_em=[data.A_em[0].view(dtype=bool)])
        else:
            # Create the boolean mask
            data.add(A_em=[data.A_em[0] == accept])

    if 'nanotimes_specs' in ph_data:
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
        data.add(nanotimes_params=nanotimes_params)

    if 'ALEX' in meas_type:
        data.add(
            D_ON = meas_specs.alex_period_spectral_ch1.read(),
            A_ON = meas_specs.alex_period_spectral_ch2.read())

    if meas_type == 'smFRET':
        data.add(ALEX=False, lifetime=False)
    elif meas_type == 'smFRET-usALEX':
        data.add(ALEX=True, lifetime=False,
                 alex_period=meas_specs.alex_period.read())
    elif meas_type == 'smFRET-nsALEX':
        data.add(ALEX=True, lifetime=True,
                 alex_period=meas_specs.laser_pulse_rate.read())


def _photon_hdf5_multich(h5data, data, ondisk=True):

    ph_times_dict = phc.hdf5.photon_data_mapping(h5data)
    nch = np.max(list(ph_times_dict.keys())) + 1

    data.add(nch=nch)
    for ich in range(data.nch):
        ph_data_name = '/photon_data%d' % ich

        if ph_data_name not in h5data:
            ph_times = np.array([], dtype='int64')
            _append_data_ch(data, 'ph_times_m', ph_times)

            a_em = np.array([], dtype=bool)
            _append_data_ch(data, 'A_em', a_em)
        else:
            ph_data = h5data._f_get_child(ph_data_name)

            _load_from_group(data, ph_data, name='timestamps',
                             dest_name='ph_times_m', allow_missing=False,
                             ich=ich, ondisk=ondisk)
            data.add(clk_p=ph_data.timestamps_specs.timestamps_unit.read())

            if 'detectors' not in ph_data:
                a_em = slice(None)
            else:
                assert 'measurement_specs' in ph_data
                meas_specs = ph_data.measurement_specs
                assert 'detectors_specs' in meas_specs
                det_specs = meas_specs.detectors_specs

                donor = det_specs.spectral_ch1.read()
                accept = det_specs.spectral_ch2.read()

                det = ph_data.detectors.read()
                if det.dtype.itemsize == 1 and donor == 0:
                    a_em = det.view(bool)
                elif det.dtype.itemsize == 1 and accept == 0:
                    a_em = np.logical_not(det.view(bool))
                else:
                    a_em = (det == accept)
                    d_em = (det == donor)
                    assert not (a_em*d_em).any()
                    assert (a_em + d_em).all()

            _append_data_ch(data, 'A_em', a_em)


def photon_hdf5(filename, ondisk=False, strict=False):
    """Load a data file saved in Photon-HDF5 format version 0.3 or higher.

    Any :class:`fretbursts.burstlib.Data` object can be saved in HDF5 format
    using :func:`fretbursts.hdf5.store` .

    For description and specs of the Photon-HDF5 format see:
    http://photon-hdf5.readthedocs.org/

    Arguments:
        ondisk (bool): if True do not load in the timestamp arrays, just
            load the array reference. Multi-spot only. Default False.

    Returns:
        A :class:`fretbursts.burstlib.Data` object containing the data.
    """
    version = phc.hdf5._check_version(filename)
    if version == '0.2':
        return hdf5(filename)

    h5data = phc.hdf5.load_photon_hdf5(filename, strict=strict)
    d = Data(fname=filename, data_file=h5data._v_file)

    for grp_name in ['setup', 'sample', 'provenance', 'identity']:
        if grp_name in h5data:
            d.add(**{grp_name:
                     phc.hdf5.dict_from_group(h5data._f_get_child(grp_name))})

    for field_name in ['comment', 'acquisition_time']:
        if field_name in h5data:
            d.add(**{field_name: h5data._f_get_child(field_name).read()})

    if _is_multich(h5data):
        _photon_hdf5_multich(h5data, d, ondisk=ondisk)
    else:
        _photon_hdf5_1ch(h5data, d, ondisk=ondisk)

    return d

def assert_valid_photon_hdf5(h5file):
    meta = dict(format_name=b'Photon-HDF5', format_version=b'0.2')
    msg = ''
    for attr, value in meta.items():
        if attr not in h5file.root._v_attrs:
            msg += ' * Error: %s not in %s\n' % (attr, h5file.root._v_attrs)
        stored_value = h5file.root._v_attrs[attr]
        if stored_value != value:
            msg += ' * Error: %s != %s\n' % (stored_value, value)
    if msg is not '':
        h5file.close()
        msg = 'Not a valid Photon-HDF5 file. \n' + msg
        raise IOError(msg)

def _is_basic_layout(h5file):
    return 'photon_data' in h5file.root

class H5Loader():

    def __init__(self, h5file, data):
        self.h5file = h5file
        self.data = data

    def load_data(self, where, name, dest_name=None, ich=None,
                  allow_missing=False, ondisk=False):
        try:
            node = self.h5file.get_node(where, name)
        except tables.NoSuchNodeError:
            if allow_missing:
                node_value = np.array([])
            else:
                self.h5file.close()
                raise IOError("Invalid file format: '%s' is missing." % name)
        else:
            node_value = node if ondisk else node.read()

        if dest_name is None:
            dest_name = hdf5_data_map.get(name, name)

        if ich is None:
            self.data.add(**{dest_name: node_value})
        else:
            if dest_name not in self.data:
                self.data.add(**{dest_name: [node_value]})
            else:
                self.data[dest_name].append(node_value)

def hdf5(fname, ondisk=False):
    """Load a data file saved in Photon-HDF5 format, version 0.2 only.

    The format version 0.2 is obsolete, please use :func:`photon_hdf5`
    for Photon-HDF5 version 0.3 or higher. If you have files in
    Photon-HDF5 format 0.2 please convert them to version 0.3 or higher.

    For description and specs of the Photon-HDF5 format see:
    http://photon-hdf5.readthedocs.org/

    Arguments:
        ondisk (bool): if True do not load in the timestamp arrays, just
            load the array reference. Multi-spot only. Default False.

    Returns:
        A :class:`fretbursts.burstlib.Data` object containing the data.
    """
    print('WARNING: You are using Photon-HDF5 version 0.2 which is '
          'obsolete. Please update you file to version 0.3 or higher.')

    # This used for reading Photon-HDF5 version 0.2.
    mandatory_root_fields = ['timestamps_unit', 'num_spots', 'num_detectors',
                             'num_spectral_ch', 'num_polariz_ch',
                             'alex', 'lifetime',]

    if not os.path.isfile(fname):
        raise IOError('File not found.')
    data_file = tables.open_file(fname, mode="r")
    assert_valid_photon_hdf5(data_file)

    d = Data(fname=fname)
    loader = H5Loader(data_file, d)

    # Load mandatory parameters
    mand_fields = [f for f in mandatory_root_fields if f != 'num_detectors']
    for field in mand_fields:
        loader.load_data('/', field)
    try:
        data_file.get_node('/', 'num_detectors')
    except tables.NoSuchNodeError:
        pprint("WARNING: No 'num_detectors' field found in the HDF5 file.\n")
    else:
        loader.load_data('/', 'num_detectors')

    if d.ALEX:
        if d.lifetime:
            loader.load_data('/', 'laser_pulse_rate', allow_missing=True)
        else:
            loader.load_data('/', 'alex_period')
        loader.load_data('/', 'alex_period_donor', allow_missing=True)
        loader.load_data('/', 'alex_period_acceptor', allow_missing=True)

    ## Load metadata groups
    metagroups = ['/setup', '/provenance', '/identity']
    for group in metagroups:
        if group in data_file:
            d.add(**{group[1:]: dict_from_group(data_file.get_node(group))})

    if '/acquisition_time' in data_file:
        d.add(acquisition_time=data_file.root.acquisition_time.read())

    if _is_basic_layout(data_file):
        ph_group = data_file.root.photon_data

    if d.lifetime:
        try:
            assert 'nanotimes' in ph_group
            assert 'nanotimes_specs' in ph_group
        except AssertionError:
            data_file.close()
            raise IOError(('The lifetime flag is True but the TCSPC '
                           'data is missing.'))

    if d.nch == 1:
        # load single-spot data from "basic layout"
        if not d.ALEX:
            mapping = {'timestamps': 'ph_times_m', 'detectors': 'A_em',
                       'nanotimes': 'nanotimes', 'particles': 'particles'}
            ich = 0  # Created a 1-element list for each field
        else:
            mapping = {'timestamps': 'ph_times_t', 'detectors': 'det_t',
                       'nanotimes': 'nanotimes_t', 'particles': 'particles_t'}
            ich = None  # don't warp the arrays in a list
        for name, dest_name in mapping.items():
            if name in ph_group:
                loader.load_data(ph_group, name, dest_name=dest_name, ich=ich)

        if 'detectors_specs' in ph_group:
            det_specs = ph_group.detectors_specs
            if 'donor' in det_specs and 'acceptor' in det_specs:
                donor = det_specs.donor.read()
                accept = det_specs.acceptor.read()
                d.add(det_donor_accept=(donor, accept))
                if not d.ALEX:
                    # We still have to convert A_em to a boolean mask
                    # because until now it is just the detector list

                    # Make sure there are at most 2 detectors
                    assert len(np.unique(d.A_em[0])) <= 2
                    if accept and not donor:
                        # In this case we can convert without a copy
                        d.add(A_em=[d.A_em[0].view(dtype=bool)])
                    else:
                        # Create the boolean mask
                        d.add(A_em=[d.A_em[0] == accept])

        if 'nanotimes_specs' in ph_group:
            nanot_specs = ph_group.nanotimes_specs
            nanotimes_params = {}
            for name in ['tcspc_unit', 'tcspc_num_bins', 'tcspc_range']:
                value = nanot_specs._f_get_child(name).read()
                nanotimes_params.update(**{name: value})
            for name in ['tau_accept_only', 'tau_donor_only',
                         'tau_fret_donor', 'inverse_fret_rate']:
                if name in nanot_specs:
                    value = nanot_specs._f_get_child(name).read()
                    nanotimes_params.update(**{name: value})
            d.add(nanotimes_params=nanotimes_params)

    else:
        # Load multi-spot data from multi-spot layout
        for ich in range(d.nch):
            ph_group_name = '/photon_data_%d' % ich
            loader.load_data(ph_group_name, 'timestamps', allow_missing=True,
                             dest_name='ph_times_m', ich=ich, ondisk=ondisk)

            if ph_group_name not in data_file:
                a_em = np.array([], dtype=bool)
            elif ph_group_name + '/detectors' not in data_file:
                a_em = slice(None)
            else:
                ph_group = data_file.root._f_get_child(ph_group_name)
                det_specs = ph_group.detectors_specs
                donor = det_specs.donor.read()
                accept = det_specs.acceptor.read()
                if ph_group.detectors.dtype == np.bool:
                    a_em = ph_group.detectors.read()
                    if not accept:
                        np.logical_not(a_em, out=a_em)
                else:
                    det = ph_group.detectors.read()
                    a_em = (det == accept)
                    d_em = (det == donor)
                    assert not (a_em*d_em).any()
                    assert (a_em + d_em).all()

            if 'A_em' not in d:
                d.add(A_em=[a_em])
            else:
                d.A_em.append(a_em)

    d.add(data_file=data_file)
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
        ## Load data from raw file and store it in a HDF5 file
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
        ## There is a HDF5 file
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

    ## Current data has only acceptor ch
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
    return ((times % period) > edges[0]) + ((times % period) < edges[1])

def _select_inner_range(times, period, edges):
    return ((times % period) > edges[0]) * ((times % period) < edges[1])

def _select_range(times, period, edges):
    return _select_inner_range(times, period, edges) if edges[0] < edges[1] \
            else _select_outer_range(times, period, edges)

def usalex(fname, leakage=0, gamma=1., header=None, start=None, stop=None,
           BT=None):
    """Load a usALEX file and return a Data() object.

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
    ph_times_t, det_t, labels = load_sm(fname, start=start, stop=stop,
                                        return_labels=True)
    print(" [DONE]\n")

    DONOR_ON = (2850, 580)
    ACCEPT_ON = (930, 2580)
    alex_period = 4000

    dx = Data(fname=fname, clk_p=12.5e-9, nch=1, leakage=leakage, gamma=gamma,
              ALEX=True, lifetime=False,
              D_ON=DONOR_ON, A_ON=ACCEPT_ON, alex_period=alex_period,
              ph_times_t=ph_times_t, det_t=det_t, det_donor_accept=(0, 1),
              ch_labels=labels)
    return dx

def usalex_apply_period(d, delete_ph_t=True, remove_d_em_a_ex=False):
    """Applies to the Data object `d` the alternation period previously set.

    Note that you first need to load the data in a variable `d` and then
    set the alternation parameters using `d.add(D_ON=..., A_ON=...)`.

    The typical pattern for loading ALEX data is the following::

        d = loader.hdf5(fname=fname)
        d.add(D_ON=(2850, 580), A_ON=(900, 2580))
        alex_plot_alternation(d)

    If the plot looks good, apply the alternation with::

        loader.alex_apply_period(d)

    Now `d` is ready for futher processing such as background estimation,
    burst search, etc...

    *See also:* :func:`alex_apply_period`.
    """
    donor_ch, accept_ch = d.det_donor_accept
    # Remove eventual ch different from donor or acceptor
    d_ch_mask_t = (d.det_t == donor_ch)
    a_ch_mask_t = (d.det_t == accept_ch)
    valid_mask = d_ch_mask_t + a_ch_mask_t
    ph_times_val = d.ph_times_t[valid_mask]
    d_ch_mask_val = d_ch_mask_t[valid_mask]
    a_ch_mask_val = a_ch_mask_t[valid_mask]
    assert (d_ch_mask_val + a_ch_mask_val).all()
    assert not (d_ch_mask_val * a_ch_mask_val).any()

    # Build masks for excitation windows
    d_ex_mask_val = _select_range(ph_times_val, d.alex_period, d.D_ON)
    a_ex_mask_val = _select_range(ph_times_val, d.alex_period, d.A_ON)
    # Safety check: each ph is either D or A ex (not both)
    assert not (d_ex_mask_val * a_ex_mask_val).any()

    mask = d_ex_mask_val + a_ex_mask_val  # Removes alternation transients

    # Assign the new ph selection mask
    ph_times = ph_times_val[mask]
    d_em = d_ch_mask_val[mask]
    a_em = a_ch_mask_val[mask]
    d_ex = d_ex_mask_val[mask]
    a_ex = a_ex_mask_val[mask]

    if remove_d_em_a_ex:
        # Removes donor-ch photons during acceptor excitation
        mask = a_em + d_em*d_ex
        assert (mask == -(a_ex*d_em)).all()

        ph_times = ph_times[mask]
        d_em = d_em[mask]
        a_em = a_em[mask]
        d_ex = d_ex[mask]
        a_ex = a_ex[mask]

    assert d_em.sum() + a_em.sum() == ph_times.size
    assert (d_em * a_em).any() == False
    assert a_ex.size == a_em.size == d_ex.size == d_em.size == ph_times.size
    print("#donor: %d  #acceptor: %d \n" % (d_em.sum(), a_em.sum()))

    d.add(ph_times_m=[ph_times],
          D_em=[d_em], A_em=[a_em], D_ex=[d_ex], A_ex=[a_ex],)

    assert d.ph_times_m[0].size == d.A_em[0].size

    if delete_ph_t:
        d.delete('ph_times_t')
        d.delete('det_t')
    return d


##
# nsALEX loader functions
#

def nsalex(fname):
    """Load a nsALEX file and return a Data() object.

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
              nanotimes_params = {'tcspc_num_bins': nanotimes_nbins},
              ph_times_t=ph_times_t, det_t=det_t, nanotimes_t=nanotimes,
              det_donor_accept=(4, 6))
    return dx

def nsalex_apply_period(d, delete_ph_t=True):
    """Applies to the Data object `d` the alternation period previously set.

    Note that you first need to load the data in a variable `d` and then
    set the alternation parameters using `d.add(D_ON=..., A_ON=...)`.

    The typical pattern for loading ALEX data is the following::

        d = loader.hdf5(fname=fname)
        d.add(D_ON=(2850, 580), A_ON=(900, 2580))
        alex_plot_alternation(d)

    If the plot looks good, apply the alternation with::

        loader.alex_apply_period(d)

    Now `d` is ready for futher processing such as background estimation,
    burst search, etc...

    *See also:* :func:`alex_apply_period`.
    """
    # Note: between boolean arrays * is equavilent to logical AND,
    #       and + is equivalent to logical OR.

    donor_ch, accept_ch = d.det_donor_accept

    # Mask for donor + acceptor detectors (discard other detectors)
    d_ch_mask_t = (d.det_t == donor_ch)
    a_ch_mask_t = (d.det_t == accept_ch)
    da_ch_mask_t = d_ch_mask_t + a_ch_mask_t

    # Masks for excitation periods
    d_ex_mask_t = (d.nanotimes_t > d.D_ON[0]) * (d.nanotimes_t < d.D_ON[1])
    a_ex_mask_t = (d.nanotimes_t > d.A_ON[0]) * (d.nanotimes_t < d.A_ON[1])
    ex_mask_t = d_ex_mask_t + a_ex_mask_t  # Select only ph during Dex or Aex

    # Total mask: D+A photons, and only during the excitation periods
    mask = da_ch_mask_t * ex_mask_t  # logical AND

    # Apply selection to timestamps and nanotimes
    ph_times = d.ph_times_t[mask]
    nanotimes = d.nanotimes_t[mask]

    # Apply selection to the emission masks
    d_em = d_ch_mask_t[mask]
    a_em = a_ch_mask_t[mask]
    assert (d_em + a_em).all()
    assert not (d_em * a_em).any()

    # Apply selection to the excitation masks
    d_ex = d_ex_mask_t[mask]
    a_ex = a_ex_mask_t[mask]
    assert (d_ex + a_ex).all()
    assert not (d_ex * a_ex).any()

    d.add(ph_times_m=[ph_times], nanotimes=nanotimes,
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

        d = loader.hdf5(fname=fname)
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
