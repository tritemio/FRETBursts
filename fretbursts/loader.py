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

import os
import numpy as np
import cPickle as pickle
import tables

from dataload.multi_ch_reader import load_data_ordered16
from dataload.smreader import load_sm
from dataload.spcreader import load_spc
from dataload.manta_reader import (load_manta_timestamps,
                                   load_xavier_manta_data,
                                   get_timestamps_detectors,
                                   #process_timestamps,
                                   process_store,
                                   load_manta_timestamps_pytables)
from utils.misc import pprint, deprecate
from burstlib import Data
from dataload.pytables_array_list import PyTablesList
from hdf5 import hdf5_data_map


def _is_valid_hdf5_phdata(h5file):
    meta = dict(format_name = 'HDF5-Ph-Data', format_version = '0.2')
    for attr, value in meta.items():
        if attr not in h5file.root._v_attrs:
            return False
        if h5file.root._v_attrs[attr] != meta[attr]:
            return False
    return True

def _is_basic_layout(h5file):
    return 'photon_data' in h5file.root

class H5Loader():

    def __init__(self, h5file, data):
        self.h5file = h5file
        self.data = data

    def load_data(self, where, name, dest_name=None):
        try:
            node = self.h5file.get_node(where, name)
        except tables.NoSuchNodeError:
            raise (IOError, "Invalid file format: '%s' is missing." % name)

        if dest_name is None:
            dest_name = hdf5_data_map[name]
        self.data.add(**{dest_name: node.read()})

def hdf5_phdata(fname):
    """Load a data file saved in HDF5-Ph-Data format.

    Any :class:`fretbursts.burstlib.Data` object can be saved in HDF5 format
    using :func:`fretbursts.hdf5.store` .
    """
    if not os.path.isfile(fname):
        raise IOError, 'File not found.'
    data_file = tables.open_file(fname, mode = "r")
    if not _is_valid_hdf5_phdata(data_file):
        raise (IOError, 'The file is not a valid HDF5-Ph-Data format.')

    # Default values for some parameters
    params = dict(leakage=0., gamma=1.)
    d = Data(fname=fname, **params)
    loader = H5Loader(data_file, d)

    # Load mandatory parameters
    mandatory_fields = ['timestamps_unit', 'number_confocal_spots', 'alex',
                        'lifetime']
    for field in mandatory_fields:
        loader.load_data('/', field)

    if d.ALEX:
        loader.load_data('/', 'alternation_period')

    if d.lifetime:
        loader.load_data('/', 'nanotime_unit')

    if _is_basic_layout(data_file):
        ph_group = data_file.root.photon_data
        loader.load_data(ph_group, 'timestamps')

        if d.lifetime:
            try:
                assert 'nanotimes' in ph_group
                assert 'nanotimes_specs' in ph_group
            except AssertionError:
                raise (IOError, ('The lifetime flag is True but the TCSPC '
                                 'data is missing.'))

        for ph_data in ['detectors', 'nanotimes', 'particles']:
            if ph_data in ph_group:
                loader.load_data(ph_group, ph_data)

        if 'detectors_specs' in ph_group:
            det_specs = ph_group.detectors_specs
            if 'donor' in det_specs and 'acceptor' in det_specs:
                donor = det_specs.donor.read()
                accept = det_specs.acceptor.read()
                d.add(det_donor_accept=(donor, accept))

        if 'nanotimes_specs' in ph_group:
            nanot_specs = ph_group.nanotimes_specs
            nanotime_params = {}
            for name in ['tcspc_bin', 'tcspc_nbins', 'tcspc_range']:
                value = nanot_specs._f_get_child(name).read()
                nanotime_params.update(**{name: value})
            for name in ['tau_accept_only', 'tau_donor_only',
                         'tau_fret_trans']:
                if name in nanot_specs:
                    value = nanot_specs._f_get_child(name).read()
                    nanotime_params.update(**{name: value})
            d.add(nanotime_params=nanotime_params)
    else:
        raise NotImplemented

    d.add(data_file=data_file)
    return d



def hdf5(fname):
    """Load a data file saved in HDF5 smFRET format.

    Any :class:`fretbursts.burstlib.Data` object can be saved in HDF5 format
    using :func:`fretbursts.hdf5.store` .
    """
    if not os.path.isfile(fname):
        raise IOError, 'File not found.'
    data_file = tables.open_file(fname, mode = "r")
    file_format = ('smFRET_format_version', '0.1')
    if file_format[0] not in data_file.root._v_attrs:
        print "WARNING: Attribute '%s' not found." % file_format[0]
    else:
        assert file_format[1] == data_file.root._v_attrs[file_format[0]]

    # Default values for optional parameters
    params = dict(leakage=0., gamma=1.)

    # Load mandatory parameters
    for field in ('clk_p', 'nch', 'ALEX'):
        if not '/' + field in data_file:
            raise(ValueError, "Filed '%s' not found" % field)
        params[field] = data_file.get_node('/', name=field).read()

    # Load optional parameter (overwriting defaults)
    for field in ('leakage', 'gamma'):
        if '/' + field in data_file:
            params[field] = data_file.get_node('/', name=field).read()

    if not params['ALEX']:
        # NOT ALEX single and multi-spot
        params['ph_times_m'] = PyTablesList(data_file,
                                            group_name='timestamps',
                                            load_array=True)
        params['A_em'] = PyTablesList(data_file,
                                      group_name='acceptor_emission',
                                      load_array=True)
    elif params['nch'] == 1:
        # Single-spot ALEX
        ph_times_t = data_file.root.timestamps_t.read()
        det_t = data_file.root.detectors_t.read()
        params.update(ph_times_t=ph_times_t, det_t=det_t)
    else:
        # Multi-spot ALEX
        raise NotImplemented

    d = Data(fname=fname, **params)

    if '/particles' in data_file:
        par = PyTablesList(data_file, group_name='particles',
                           load_array=True)
        d.add(par=par)

    if '/nanotime' in data_file:
        nanotime = data_file.get_node('/nanotime').read()
        d.add(nanotime=nanotime)

    if '/nanotime_params' in data_file:
        nanotime_params = dict()
        nanot_group = data_file.get_node('/nanotime_params')
        for node in nanot_group._f_list_nodes():
            nanotime_params[node.name] = node.read()
        d.add(nanotime_params=nanotime_params)

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
    dx = Data(fname=fname, clk_p=12.5e-9, nch=8, leakage=leakage, gamma=gamma)
    ph_times_m, A_em, ph_times_det = load_data_ordered16(fname=fname,
            n_bytes_to_read=bytes_to_read, swap_D_A=swap_D_A)
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
        print 'WARNING: Big-FIFO full, flags saved in Data()'
        dx.add(big_fifo=big_fifo)
    if ch_fifo_full:
        print 'WARNING: CH-FIFO full, flags saved in Data()'
        dx.add(ch_fifo=ch_fifo)
    return dx

def multispot48(fname, leakage=0, gamma=1., reprocess=False,
                i_start=0, i_stop=None, debug=False):
    """Load a 48-ch multispot file and return a Data() object.
    """
    import tables
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
        ph_times_m, big_fifo, ch_fifo = process_store(timestamps, det,
                        out_fname=fname_h5, fifo_flag=True, debug=False)
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
        print 'WARNING: Big-FIFO full, flags saved in Data()'
        dx.add(big_fifo=big_fifo)
    if ch_fifo_full:
        print 'WARNING: CH-FIFO full, flags saved in Data()'
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

def usalex(fname, leakage=0, gamma=1., header=166, bytes_to_read=-1, BT=None):
    """Load a usALEX file and return a Data() object.

    This function returns a Data() object to which you need to apply
    an alternation selection before performing further analysis (background
    estimation, burst search, etc.).

    The pattern to load usALEX data is the following::

        d = loader.usalex(fname=fname)
        d.add(D_ON=(2850, 580), A_ON=(900, 2580), alex_period=4000)
        plot_alternation_hist(d)

    If the plot looks good apply the alternation with::

        loader.usalex_apply_period(d)

    Now `d` is ready for futher processing such as background estimation,
    burst search, etc...
    """
    if BT is not None:
        print 'WARNING: `BT` argument is deprecated, use `leakage` instead.'
        leakage = BT
    print " - Loading '%s' ... " % fname
    ph_times_t, det_t = load_sm(fname, header=header)
    print " [DONE]\n"

    DONOR_ON = (2850, 580)
    ACCEPT_ON = (930, 2580)
    alex_period = 4000

    dx = Data(fname=fname, clk_p=12.5e-9, nch=1, leakage=leakage, gamma=gamma,
              ALEX=True,
              D_ON=DONOR_ON, A_ON=ACCEPT_ON, alex_period=alex_period,
              ph_times_t=ph_times_t, det_t=det_t, det_donor_accept=(0, 1),
              )
    return dx

def usalex_apply_period(d, delete_ph_t=True, remove_d_em_a_ex=False):
    """Applies to the Data object `d` the alternation period previously set.

    Note that you first need to load the data with :func:`usalex` and then
    to set the alternation parameters using `d.add()`.

    The pattern to load usALEX data is the following::

        d = loader.usalex(fname=fname)
        d.add(D_ON=(2850, 580), A_ON=(900, 2580), alex_period=4000)
        plot_alternation_hist(d)

    If the plot looks good apply the alternation with::

        loader.usalex_apply_period(d)

    Now `d` is ready for futher processing such as background estimation,
    burst search, etc...
    """
    donor_ch, accept_ch  = d.det_donor_accept
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
    print "#donor: %d  #acceptor: %d \n" % (d_em.sum(), a_em.sum())

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

def nsalex(fname, leakage=0, gamma=1.):
    """Load a nsALEX file and return a Data() object.

    This function returns a Data() object to which you need to apply
    an alternation selection before performing further analysis (background
    estimation, burst search, etc.).

    The pattern to load nsALEX data is the following::

        d = loader.nsalex(fname=fname)
        d.add(D_ON=(2850, 580), A_ON=(900, 2580))
        nsalex_plot_alternation(d)

    If the plot looks good apply the alternation with::

        loader.nsalex_apply_period(d)

    Now `d` is ready for futher processing such as background estimation,
    burst search, etc...
    """
    ph_times_t, det_t, nanotime = load_spc(fname)

    DONOR_ON = (10, 1500)
    ACCEPT_ON = (2000, 3500)
    nanotime_fsr = 4095

    dx = Data(fname=fname, clk_p=50e-9, nch=1, ALEX=True, nsALEX=True,
              D_ON=DONOR_ON, A_ON=ACCEPT_ON,
              nanotime_fsr=nanotime_fsr,
              ph_times_t=ph_times_t, det_t=det_t, nanotime_t=nanotime,
              det_donor_accept=(4, 6),
              )
    return dx

def nsalex_apply_period(d, delete_ph_t=True):
    """Applies to the Data object `d` the alternation period previously set.

    Note that you first need to load the data with :func:`nsalex` and then
    to set the alternation parameters using `d.add()`.

    The pattern to load nsALEX data is the following::

        d = loader.nsalex(fname=fname)
        d.add(D_ON=(2850, 580), A_ON=(900, 2580))
        nsalex_plot_alternation(d)

    If the plot looks good apply the alternation with::

        loader.nsalex_apply_period(d)

    Now `d` is ready for futher processing such as background estimation,
    burst search, etc...
    """
    # Note: between boolean arrays * is equavilent to logical AND,
    #       and + is equivalent to logical OR.

    donor_ch, accept_ch  = d.det_donor_accept

    # Mask for donor + acceptor detectors (discard other detectors)
    d_ch_mask_t = (d.det_t == donor_ch)
    a_ch_mask_t = (d.det_t == accept_ch)
    da_ch_mask_t = d_ch_mask_t + a_ch_mask_t

    # Masks for excitation periods
    d_ex_mask_t = (d.nanotime_t > d.D_ON[0]) * (d.nanotime_t < d.D_ON[1])
    a_ex_mask_t = (d.nanotime_t > d.A_ON[0]) * (d.nanotime_t < d.A_ON[1])
    ex_mask_t = d_ex_mask_t + a_ex_mask_t  # Select only ph during Dex or Aex

    # Total mask: D+A photons, and only during the excitation periods
    mask = da_ch_mask_t * ex_mask_t  # logical AND

    # Apply selection to timestamps and nanotime
    ph_times = d.ph_times_t[mask]
    nanotime = d.nanotime_t[mask]

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

    d.add(ph_times_m=[ph_times], nanotime=nanotime,
          D_em=[d_em], A_em=[a_em], D_ex=[d_ex], A_ex=[a_ex],)

    if delete_ph_t:
        d.delete('ph_times_t')
        d.delete('det_t')
        d.delete('nanotime_t')

