#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
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
import tables

from .utils.misc import pprint, deprecate
from .burstlib import Data
from .hdf5 import hdf5_data_map

from phconvert.hdf5 import dict_from_group


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
        :class:`fretbursts.burstlib.Data` object containing the data.
    """
    print('WARNING: You are using Photon-HDF5 version 0.2 which is '
          'obsolete. Please update you file to version 0.3 or higher.')

    # This used for reading Photon-HDF5 version 0.2.
    mandatory_root_fields = ['timestamps_unit', 'num_spots', 'num_detectors',
                             'num_spectral_ch', 'num_polariz_ch',
                             'alex', 'lifetime']

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
                    assert not (a_em * d_em).any()
                    assert (a_em + d_em).all()

            if 'A_em' not in d:
                d.add(A_em=[a_em])
            else:
                d.A_em.append(a_em)

    d.add(data_file=data_file)
    return d
