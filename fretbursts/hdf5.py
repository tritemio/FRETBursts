#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
Module to store a Data object to disk in HDF5 format.
"""

import os
import tables
from dataload.pytables_array_list import PyTablesList


def store(d, compression=dict(complevel=6, complib='blosc')):
    """
    Saves the `Data` object `d` in an HDF5 file using pytables.

    The file name is d.fname, replacing the extension with '.hdf5'.

    **HDF5 file structure**

    Attributes on the root node:
        - smFRET_format_title: a string description for the file format
        - smFRET_format_version: file-format version string ('0.1')

    Compulsory parameters:
        - /nch: number of excitation spots
        - /clk_p: clock period for the timestamps
        - /ALEX: (bool) 1 or True if ALEX, else 0 or False

    Optional parameters:
        - /BT
        - /gamma
        - /nanotime
        - /nanotime_params (group)

    Saved timestamps if ALEX:
        - ph_times_t -> /timestamps_t
        - det_t      -> /detectors_t

    Saved timestamps if NOT ALEX:
        - ph_times_m -> /timestamps/ts_0 , ts_1, ...
        - A_em       -> /acceptor_emssion/a_em_0, a_em1, ...

    Additional optional parameters (TODO):
        - /orig_data_file (string)

    """
    comp_filter = tables.Filters(**compression)

    basename, extension = os.path.splitext(d.fname)
    h5_fname = basename + '.hdf5'
    data_file = tables.open_file(h5_fname, mode = "w",
                                 title = "Confocal smFRET data")
    # Save the metadata
    smFRET_format_title = ('HDF5-based file format for confocal single-'
                           'molecule FRET data')
    smFRET_format_version = '0.1'
    data_file.root._v_attrs.smFRET_format_title = smFRET_format_title
    data_file.root._v_attrs.smFRET_format_version = smFRET_format_version

    # Save the compulsory parameters
    data_file.create_array('/', 'clk_p', obj=d.clk_p,
                           title='Clock period for the timestamps')

    data_file.create_array('/', 'nch', obj=d.nch,
                           title='Number of smFRET excitation spots')

    data_file.create_array('/', 'ALEX', obj=d.ALEX,
                           title='True if ALternated EXcitation was used.')

    # Save the timestamps
    if d.ALEX:
        assert 'ph_times_t' in d
        assert 'det_t' in d
        # ALEX case: save all the timestamps before alternation selection
        data_file.create_carray('/', 'timestamps_t', obj=d.ph_times_t,
                               title='Array of all the timestamps',
                               filters=comp_filter,
                               )
        data_file.create_carray('/', 'detectors_t', obj=d.det_t,
                               title=('Array of detector number foe each '
                                      'timestamps'),
                               filters=comp_filter,
                               )
    else:
        # Non-ALEX case: save directly the list of timestamps per-ch
        d.ts_list = PyTablesList(data_file, group_name='timestamps',
                        group_descr='Photon timestamp arrays', prefix='ts_',
                        compression=compression)
        for ph in d.ph_times_m:
            d.ts_list.append(ph)

        d.a_em_list = PyTablesList(data_file, group_name='acceptor_emission',
                        group_descr='Boolean masks for acceptor emission',
                        prefix='a_em_', compression=compression)
        for aem in d.A_em:
            d.a_em_list.append(aem)

    # Optional parameters
    if 'nanotime' in d:
        data_file.create_array('/', 'nanotime', obj=d.nanotime,
                               title=('Photon arrival time (with sub-ns '
                                      'resolution) respect to a laser sync.'))
    if 'BT' in d:
        data_file.create_array('/', 'BT', obj=d.BT,
                               title=('Fraction of donor emission detected by '
                                      'the acceptor channel'))
    if 'gamma' in d:
        data_file.create_array('/', 'gamma', obj=d.gamma,
                               title='smFRET Gamma-factor')

    if 'par' in d:
        # Array of particles, used for simulated data
        d.par_list = PyTablesList(data_file, group_name='particles',
                    group_descr='Particle No for each emitted timestamp',
                    prefix='par_')
        for par in d.par:
            d.par_list.append(par)

    if 'nanotime_params' in d:
        # Parameters for the TCSPC configuration
        data_file.create_group('/', 'nanotime_params',
                               title='Parameters of TCSPC and lifetime data')

        for key, val in d.nanotime_params.iteritems():
            if type(val) is tuple:
                obj, title = val
            else:
                obj, title = val, ''
            data_file.create_array('/nanotime_params', key, obj=obj,
                                   title=title)

    data_file.flush()
    d.add(data_file=data_file)

