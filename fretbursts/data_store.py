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


def store(d):
    """
    Compulsory parameters:
    /nch
    /clk_p
    /ALEX

    Optional parameters:
    /BT
    /gamma
    /nanotime
    /nanotime_params (group)

    Timestamps and detectors:
    ph_times_m -> '/timestamps/ts_0' , ...
    A_em       -> '/acceptor_emssion/a_em_0', ...

    Additional paramenters (TODO):
    /orig_data_file (1024 bytes string)

    """
    basename, extension = os.path.splitext(d.fname)
    h5_fname = basename + '.hdf5'
    data_file = tables.open_file(h5_fname, mode = "w",
                                 title = "Confocal smFRET data")

    data_file.create_array('/', 'clk_p', obj=d.clk_p,
                           title='Clock period for the timestamps')

    data_file.create_array('/', 'nch', obj=d.nch,
                           title='Number of smFRET excitation spots')

    data_file.create_array('/', 'ALEX', obj=d.ALEX,
                           title='True if ALternated EXcitation was used.')

    d.ts_list = PyTablesList(data_file, group_name='timestamps',
                    group_descr='Photon timestamp arrays', prefix='ts_')
    for ph in d.ph_times_m:
        d.ts_list.append(ph)

    d.a_em_list = PyTablesList(data_file, group_name='acceptor_emission',
                    group_descr='Boolean masks for acceptor emission',
                    prefix='a_em_')
    for aem in d.A_em:
        d.a_em_list.append(aem)

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
        d.par_list = PyTablesList(data_file, group_name='particles',
                    group_descr='Particle No for each emitted timestamp',
                    prefix='par_')
        for par in d.par:
            d.par_list.append(par)

    if 'nanotime_params' in d:
        data_file.create_group('/', 'nanotime_params',
                               title='Parameters of TCSPC and lifetime data')

        for key, val in d.nanotime_params.iteritems():
            if type(val) is tuple:
                obj, title = val
            else:
                obj, title = val, ''
            data_file.create_array('/nanotime_params', key, obj=obj,
                                   title=title)

    #TODO: save also fname
    data_file.flush()

    d.data_file = data_file

