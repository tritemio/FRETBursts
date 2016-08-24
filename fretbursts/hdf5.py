#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014-2016 The Regents of the University of California,
#               Antonino Ingargiola <tritemio@gmail.com>
#
"""
This module contains a function to store :class:`fretbursts.burstlib.Data`
objects to disk in **Photon-HDF5** format.

Utility functions to print the HDF5 file structure and data-attributes are
also provided.
"""

from __future__ import print_function, absolute_import
from builtins import range, zip

import phconvert as phc


hdf5_data_map = dict(
    filename='fname',
    timestamps_unit='clk_p',
    num_spots='nch',
    alex='ALEX',
    #lifetime
    #alex_period
    timestamps='ph_times_t',
    detectors='det_t',
    #nanotimes='nanotime',
    #particles='par',
    alex_period_donor='D_ON',
    alex_period_acceptor='A_ON',
)

hdf5_data_map_r = {v: k for k, v in hdf5_data_map.items()}

def store(d, compression=dict(complevel=6, complib='zlib'), h5_fname=None,
          verbose=True, num_spectral_ch=2):
    """
    Saves the `Data` object `d` in the Photon-HDF5 format.

    As a side effect the `d` object is modified by adding the attribute
    `data_file` that contains a reference to the pytables file.

    Arguments:
        d (Data object): the Data object containing the smFRET measurement.
        compression (dict): a dictionary containing the compression type
            and level. Passed to pytables `tables.Filters()`.
        h5_fname (string or None): if not None, contains the file name
            to be used for the HDF5 file. If None, the file name is generated
            from `d.fname`, by replacing the original extension with '.hdf5'.
        verbose (bool): if True prints the name of the saved file.

    For description and specs of the Photon-HDF5 format see:
    http://photon-hdf5.readthedocs.org/
    """
    print('DEPRECATED: This function saves the Photon-HDF5 0.2 format '
          'which is deprecated. \n            Please use '
          '`save_photon_hdf5()` in `phconvert` to save version >=0.3.')
    #comp_filter = tables.Filters(**compression)
    if 'lifetime' not in d:
        # Test on different fields for ALEX and non-ALEX
        d.add(lifetime=('nanotimes_t' in d) or ('nanotimes' in d))

    # Add default values for missing mandatory fields
    if 'num_spectral_ch' not in d:
        d.add(num_spectral_ch=num_spectral_ch)
    if 'num_polariz_ch' not in d:
        d.add(num_polariz_ch=1)
    if 'num_detectors' not in d:
        d.add(num_detectors=d.nch*d.num_spectral_ch)

    data = {hdf5_data_map_r.get(k, k): v for k, v in d.items()}

    if d.lifetime:
        data.update(d.nanotimes_params)

    if d.nch == 1:
        if d.ALEX:
            data['donor'], data['acceptor'] = d.det_donor_accept
            pass
        else:
            data['timestamps'] = data.pop('ph_times_m')[0]
            data['detectors'] = data.pop('A_em')[0]
            data['donor'], data['acceptor'] = 1, 0
            if 'par' in d:
                data['particles'] = data.pop('par')[0]
        phc.hdf5.photon_hdf5(data, compression=compression)
    else:
        if d.ALEX:
            raise NotImplementedError
        else:
            data['donor'], data['acceptor'] = 1, 0
            det = [None if isinstance(a_em, slice) else a_em.view('uint8')
                   for a_em in d.A_em]
            phc.hdf5.photon_hdf5(data, compression=compression,
                                 iter_timestamps=d.iter_ph_times(),
                                 iter_detectors=iter(det))
    d.add(data_file=data['data_file'])
