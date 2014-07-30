#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
This module contains a function to store :class:`fretbursts.burstlib.Data`
objects to disk in **HDF5-smFRET** format.

Utility functions to print the HDF5 file structure and data-attributes are
also provided.
"""

import os
import tables
from dataload.pytables_array_list import PyTablesList


_hdf5_smfret_meta = dict(
    clk_p = 'Duration in seconds of 1 timestamp unit (clock period).',
    nch = 'Number of confocal excitation spots',
    ALEX = 'If True the file contains ALternated EXcitation data.',

    # ALEX
    timestamps_t = 'Array of all the timestamps (before alternation selection)',
    detectors_t = 'Array of detector numbers for each timestamp',

    # Non-ALEX
    timestamps = ('Donor + Acceptor timestamp arrays. Contains a series of '
                  'timestamp arrays, one for each spot.'),
    acceptor_emission = ('Contains a list of bool arrays, one for each spot. '
                         'Each bool element indicates whether the '
                         'corresponding timestamp has been detected in the '
                         'acceptor channel (True) or in the donor channel '
                         '(False).'),

    # Lifetime
    nanotime = 'TCSPC photon arrival time (nanotime)',
    nanotime_params = 'TCSPC hardware and lifetime data parameters',

    # Simulation
    particles = ('Particle label (number) for each timestamp.'),

    # Optional parameters
    BT = ('Bleed-through or leakage coefficient (donor emission in the '
          'acceptor ch)'),
    gamma = 'Gamma factor',
)

_hdf5_map = {key: key for key in _hdf5_smfret_meta.keys()}
_hdf5_map['timestamps_t'] = 'ph_times_t'
_hdf5_map['detectors_t'] = 'det_t'
_hdf5_map['particles'] = 'par'


def store(d, compression=dict(complevel=6, complib='zlib')):
    """
    Saves the `Data` object `d` in an HDF5 file using pytables.

    The file name is obtained from d.fname, by replacing the extension
    with '.hdf5'.

    **HDF5-smFRET file structure**

    Attributes on the root node:
        - smFRET_format_title: a string description for the file format
        - smFRET_format_version: file-format version string ('0.1')

    Mandatory parameters:
        - /nch: number of confocal excitation spots
        - /clk_p: (float) duration in seconds of 1 timestamp unit (clock period)
        - /ALEX: (bool) If True the file contains ALternated EXcitation data

    Optional parameters:
        - /BT: Bleed-through (or leakage) coefficient
        - /gamma: gamma factor
        - /nanotime: TCSPC photon arrival time (nanotime)
        - /nanotime_params (group): TCSPC hardware and lifetime data parameters
        - /particles: particle label (number) for each timestamp.

    Saved timestamps if ALEX:
        - /timestamps_t: Array of all the timestamps (before alternation
          selection. Maps to `Data` attribute `ph_times_t`.
        - /detectors_t: Array of detector numbers for each timestamp.
          Maps to `Data` attribute `det_t`.

    Saved timestamps if NOT ALEX:
        - /timestamps: (group) Donor + Acceptor timestamp arrays. Contains a
          series of timestamp arrays, one for each spot.
          Arrays names: /timestamps/ts_0 , ts_1, etc. Maps to `Data` attribute
          `ph_times_m`.
        - /acceptor_emission: (group): Contains a list of bool arrays, one for
          each spot. Each bool element indicates whether the corresponding
          timestamp has been detected in the acceptor channel (True) or in
          the donor channel (False).
          Arrays names /acceptor_emssion/a_em_0, a_em1, etc. Maps to `Data`
          attribute `A_em`.

    Additional optional parameters (TODO):
        - /orig_data_file (string)

    """
    comp_filter = tables.Filters(**compression)

    basename, extension = os.path.splitext(d.fname)
    h5_fname = basename + '.hdf5'
    data_file = tables.open_file(h5_fname, mode = "w",
                                 title = "Confocal smFRET data")

    # Save the metadata for the file and the root node
    smFRET_format_title = ('HDF5-based file format for confocal single-'
                           'molecule FRET data')
    smFRET_format_version = '0.1'
    data_file.root._v_attrs.smFRET_format_title = smFRET_format_title
    data_file.root._v_attrs.smFRET_format_version = smFRET_format_version

    # Save the mandatory parameters
    mandatory_fields = ['clk_p', 'nch', 'ALEX']
    for field in mandatory_fields:
        data_file.create_array('/', field, obj=d[field],
                               title=_hdf5_smfret_meta[field])

    # Save the timestamps
    if d.ALEX:
        alex_fields = ['timestamps_t', 'detectors_t']

        for field in alex_fields:
            assert _hdf5_map[field] in d
            data_file.create_carray('/', field, obj=d[_hdf5_map[field]],
                                    title=_hdf5_smfret_meta[field],
                                    filters=comp_filter)
    else:
        # Non-ALEX case: save directly the list of timestamps per-ch
        d.ts_list = PyTablesList(data_file, group_name='timestamps',
                        group_descr=_hdf5_smfret_meta['timestamps'],
                        prefix='ts_', compression=compression)
        for ph in d.ph_times_m:
            d.ts_list.append(ph)

        d.a_em_list = PyTablesList(data_file, group_name='acceptor_emission',
                        group_descr=_hdf5_smfret_meta['acceptor_emission'],
                        prefix='a_em_', compression=compression)
        for aem in d.A_em:
            d.a_em_list.append(aem)

    # Optional parameters
    optional_fields = ['nanotime', 'BT', 'gamma']
    for field in optional_fields:
        if _hdf5_map[field] in d:
            data_file.create_array('/', field, obj=d[_hdf5_map[field]],
                                   title=_hdf5_smfret_meta[field])

    if 'par' in d:
        # Array of particles, used for simulated data
        d.par_list = PyTablesList(data_file, group_name='particles',
                                  group_descr=_hdf5_smfret_meta['particles'],
                                  prefix='par_', compression=compression)
        for par in d.par:
            d.par_list.append(par)

    if 'nanotime_params' in d:
        # Parameters for the TCSPC configuration
        data_file.create_group('/', 'nanotime_params',
                               title=_hdf5_smfret_meta['nanotime_params'])

        for key, val in d.nanotime_params.iteritems():
            if type(val) is tuple:
                obj, title = val
            else:
                obj, title = val, ''
            data_file.create_array('/nanotime_params', key, obj=obj,
                                   title=title)

    data_file.flush()
    d.add(data_file=data_file)


def print_attrs(data_file, node_name='/', which='user'):
    """Print the HDF5 attributes for `node_name`.

    Parameters:
        data_file (pytables HDF5 file object): the data file to print
        node_name (string): name of the path inside the file to be printed.
            Can be either a group or a leaf-node. Default: '/', the root node.
        which (string): Valid values are 'user' for user-defined attributes,
            'sys' for pytables-specific attributes and 'all' to print both
            groups of attributes. Default 'user'.
    """
    node = data_file.get_node(node_name)
    print 'List of attributes for:\n  %s\n' % node
    for attr in node._v_attrs._f_list():
        print '\t%s' % attr
        print "\t    %s" % repr(node._v_attrs[attr])

def print_children(data_file, group='/'):
    """Print all the sub-groups in `group` and leaf-nodes children of `group`.

    Parameters:
        data_file (pytables HDF5 file object): the data file to print
        group (string): path name of the group to be printed.
            Default: '/', the root node.
    """
    base = data_file.get_node(group)
    print 'Groups in:\n  %s\n' % base

    for node in base._f_walk_groups():
        if node is not base:
            print '    %s' % node

    print '\nLeaf-nodes in %s:' % group
    for node in base._v_leaves.itervalues():
        print '\t%s %s' % (node.name, node.shape)
        if len(node.title) > 0:
            print '\t    %s' % node.title
