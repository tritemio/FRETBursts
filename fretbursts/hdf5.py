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


# Metadata for the HDF5 root node
_format_meta = dict(
    format_name = 'HDF5-Ph-Data',
    format_title = 'HDF5-based format for time-series of photon data.',
    format_version = '0.2'
    )

# Metadata for different fields (arrays) in the HDF5 format
_fields_meta = dict(
    # Global data
    timestamps_unit = 'Time in seconds of 1-unit increment in timestamps.',
    num_spots = 'Number of excitation or detection spots',
    alex = 'If True the file contains ALternated EXcitation data.',
    lifetime = 'If True the data contains nanotimes from TCSPC hardware',
    alternation_period = ('The duration of the excitation alternation using '
                          'the same units as the timestamps.'),
    nanotimes_unit = 'TCSPC time bin duration in seconds.',

    # Photon-data
    photon_data = ('Group containing arrays of photon-data (one element per '
                   'photon'),
    timestamps = 'Array of photon timestamps',
    detectors = 'Array of detector numbers for each timestamp',
    nanotimes = 'TCSPC photon arrival time (nanotimes)',
    particles = ('Particle label (number) for each timestamp.'),

    detectors_specs = 'Group for detector-specific data.',
    donor = 'Detectors for the donor spectral range',
    acceptor = 'Detectors for the acceptor spectral range',
    polariz_paral = 'Detectors for polarization parallel to excitation',
    polariz_perp = 'Detectors for polarization perpendicular to excitation',

    nanotimes_specs =  'Group for nanotime-specific data.',
    tcspc_bin = 'TCSPC time bin duration in seconds.',
    tcspc_nbins = 'Number of TCSPC bins.',
    tcspc_range = 'TCSPC full-scale range in seconds.',
)

hdf5_data_map = {key: key for key in _fields_meta.keys()}
hdf5_data_map.update(
            timestamps_unit = 'clk_p',
            num_spots = 'nch',
            alex = 'ALEX',
            #lifetime
            alternation_period = 'alex_period',
            #nanotimes_unit
            timestamps = 'ph_times_t',
            detectors = 'det_t',
            nanotimes = 'nanotime',
            particles = 'par',
            )

class H5Writer():

    def __init__(self, h5file, data, comp_filter):
        self.h5file = h5file
        self.data = data
        self.comp_filter = comp_filter

    def _add_data(self, where, name, func, obj=None, **kwargs):
        if obj is None:
            assert hdf5_data_map[name] in self.data
            obj = self.data[hdf5_data_map[name]]

        func(where, name, obj=obj,
             title=_fields_meta[name],
             **kwargs)

    def add_carray(self, where, name, obj=None):
        self._add_data(where, name, self.h5file.create_carray, obj=obj,
                       filters=self.comp_filter)

    def add_array(self, where, name, obj=None):
        self._add_data(where, name, self.h5file.create_array, obj=obj)

    def add_group(self, where, name, metakey=None):
        if metakey is None:
            metakey = name
        return self.h5file.create_group(where, name,
                                        title=_fields_meta[metakey])



def store(d, compression=dict(complevel=6, complib='zlib'), h5_fname=None):
    """
    Saves the `Data` object `d` in HDF5 format.

    The file name is `h5_fname` or generated from d.fname, by replacing
    the original extension with '.hdf5'.

    **HDF5-smFRET file structure**

    For complete specs of this specific HDF5 format see:
    https://github.com/tritemio/FRETBursts/wiki/HDF5-smFRET-format---Draft-0.2-(WIP)

    """
    comp_filter = tables.Filters(**compression)
    if 'lifetime' not in d:
        d.add(lifetime='nanotime' in d)

    if h5_fname is None:
        basename, extension = os.path.splitext(d.fname)
        h5_fname = basename + '.hdf5'
    data_file = tables.open_file(h5_fname, mode = "w",
                                 title = "Confocal smFRET data")
    writer = H5Writer(data_file, d, comp_filter)

    ## Save the root-node metadata
    for name, value in _format_meta.items():
        data_file.root._f_setattr(name, value)

    ## Save the mandatory parameters
    mandatory_fields = ['timestamps_unit', 'num_spots', 'alex', 'lifetime']
    for field in mandatory_fields:
        writer.add_array('/', field)

    if d.ALEX:
        writer.add_array('/', 'alternation_period')

    ## Save the photon-data
    if d.nch == 1:
        # Single-spot: using "basic layout"
        ph_group = writer.add_group('/', 'photon_data')

        for field in ['timestamps', 'detectors']:
            writer.add_carray(ph_group, field)

        # If present save detector_specs
        if 'det_donor_accept' in d:
            donor, accept = d.det_donor_accept
            det_group = writer.add_group(ph_group, 'detectors_specs')
            writer.add_array(det_group, 'donor', obj=donor)
            writer.add_array(det_group, 'acceptor', obj=accept)

        # If present save nanotime data
        if d.lifetime:
            writer.add_carray(ph_group, 'nanotime')
            nt_group = writer.add_group(ph_group, 'nanotime_specs')
            for spec in ['tcspc_bin', 'tcspc_nbins', 'tcspc_range']:
                writer.add(nt_group, spec, obj=d.nanotime_params[spec])

        if 'par' in d:
            writer.add_carray(ph_group, 'particles')

    else:
        # Multi-spot: using "multi-spot layout"
        for ich, ph in enumerate(d.iter_ph_times()):
            ch_group = writer.add_group('/', 'photon_data_%d' % ich,
                                        metakey='photon_data')

            writer.add_carray(ch_group, 'timestamps', obj=ph)

            # If A_em[ich] is a slice we have a single color so we don't
            # save the detector (there is only one detector per channel).
            a_em = d.A_em[ich]
            if type(a_em) is not slice:
                writer.add_carray(ch_group, 'detectors', obj=a_em)
                # Detector specs
                det_group = writer.add_group(ch_group, 'detectors_specs')
                writer.add_array(det_group, 'donor', obj=False)
                writer.add_array(det_group, 'acceptor', obj=True)

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
