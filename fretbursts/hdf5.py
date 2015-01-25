#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
This module contains a function to store :class:`fretbursts.burstlib.Data`
objects to disk in **Photon-HDF5** format.

Utility functions to print the HDF5 file structure and data-attributes are
also provided.
"""

from __future__ import print_function, absolute_import
import os
import time
import tables
from collections import OrderedDict

from .utils.misc import pprint


# Metadata for the HDF5 root node
_root_attributes = OrderedDict([
    ('format_name', 'Photon-HDF5'),
    ('format_title', 'HDF5-based format for time-series of photon data.'),
    ('format_version', '0.2'),
    ('format_url', 'http://photon-hdf5.readthedocs.org/'),
])

# Metadata for different fields (arrays) in the HDF5 format
fields_descr = OrderedDict([
    # Root parameters
    ('num_spots', 'Number of excitation or detection spots.'),
    ('num_spectral_ch', ('Number of different spectral bands in the detection '
                         'channels (i.e. 2 for 2-colors smFRET).')),
    ('num_polariz_ch', ('Number of different polarization in the detection '
                        'channels. The value is 1 if no polarization selection '
                        'is performed and 2 if two orthogonal polarizations '
                        'are recorded.')),
    ('num_detectors', ('Total number of detector pixels used in the '
                       'measurement.')),
    ('measurement_duration', 'Measurement duration in seconds.'),
    ('lifetime', ('If True (or 1) the data contains nanotimes from TCSPC '
                  'hardware')),
    ('alex', 'If True (or 1) the file contains ALternated EXcitation data.'),
    ('alex_period', ('The duration of the us-ALEX excitation alternation '
                     'in the same units as the timestamps.')),
    ('laser_pulse_rate', 'The laser(s) pulse rate in Hertz.'),
    ('alex_period_donor', ('Start and stop values identifying the donor '
                           'emission period.')),
    ('alex_period_acceptor', ('Start and stop values identifying the acceptor '
                              'emission period.')),
    ('timestamps_unit', 'Time in seconds of 1-unit increment in timestamps.'),

    # Photon-data
    ('photon_data', ('Group containing arrays of photon-data (one element per '
                     'photon)')),
    ('timestamps', 'Array of photon timestamps'),

    ('detectors', 'Array of detector numbers for each timestamp'),
    ('detectors_specs', 'Group for detector-specific data.'),
    ('donor', 'Detectors for the donor spectral range'),
    ('acceptor', 'Detectors for the acceptor spectral range'),
    ('polarization1', ('Detectors ID for the "polarization1". By default is '
                       'the polarization parallel to the excitation, '
                       'unless specified differently in the "/setup_specs".')),
    ('polarization2', ('Detectors ID for the "polarization2". By default is '
                       'the polarization perpendicular to the excitation, '
                       'unless specified differently in the "/setup_specs".')),

    ('nanotimes', 'TCSPC photon arrival time (nanotimes)'),
    ('nanotimes_specs', 'Group for nanotime-specific data.'),
    ('tcspc_unit', 'TCSPC time bin duration in seconds (nanotimes unit).'),
    ('tcspc_num_bins', 'Number of TCSPC bins.'),
    ('tcspc_range', 'TCSPC full-scale range in seconds.'),
    ('tau_accept_only', 'Intrinsic Acceptor lifetime (seconds).'),
    ('tau_donor_only', 'Intrinsic Donor lifetime (seconds).'),
    ('tau_fret_donor', 'Donor lifetime in presence of Acceptor (seconds).'),
    ('inverse_fret_rate', ('FRET energy transfer lifetime (seconds). Inverse '
                           'of the rate of D*A -> DA*.')),

    ('particles', 'Particle label (integer) for each timestamp.'),

    ## Setup group
    ('setup', 'Information about the experimental setup.'),
    ('excitation_wavelengths', 'Array of excitation wavelengths (meters).'),
    ('excitation_powers', ('Array of excitation powers (in the same order as '
                           'excitation_wavelengths). Units: Watts.')),
    ('excitation_polarizations', ('Polarization angle (in degrees), one for '
                                  'each laser.')),
    ('detection_polarization1', ('Polarization angle (in degrees) for '
                                 '"polarization1".')),
    ('detection_polarization2', ('Polarization angle (in degrees) for '
                                 '"polarization2".')),

    ## Provenance group
    ('provenance', 'Information about the original data file.'),
    ('filename', 'Original file name.'),
    ('full_filename', 'Original full file name, including the folder.'),
    ('creation_time', 'Original file creation time.'),
    ('modification_time', 'Original file time of last modification.'),
    ])


hdf5_data_map = {key: key for key in fields_descr.keys()}
hdf5_data_map.update(
    timestamps_unit='clk_p',
    num_spots='nch',
    alex='ALEX',
    #lifetime
    #alex_period
    #nanotimes_unit
    timestamps='ph_times_t',
    detectors='det_t',
    #nanotimes='nanotime',
    #particles='par',
    alex_period_donor='D_ON',
    alex_period_acceptor='A_ON',
)

mandatory_root_fields = ['timestamps_unit', 'num_spots', 'num_detectors',
                         'num_spectral_ch', 'num_polariz_ch',
                         'alex', 'lifetime',]

setup_fields = ['excitation_wavelengths', 'excitation_powers',
                'excitation_polarizations', 'detection_polarization1',
                'detection_polarization2']

provenance_fields = ['filename', 'full_filename', 'creation_time',
                     'modification_time']

class H5Writer(object):
    """Helper class for writing fields of a Data() object in HDF5.

    It uses the global mapping and meta-data dictionaries defined above
    to retrive the field content and the description from the field name.
    """

    def __init__(self, h5file, data, comp_filter):
        self.h5file = h5file
        self.data = data
        self.comp_filter = comp_filter

    def _add_data(self, where, name, func, obj=None, **kwargs):
        if obj is None:
            assert hdf5_data_map[name] in self.data
            obj = self.data[hdf5_data_map[name]]

        func(where, name, obj=obj,
             title=fields_descr[name],
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
                                        title=fields_descr[metakey])


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
    comp_filter = tables.Filters(**compression)
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

    if h5_fname is None:
        basename, extension = os.path.splitext(d.fname)
        if compression['complib'] == 'blosc':
            basename += '_blosc'
        h5_fname = basename + '.hdf5'

    orig_file_metadata = {}
    if os.path.isfile(d.fname):
        orig_file_metadata = get_file_metadata(d.fname)
    else:
        pprint("WARNING: Could locate original file '%s'\n" % d.fname)
        pprint("         Provenance info not saved.\n")

    if os.path.isfile(h5_fname):
        basename, extension = os.path.splitext(h5_fname)
        h5_fname = basename + '_new_copy.hdf5'

    pprint('Saving: %s' % h5_fname, not verbose)
    data_file = tables.open_file(h5_fname, mode="w",
                                 title="Confocal smFRET data")
    d.add(data_file=data_file)
    writer = H5Writer(data_file, d, comp_filter)

    ## Save the root-node metadata
    for name, value in _root_attributes.items():
        data_file.root._f_setattr(name, value)

    ## Save the mandatory parameters
    for field in mandatory_root_fields:
        writer.add_array('/', field)

    if d.ALEX:
        if d.lifetime:
            writer.add_array('/', 'laser_pulse_rate')
        else:
            writer.add_array('/', 'alex_period')
        writer.add_array('/', 'alex_period_donor')
        writer.add_array('/', 'alex_period_acceptor')

    ## Add provenance metadata
    prov_group = writer.add_group('/', 'provenance')
    for field, value in orig_file_metadata.items():
        assert field in provenance_fields
        writer.add_array(prov_group, field, obj=value.encode('latin-1'))

    ## Add setup info, if present in d
    setup_group = writer.add_group('/', 'setup')
    for field in setup_fields:
        if field in d:
            writer.add_array(setup_group, field)

    ## Save the photon-data
    if d.nch == 1:
        # Single-spot: using "basic layout"
        ph_group = writer.add_group('/', 'photon_data')

        if d.ALEX:
            for field in ['timestamps', 'detectors']:
                writer.add_carray(ph_group, field)
            donor, accept = d.det_donor_accept
        else:
            writer.add_carray(ph_group, 'timestamps', obj=d.ph_times_m[0])
            writer.add_carray(ph_group, 'detectors',
                              obj=d.A_em[0].view(dtype='uint8'))
            donor, accept = 0, 1

        det_group = writer.add_group(ph_group, 'detectors_specs')
        writer.add_array(det_group, 'donor', obj=donor)
        writer.add_array(det_group, 'acceptor', obj=accept)

        # If present save nanotime data
        if d.lifetime:
            if d.ALEX:
                writer.add_carray(ph_group, 'nanotimes', obj=d.nanotimes_t)
            else:
                writer.add_carray(ph_group, 'nanotimes')
            nt_group = writer.add_group(ph_group, 'nanotimes_specs')

            # Mandatory specs
            nanotimes_specs = ['tcspc_bin', 'tcspc_nbins', 'tcspc_range']
            for spec in nanotimes_specs:
                writer.add_array(nt_group, spec, obj=d.nanotimes_params[spec])

            # Optional specs
            nanotimes_specs = ['tau_accept_only', 'tau_donor_only',
                               'tau_fret_donor', 'inverse_fret_rate']
            for spec in nanotimes_specs:
                if spec in d.nanotimes_params:
                    writer.add_array(nt_group, spec,
                                     obj=d.nanotimes_params[spec])

        if 'par' in d:
            writer.add_carray(ph_group, 'particles', obj=d.par[0])

    else:
        # Multi-spot: using "multi-spot layout"
        for ich, ph in enumerate(d.iter_ph_times()):
            if ph.size == 0: continue
            ch_group = writer.add_group('/', 'photon_data_%d' % ich,
                                        metakey='photon_data')

            writer.add_carray(ch_group, 'timestamps', obj=ph)

            # If A_em[ich] is a slice we have a single color so we don't
            # save the detector (there is only one detector per channel).
            a_em = d.A_em[ich]
            if not isinstance(a_em, slice):
                writer.add_carray(ch_group, 'detectors',
                                  obj=a_em.view(dtype='uint8'))
                # Detector specs
                det_group = writer.add_group(ch_group, 'detectors_specs')
                writer.add_array(det_group, 'donor', obj=False)
                writer.add_array(det_group, 'acceptor', obj=True)

    data_file.flush()


def get_file_metadata(fname):
    """Return a dict with file metadata.
    """
    assert os.path.isfile(fname)

    full_filename = os.path.abspath(fname)
    filename = os.path.basename(full_filename)

    # Creation and modification time (but not exactly on *NIX)
    # see https://docs.python.org/2/library/os.path.html#os.path.getctime)
    ctime = time.localtime(os.path.getctime(full_filename))
    mtime = time.localtime(os.path.getmtime(full_filename))

    ctime_str = time.strftime("%Y-%m-%d %H:%M:%S", ctime)
    mtime_str = time.strftime("%Y-%m-%d %H:%M:%S", mtime)

    metadata = dict(filename=filename, full_filename=full_filename,
                    creation_time=ctime_str, modification_time=mtime_str)
    return metadata


def dict_from_group(group):
    """Return a dict with the content of a PyTables group."""
    return {node.name: node.read() for node in group}

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
    print('List of attributes for:\n  %s\n' % node)
    for attr in node._v_attrs._f_list(which):
        print('\t%s' % attr)
        print('\t    %s' % repr(node._v_attrs[attr]))

def print_children(data_file, group='/'):
    """Print all the sub-groups in `group` and leaf-nodes children of `group`.

    Parameters:
        data_file (pytables HDF5 file object): the data file to print
        group (string): path name of the group to be printed.
            Default: '/', the root node.
    """
    base = data_file.get_node(group)
    print('Groups in:\n  %s\n' % base)

    for node in base._f_walk_groups():
        if node is not base:
            print('    %s' % node)

    print('\nLeaf-nodes in %s:' % group)
    for node in base._v_leaves.itervalues():
        info = node.shape
        if len(info) == 0:
            info = node.read()
        print('\t%s, %s' % (node.name, info))
        if len(node.title) > 0:
            print('\t    %s' % node.title)
