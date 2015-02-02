#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
This module implements a list of arrays stored into a file with pytables.

The list is created empty (if the file does not exist) and must be populated
with the `append()` method.

If the file-name exists the list is populated with arrays stored
in the file.

Each list element is a reference to a pytable array. To read the array in
memory use the slicing notation (like pytable_array[:]).
"""

from __future__ import print_function
from builtins import range, zip

import os
import tables

_default_compression = dict(complevel=6, complib='blosc')


class PyTablesList(list):
    def __init__(self, file, overwrite=False, parent_node='/',
                 group_name='array_list', group_descr='List of arrays',
                 prefix='data', compression=_default_compression,
                 load_array=False):
        """List of arrays stored in a pytables file.

        The list is inizialized empty and populated with `.append()`.

        Arguments:
            load_array (bool): if True, read the data and put numpy arrays
                in the list. If False, put only pytable arrays.

        `group_descr`, `prefix`, `compression` are only used if a new group is
        created (for example for a new file).
        """
        super(PyTablesList, self).__init__()
        self.parent_node = parent_node
        self.group_name = group_name
        self.load_array = load_array

        # Ignored if group exist
        self.size = 0
        self.prefix = prefix
        self.compression = compression

        ## Retrive the file reference file
        if type(file) is tables.file.File:
            self.data_file = file
        elif os.path.exists(file) and not overwrite:
            self.data_file = tables.open_file(file, mode = "a")
        else:
            self.data_file = tables.open_file(file, mode = "w",
                               title = "Container for lists of arrays")

        ## Create the group if not existent
        if group_name not in self.data_file.get_node(parent_node):
            self.data_file.create_group(parent_node, group_name,
                                        title=group_descr)
        self.group = self.data_file.get_node(parent_node, group_name)

        if 'size' in self.group._v_attrs:
            ## If the group was already present read the data
            self.size = self.group._v_attrs.size
            self.prefix = self.group._v_attrs.prefix
            for i in range(self.group._v_attrs.size):
                array_ = self.group._f_get_child(self.get_name(i))
                if self.load_array:
                    array_ = array_[:]
                super(PyTablesList, self).append(array_)
        else:
            ## If a new group save some metadata
            self.group._v_attrs.size = self.size
            self.group._v_attrs.prefix = self.prefix
            self.group._v_attrs.load_array = self.load_array

    def get_name(self, i=None):
        if i is None:
            i = self.size
        return self.prefix + str(i)

    def append(self, ndarray):
        name = self.get_name()
        comp_filter = tables.Filters(**self.compression)
        tarray = self.data_file.create_carray(self.group, name, obj=ndarray,
                                             filters=comp_filter)
        self.data_file.flush()
        super(PyTablesList, self).append(tarray)
        #print(self.prefix+str(self.size), ndarray)
        self.size += 1
        self.group._v_attrs.size = self.size

    def get_array_list(self):
        return [array_[:] for array_ in self]