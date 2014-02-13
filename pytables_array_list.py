# -*- coding: utf-8 -*-
"""
This module implements a list of arrays stored to a file with pytables.

The list is created empty (if the file does not exist) and must be populated 
with the `append()` method.

TODO: If the file-name exist the list is populated with arrays stored 
in the file.

Each list element is a referent to a pytable array. To read the array in
memory use the slicing notation (like pytable_array[:]).
"""

import os
import tables
import numpy as np

_default_compression = dict(complevel=6, complib='blosc')

class PyTablesList(list):
    def __init__(self, fname, overwrite=False, parent_node='/', 
                 group_name='array_list', group_descr='List of arrays',  
                 prefix='data', compression=_default_compression):
        """List of arrays stored in pytables file.
        
        The list is inizialized empty and populated with `.append()`.
        """
        super(PyTablesList, self).__init__()
        self.prefix = prefix
        self.compression = compression
        self.size = 0
        
        if not overwrite and os.path.exists(fname):
            self.data_file = tables.open_file(fname, mode = "a")
        else:
            self.data_file = tables.open_file(fname, mode = "w",
                               title = "Container for lists of arrays")
        if group_name not in self.data_file.get_node(parent_node):
            # Create the group
            self.data_file.create_group(parent_node, group_name, 
                                        title=group_descr)
        self.group = self.data_file.get_node(parent_node, group_name)
        ##TODO: save prefix, compression and size inside the group
        
    def append(self, ndarray):
        name = self.prefix + str(self.size)
        comp_filter = tables.Filters(**self.compression)
        tarray = self.data_file.create_carray(self.group, name, obj=ndarray,
                                             filters=comp_filter)
        self.data_file.flush()
        super(PyTablesList, self).append(tarray)
        print self.prefix+str(self.size), ndarray
        self.size += 1
        