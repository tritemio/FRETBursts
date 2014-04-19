#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
In this module you can define a common paths for data folder
"""

# Variable: DATA_DIR
# Path (usually) used as base dir for file names when loading data.
# You can put this path in an environment variable FRETBURSTS_DATA_DIR
# that, if defined, has the precedence over the variables set here
_DATA_DIR_WIN = '../../../data/'
_DATA_DIR_LINUX = '../data/'

import os

if os.name == 'posix':
    # Linux or Mac
    data_dir = _DATA_DIR_LINUX
elif os.name == 'nt':
    # Windows
    data_dir = _DATA_DIR_WIN
else:
    raise OSError ("Operating system not recognized (%s)." % os.name)

# Overwrites data_dir if FRETBURSTS_DATA_DIR is defined
if 'FRETBURSTS_DATA_DIR' in os.environ:
    data_dir = os.environ['FRETBURSTS_DATA_DIR']

data_dir = os.path.abspath(data_dir) + '/'

# Check that all the dir names end with '/'
for dir_name in [data_dir]:
    assert dir_name.endswith('/')

