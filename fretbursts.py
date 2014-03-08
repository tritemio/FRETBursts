#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
This is the main file that loads the FRETBursts software, except for the plot 
functions.

This file imports two modules:

- `burstlib.py`
- `loaders.py`
"""

from fretbursts_path_def import data_dir
from burstlib import *
from loaders import *

#ip = get_ipython()
#ip.magic("run -i burstlib.py")
#ip.magic("run -i loaders.py")

#ip.magic("run -i burst_selection.py")
#ip.magic("run -i burstsearch.py")
#ip.magic("run -i background.py")