# -*- coding: utf-8 -*-
"""
Main module to run the fretbursts software.
"""

from fretbursts_path_def import data_dir
#from burstlib import *
#from loaders import *

ip = get_ipython()
ip.magic("run -i burstlib.py")
ip.magic("run -i loaders.py")

#ip.magic("run -i burst_selection.py")
#ip.magic("run -i burstsearch.py")
#ip.magic("run -i background.py")