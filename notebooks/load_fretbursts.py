#%%writefile load_fretbursts.py
"""
Helper script to load the FRETBursts from a IPython Notebook.

Run this file from a notebook as follows:

    %run load_fretbursts

"""

import sys
import os

# Import some useful functions for the ipyton notebook
from IPython.display import display, Math, clear_output

# Process the command-line arguments
no_gui = len(sys.argv) > 1 and sys.argv[1] == '--nogui'

## Find FRETBursts sources folder
if os.name == 'posix':
    # Linux or Mac
    HOME = os.environ['HOME'] + '/'
elif os.name == 'nt':
    # Windows
    HOME = os.environ['HOMEPATH'] + '/'
else:
    raise OSError ("Operating system not recognized (%s)." % os.name)

config_file_name = '.fretbursts'
with open(HOME + config_file_name) as f:
    FRETBURSTS_DIR = f.read().strip()

## Save current dir as NOTEBOOK_DIR
if not 'NOTEBOOK_DIR' in globals():
    NOTEBOOK_DIR = os.path.abspath('.')

## Change to FRETBursts folder and import the software
os.chdir(FRETBURSTS_DIR)

from fretbursts import *
from fretbursts.utils import git

git.print_summary('FRETBursts')
citation()

os.chdir(NOTEBOOK_DIR)

## Enable inline plots and QT windows.
# Workaround for open-file dialog in IPython Notebook
# see https://github.com/ipython/ipython/issues/5798
ip = get_ipython()
ip.enable_matplotlib("inline")
if not no_gui:
    ip.enable_gui("qt")

