#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
`load_fretbursts.py` is a small helper script used to import FRETBursts
in an ipython notebook when you want to load fretbursts from the sources
folder (not from system installation).

The script is typically ran at the beginning of a notebook as::

    %run load_fretbursts --nogui --source

See also http://fretbursts.readthedocs.org/en/latest/files_description.html#load-fretbursts-py
"""

import sys
import os

# Import some useful functions for the ipyton notebook
from IPython.display import display, Math, clear_output

# Process the command-line arguments
enable_qt_gui = not (len(sys.argv) > 1 and '--nogui' in sys.argv[1:])
enable_mpl = not (len(sys.argv) > 1 and '--nompl' in sys.argv[1:])
load_mpl_style = not (len(sys.argv) > 1 and '--nostyle' in sys.argv[1:])
load_from_source = len(sys.argv) > 1 and '--source' in sys.argv[1:]

## Enable inline plots and QT windows.
# Workaround for open-file dialog in IPython Notebook
# see https://github.com/ipython/ipython/issues/5798
ip = get_ipython()
if enable_mpl:
    ip.enable_matplotlib("inline")
if enable_qt_gui:
    ip.enable_gui("qt")

## Save current dir as NOTEBOOK_DIR
if not 'NOTEBOOK_DIR' in globals():
    NOTEBOOK_DIR = os.path.abspath('.')

try:
    if load_from_source:
        raise ImportError
    ## Try to import FRETBursts
    import fretbursts as fb
    system_install = True
except ImportError:
    ## Try FRETBursts in-place execution

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

    ## Change to FRETBursts folder and import the software
    os.chdir(FRETBURSTS_DIR)
    import fretbursts as fb
    system_install = False

finally:
    from fretbursts import *
    from fretbursts.utils import git

if load_mpl_style:
    import fretbursts.style

if not system_install:
    git.print_summary('FRETBursts')
    os.chdir(NOTEBOOK_DIR)
else:
    print("\nFRETBursts version: %s" % fb.__version__)

citation()



