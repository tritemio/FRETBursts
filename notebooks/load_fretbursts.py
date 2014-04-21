#%%writefile load_fretbursts.py
"""
Helper script to load the FRETBursts from a IPython Notebook.

Run this file from a notebook as follows:

    run load_fretbursts

"""

import os
from IPython.display import display, Math, clear_output

HOME = os.environ['HOME'] if 'HOME' in os.environ else ''


# Modify these to point to your FRETBursts source folder
# or set an environment variable FRETBURSTS_DIR containing the path
# (the enviroment variable, if set, has the precedence).
FRETBURSTS_DIR_WIN = r"C:\Data\Antonio\software\src\fretbursts"
FRETBURSTS_DIR_POSIX = HOME + "/src/fretbursts"


if 'FRETBURSTS_DIR' in os.environ:
    FRETBURSTS_DIR = os.environ['FRETBURSTS_DIR']
elif os.name == 'posix':
    # Runnning Mac OSX or Linux
    FRETBURSTS_DIR = FRETBURSTS_DIR_POSIX
elif os.name == 'nt':
    # Running Windows
    FRETBURSTS_DIR = FRETBURSTS_DIR_WIN
else:
    raise OSError ("Operating system not recognized (%s)." % os.name)

ip = get_ipython()

# Save current dir as NOTEBOOK_DIR
if not 'NOTEBOOK_DIR' in globals():
    NOTEBOOK_DIR = ip.magic('%pwd')

ip.magic('%matplotlib inline')
ip.magic('%cd "$FRETBURSTS_DIR"')

from fretbursts import *
from fretbursts.utils.gui import gui_fname
from fretbursts.utils import git

git.print_summary('FRETBursts')
