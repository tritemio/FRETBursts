#%%writefile load_fretbursts.py
"""
Helper script to load the FretBurst software.

This should be run from inside a notebook.
"""

# Modify these to point to your FretBursts folder
# or set an environment variable FRETBURSTS_DIR containing the path
# (the enviroment variable, if set, has the precedence).
FRETBURSTS_DIR_WIN = r"C:\Data\Antonio\software\src\fretbursts"
FRETBURSTS_DIR_LINUX = r"/home/user/src/fretbursts"


import os
from subprocess import check_output        
from IPython.display import display, Math, clear_output
from glob import glob                      # helps finding files


if 'FRETBURSTS_DIR' in os.environ:
    BURST_DIR = os.environ['FRETBURSTS_DIR']
elif os.name == 'posix':
    BURST_DIR = FRETBURSTS_DIR_LINUX
elif os.name == 'nt':
    BURST_DIR = FRETBURSTS_DIR_WIN
else:
    raise OSError ("Operating system not recognized (%s)." % os.name)

ip = get_ipython()

# Save current dir as NOTEBOOK_DIR
if not 'NOTEBOOK_DIR' in globals():
    NOTEBOOK_DIR = ip.magic('%pwd')

ip.magic('%matplotlib inline')

ip.magic('%cd "$BURST_DIR"')
ip.magic('%run -i fretbursts.py')
ip.magic('%run -i burst_plot.py')
ip.magic('%run -i style')

from utils.gui import gui_fname
from utils import git

# If git is available, check fretbursts version
if not git.git_path_valid():
    print('\nSoftware revision unknown (git not found).')
else:
    last_commit = git.get_last_commit_line()
    print('\nCurrent software revision:\n {}\n'.format(last_commit))
    if not git.check_clean_status():
        print('\nWARNING -> Uncommitted changes:')
        print(git.get_status())

