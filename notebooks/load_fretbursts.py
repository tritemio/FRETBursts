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
ip.magic('%run -i fretbursts.py')
ip.magic('%run -i burst_plot.py')
ip.magic('%run -i style')

from utils.gui import gui_fname
from utils import git

# If git is available, check fretbursts version
if not git.git_path_valid():
    print('\nFRETBursts revision unknown (git not found).')
else:
    last_commit = git.get_last_commit_line()
    print('\nCurrent FRETBursts revision:\n {}\n'.format(last_commit))
    if not git.check_clean_status():
        print('\nWARNING -> Uncommitted changes:')
        print(git.get_status())

