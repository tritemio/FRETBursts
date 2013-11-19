"""
Helper script to load the FretBurst software.

This should be run from inside a notebook.
"""

#%%writefile load_fretburst.py

from subprocess import check_output
import os
if os.name == 'posix':
    BURST_DIR = r"/home/anto/Documents/ucla/src/burst/sources/smfretbursts"
elif os.name == 'nt':
    BURST_DIR = r"C:\Data\Antonio\software\src\fretburst"

ip = get_ipython()
if not 'NOTEBOOK_DIR' in globals():
    NOTEBOOK_DIR = ip.magic('%pwd')
    
from IPython.display import clear_output
from glob import glob

#ip.magic('%pylab inline')
ip.magic('%matplotlib inline')
ip.magic('%cd "$BURST_DIR"')

ip.magic('%run -i burstlib.py')
ip.magic('%run -i burst_plot.py')
ip.magic('%run -i style')
from utils.gui import gui_fname

if not git.git_path_valid():
    print('Software revision unknown (git not found).')
else:
    last_commit = git.get_last_commit()
    print('Current software revision: {}'.format(last_commit))
    if not git.check_clean_status():
        print('\nWARNING -> Uncommitted changes:')
        print(git.get_status())

