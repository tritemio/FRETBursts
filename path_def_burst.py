"""
This module defines paths for common folders
"""

# WARNING: Always use trailing slash "/" after a folder name.

import os


# This allows different paths between linux and windows OS
if os.name == 'posix':
    # Linux or Mac
    data_dir = '../data/'
    GIT_PATH = 'git'  # On *nix assumes that git is in the PATH
elif os.name == 'nt':
    # Windows
    data_dir = '../../../data/'
    GIT_PATH = r'C:\Users\temp\AppData\Local\Atlassian\SourceTree\git_local\bin\git.exe'
else:
    raise OSError ("Operating system not recognized (%s)." % os.name)

alex_data_dir = data_dir+'/alex/'
nsalex_data_dir = data_dir+'/nsAlex/'

fig_dir = '../figure/'
log_dir = '../log/'
cache_dir = '../cache/'

# Return a path with only the last subfolder (i.e. date for measurements)
def shorten_fname(f): 
    return '/'.join(f.split('/')[-2:])

# Check that all the dir names end with '/'
for dir_name in [data_dir, alex_data_dir, fig_dir, log_dir, cache_dir]:
    assert dir_name.endswith('/')

