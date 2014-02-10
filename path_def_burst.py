"""
In this module you can define common paths (data folder, git path, etc...).

WARNING: Always use trailing slash "/" after a folder name.
"""

# Variable: DATA_DIR
# Usually prepended to file names when loading data.
# You can put this path in an environment variable FRETBURSTS_DATA_DIR
# that, if defined, has the precedence on the variable set here
_DATA_DIR_WIN = '../../../data/'
_DATA_DIR_LINUX = '../data/'

import os

if os.name == 'posix':
    # Linux or Mac
    _DEFAULT_GIT_PATH_LINUX = 'git'   # On *nix assumes that git is in the PATH
    data_dir = _DATA_DIR_LINUX
    GIT_PATH = _DEFAULT_GIT_PATH_LINUX
elif os.name == 'nt':
    # Windows
    _DEFAULT_GIT_PATH_WIN = (os.environ['homepath'] + \
                r'\AppData\Local\Atlassian\SourceTree\git_local\bin\git.exe')
    data_dir = _DATA_DIR_WIN
    GIT_PATH = _DEFAULT_GIT_PATH_WIN
else:
    raise OSError ("Operating system not recognized (%s)." % os.name)

# Variable: GIT_PATH
# Git is used to check the software revision during execution
# Here you can set the full path to a git executable
#GIT_PATH = r'C:\git.exe'


# Overwrites data_dir if FRETBURSTS_DATA_DIR is defined
if 'FRETBURSTS_DATA_DIR' in os.environ:
    data_dir = os.environ['FRETBURSTS_DATA_DIR']


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

