#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
Misc utility functions
"""

import os
import sys
import numpy as np


def clk_to_s(t_ck, clk_p=12.5*1e-9):
    """Convert clock cycles to seconds."""
    return t_ck*clk_p

def pprint(s, mute=False):
    """Print immediately, even if inside a busy loop."""
    if mute: return
    sys.stdout.write(s)
    sys.stdout.flush()

def deprecate(function, old_name, new_name):
    def deprecated_function(*args, **kwargs):
        pprint("Function name %s is deprecated, use %s instead.\n" %\
                (old_name, new_name))
        res = function(*args, **kwargs)
        return res
    return deprecated_function

def shorten_fname(f):
    """Return a path with only the last subfolder (i.e. measurement date)."""
    return '/'.join(f.split('/')[-2:])

def binning(times, bin_width_ms=1, max_num_bins=1e5, clk_p=12.5e-9):
    """Return the binned histogram of array times."""
    bin_width_clk = (bin_width_ms*1e-3)/clk_p
    num_bins = min(times.max()/bin_width_clk, max_num_bins)
    h = np.histogram(times[times<(num_bins*bin_width_clk)], bins=num_bins)
    return h


def mkdir_p(path):
    """Create the path if not existent, otherwise do nothing.
    If path exists and is not a dir raise an exception.
    """
    import errno
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def download(url, save_dir='./'):
    """Download a file from `url` streaming it to disk.

    The file name is taken from `url` and left unchanged.
    The destination dir can be set using `save_dir`
    (Default: current dir).
    """
    from IPython.display import clear_output
    import requests

    fname = url.split('/')[-1]
    path = '/'.join([os.path.abspath(save_dir), fname])
    mkdir_p(save_dir)

    r = requests.get(url, stream=True)
    chunk_size = 2**16
    downloaded_size = 0
    if r.status_code == 200:
        with open(path, 'wb') as f:
            for chunk in r.iter_content(chunk_size):
                f.write(chunk)
                downloaded_size += chunk_size/(2.**20)
                clear_output()
                print 'Saving %s ...' % path
                print 'Downloaded %5.1f MB' % downloaded_size
                sys.stdout.flush()
    else:
        print 'URL not found.'
