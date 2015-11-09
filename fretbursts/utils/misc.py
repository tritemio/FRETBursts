#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
Misc utility functions
"""

from __future__ import division, print_function
import os
import sys
import numpy as np


def _is_list_of_arrays(obj):
    return isinstance(obj, list) and np.all([isinstance(v, np.ndarray)
                                             for v in obj])

class HistData(object):
    """Stores histogram counts and bins and provides derived fields.

    Attributes:
        counts (array, ints): array of counts in each bin
        bins (array): array of bin edges. Size is size(counts) + 1.
        bincenters (array): array of bin  centers. Size is size(counts).
        pdf (array, floats): array of normalized counts (aka PDF)
    """
    def __init__(self, counts, bins):
        self.counts = counts
        self.bins = bins
        self.binwidth = bins[1] - bins[0]

    @property
    def bincenters(self):
        if not hasattr(self, '_bincenters'):
            self._bincenters = self.bins[:-1] + 0.5*self.binwidth
        return self._bincenters

    @property
    def pdf(self):
        if not hasattr(self, '_pdf'):
            self._pdf = np.array(self.counts, dtype=np.float)
            self._pdf /= (self.counts.sum() * self.binwidth)
        return self._pdf


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
    If `path` exists, and is not a dir, raise an exception.
    """
    import errno
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise exc

def download_file(url, save_dir='./'):
    """Download a file from `url` saving it to disk.

    The file name is taken from `url` and left unchanged.
    The destination dir can be set using `save_dir`
    (Default: the current dir).
    """
    ## Check if local path already exist
    fname = url.split('/')[-1]
    print('URL:  %s' % url)
    print('File: %s\n ' % fname)

    path = '/'.join([os.path.abspath(save_dir), fname])
    if os.path.exists(path):
        print('File already on disk: %s \nDelete it to re-download.' % path)
        return

    from future.standard_library import install_aliases
    install_aliases()
    from urllib.request import urlopen, urlretrieve
    from urllib.error import HTTPError, URLError

    ## Check if the URL is valid
    try:
        urlopen(url)
    except URLError as e:
        print('Wrong URL or no connection.\n\nError:\n%s\n' % e)
    except HTTPError:
        print('URL not found: ' + url)
        return

    ## Donwload the file
    def _report(blocknr, blocksize, size):
        current = blocknr*blocksize/2**20
        sys.stdout.write(
            "\rDownloaded {0:4.1f} / {1:4.1f} MB".format(current, size/2**20))
    mkdir_p(save_dir)
    urlretrieve(url, path, _report)

