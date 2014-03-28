#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
Misc utility functions
"""

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

def shorten_fname(f):
    """Return a path with only the last subfolder (i.e. measurement date)."""
    return '/'.join(f.split('/')[-2:])

def binning(times, bin_width_ms=1, max_num_bins=1e5, clk_p=12.5e-9):
    """Return the binned histogram of array times."""
    bin_width_clk = (bin_width_ms*1e-3)/clk_p
    num_bins = min(times.max()/bin_width_clk, max_num_bins)
    h = np.histogram(times[times<(num_bins*bin_width_clk)], bins=num_bins)
    return h
