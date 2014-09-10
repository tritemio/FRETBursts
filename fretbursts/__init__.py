#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
 FRETBursts - A single-molecule FRET burst analysis toolkit.

 Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>

 Import this package as

     from fretbursts import *

 to load the common namespace used in the example notebooks.
"""

__version__ = '0.4.dev'

import warnings

try:
    import pandas
except ImportError:
    has_pandas = False
    warnings.warn((' - Cannot import pandas. Some functionality will not be '
                   'avalable.'))
else:
    has_pandas = True

try:
    import matplotlib
except ImportError:
    has_matplotlib = False
    warnings.warn((' - Cannot import matplotlib. Plotting will not be '
                   'avalable.'))
else:
    has_matplotlib = True

try:
    import lmfit
except ImportError:
    has_lmfit = False
    warnings.warn((' - Cannot import lmfit. Some fitting functionalities '
                   ' will not be avalable.'))
else:
    has_lmfit = True


__all__numpy = ["np", "r_", "zeros"]

__all__matplotlib = [
        # Library modules and functions
        "plt", "rcParams", "matplotlib", "plot", "hist",
        "grid", "xlim", "ylim", "gca", "gcf",]

__all_local_names = [
        # Local modules
        "loader", "select_bursts", "bl", "bg", "bpl", "bext", "bg_cache",
        "hdf5", "fretmath", "mfit",

        # Classes, functions, variables
        "Data", "Sel", "Sel_mask", "Sel_mask_apply", "gui_fname", "Ph_sel",
        "download_file",

        # Standalone plots or plots as a function of ch
        "mch_plot_bg", "plot_alternation_hist",

        # Plots types used for 1ch of multi-ch plots through `dplot`
        "timetrace", "timetrace_da", "ratetrace", "ratetrace_da",
        "timetrace_alex", "timetrace_fret",
        "timetrace_bg",
        "hist_width", "hist_size", "hist_size_all", "hist_fret",
        "hist2d_alex", "hist_S", "hist_sbr", "hist_asymmetry",
        "hist_bg_fit_single", "hist_bg_fit", "hist_ph_delays", "hist_mdelays",
        "hist_mrates", "hist_rate_in_burst", "hist_burst_delays",
        "scatter_width_size", "scatter_rate_da", "scatter_fret_size",
        "scatter_fret_nd_na", "scatter_fret_width", "scatter_da",
        "scatter_naa_nt", "scatter_alex",

        # Wrapper functions that create a plot for each channel
        "dplot", "dplot_48ch", "dplot_8ch", "dplot_1ch",
        ]

__all__ = __all__numpy + __all_local_names

if has_matplotlib:
    __all__ += __all__matplotlib

import numpy as np
from numpy import r_, zeros

if has_matplotlib:
    from matplotlib import rcParams
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import plot, hist, grid, xlim, ylim, gca, gcf
    import style

import background as bg
import burstlib as bl
from .burstlib import Data, Sel, Sel_mask, Sel_mask_apply
import burstlib_ext as bext
from .ph_sel import Ph_sel
import fretmath
import mfit

import burst_plot as bpl
from burst_plot import (
        # Standalone plots as a function of ch
        mch_plot_bg, plot_alternation_hist,

        # Plots types used for 1ch of multi-ch plots through `dplot`
        timetrace, timetrace_da, ratetrace, ratetrace_da,
        timetrace_alex, timetrace_fret,
        timetrace_bg,
        hist_width, hist_size, hist_size_all, hist_fret,
        hist2d_alex, hist_S, hist_sbr, hist_asymmetry,
        hist_bg_fit_single, hist_bg_fit, hist_ph_delays, hist_mdelays,
        hist_mrates, hist_rate_in_burst, hist_burst_delays,
        scatter_width_size, scatter_rate_da, scatter_fret_size,
        scatter_fret_nd_na, scatter_fret_width, scatter_da,
        scatter_naa_nt, scatter_alex,

        # Wrapper functions that create a plot for each channel
        dplot, dplot_48ch, dplot_8ch, dplot_1ch,
        )
from utils.gui import gui_fname
from utils.misc import download_file
