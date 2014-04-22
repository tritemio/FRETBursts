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

__all__ = [
        # Library modules and functions
        "np", "r_", "zeros", "plt", "rcParams", "plot", "hist", "grid",
        "xlim", "ylim", "gca", "gcf",

        # Local modules
        "loader", "select_bursts", "bl", "bg", "bpl", "bext",

        # Classes, functions, variables
        "data_dir", "Data", "Sel", "Sel_mask", "Sel_mask_apply",

        # Standalone plots or plots as a function of ch
        "mch_plot_bg", "plot_alternation_hist",

        # Plots types used for 1ch of multi-ch plots through `dplot`
        "timetrace", "timetrace_da", "ratetrace", "ratetrace_da",
        "timetrace_alex", "timetrace_fret",
        "hist_width", "hist_size", "hist_fret", "kde_fret", "hist_fret_kde",
        "hist2d_alex", "hist_S", "hist_sbr",
        "hist_bg_fit_single", "hist_bg_fit", "hist_ph_delays", "hist_mdelays",
        "hist_mrates", "hist_rate_in_burst", "hist_burst_delays",
        "scatter_width_size", "scatter_rate_da", "scatter_fret_size",
        "scatter_fret_nd_na", "scatter_fret_width", "scatter_da",
        "scatter_naa_nt", "scatter_alex",

        # Wrapper functions that create a plot for each channel
        "dplot", "dplot_48ch", "dplot_8ch", "dplot_1ch",
        ]

import numpy as np
from numpy import r_, zeros
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, hist, grid, xlim, ylim, gca, gcf

from .path_def import data_dir
import background as bg
import burstlib as bl
from .burstlib import Data, Sel, Sel_mask, Sel_mask_apply
import burstlib_ext as bext

import burst_plot as bpl
from burst_plot import (
        # Standalone plots as a function of ch
        mch_plot_bg, plot_alternation_hist,

        # Plots types used for 1ch of multi-ch plots through `dplot`
        timetrace, timetrace_da, ratetrace, ratetrace_da,
        timetrace_alex, timetrace_fret,
        hist_width, hist_size, hist_fret, kde_fret, hist_fret_kde,
        hist2d_alex, hist_S, hist_sbr,
        hist_bg_fit_single, hist_bg_fit, hist_ph_delays, hist_mdelays,
        hist_mrates, hist_rate_in_burst, hist_burst_delays,
        scatter_width_size, scatter_rate_da, scatter_fret_size,
        scatter_fret_nd_na, scatter_fret_width, scatter_da,
        scatter_naa_nt, scatter_alex,

        # Wrapper functions that create a plot for each channel
        dplot, dplot_48ch, dplot_8ch, dplot_1ch,
        )

import style
