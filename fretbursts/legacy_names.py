#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
 FRETBursts - A single-molecule FRET burst analysis toolkit.

 Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>

 Load this module as

     from legacy_names import *

 to load all the old names used before the pakage refactoring.
"""

__all__ = [
        # Library modules and functions
        "np", "r_", "zeros", "arange", "size", "SS",
        "plt", "rcParams", "plot", "hist", "xlabel", "ylabel", "grid",
        "title", "legend", "gca", "gcf", "xlim", "ylim",

        # Local modules
        "loader", "select_bursts", "bl", "bg", "bpl", "bext",

        # Classes, functions, variables
        "data_dir", "Data", "Sel", "Sel_mask", "Sel_mask_apply",
        "gamma_correct_E",
        "load_multispot8", # <- Deprecated function (only for comatibility)

        # Generic fit functions
        "gaussian_fit_hist",
        "gaussian_fit_cdf",
        "two_gaussian_fit_hist",
        "two_gaussian_fit_hist_min",
        "two_gaussian_fit_hist_min_ab",
        "two_gaussian_fit_EM",
        "two_gauss_mix_pdf",
        "two_gauss_mix_ab",

        # Standalone plots or plots as a function of ch
        "bsavefig", "mch_plot_bg", "plot_alternation_hist",

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
from numpy import r_, zeros, arange, size
import scipy.stats as SS
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.pyplot import (plot, hist, xlabel, ylabel, grid, title, legend,
                               gca, gcf)
from .path_def import data_dir
import background as bg
import burstlib as bl
from .burstlib import Data, Sel, Sel_mask, Sel_mask_apply, gamma_correct_E
import burstlib_ext as bext
from loader import load_multispot8
from fit.gaussian_fitting import (gaussian_fit_hist,
                                  gaussian_fit_cdf,
                                  two_gaussian_fit_hist,
                                  two_gaussian_fit_hist_min,
                                  two_gaussian_fit_hist_min_ab,
                                  two_gaussian_fit_EM,
                                  two_gauss_mix_pdf,
                                  two_gauss_mix_ab,)
import burst_plot as bpl
from burst_plot import (
        # Standalone plots as a function of ch
        bsavefig, mch_plot_bg, plot_alternation_hist,

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
