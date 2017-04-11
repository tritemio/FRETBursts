#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014-2016 The Regents of the University of California,
#               Antonino Ingargiola <tritemio@gmail.com>
#

from __future__ import print_function, absolute_import

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


## Citation information
_CITATION = """
   FRETBursts: An Open Source Toolkit for Analysis of Freely-Diffusing Single-Molecule FRET
   Ingargiola et al. (2016). http://dx.doi.org/10.1371/journal.pone.0160716 """

_INFO_CITATION = (' You are running FRETBursts (version {}).\n\n'
                  ' If you use this software please cite the following'
                  ' paper:\n{}\n\n').format(__version__, _CITATION)

def citation(bar=True):
    cit = _INFO_CITATION
    if bar:
        cit = ('-' * 62) + '\n' + _INFO_CITATION +  ('-' * 62)
    print(cit)


import warnings

try:
    import pandas
except ImportError:
    has_pandas = False
    warnings.warn((' - Cannot import pandas. Some functionality will not be '
                   'available.'))
else:
    has_pandas = True

try:
    import matplotlib
except ImportError:
    has_matplotlib = False
    warnings.warn((' - Cannot import matplotlib. Plotting will not be '
                   'available.'))
else:
    has_matplotlib = True

try:
    try:
        from PyQt5 import QtWidgets, QtCore
        QtGui = QtWidgets
    except ImportError:
        try:
            from PyQt4 import QtGui, QtCore
        except ImportError:
            from PySide import QtGui, QtCore
except ImportError:
    has_qt = False
    # This catches ImportError or other errors due to broken QT installation
    warnings.warn((' - Cannot import QT, custom GUI widgets disabled.'))
else:
    has_qt = True


try:
    import lmfit
except ImportError:
    has_lmfit = False
    warnings.warn((' - Cannot import lmfit. Some fitting functionalities '
                   ' will not be available.'))
else:
    has_lmfit = True


__all__numpy = ["np"]

__all__matplotlib = [
        # Library modules and functions
        "plt", "rcParams", "matplotlib", "plot", "hist",
        "grid", "xlim", "ylim", "gca", "gcf",]

__all_local_names = [
        # Local modules
        "loader", "select_bursts", "bl", "bg", "bpl", "bext", "bg_cache",
        "hdf5", "fretmath", "mfit", "citation", "git",

        # Classes, functions, variables
        "Data", "Sel", "Ph_sel",
        "download_file", "init_notebook",

        # Standalone plots or plots as a function of ch
        "mch_plot_bg", "plot_alternation_hist", "alex_jointplot",

        # Plots types used for 1ch of multi-ch plots through `dplot`
        "timetrace", "timetrace_single", "ratetrace", "ratetrace_single",
        "timetrace_fret", "timetrace_bg",
        "hist_width", "hist_size", "hist_size_all", "hist_brightness",
        "hist_fret", "hist_burst_data",
        "hist2d_alex", "hist_S", "hist_sbr", "hist_asymmetry",
        "hist_interphoton_single", "hist_interphoton",
        "hist_bg_single", "hist_bg", "hist_ph_delays", "hist_mdelays",
        "hist_mrates", "hist_burst_phrate", "hist_burst_delays",
        "scatter_width_size", "scatter_rate_da", "scatter_fret_size",
        "scatter_fret_nd_na", "scatter_fret_width", "scatter_da",
        "scatter_naa_nt", "scatter_alex", "hexbin_alex",

        # Wrapper functions that create a plot for each channel
        "dplot", "dplot_48ch", "dplot_8ch", "dplot_1ch",
        ]

__all__ = __all__numpy + __all_local_names

import numpy as np

if has_qt:
    __all__ += ['OpenFileDialog']
    from .utils.gui import OpenFileDialog

if has_matplotlib:
    __all__ += __all__matplotlib
    from matplotlib import rcParams
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import plot, hist, grid, xlim, ylim, gca, gcf

# Import plain module names
from . import loader, hdf5, select_bursts, fretmath

# Import modules with custom names
from . import background as bg
from . import burstlib as bl

# Import objects
from .burstlib import Data, Sel
from .ph_sel import Ph_sel

if has_pandas and has_lmfit:
    from . import burstlib_ext as bext

if has_matplotlib and has_pandas and has_lmfit:
    from . import mfit
    from . import burst_plot as bpl
    from .burst_plot import (
            # Standalone plots as a function of ch
            mch_plot_bg, plot_alternation_hist, alex_jointplot,

            # Single-ch plots used in multi-ch plots through `dplot`
            timetrace, timetrace_single, ratetrace, ratetrace_single,
            timetrace_fret, timetrace_bg,
            hist_width, hist_size, hist_size_all, hist_brightness,
            hist_fret, hist_burst_data,
            hist2d_alex, hist_S, hist_sbr, hist_asymmetry,
            hist_interphoton_single, hist_interphoton,
            hist_bg_single, hist_bg, hist_ph_delays, hist_mdelays,
            hist_mrates, hist_burst_phrate, hist_burst_delays,
            scatter_width_size, scatter_rate_da, scatter_fret_size,
            scatter_fret_nd_na, scatter_fret_width, scatter_da,
            scatter_naa_nt, scatter_alex, hexbin_alex,

            # Wrapper functions that create a plot for each channel
            dplot, dplot_48ch, dplot_8ch, dplot_1ch,
            )

from .utils.misc import download_file
from .utils import git


def init_notebook(fs=13, savefig_dpi=65, seaborn_style='darkgrid',
                  mpl_backend='inline'):
    """
    Set the default plot style for inline plots using the seaborn library.

    This function must be called from an ipython notebook or
    ipython QT console.

    Arguments:
        fs (int): base font size for text labels (not for title)
        savefig_dpi (int): this value determines the figure size in
            the notebook. It is assigned to
            matplotlib.rcParams['savefig.dpi']

    Returns:
        The imported seaborn library. By saving the return value you
        don't need to import seaborn again.
    """
    if mpl_backend is not None:
        ip = get_ipython()
        ip.enable_matplotlib(mpl_backend)

    import seaborn as sns

    rc={'font.size': fs, 'axes.labelsize': fs, 'legend.fontsize': fs,
        'axes.titlesize': fs*1.1,
        'xtick.labelsize': fs, 'ytick.labelsize': fs,
        'savefig.dpi': savefig_dpi,
        'font.sans-serif': ['Arial', 'Liberation Sans'],
    }
    sns.set(rc=rc)
    blue = '#0055d4'
    green = '#2ca02c'
    color_brewer = sns.color_palette("Set1", 9)
    colors = np.array(color_brewer)[(1,0,2,3,4,8,6,7), :]
    colors = list(colors)
    colors[:3] = (blue, colors[1], green)
    sns.set_palette(colors, 8)
    sns.colors = colors
    sns.set_style(seaborn_style)
    return sns

citation()
