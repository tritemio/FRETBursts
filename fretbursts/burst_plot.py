# encoding: utf-8
#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2013-2016 The Regents of the University of California,
#               Antonino Ingargiola <tritemio@gmail.com>
#
"""
This module defines all the plotting functions for the
:class:`fretbursts.burstlib.Data` object.

The main plot function is `dplot()` that takes, as parameters, a `Data()`
object and a 1-ch-plot-function and creates a subplot for each channel.

The 1-ch plot functions are usually called through `dplot` but can also be
called directly to make a single channel plot.

The 1-ch plot functions names all start with the plot type (`timetrace`,
`ratetrace`, `hist` or `scatter`).

**Example 1** - Plot the timetrace for all ch::

    dplot(d, timetrace, scroll=True)

**Example 2** - Plot a FRET histogramm for each ch with a fit overlay::

    dplot(d, hist_fret, show_model=True)

For more examples refer to
`FRETBurst notebooks <http://nbviewer.ipython.org/github/tritemio/FRETBursts_notebooks/tree/master/notebooks/>`_.

"""

from __future__ import division, print_function, absolute_import
from builtins import range

import warnings

# Numeric imports
import numpy as np
from numpy import arange, r_
from matplotlib.mlab import normpdf
from scipy.stats import erlang
from scipy.interpolate import UnivariateSpline

# Graphics imports
import matplotlib.pyplot as plt
from matplotlib.pyplot import (plot, hist, xlabel, ylabel, grid, title, legend,
                               gca, gcf)
from matplotlib.patches import Rectangle, Ellipse
from matplotlib.collections import PatchCollection, PolyCollection
import seaborn as sns

# Local imports
from .ph_sel import Ph_sel
from . import burstlib as bl
from .phtools import phrates
from . import burstlib_ext as bext
from . import background as bg
from .utils.misc import pprint, HistData, _is_list_of_arrays
from .scroll_gui import ScrollingToolQT
from . import gui_selection as gs


params = {
    'font.size': 12,
    'legend.fontsize': 11,
}
plt.rcParams.update(params)


##
# Globals
#
blue = '#0055d4'
green = '#2ca02c'
red = '#e74c3c'  # '#E41A1C'
purple = '#9b59b6'

_ph_sel_color_dict = {Ph_sel('all'): blue, Ph_sel(Dex='Dem'): green,
                      Ph_sel(Dex='Aem'): red, Ph_sel(Aex='Aem'): purple,
                      Ph_sel(Aex='Dem'): 'c', }
_ph_sel_label_dict = {Ph_sel('all'): 'All-ph', Ph_sel(Dex='Dem'): 'DexDem',
                      Ph_sel(Dex='Aem'): 'DexAem', Ph_sel(Aex='Aem'): 'AexAem',
                      Ph_sel(Aex='Dem'): 'AexDem'}

# Global store for plot status
_plot_status = {}

# Global store for GUI handlers
gui_status = {'first_plot_in_figure': True}


##
#  Utility functions
#
def _normalize_kwargs(kwargs, kind='patch'):
    """Convert matplotlib keywords from short to long form."""
    if kwargs is None:
        return {}

    if kind == 'line2d':
        long_names = dict(c='color', ls='linestyle', lw='linewidth',
                          mec='markeredgecolor', mew='markeredgewidth',
                          mfc='markerfacecolor', ms='markersize',)
    elif kind == 'patch':
        long_names = dict(c='color', ls='linestyle', lw='linewidth',
                          ec='edgecolor', fc='facecolor',)
    for short_name in long_names:
        if short_name in kwargs:
            kwargs[long_names[short_name]] = kwargs.pop(short_name)
    return kwargs

def bsavefig(d, s):
    """Save current figure with name in `d`, appending the string `s`."""
    plt.savefig(d.Name() + s)

##
#  Multi-channel plot functions
#

def mch_plot_bg(d, **kwargs):
    """Plot background vs channel for DA, D and A photons."""
    bg = d.bg_from(Ph_sel('all'))
    bg_dd = d.bg_from(Ph_sel(Dex='Dem'))
    bg_ad = d.bg_from(Ph_sel(Dex='Aem'))
    plot(r_[1:d.nch+1], [b.mean()*1e-3 for b in bg], lw=2, color=blue,
         label=' T', **kwargs)
    plot(r_[1:d.nch+1], [b.mean()*1e-3 for b in bg_dd], color=green, lw=2,
         label=' D', **kwargs)
    plot(r_[1:d.nch+1], [b.mean()*1e-3 for b in bg_ad], color=red, lw=2,
         label=' A', **kwargs)
    xlabel("CH"); ylabel("kcps"); grid(True); legend(loc='best')
    title(d.name)

def mch_plot_bg_ratio(d):
    """Plot ratio of A over D background vs channel."""
    bg_dd = d.bg_from(Ph_sel(Dex='Dem'))
    bg_ad = d.bg_from(Ph_sel(Dex='Aem'))
    plot(r_[1:d.nch+1],
         [ba.mean()/bd.mean() for bd, ba in zip(bg_dd, bg_ad)],
         color=green, lw=2, label='A/D')
    xlabel("CH"); ylabel("BG Ratio A/D"); grid(True)
    title("BG Ratio A/D "+d.name)

def mch_plot_bsize(d):
    """Plot mean burst size vs channel."""
    CH = np.arange(1, d.nch+1)
    plot(CH, [b.mean() for b in d.nt], color=blue, lw=2, label=' T')
    plot(CH, [b.mean() for b in d.nd], color=green, lw=2, label=' D')
    plot(CH, [b.mean() for b in d.na], color=red, lw=2, label=' A')
    xlabel("CH"); ylabel("Mean burst size")
    grid(True)
    legend(loc='best')
    title(d.name)


##
#  ALEX alternation period plots
#
def plot_alternation_hist(d, bins=None, ax=None, **kwargs):
    """Plot the ALEX alternation histogram for the variable `d`.

    This function works both for us-ALEX and ns-ALEX data.

    This function must be called on ALEX data **before** calling
    :func:`fretbursts.loader.alex_apply_period`.
    """
    assert d.ALEX
    if d.lifetime:
        plot_alternation = plot_alternation_hist_nsalex
    else:
        plot_alternation = plot_alternation_hist_usalex
    plot_alternation(d, bins=bins, ax=ax, **kwargs)

def plot_alternation_hist_usalex(d, bins=None, ax=None, ich=0,
                                 hist_style={}, span_style={}):
    """Plot the us-ALEX alternation histogram for the variable `d`.

    This function must be called on us-ALEX data **before** calling
    :func:`fretbursts.loader.alex_apply_period`.
    """
    if ax is None:
        _, ax = plt.subplots()

    if bins is None:
        bins = 100

    D_ON, A_ON = d._D_ON_multich[ich], d._A_ON_multich[ich]
    d_ch, a_ch = d._det_donor_accept_multich[ich]
    offset = d.get('offset', 0)
    ph_times_t, det_t, period = d.ph_times_t[ich], d.det_t[ich], d.alex_period
    d_em_t = (det_t == d_ch)
    hist_style_ = dict(bins=bins, histtype='step', lw=2, alpha=0.9, zorder=2)
    hist_style_.update(hist_style)

    span_style_ = dict(alpha=0.2, zorder=1)
    span_style_.update(span_style)

    D_label = 'Donor: %d-%d' % (D_ON[0], D_ON[1])
    A_label = 'Accept: %d-%d' % (A_ON[0], A_ON[1])

    ax.hist((ph_times_t[d_em_t] - offset) % period, color=green, label=D_label,
            **hist_style_)
    ax.hist((ph_times_t[~d_em_t] - offset) % period, color=red, label=A_label,
            **hist_style_)
    ax.set_xlabel('Timestamp MODULO Alternation period')

    if D_ON[0] < D_ON[1]:
        ax.axvspan(D_ON[0], D_ON[1], color=green, **span_style_)
    else:
        ax.axvspan(0, D_ON[1], color=green, **span_style_)
        ax.axvspan(D_ON[0], period, color=green, **span_style_)

    if A_ON[0] < A_ON[1]:
        ax.axvspan(A_ON[0], A_ON[1], color=red, **span_style_)
    else:
        ax.axvspan(0, A_ON[1], color=red, **span_style_)
        ax.axvspan(A_ON[0], period, color=red, **span_style_)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)

def plot_alternation_hist_nsalex(d, bins=None, ax=None, ich=0,
                                 hist_style={}, span_style={}):
    """Plot the ns-ALEX alternation histogram for the variable `d`.

    This function must be called on ns-ALEX data **before** calling
    :func:`fretbursts.loader.alex_apply_period`.
    """
    if ax is None:
        _, ax = plt.subplots()

    if bins is None:
        bins = np.arange(d.nanotimes_params[ich]['tcspc_num_bins'])

    D_ON_multi, A_ON_multi = d._D_ON_multich[ich], d._A_ON_multich[ich]
    D_ON = [(D_ON_multi[i], D_ON_multi[i+1]) for i in range(0, len(D_ON_multi), 2)]
    A_ON = [(A_ON_multi[i], A_ON_multi[i+1]) for i in range(0, len(A_ON_multi), 2)]

    d_ch, a_ch = d._det_donor_accept_multich[ich]
    hist_style_ = dict(bins=bins, histtype='step', lw=1.3, alpha=0.9, zorder=2)
    hist_style_.update(hist_style)

    span_style_ = dict(alpha=0.2, zorder=1)
    span_style_.update(span_style)

    D_label = 'Donor: '
    for d_on in D_ON:
        D_label += '%d-%d' % (d_on[0], d_on[1])
    A_label = 'Accept: '
    for a_on in A_ON:
        A_label += '%d-%d' % (a_on[0], a_on[1])

    nanotimes_d = d.nanotimes_t[ich][d.det_t[ich] == d_ch]
    nanotimes_a = d.nanotimes_t[ich][d.det_t[ich] == a_ch]

    ax.hist(nanotimes_d, label=D_label, color=green, **hist_style_)
    ax.hist(nanotimes_a, label=A_label, color=red, **hist_style_)
    ax.set_xlabel('Nanotime bin')
    ax.set_yscale('log')
    for d_on in D_ON:
        ax.axvspan(d_on[0], d_on[1], color=green, **span_style_)
    for a_on in A_ON:
        ax.axvspan(a_on[0], a_on[1], color=red, **span_style_)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##  Multi-channel plots
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

##
#  Timetrace plots
#
def _plot_bursts(d, i, t_max_clk, pmax=1e3, pmin=0, color="#999999"):
    """Highlights bursts in a timetrace plot."""
    b = d.mburst[i]
    if b.num_bursts == 0:
        return
    pprint("CH %d burst..." % (i + 1))
    bs = b[b.start < t_max_clk]
    start = bs.start * d.clk_p
    end = bs.stop * d.clk_p
    R = []
    width = end - start
    ax = gca()
    for s, w in zip(start, width):
        r = Rectangle(xy=(s, pmin), height=pmax - pmin, width=w)
        r.set_clip_box(ax.bbox)
        r.set_zorder(0)
        R.append(r)
    ax.add_artist(PatchCollection(R, lw=0, color=color))
    pprint("[DONE]\n")

def _plot_rate_th(d, i, F, ph_sel, invert=False, scale=1,
                  plot_style_={}, rate_th_style={}):
    """Plots background_rate*F as a function of time.

    `plot_style_` is the style of a timetrace/ratetrace plot used as starting
    style. Linestyle and label are changed. Finally, `rate_th_style` is
    applied and can override any style property.

    If rate_th_style_['label'] is 'auto' the label is generated from
    plot_style_['label'] and F.
    """
    if F is None:
        F = d.F if F in d else 6

    rate_th_style_ = dict(plot_style_)
    rate_th_style_.update(linestyle='--', label='auto')
    rate_th_style_.update(_normalize_kwargs(rate_th_style, kind='line2d'))
    if rate_th_style_['label'] is 'auto':
        rate_th_style_['label'] = 'bg_rate*%d %s' % \
                                  (F, plot_style_['label'])
    x_rate = np.hstack(d.Ph_p[i]) * d.clk_p
    y_rate = F * np.hstack([(rate, rate) for rate in d.bg_from(ph_sel)[i]])
    y_rate *= scale
    if invert:
        y_rate *= -1
    plot(x_rate, y_rate, **rate_th_style_)

def _gui_timetrace_burst_sel(d, fig, ax):
    """Add GUI burst selector via mouse click to the current plot."""
    global gui_status
    if gui_status['first_plot_in_figure']:
        gui_status['burst_sel'] = gs.MultiAxPointSelection(fig, ax, d)
    else:
        gui_status['burst_sel'].ax_list.append(ax)

def _gui_timetrace_scroll(fig):
    """Add GUI to scroll a timetrace wi a slider."""
    global gui_status
    if gui_status['first_plot_in_figure']:
        gui_status['scroll_gui'] = ScrollingToolQT(fig)

def timetrace_single(d, i=0, binwidth=1e-3, bins=None, tmin=0, tmax=200,
                     ph_sel=Ph_sel('all'), invert=False, bursts=False,
                     burst_picker=True, scroll=False, cache_bins=True,
                     plot_style=None, show_rate_th=True, F=None,
                     rate_th_style={}, set_ax_limits=True,
                     burst_color='#BBBBBB'):
    """Plot the timetrace (histogram) of timestamps for a photon selection.

    See :func:`timetrace` to plot multiple photon selections (i.e.
    Donor and Acceptor photons) in one step.
    """
    if tmax is None or tmax < 0 or tmax > d.time_max:
        tmax = d.time_max

    def _get_cache():
        return (timetrace_single.bins, timetrace_single.x,
                timetrace_single.binwidth,
                timetrace_single.tmin, timetrace_single.tmax)

    def _set_cache(bins, x, binwidth, tmin, tmax):
        cache = dict(bins=bins, x=x, binwidth=binwidth, tmin=tmin, tmax=tmax)
        for name, value in cache.items():
            setattr(timetrace_single, name, value)

    def _del_cache():
        names = ['bins', 'x', 'binwidth', 'tmin', 'tmax']
        for name in names:
            delattr(timetrace_single, name)

    def _has_cache():
        return hasattr(timetrace_single, 'bins')

    def _has_cache_for(binwidth, tmin, tmax):
        if _has_cache():
            return (binwidth, tmin, tmax) == _get_cache()[2:]
        return False

    # If cache_bins is False delete any previously saved attribute
    if not cache_bins and _has_cache:
        _del_cache()

    tmin_clk, tmax_clk = tmin / d.clk_p, tmax / d.clk_p
    binwidth_clk = binwidth / d.clk_p

    # If bins is not passed try to use the cached one
    if bins is None:
        if cache_bins and _has_cache_for(binwidth, tmin, tmax):
            bins, x = timetrace_single.bins, timetrace_single.x
        else:
            bins = np.arange(tmin_clk, tmax_clk + 1, binwidth_clk)
            x = bins[:-1] * d.clk_p + 0.5 * binwidth
            if cache_bins:
                _set_cache(bins, x, binwidth, tmin, tmax)

    # Compute histogram
    ph_times = d.get_ph_times(i, ph_sel=ph_sel)
    timetrace, _ = np.histogram(ph_times, bins=bins)
    if invert:
        timetrace *= -1

    # Plot bursts
    if bursts:
        t_max_clk = int(tmax / d.clk_p)
        print(burst_color)
        _plot_bursts(d, i, t_max_clk, pmax=500, pmin=-500, color=burst_color)

    # Plot timetrace
    plot_style_ = dict(linestyle='-', linewidth=1.2, marker=None)
    if ph_sel in _ph_sel_color_dict:
        plot_style_['color'] = _ph_sel_color_dict[ph_sel]
        plot_style_['label'] = _ph_sel_label_dict[ph_sel]
    else:
        plot_style_['label'] = str(ph_sel)
    plot_style_.update(_normalize_kwargs(plot_style, kind='line2d'))
    plot(x, timetrace, **plot_style_)

    # Plot burst-search rate-threshold
    if show_rate_th and 'bg' in d:
        _plot_rate_th(d, i, F=F, ph_sel=ph_sel, invert=invert,
                      scale=binwidth, plot_style_=plot_style_,
                      rate_th_style=rate_th_style)

    plt.xlabel('Time (s)')
    plt.ylabel('# ph')
    if burst_picker and 'mburst' in d:
        _gui_timetrace_burst_sel(d, gcf(), gca())
    if scroll:
        _gui_timetrace_scroll(gcf())

    if set_ax_limits:
        plt.xlim(tmin, tmin + 1)
        if not invert:
            plt.ylim(ymax=100)
        else:
            plt.ylim(ymin=-100)
        _plot_status['timetrace_single'] = {'autoscale': False}


def timetrace(d, i=0, binwidth=1e-3, bins=None, tmin=0, tmax=200,
              bursts=False, burst_picker=True, scroll=False,
              show_rate_th=True, F=None, rate_th_style={'label': None},
              show_aa=True, legend=False, set_ax_limits=True,
              burst_color='#BBBBBB', plot_style=None,
              #dd_plot_style={}, ad_plot_style={}, aa_plot_style={}
             ):
    """Plot the timetraces (histogram) of photon timestamps.

    Arguments:
        d (Data object): the measurement's data to plot.
        i (int): the channel to plot. Default 0.
        binwidth (float): the bin width (seconds) of the timetrace histogram.
        bins (array or None): If not None, defines the bin edges for the
            timetrace (overriding `binwidth`). If None, `binwidth` is use
            to generate uniform bins.
        tmin, tmax (float): min and max time (seconds) to include in the
            timetrace. Note that a long time range and a small `binwidth`
            can require a significant amount of memory.
        bursts (bool): if True, plot the burst start-stop times.
        burst_picker (bool): if True, enable the ability to click on bursts
            to obtain burst info. This function requires the matplotlib's QT
            backend.
        scroll (bool): if True, activate a scrolling bar to quickly scroll
            through the timetrace. This function requires the matplotlib's QT
            backend.
        show_rate_th (bool): if True, plot the burst search threshold rate.
        F (bool): if `show_rate` is True, show a rate `F` times larger
            than the background rate.
        rate_th_style (dict): matplotlib style for the rate line.
        show_aa (bool): if True, plot a timetrace for the AexAem photons.
            If False (default), plot timetraces only for DexDem and DexAem
            streams.
        legend (bool): whether to show the legend or not.
        set_ax_limits (bool): if True, set the xlim to zoom on a small portion
            of timetrace. If False, do not set the xlim, display the full
            timetrace.
        burst_color (string): string containing the the HEX RGB color to use
            to highlight the burst regions.
        plot_style (dict): matplotlib's style for the timetrace lines.
    """
    # Plot bursts
    if bursts:
        t_max_clk = int(tmax / d.clk_p)
        _plot_bursts(d, i, t_max_clk, pmax=500, pmin=-500, color=burst_color)

    # Plot multiple timetraces
    ph_sel_list = [Ph_sel(Dex='Dem'), Ph_sel(Dex='Aem')]
    invert_list = [False, True]
    burst_picker_list = [burst_picker, False]
    scroll_list = [scroll, False]
    if d.ALEX and show_aa:
        ph_sel_list.append(Ph_sel(Aex='Aem'))
        invert_list.append(True)
        burst_picker_list.append(False)
        scroll_list.append(False)

    for ix, (ph_sel, invert) in enumerate(zip(ph_sel_list, invert_list)):
        if not bl.mask_empty(d.get_ph_mask(i, ph_sel=ph_sel)):
            timetrace_single(
                d, i, binwidth=binwidth, bins=bins, tmin=tmin,
                tmax=tmax, ph_sel=ph_sel, invert=invert, bursts=False,
                burst_picker=burst_picker_list[ix],
                scroll=scroll_list[ix], cache_bins=True,
                show_rate_th=show_rate_th, F=F,
                rate_th_style=rate_th_style, set_ax_limits=set_ax_limits,
                plot_style=plot_style)
    if legend:
        plt.legend(loc='best', fancybox=True)


def ratetrace_single(d, i=0, m=None, max_num_ph=1e6, tmin=0, tmax=200,
                     ph_sel=Ph_sel('all'), invert=False, bursts=False,
                     burst_picker=True, scroll=False, plot_style={},
                     show_rate_th=True, F=None, rate_th_style={},
                     set_ax_limits=True, burst_color='#BBBBBB'):
    """Plot the ratetrace of timestamps for a photon selection.

    See :func:`ratetrace` to plot multiple photon selections (i.e.
    Donor and Acceptor photons) in one step.
    """
    if tmax is None or tmax < 0:
        tmax = d.time_max

    if m is None:
        m = d.m if m in d else 10

    # Plot bursts
    if bursts:
        t_max_clk = int(tmax/d.clk_p)
        _plot_bursts(d, i, t_max_clk, pmax=500, pmin=-500, color=burst_color)

    # Compute ratetrace
    tmin_clk, tmax_clk = tmin/d.clk_p, tmax/d.clk_p
    ph_times = d.get_ph_times(i, ph_sel=ph_sel)
    iph1 = np.searchsorted(ph_times, tmin_clk)
    iph2 = np.searchsorted(ph_times, tmax_clk)
    if iph2 - iph1 > max_num_ph:
        iph2 = iph1 + max_num_ph
        tmax = ph_times[iph2]*d.clk_p
        warnings.warn(('Max number of photons reached in ratetrace_single().'
                       '\n    tmax is reduced to %ds. To plot a wider '
                       'time range increase `max_num_ph`.') % tmax,
                      UserWarning)
    ph_times = ph_times[iph1:iph2]
    rates = 1e-3*phrates.mtuple_rates(ph_times, m)/d.clk_p
    if invert:
        rates *= -1
    times = phrates.mtuple_rates_t(ph_times, m)*d.clk_p

    # Plot ratetrace
    plot_style_ = dict(linestyle='-', linewidth=1.2, marker=None)
    if ph_sel in _ph_sel_color_dict:
        plot_style_['color'] = _ph_sel_color_dict[ph_sel]
        plot_style_['label'] = _ph_sel_label_dict[ph_sel]
    plot_style_.update(_normalize_kwargs(plot_style, kind='line2d'))
    plot(times, rates, **plot_style_)

    # Plot burst-search rate-threshold
    if show_rate_th and 'bg' in d:
        _plot_rate_th(d, i, F=F, scale=1e-3, ph_sel=ph_sel, invert=invert,
                      plot_style_=plot_style_, rate_th_style=rate_th_style)

    xlabel('Time (s)'); ylabel('Rate (kcps)')
    if burst_picker:
        _gui_timetrace_burst_sel(d, gcf(), gca())
    if scroll:
        _gui_timetrace_scroll(gcf())

    if set_ax_limits:
        plt.xlim(tmin, tmin + 1)
        if not invert:
            plt.ylim(ymax=100)
        else:
            plt.ylim(ymin=-100)
        _plot_status['ratetrace_single'] = {'autoscale': False}


def ratetrace(d, i=0, m=None, max_num_ph=1e6, tmin=0, tmax=200,
              bursts=False, burst_picker=True, scroll=False,
              show_rate_th=True, F=None, rate_th_style={'label': None},
              show_aa=True, legend=False, set_ax_limits=True,
              #dd_plot_style={}, ad_plot_style={}, aa_plot_style={}
              burst_color='#BBBBBB'):
    """Plot the rate timetraces of photon timestamps.

    Arguments:
        d (Data object): the measurement's data to plot.
        i (int): the channel to plot. Default 0.
        max_num_ph (int): Clip the rate timetrace after the
            max number of photons `max_num_ph` is reached.
        tmin, tmax (float): min and max time (seconds) to include in the
            timetrace. Note that a long time range and a small `binwidth`
            can require a significant amount of memory.
        bursts (bool): if True, plot the burst start-stop times.
        burst_picker (bool): if True, enable the ability to click on bursts
            to obtain burst info. This function requires the matplotlib's QT
            backend.
        scroll (bool): if True, activate a scrolling bar to quickly scroll
            through the timetrace. This function requires the matplotlib's QT
            backend.
        show_rate_th (bool): if True, plot the burst search threshold rate.
        F (bool): if `show_rate` is True, show a rate `F` times larger
            than the background rate.
        rate_th_style (dict): matplotlib style for the rate line.
        show_aa (bool): if True, plot a timetrace for the AexAem photons.
            If False (default), plot timetraces only for DexDem and DexAem
            streams.
        legend (bool): whether to show the legend or not.
        set_ax_limits (bool): if True, set the xlim to zoom on a small portion
            of timetrace. If False, do not set the xlim, display the full
            timetrace.
        burst_color (string): string containing the the HEX RGB color to use
            to highlight the burst regions.
    """

    # Plot bursts
    if bursts:
        t_max_clk = int(tmax/d.clk_p)
        _plot_bursts(d, i, t_max_clk, pmax=500, pmin=-500, color=burst_color)

    # Plot multiple timetraces
    ph_sel_list = [Ph_sel(Dex='Dem'), Ph_sel(Dex='Aem')]
    invert_list = [False, True]
    burst_picker_list = [burst_picker, False]
    scroll_list = [scroll, False]
    if d.ALEX and show_aa:
        ph_sel_list.append(Ph_sel(Aex='Aem'))
        invert_list.append(True)
        burst_picker_list.append(False)
        scroll_list.append(False)

    for ix, (ph_sel, invert) in enumerate(zip(ph_sel_list, invert_list)):
        if not bl.mask_empty(d.get_ph_mask(i, ph_sel=ph_sel)):
            ratetrace_single(
                d, i, m=m, max_num_ph=max_num_ph, tmin=tmin,
                tmax=tmax, ph_sel=ph_sel, invert=invert, bursts=False,
                burst_picker=burst_picker_list[ix],
                scroll=scroll_list[ix],
                show_rate_th=show_rate_th, F=F,
                rate_th_style=rate_th_style, set_ax_limits=set_ax_limits)
    if legend:
        plt.legend(loc='best', fancybox=True)


def sort_burst_sizes(sizes, levels=np.arange(1, 102, 20)):
    """Return a list of masks that split `sizes` in levels.
    Used by timetrace_fret to select burst based on size groups.
    """
    masks = []
    for level1, level2 in zip(levels[:-1], levels[1:]):
        masks.append((sizes >= level1)*(sizes < level2))
    masks.append(sizes >= level2)
    return masks

def timetrace_fret(d, i=0, gamma=1., **kwargs):
    """Timetrace of burst FRET vs time. Uses `plot`."""
    b = d.mburst[i]
    bsizes = bl.select_bursts.get_burst_size(d, ich=i, gamma=gamma)

    style_kwargs = dict(marker='o', mew=0.5, color=blue, mec='grey',
                        alpha=0.4, ls='')
    style_kwargs.update(**kwargs)

    t, E = b.start*d.clk_p, d.E[i]
    levels = sort_burst_sizes(bsizes)
    for ilev, level in enumerate(levels):
        plt.plot(t[level], E[level], ms=np.sqrt((ilev+1)*15),
                 **style_kwargs)
    plt.plot(b.start*d.clk_p, d.E[i], '-k', alpha=0.1, lw=1)
    xlabel('Time (s)'); ylabel('E')
    _gui_timetrace_burst_sel(d, i, timetrace_fret, gcf(), gca())

def timetrace_fret_scatter(d, i=0, gamma=1., **kwargs):
    """Timetrace of burst FRET vs time. Uses `scatter` (slow)."""
    b = d.mburst[i]
    bsizes = d.burst_sizes_ich(d, ich=i, gamma=gamma)

    style_kwargs = dict(s=bsizes, marker='o', alpha=0.5)
    style_kwargs.update(**kwargs)
    plt.scatter(b.start*d.clk_p, d.E[i], **style_kwargs)
    xlabel('Time (s)'); ylabel('E')


def timetrace_bg(d, i=0, nolegend=False, ncol=2, plot_style={}):
    """Timetrace of background rates."""
    bg = d.bg_from(Ph_sel('all'))
    bg_dd = d.bg_from(Ph_sel(Dex='Dem'))
    bg_ad = d.bg_from(Ph_sel(Dex='Aem'))
    t = arange(bg[i].size) * d.bg_time_s
    plot_style_ = dict(linewidth=2, marker='o', markersize=6)
    plot_style_.update(_normalize_kwargs(plot_style, kind='line2d'))
    label = "T: %d cps" % d.bg_mean[Ph_sel('all')][i]
    plot(t, 1e-3 * bg[i], color='k', label=label, **plot_style_)
    label = "DD: %d cps" % d.bg_mean[Ph_sel(Dex='Dem')][i]
    plot(t, 1e-3 * bg_dd[i], color=green, label=label, **plot_style_)
    label = "AD: %d cps" % d.bg_mean[Ph_sel(Dex='Aem')][i]
    plot(t, 1e-3 * bg_ad[i], color=red, label=label, **plot_style_)
    if d.ALEX:
        bg_aa = d.bg_from(Ph_sel(Aex='Aem'))
        label = "AA: %d cps" % d.bg_mean[Ph_sel(Aex='Aem')][i]
        plot(t, 1e-3 * bg_aa[i], label=label, color=purple, **plot_style_)
    if not nolegend:
        legend(loc='best', frameon=False, ncol=ncol)
    xlabel("Time (s)"); ylabel("BG rate (kcps)"); grid(True)
    plt.ylim(ymin=0)


def timetrace_b_rate(d, i=0):
    """Timetrace of bursts-per-second in each period."""
    t = arange(d.bg[i].size)*d.bg_time_s
    b_rate = r_[[(d.bp[i] == p).sum() for p in range(d.bp[i].max()+1)]]
    b_rate /= d.bg_time_s
    if t.size == b_rate.size+1:
        t = t[:-1] # assuming last period without bursts
    else:
        assert t.size == b_rate.size
    plot(t, b_rate, lw=2, label="CH%d" % (i+1))
    legend(loc='best', fancybox=True, frameon=False, ncol=3)
    xlabel("Time (s)"); ylabel("Burst per second"); grid(True)
    plt.ylim(ymin=0)

def time_ph(d, i=0, num_ph=1e4, ph_istart=0):
    """Plot 'num_ph' ph starting at 'ph_istart' marking burst start/end.
    TODO: Update to use the new matplotlib eventplot.
    """
    b = d.mburst[i]
    SLICE = slice(ph_istart, ph_istart+num_ph)
    ph_d = d.ph_times_m[i][SLICE][~d.A_em[i][SLICE]]
    ph_a = d.ph_times_m[i][SLICE][d.A_em[i][SLICE]]

    BSLICE = (b.stop < ph_a[-1])
    start, end = b[BSLICE].start, b[BSLICE].stop

    u = d.clk_p # time scale
    plt.vlines(ph_d*u, 0, 1, color='k', alpha=0.02)
    plt.vlines(ph_a*u, 0, 1, color='k', alpha=0.02)
    plt.vlines(start*u, -0.5, 1.5, lw=3, color=green, alpha=0.5)
    plt.vlines(end*u, -0.5, 1.5, lw=3, color=red, alpha=0.5)
    xlabel("Time (s)")


##
#  Histogram plots
#
def _bins_array(bins):
    """When `bins` is a 3-element sequence returns an array of bin edges.
    Otherwise returns the `bins` unchaged.
    """
    if np.size(bins) == 3:
        bins = np.arange(*bins)
    return bins


def _hist_burst_taildist(data, bins, pdf,  weights=None, yscale='log',
                         color=None, label=None, plot_style=None):
    hist = HistData(*np.histogram(data, bins=_bins_array(bins),
                                  weights=weights))
    ydata = hist.pdf if pdf else hist.counts

    default_plot_style = dict(marker='o')
    if plot_style is None:
        plot_style = {}
    if color is not None:
        plot_style['color'] = color
    if label is not None:
        plot_style['label'] = label
    default_plot_style.update(_normalize_kwargs(plot_style, kind='line2d'))
    plt.plot(hist.bincenters, ydata, **default_plot_style)
    plt.yscale(yscale)
    if pdf:
        plt.ylabel('PDF')
    else:
        plt.ylabel('# Bursts')


def hist_width(d, i=0, bins=(0, 10, 0.025), pdf=True, weights=None,
               yscale='log', color=None, plot_style=None):
    """Plot histogram of burst durations.

    Parameters:
        d (Data): Data object
        i (int): channel index
        bins (array or None): array of bin edges. If len(bins) == 3
            then is interpreted as (start, stop, step) values.
        pdf (bool): if True, normalize the histogram to obtain a PDF.
        color (string or tuple or None): matplotlib color used for the plot.
        yscale (string): 'log' or 'linear', sets the plot y scale.
        plot_style (dict): dict of matplotlib line style passed to `plot`.
    """
    weights = weights[i] if weights is not None else None
    burst_widths = d.mburst[i].width*d.clk_p*1e3

    _hist_burst_taildist(burst_widths, bins, pdf, weights=weights,
                         yscale=yscale, color=color, plot_style=plot_style)
    plt.xlabel('Burst width (ms)')
    plt.xlim(xmin=0)


def hist_brightness(d, i=0, bins=(0, 60, 1), pdf=True, weights=None,
                    yscale='log', gamma=1, add_naa=False, beta=1.,
                    donor_ref=True, label_prefix=None,
                    color=None, plot_style=None):
    """Plot histogram of burst brightness, i.e. burst size / duration.

    Parameters:
        d (Data): Data object
        i (int): channel index
        bins (array or None): array of bin edges. If len(bins) == 3
            then is interpreted as (start, stop, step) values.
        gamma, beta (floats): factors used to compute the corrected burst
            size. See :meth:`fretbursts.burstlib.Data.burst_sizes_ich`.
        add_naa (bool): if True, include `naa` to the total burst size.
        donor_ref (bool): convention used for corrected burst size computation.
            See :meth:`fretbursts.burstlib.Data.burst_sizes_ich` for details.
        label_prefix (string or None): a custom prefix for the legend label.
        color (string or tuple or None): matplotlib color used for the plot.
        pdf (bool): if True, normalize the histogram to obtain a PDF.
        yscale (string): 'log' or 'linear', sets the plot y scale.
        plot_style (dict): dict of matplotlib line style passed to `plot`.
    """
    weights = weights[i] if weights is not None else None
    if plot_style is None:
        plot_style = {}

    burst_widths = d.mburst[i].width*d.clk_p*1e3
    sizes = d.burst_sizes_ich(ich=i, gamma=gamma, beta=beta, add_naa=add_naa,
                              donor_ref=donor_ref)
    brightness = sizes / burst_widths
    label = 'nd + na/g' if donor_ref else 'g*nd + na'
    if add_naa:
        label += " + naa/(b*g)" if donor_ref else ' + naa/b'
    label = '(' + label + ') / w'
    if label_prefix is not None:
        label = label_prefix + ' ' + label

    # Use default label (with optional prefix) only if not explicitly
    # specified in `plot_style`
    if 'label' not in plot_style:
        plot_style['label'] = label

    _hist_burst_taildist(brightness, bins, pdf, weights=weights,
                         yscale=yscale, color=color, plot_style=plot_style)
    plt.xlabel('Burst brightness (kHz)')
    plt.legend(loc='best')


def hist_size(d, i=0, which='all', bins=(0, 600, 4), pdf=False, weights=None,
              yscale='log', gamma=1, add_naa=False, beta=1, donor_ref=True,
              label_prefix=None, legend=True, color=None, plot_style=None):
    """Plot histogram of burst sizes.

    Parameters:
        d (Data): Data object
        i (int): channel index
        bins (array or None): array of bin edges. If len(bins) == 3
            then is interpreted as (start, stop, step) values.
        which (string): what photons to include in "size". Valid values are
            'all', 'nd', 'na', 'naa'. When 'all', sizes are computed with
            `d.burst_sizes()` (by default nd + na); 'nd', 'na', 'naa' get
            counts from `d.nd`, `d.na`, `d.naa` (respectively Dex-Dem,
            Dex-Aem, Aex-Aem).
        gamma, beta (floats): factors used to compute the corrected burst
            size. Ignored when `which` != 'all'.
            See :meth:`fretbursts.burstlib.Data.burst_sizes_ich`.
        add_naa (bool): if True, include `naa` to the total burst size.
        donor_ref (bool): convention used for corrected burst size computation.
            See :meth:`fretbursts.burstlib.Data.burst_sizes_ich` for details.
        label_prefix (string or None): a custom prefix for the legend label.
        color (string or tuple or None): matplotlib color used for the plot.
        pdf (bool): if True, normalize the histogram to obtain a PDF.
        yscale (string): 'log' or 'linear', sets the plot y scale.
        legend (bool): if True add legend to plot
        plot_style (dict): dict of matplotlib line style passed to `plot`.
    """
    weights = weights[i] if weights is not None else None
    if plot_style is None:
        plot_style = {}

    which_dict = {'all': 'k', 'nd': green, 'na': red, 'naa': purple}
    assert which in which_dict
    if which == 'all':
        sizes = d.burst_sizes_ich(ich=i, gamma=gamma, add_naa=add_naa,
                                  beta=beta, donor_ref=donor_ref)
        label = 'nd + na/g' if donor_ref else 'g*nd + na'
        if add_naa:
            label += " + naa/(b*g)" if donor_ref else ' + naa/b'
    else:
        sizes = d[which][i]
        label = which

    # Use default label (with optional prefix) only if not explicitly
    # specified in `plot_style`
    if 'label' not in plot_style:
        if label_prefix is not None:
            label = label_prefix + ' ' + label
        plot_style['label'] = label

    # Use default color only if not specified in `color` or `plot_style`
    if color is None and 'color' not in plot_style:
        plot_style['color'] = which_dict[which]
    elif color is not None:
        plot_style['color'] = color

    _hist_burst_taildist(sizes, bins, pdf, weights=weights, yscale=yscale,
                         plot_style=plot_style)
    plt.xlabel('Burst size')
    if legend:
        plt.legend(loc='best')

def hist_size_all(d, i=0, **kwargs):
    """Plot burst sizes for all the combinations of photons.

    Calls :func:`hist_size` multiple times with different `which` parameters.
    """
    fields = ['all', 'nd', 'na']
    if d.ALEX:
        fields.append('naa')
    for which in fields:
        hist_size(d, i, which=which, **kwargs)


def _get_fit_E_text(d, pylab=True):
    """Return a formatted string for fitted E."""
    delta = (d.E_fit.max()-d.E_fit.min())*100
    fit_text = r'\langle{E}_{fit}\rangle = %.3f \qquad ' % d.E_fit.mean()
    fit_text += r'\Delta E_{fit} = %.2f \%%' % delta
    if pylab: fit_text = r'$'+fit_text+r'$'
    return fit_text

def _fitted_E_plot(d, i=0, F=1, no_E=False, ax=None, show_model=True,
                   verbose=False, two_gauss_model=False, lw=2.5, color='k',
                   alpha=0.5, fillcolor=None):
    """Plot a fitted model overlay on a FRET histogram."""
    if ax is None:
        ax2 = gca()
    else:
        ax2 = plt.twinx(ax=ax)
        ax2.grid(False)

    if d.fit_E_curve and show_model:
        x = r_[-0.2:1.21:0.002]
        y = d.fit_E_model(x, d.fit_E_res[i, :])
        scale = F*d.fit_E_model_F[i]
        if two_gauss_model:
            assert d.fit_E_res.shape[1] > 2
            if d.fit_E_res.shape[1] == 5:
                m1, s1, m2, s2, a1 = d.fit_E_res[i, :]
                a2 = (1-a1)
            elif d.fit_E_res.shape[1] == 6:
                m1, s1, a1, m2, s2, a2 = d.fit_E_res[i, :]
            y1 = a1*normpdf(x, m1, s1)
            y2 = a2*normpdf(x, m2, s2)
            ax2.plot(x, scale*y1, ls='--', lw=lw, alpha=alpha, color=color)
            ax2.plot(x, scale*y2, ls='--', lw=lw, alpha=alpha, color=color)
        if fillcolor == None:
            ax2.plot(x, scale*y, lw=lw, alpha=alpha, color=color)
        else:
            ax2.fill_between(x, scale*y, lw=lw, alpha=alpha, edgecolor=color,
                             facecolor=fillcolor, zorder=10)
        if verbose:
            print('Fit Integral:', np.trapz(scale*y, x))

    ax2.axvline(d.E_fit[i], lw=3, color=red, ls='--', alpha=0.6)
    xtext = 0.6 if d.E_fit[i] < 0.6 else 0.2
    if d.nch > 1 and not no_E:
        ax2.text(xtext, 0.81, "CH%d: $E_{fit} = %.3f$" % \
                     (i+1, d.E_fit[i]),
                 transform=gca().transAxes, fontsize=16,
                 bbox=dict(boxstyle='round', facecolor='#dedede', alpha=0.5))

def hist_burst_data(
        d, i=0, data_name='E', ax=None, binwidth=0.03, bins=None,
        vertical=False, pdf=False, hist_style='bar',
        weights=None, gamma=1., add_naa=False,            # weights args
        show_fit_stats=False, show_fit_value=False, fit_from='kde',
        show_kde=False, bandwidth=0.03, show_kde_peak=False,  # kde args
        show_model=False, show_model_peaks=True,
        hist_bar_style=None, hist_plot_style=None, model_plot_style=None,
        kde_plot_style=None, verbose=False):
    """Plot burst_data (i.e. E, S, etc...) histogram and KDE.

    This a generic function to plot histograms for any burst data.
    In particular this function is called by :func:`hist_fret` and
    :func:`hist_S` to make E and S histograms respectively.

    Histograms and KDE can be plotted on any `Data` variable after
    burst search. To show a model, a model must be fitted first by calling
    `d.E_fitter.fit_histogram()`. To show the KDE peaks position, they
    must be computed first with `d.E_fitter.find_kde_max()`.

    The arguments are shown below grouped in logical sections.

    **Generic arguments**

    Args:
        data_name (string): name of the burst data (i.e. 'E' or 'S')
        ax (None or matplotlib axis): optional axis instance to plot in.
        vertical (bool): if True the x axis is oriented vertically.
        verbose (bool): if False, suppress any printed output.

    **Histogram arguments**: control the histogram appearance

    Args:
        hist_style (string): if 'bar' use a classical bar histogram,
            otherwise do a normal line plot of bin counts vs bin centers
        bins (None or array): if None the bins are computed according to
            `binwidth`. If not None contains the arrays of bin edges
            and overrides `binwidth`.
        binwidth (float): bin width for the histogram.
        pdf (bool): if True, normalize the histogram to obtain a PDF.
        hist_bar_style (dict): style dict for the histogram when
            `hist_style == 'bar'`.
        hist_plot_style (dict): style dict for the histogram when
            `hist_style != 'bar'`.

    **Model arguments**: control the model plot

    Args:
        show_model (bool): if True shows the model fitted to the histogram
        model (lmfit.Model object or None): lmfit Model used for histogram
            fitting. If None the histogram is not fitted.
        show_model_peaks (bool): if True marks the position of model peaks
        model_plot_style (dict): style dict for the model plot

    **KDE arguments**: control the KDE plot

    Args:
        show_kde (bool): if True shows the KDE curve
        show_kde_peak (bool): if True marks the position of the KDE peak
        bandwidth (float or None): bandwidth used to compute the KDE
            If None the KDE is not computed.
        kde_plot_style (dict): style dict for the KDE curve

    **Weights arguments** (weights are used to weight bursts according to
    their size, affecting histograms and KDEs).

    Args:
        weights (string or None): kind of burst-size weights.
            See :func:`fretbursts.fret_fit.get_weights`.
        gamma (float): gamma factor passed to `get_weights()`.
        add_naa (bool): if True adds `naa` to the burst size.

    **Fit text arguments**: control how to print annotation with
    fit information.

    Args:
        fit_from (string): determines how to obtain the fit value. If 'kde'
            the fit value is the KDE peak. Otherwise it must be the name
            of a model parameter that will be used as fit value.
        show_fit_value (bool): if True annotate the plot with fit value.
        show_fit_stats (bool): if True annotate the figure with mean fit
            value and max deviation across the channels (for multi-spot).
    """

    assert data_name in d
    fitter_name = data_name + '_fitter'

    if ax is None:
        ax = gca()
    pline = ax.axhline if vertical else ax.axvline
    bar = ax.barh if vertical else ax.bar
    xlabel, ylabel = ax.set_xlabel, ax.set_ylabel
    xlim, ylim = ax.set_xlim, ax.set_ylim
    if vertical:
        xlabel, ylabel = ylabel, xlabel
        xlim, ylim = ylim, xlim
    weights_tuple = (weights, float(gamma), add_naa)
    if not hasattr(d, fitter_name) or _is_list_of_arrays(weights) \
            or getattr(d, data_name+'_weights') != weights_tuple:
        if hasattr(d, fitter_name):
            print(' - Overwriting the old %s object with the new weights.' %\
                    fitter_name)
            if verbose:
                print('   Old weights:', getattr(d, data_name+'_weights'))
                print('   New weights:', weights_tuple)
        bext.bursts_fitter(d, burst_data=data_name, weights=weights,
                           gamma=gamma, add_naa=add_naa)

    # fitter_name is only an attribute of Data, not a key in the dictionary
    fitter = getattr(d, fitter_name)
    fitter.histogram(binwidth=binwidth, bins=bins, verbose=verbose)
    if pdf:
        ylabel('PDF')
        hist_vals = fitter.hist_pdf[i]
    else:
        ylabel('# Bursts')
        hist_vals = fitter.hist_counts[i]
    xlabel(data_name)
    if data_name in ['E', 'S']:
        xlim(-0.19, 1.19)

    hist_bar_style_ = dict(facecolor='#80b3ff', edgecolor='#5f8dd3',
                           linewidth=1.5, alpha=0.7, label='E Histogram')
    hist_bar_style_.update(**_normalize_kwargs(hist_bar_style))

    hist_plot_style_ = dict(linestyle='-', marker='o', markersize=6,
                            linewidth=2, alpha=0.6, label='E Histogram')
    hist_plot_style_.update(_normalize_kwargs(hist_plot_style,
                                              kind='line2d'))
    if hist_style == 'bar':
        bar(fitter.hist_bins[:-1], hist_vals, fitter.hist_binwidth,
            **hist_bar_style_)
    else:
        if vertical:
            ax.plot(hist_vals, fitter.hist_axis, **hist_plot_style_)
        else:
            ax.plot(fitter.hist_axis, hist_vals, **hist_plot_style_)

    if show_model or show_kde:
        if pdf:
            scale = 1
        else:
            scale = fitter.hist_binwidth * d.num_bursts[i]

    if show_model:
        model_plot_style_ = dict(color='k', alpha=0.8, label='Model')
        model_plot_style_.update(_normalize_kwargs(model_plot_style,
                                                   kind='line2d'))
        fit_res = fitter.fit_res[i]
        x = fitter.x_axis
        y = fit_res.model.eval(x=x, **fit_res.values)
        xx, yy = (y, x) if vertical else (x, y)
        ax.plot(xx, yy, **model_plot_style_)
        if fit_res.model.components is not None:
            for component in fit_res.model.components:
                model_plot_style_.update(ls='--', label='Model component')
                y = component.eval(x=x, **fit_res.values)
                xx, yy = (y, x) if vertical else (x, y)
                ax.plot(xx, yy, **model_plot_style_)
        if show_model_peaks:
            for param in fitter.params:
                if param.endswith('center'):
                    pline(fitter.params[param][i], ls='--', color=red)
    if show_kde:
        x = fitter.x_axis
        fitter.calc_kde(bandwidth=bandwidth)
        kde_plot_style_ = dict(linewidth=1.5, color='k', alpha=0.8,
                               label='KDE')
        kde_plot_style_.update(_normalize_kwargs(kde_plot_style,
                                                 kind='line2d'))
        y = scale * fitter.kde[i](x)
        xx, yy = (y, x) if vertical else (x, y)
        ax.plot(xx, yy, **kde_plot_style_)
    if show_kde_peak:
        pline(fitter.kde_max_pos[i], ls='--', color='orange')

    if show_fit_value or show_fit_stats:
        if fit_from == 'kde':
            fit_arr = fitter.kde_max_pos
        else:
            assert fit_from in fitter.params
            fit_arr = fitter.params[fit_from]

        if i == 0:
            if show_fit_stats:
                plt.figtext(0.4, 0.01, _get_fit_text_stats(fit_arr),
                            fontsize=16)
        if show_fit_value:
            _plot_fit_text_ch(fit_arr, i, ax=ax)

def hist_fret(
        d, i=0, ax=None, binwidth=0.03, bins=None, pdf=True,
        hist_style='bar',
        weights=None, gamma=1., add_naa=False,            # weights args
        show_fit_stats=False, show_fit_value=False, fit_from='kde',
        show_kde=False, bandwidth=0.03, show_kde_peak=False,  # kde args
        show_model=False, show_model_peaks=True,
        hist_bar_style=None, hist_plot_style=None, model_plot_style=None,
        kde_plot_style=None, verbose=False):
    """Plot FRET histogram and KDE.

    The most used argument is `binwidth` that sets the histogram bin width.

    For detailed documentation see :func:`hist_burst_data`.
    """

    hist_burst_data(
        d, i, data_name='E', ax=ax, binwidth=binwidth, bins=bins,
        pdf=pdf, weights=weights, gamma=gamma, add_naa=add_naa,
        hist_style=hist_style, show_fit_stats=show_fit_stats,
        show_fit_value=show_fit_value, fit_from=fit_from,
        show_kde=show_kde, bandwidth=bandwidth,
        show_kde_peak=show_kde_peak,  # kde args
        show_model=show_model, show_model_peaks=show_model_peaks,
        hist_bar_style=hist_bar_style, hist_plot_style=hist_plot_style,
        model_plot_style=model_plot_style, kde_plot_style=kde_plot_style,
        verbose=verbose)

def hist_S(
        d, i=0, ax=None, binwidth=0.03, bins=None, pdf=True,
        hist_style='bar',
        weights=None, gamma=1., add_naa=False,                # weights args
        show_fit_stats=False, show_fit_value=False, fit_from='kde',
        show_kde=False, bandwidth=0.03, show_kde_peak=False,  # kde args
        show_model=False, show_model_peaks=True,
        hist_bar_style=None, hist_plot_style=None, model_plot_style=None,
        kde_plot_style=None, verbose=False):
    """Plot S histogram and KDE.

    The most used argument is `binwidth` that sets the histogram bin width.

    For detailed documentation see :func:`hist_burst_data`.    """

    hist_burst_data(
        d, i, data_name='S', ax=ax, binwidth=binwidth, bins=bins,
        pdf=pdf, weights=weights, gamma=gamma, add_naa=add_naa,
        hist_style=hist_style, show_fit_stats=show_fit_stats,
        show_fit_value=show_fit_value, fit_from=fit_from,
        show_kde=show_kde, bandwidth=bandwidth,
        show_kde_peak=show_kde_peak,  # kde args
        show_model=show_model, show_model_peaks=show_model_peaks,
        hist_bar_style=hist_bar_style, hist_plot_style=hist_plot_style,
        model_plot_style=model_plot_style, kde_plot_style=kde_plot_style,
        verbose=verbose)

def _get_fit_text_stats(fit_arr, pylab=True):
    """Return a formatted string for mean E and max delta-E."""
    delta = (fit_arr.max() - fit_arr.min())*100
    fit_text = r'\langle{E}_{fit}\rangle = %.3f \qquad ' % fit_arr.mean()
    fit_text += r'\Delta E_{fit} = %.2f \%%' % delta
    if pylab: fit_text = r'$'+fit_text+r'$'
    return fit_text

def _plot_fit_text_ch(
        fit_arr, ich, fmt_str="CH%d: $E_{fit} = %.3f$", ax=None,
        bbox=dict(boxstyle='round', facecolor='#dedede', alpha=0.5),
        xtext_low=0.2, xtext_high=0.6, fontsize=16):
    """Plot a text box with ch and fit value."""
    if ax is None: ax = gca()
    xtext = xtext_high if fit_arr[ich] < xtext_high else xtext_low
    ax.text(xtext, 0.81, fmt_str % (ich+1, fit_arr[ich]),
            transform=ax.transAxes, fontsize=fontsize, bbox=bbox)


def hist2d_alex(d, i=0, vmin=2, vmax=0, binwidth=0.05, S_max_norm=0.8,
                interp='bicubic', cmap='hot', under_color='white',
                over_color='white', scatter=True, scatter_ms=3,
                scatter_color='orange', scatter_alpha=0.2, gui_sel=False,
                cbar_ax=None, grid_color='#D0D0D0'):
    """Plot 2-D E-S ALEX histogram with a scatterplot overlay.
    """
    ax = plt.gca()
    d._calc_alex_hist(binwidth)
    ES_hist, E_bins, S_bins, S_ax = d.ES_hist[i], d.E_bins, d.S_bins, d.S_ax

    colormap = plt.get_cmap(cmap)
    # Heuristic for colormap range
    if vmax <= vmin:
        S_range = (S_ax < S_max_norm)
        vmax = ES_hist[:, S_range].max()
        if vmax <= vmin: vmax = 10*vmin

    if scatter:
        ax.plot(d.E[i], d.S[i], 'o', mew=0, ms=scatter_ms,
                alpha=scatter_alpha, color=scatter_color)
    im = ax.imshow(ES_hist[:, ::-1].T, interpolation=interp,
                   extent=(E_bins[0], E_bins[-1], S_bins[0], S_bins[-1]),
                   vmin=vmin, vmax=vmax, cmap=colormap)
    im.cmap.set_under(under_color)
    im.cmap.set_over(over_color)
    if cbar_ax is None:
        gcf().colorbar(im)
    else:
        cbar_ax.colorbar(im)
    ax.set_xlim(-0.2, 1.2)
    ax.set_ylim(-0.2, 1.2)
    ax.set_xlabel('E')
    ax.set_ylabel('S')
    ax.grid(color=grid_color)
    if gui_sel:
        # the selection object must be saved (otherwise will be destroyed)
        hist2d_alex.gui_sel = gs.rectSelection(gcf(), gca())

def hexbin_alex(d, i=0, figsize=(6, 5), vmin=0, vmax_fret=True, gridsize=40,
                **hexbin_kwargs):
    """Plot an hexbin 2D histogram for E-S.
    """
    fig = plt.gcf()
    fig.set_size_inches(*figsize)

    hexbin_kwargs_ = dict(extent=(-0.2, 1.2, -0.2, 1.2), gridsize=gridsize,
                          mincnt=1, cmap='GnBu_r')
    if hexbin_kwargs is not None:
        hexbin_kwargs_.update(hexbin_kwargs)
    poly = plt.hexbin(d.E[i], d.S[i], **hexbin_kwargs_)
    vmax = _alex_hexbin_vmax(poly, vmax_fret=vmax_fret)
    poly.set_clim(vmin, vmax)
    plt.xlabel('E')
    plt.ylabel('S')
    plt.colorbar()

def plot_ES_selection(ax, E1, E2, S1, S2, rect=True, **kwargs):
    """Plot an overlay ROI on top of an E-S plot (i.e. ALEX histogram).

    This function plots a rectangle and inscribed ellipsis with x-axis limits
    (E1, E2) and y-axis limits (S1, S2).

    Note that, a dict with keys (E1, E2, S1, S2, rect) can be also passed to
    :func:`fretbursts.select_bursts.ES` to apply a selection.

    Parameters:
        ax (matplotlib axis): the axis where the rectangle is plotted.
            Typically you pass the axis of a previous E-S scatter plot
            or histogram.
        E1, E2, S1, S2 (floats): limits for E and S (X and Y axis respectively)
            used to plot the rectangle.
        rect (bool): if True, the rectangle is highlighted and the ellipsis is
            grey. The color are swapped otherwise.
        **kwargs: other keywords passed to both matplotlib's `Rectangle`
            and `Ellipse`.

    See also:
        For selecting bursts according to (`E1`, `E2`, `S1`, `S2`, `rect`) see:

        - :func:`fretbursts.select_bursts.ES`
    """
    if rect:
        rect_color, ellips_color = blue, 'gray'
    else:
        rect_color, ellips_color = 'gray', blue
    patch_style = dict(fill=False, lw=1.5, alpha=0.5)
    patch_style.update(**kwargs)
    rect = Rectangle(xy=(E1, S1), height=(S2 - S1), width=(E2 - E1),
                     color=rect_color, **patch_style)
    ellips = Ellipse(xy=(0.5*(E1 +  E2), 0.5*(S1 + S2)), height=(S2 - S1),
                     width=(E2 - E1), color=ellips_color, **patch_style)
    ax.add_patch(rect)
    ax.add_patch(ellips)
    return rect, ellips

def get_ES_range():
    """Get the range of ES histogram selected via GUI.

    Prints E1, E2, S1, S2 and return a dict containig these values.
    """
    sel = None
    if hasattr(hist2d_alex.gui_sel, 'selection'):
        sel = hist2d_alex.gui_sel.selection
        print('E1={E1:.3}, E2={E2:.3}, S1={S1:.3}, S2={S2:.3}'.format(**sel))
    return sel


def hist_interphoton_single(d, i=0, binwidth=1e-4, tmax=None, bins=None,
                            ph_sel=Ph_sel('all'), period=None,
                            yscale='log', xscale='linear', xunit='ms',
                            plot_style=None):
    """Plot histogram of interphoton delays for a single photon streams.

    Arguments:
        d (Data object): the input data.
        i (int): the channel for which the plot must be done. Default is 0.
            For single-spot data the only valid value is 0.
        binwidth (float): histogram bin width in seconds.
        tmax (float or None): max timestamp delay in the histogram (seconds).
            If None (default), uses the the max timestamp delay in the stream.
            If not None, the plotted histogram may be further trimmed to
            the smallest delay with counts > 0 if this delay happens to be
            smaller than `tmax`.
        bins (array or None): specifies the bin edged (in seconds). When
            `bins` is not None then the arguments `binwidth` and `tmax`
            are ignored. When `bins` is None, the bin edges are computed
            from the `binwidth` and `tmax` arguments.
        ph_sel (Ph_sel object): photon stream for which plotting the histogram
        period (int): the background period to use for plotting the histogram.
            The background period is a time-slice of the measurement from which
            timestamps are taken. If `period` is None (default) the
            time-windows is the full measurement.
        yscale (string): scale for the y-axis. Valid values include 'log' and
            'linear'. Default 'log'.
        xscale (string): scale for the x-axis. Valid values include 'log' and
            'linear'. Default 'linear'.
        xunit (string): unit used for the x-axis. Valid values are 's', 'ms',
            'us', 'ns'. Default 'ms'.
        plot_style (dict): keyword arguments to be passed to matplotlib's
            `plot` function. Used to customize the plot style.
    """
    unit_dict = {'s': 1, 'ms': 1e3, 'us': 1e6, 'ns': 1e9}
    assert xunit in unit_dict
    scalex = unit_dict[xunit]

    # Compute interphoton delays
    if period is None:
        ph_times = d.get_ph_times(ich=i, ph_sel=ph_sel)
    else:
        ph_times = d.get_ph_times_period(ich=i, period=period, ph_sel=ph_sel)
    delta_ph_t = np.diff(ph_times) * d.clk_p
    if tmax is None:
        tmax = delta_ph_t.max()

    # Compute bin edges if not passed in
    if bins is None:
        # Shift by half clk_p to avoid "beatings" in the distribution
        # due to floating point inaccuracies.
        bins = np.arange(0, tmax + binwidth, binwidth) - 0.5 * d.clk_p
    else:
        warnings.warn('Using `bins` and ignoring `tmax` and `binwidth`.')
    t_ax = bins[:-1] + 0.5 * binwidth

    # Compute interphoton histogram
    counts, _ = np.histogram(delta_ph_t, bins=bins)

    # Max index with counts > 0
    n_trim = np.trim_zeros(counts).size + 1

    # Plot histograms
    plot_style_ = dict(marker='o', markersize=5, linestyle='none', alpha=0.6)
    if ph_sel in _ph_sel_color_dict:
        plot_style_['color'] = _ph_sel_color_dict[ph_sel]
        plot_style_['label'] = _ph_sel_label_dict[ph_sel]
    plot_style_.update(_normalize_kwargs(plot_style, kind='line2d'))
    plot(t_ax[:n_trim] * scalex, counts[:n_trim], **plot_style_)

    if yscale == 'log':
        gca().set_yscale(yscale)
        plt.ylim(1)
        _plot_status['hist_interphoton_single'] = {'autoscale': False}
    if xscale == 'log':
        gca().set_xscale(yscale)
        plt.xlim(0.5 * binwidth)
        _plot_status['hist_interphoton_single'] = {'autoscale': False}
    plt.xlabel('Inter-photon delays (%s)' % xunit.replace('us', 'μs'))
    plt.ylabel('# Delays')
    # Return interal variables so that other functions can extend the plot
    return dict(counts=counts, n_trim=n_trim, plot_style_=plot_style_,
                t_ax=t_ax, scalex=scalex)


def hist_interphoton(d, i=0, binwidth=1e-4, tmax=None, bins=None, period=None,
                     yscale='log', xscale='linear', xunit='ms', plot_style=None,
                     show_da=False, legend=True):
    """Plot histogram of photon interval for different photon streams.

    Arguments:
        d (Data object): the input data.
        i (int): the channel for which the plot must be done. Default is 0.
            For single-spot data the only valid value is 0.
        binwidth (float): histogram bin width in seconds.
        tmax (float or None): max timestamp delay in the histogram (seconds).
            If None (default), uses the the max timestamp delay in the stream.
            If not None, the plotted histogram may be further trimmed to
            the smallest delay with counts > 0 if this delay happens to be
            smaller than `tmax`.
        bins (array or None): specifies the bin edged (in seconds). When
            `bins` is not None then the arguments `binwidth` and `tmax`
            are ignored. When `bins` is None, the bin edges are computed
            from the `binwidth` and `tmax` arguments.
        period (int): the background period to use for plotting the histogram.
            The background period is a time-slice of the measurement from which
            timestamps are taken. If `period` is None (default) the
            time-windows is the full measurement.
        yscale (string): scale for the y-axis. Valid values include 'log' and
            'linear'. Default 'log'.
        xscale (string): scale for the x-axis. Valid values include 'log' and
            'linear'. Default 'linear'.
        xunit (string): unit used for the x-axis. Valid values are 's', 'ms',
            'us', 'ns'. Default 'ms'.
        plot_style (dict): keyword arguments to be passed to matplotlib's
            `plot` function. Used to customize the plot style.
        show_da (bool): If False (default) do not plot the AexDem photon stream.
            Ignored when the measurement is not ALEX.
        legend (bool): If True (default) plot a legend.
    """
    # Plot multiple timetraces
    ph_sel_list = [Ph_sel('all'), Ph_sel(Dex='Dem'), Ph_sel(Dex='Aem')]
    if d.ALEX:
        ph_sel_list.append(Ph_sel(Aex='Aem'))
        if show_da:
            ph_sel_list.append(Ph_sel(Aex='Dem'))

    for ix, ph_sel in enumerate(ph_sel_list):
        if not bl.mask_empty(d.get_ph_mask(i, ph_sel=ph_sel)):
            hist_interphoton_single(d, i=i, binwidth=binwidth, tmax=tmax,
                                    bins=bins, period=period, ph_sel=ph_sel,
                                    yscale=yscale, xscale=xscale, xunit=xunit,
                                    plot_style=plot_style)
    if legend:
        plt.legend(loc='best', fancybox=True)

    if yscale == 'log' or xscale == 'log':
        _plot_status['hist_interphoton'] = {'autoscale': False}


def hist_bg_single(d, i=0, binwidth=1e-4, tmax=0.01, bins=None,
                   ph_sel=Ph_sel('all'), period=0,
                   yscale='log', xscale='linear', xunit='ms', plot_style=None,
                   show_fit=True, fit_style=None, manual_rate=None):
    """Plot histogram of photon interval for a single photon streams.

    Optionally plots the fitted background as an exponential curve.
    Most arguments are described in :func:`hist_interphoton_single`.
    In the following we document only the additional arguments.

    Arguments:
        show_fit (bool): If True shows the fitted background rate as an
            exponential distribution.
        manual_rate (float or None): When not None use this value as background
            rate (ignoring the value saved in Data).
        fit_style (dict): arguments passed to matplotlib's `plot` for
            for plotting the exponential curve.

    For a description of all the other arguments see
    :func:`hist_interphoton_single`.
    """
    hist = hist_interphoton_single(d, i=i, binwidth=binwidth, tmax=tmax,
                                   bins=bins, ph_sel=ph_sel, period=period,
                                   yscale=yscale, xscale=xscale, xunit=xunit,
                                   plot_style=None)

    if show_fit or manual_rate is not None:
        # Compute the fit function
        if manual_rate is not None:
            bg_rate = manual_rate
        else:
            bg_rate = d.bg_from(ph_sel)[i][period]
        tau_th = bg_rate * 1e-6

        i_tau_th = np.searchsorted(hist['t_ax'], tau_th)
        counts_integral = hist['counts'][i_tau_th:].sum()
        y_fit = np.exp(- hist['t_ax'] * bg_rate)
        y_fit *= counts_integral / y_fit[i_tau_th:].sum()

        # Plot
        fit_style_ = dict(hist['plot_style_'])
        fit_style_.update(linestyle='-', marker='', label='auto')
        fit_style_.update(_normalize_kwargs(fit_style, kind='line2d'))
        if fit_style_['label'] is 'auto':
            plt_label = hist['plot_style_'].get('label', None)
            label = str(ph_sel) if plt_label is None else plt_label
            fit_style_['label'] = '%s, %.2f kcps' % (label, bg_rate * 1e-3)
        n_trim = hist['n_trim']
        plot(hist['t_ax'][:n_trim] * hist['scalex'], y_fit[:n_trim],
             **fit_style_)


def hist_bg(d, i=0, binwidth=1e-4, tmax=0.01, bins=None, period=0,
            yscale='log', xscale='linear', xunit='ms', plot_style=None,
            show_da=False, legend=True, show_fit=True, fit_style=None):
    """Plot histogram of photon interval for different photon streams.

    Optionally plots the fitted background.
    Most arguments are described in :func:`hist_interphoton`.
    In the following we document only the additional arguments.

    Arguments:
        show_fit (bool): If True shows the fitted background rate as an
            exponential distribution.
        fit_style (dict): arguments passed to matplotlib's `plot` for
            for plotting the exponential curve.

    For a description of all the other arguments see :func:`hist_interphoton`.
    """
    # Plot multiple timetraces
    ph_sel_list = [Ph_sel('all'), Ph_sel(Dex='Dem'), Ph_sel(Dex='Aem')]
    if d.ALEX:
        ph_sel_list.append(Ph_sel(Aex='Aem'))
        if show_da:
            ph_sel_list.append(Ph_sel(Aex='Dem'))

    for ix, ph_sel in enumerate(ph_sel_list):
        if not bl.mask_empty(d.get_ph_mask(i, ph_sel=ph_sel)):
            hist_bg_single(d, i=i, period=period, binwidth=binwidth,
                           bins=bins, tmax=tmax, ph_sel=ph_sel, xunit=xunit,
                           show_fit=show_fit, yscale=yscale, xscale=xscale,
                           plot_style=plot_style, fit_style=fit_style)
    if legend:
        plt.legend(loc='best', fancybox=True)

    if yscale == 'log' or xscale == 'log':
        _plot_status['hist_bg'] = {'autoscale': False}


def hist_ph_delays(
        d, i=0, time_min_s=0, time_max_s=30, bin_width_us=10, mask=None,
        yscale='log', hfit_bin_ms=1, efit_tail_min_us=1000, **kwargs):
    """Histog. of ph delays and comparison with 3 BG fitting functions.
    """
    ph = d.ph_times_m[i].copy()
    if mask is not None: ph = ph[mask[i]]
    ph = ph[(ph < time_max_s/d.clk_p)*(ph > time_min_s/d.clk_p)]
    dph = np.diff(ph)*d.clk_p
    H = hist(dph*1e6, bins=r_[0:1200:bin_width_us], histtype='step', **kwargs)
    gca().set_yscale('log')
    xlabel(u'Ph delay time (μs)'); ylabel("# Ph")

    efun = lambda t, r: np.exp(-r*t)*r
    re = bg.exp_fit(ph, tail_min_us=efit_tail_min_us)
    rg = bg.exp_hist_fit(ph, tail_min_us=efit_tail_min_us, binw=hfit_bin_ms*1e3)
    rc = bg.exp_cdf_fit(ph, tail_min_us=efit_tail_min_us)
    t = r_[0:1200]*1e-6
    F = 1 if 'normed' in kwargs else H[0].sum()*(bin_width_us)
    plot(t*1e6, 0.65*F*efun(t, rc)*1e-6, lw=3, alpha=0.5, color=purple,
         label="%d cps - Exp CDF (tail_min_p=%.2f)" % (rc, efit_tail_min_us))
    plot(t*1e6, 0.65*F*efun(t, re)*1e-6, lw=3, alpha=0.5, color=red,
         label="%d cps - Exp ML (tail_min_p=%.2f)" % (re, efit_tail_min_us))
    plot(t*1e6, 0.68*F*efun(t, rg)*1e-6, lw=3, alpha=0.5, color=green,
         label=u"%d cps - Hist (bin_ms=%d) [Δ=%d%%]" % (hfit_bin_ms, rg,
                                                        100*(rg-re)/re))
    plt.legend(loc='best', fancybox=True)

def hist_mdelays(d, i=0, m=10, bins_s=(0, 10, 0.02), period=0,
                 hold=False, bg_ppf=0.01, ph_sel=Ph_sel('all'), spline=True,
                 s=1., bg_fit=True, bg_F=0.8):
    """Histogram of m-photons delays (all-ph vs in-burst ph).
    """
    ax = gca()
    if not hold:
        #ax.clear()
        for _ind in range(len(ax.lines)): ax.lines.pop()

    results = bext.calc_mdelays_hist(
        d=d, ich=i, m=m, period=period, bins_s=bins_s,
        ph_sel=ph_sel, bursts=True, bg_fit=bg_fit, bg_F=bg_F)
    bin_x, histog_y = results[:2]
    bg_dist = results[2]
    rate_ch_kcps = 1./bg_dist.kwds['scale']  # extract the rate
    if bg_fit:
        a, rate_kcps = results[3:5]

    mdelays_hist_y = histog_y[0]
    mdelays_b_hist_y = histog_y[1]

    # Center of mass (COM)
    binw = bins_s[2]
    com = np.sum(bin_x*mdelays_hist_y)*binw
    com_b = np.sum(bin_x*mdelays_b_hist_y)*binw
    #print(com, com_b)

    # Compute a spline smoothing of the PDF
    mdelays_spline = UnivariateSpline(bin_x, mdelays_hist_y, s=s*com)
    mdelays_b_spline = UnivariateSpline(bin_x, mdelays_b_hist_y, s=s*com_b)
    mdelays_spline_y = mdelays_spline(bin_x)
    mdelays_b_spline_y = mdelays_b_spline(bin_x)
    if spline:
        mdelays_pdf_y = mdelays_spline_y
        mdelays_b_pdf_y = mdelays_b_spline_y
    else:
        mdelays_pdf_y = mdelays_hist_y
        mdelays_b_pdf_y = mdelays_b_hist_y

    # Thresholds and integrals
    max_delay_th_P = bg_dist.ppf(bg_ppf)
    max_delay_th_F = m/rate_ch_kcps/d.F

    burst_domain = bin_x < max_delay_th_F
    burst_integral = np.trapz(x=bin_x[burst_domain],
                              y=mdelays_hist_y[burst_domain])

    title("I = %.1f %%" % (burst_integral*100), fontsize='small')
    #text(0.8,0.8,"I = %.1f %%" % (integr*100), transform = gca().transAxes)

    ## MDelays plot
    plot(bin_x, mdelays_pdf_y, lw=2, color=blue, alpha=0.5,
         label="Delays dist.")
    plot(bin_x, mdelays_b_pdf_y, lw=2, color=red, alpha=0.5,
         label="Delays dist. (in burst)")
    plt.axvline(max_delay_th_P, color='k',
                label="BG ML dist. @ %.1f%%" % (bg_ppf*100))
    plt.axvline(max_delay_th_F, color=purple,
                label="BS threshold (F=%d)" % d.F)

    ## Bg distribution plots
    bg_dist_y = bg_dist.pdf(bin_x)
    ibin_x_bg_mean = np.abs(bin_x - bg_dist.mean()).argmin()
    bg_dist_y *= mdelays_pdf_y[ibin_x_bg_mean]/bg_dist_y[ibin_x_bg_mean]
    plot(bin_x, bg_dist_y, '--k', alpha=1.,
         label='BG ML dist.')
    plt.axvline(bg_dist.mean(), color='k', ls='--', label="BG mean")
    if bg_fit:
        bg_y = a*erlang.pdf(bin_x, a=m, scale=1./rate_kcps)
        plot(bin_x, bg_y, '--k', alpha=1.)
    plt.legend(ncol=2, frameon=False)
    xlabel("Time (ms)")


def hist_mrates(d, i=0, m=10, bins=(0, 4000, 100), yscale='log', pdf=False,
                dense=True, plot_style=None):
    """Histogram of m-photons rates. See also :func:`hist_mdelays`.
    """
    ph = d.get_ph_times(ich=i)
    if dense:
        ph_mrates = 1.*m/((ph[m-1:]-ph[:ph.size-m+1])*d.clk_p*1e3)
    else:
        ph_mrates = 1.*m/(np.diff(ph[::m])*d.clk_p*1e3)

    hist = HistData(*np.histogram(ph_mrates, bins=_bins_array(bins)))
    ydata = hist.pdf if pdf else hist.counts
    plot_style_ = dict(marker='o')
    plot_style_.update(_normalize_kwargs(plot_style, kind='line2d'))
    plot(hist.bincenters, ydata, **plot_style_)
    gca().set_yscale(yscale)
    xlabel("Rates (kcps)")

## Bursts stats
def hist_sbr(d, i=0, bins=(0, 30, 1), pdf=True, weights=None, color=None,
             plot_style=None):
    """Histogram of per-burst Signal-to-Background Ratio (SBR).
    """
    weights = weights[i] if weights is not None else None
    if 'sbr' not in d:
        d.calc_sbr()
    _hist_burst_taildist(d.sbr[i], bins, pdf, weights=weights, color=color,
                         plot_style=plot_style)
    plt.xlabel('SBR')


def hist_burst_phrate(d, i=0, bins=(0, 1000, 20), pdf=True, weights=None,
                      color=None, plot_style=None):
    """Histogram of max photon rate in each burst.
    """
    weights = weights[i] if weights is not None else None
    if 'max_rate' not in d:
        d.calc_max_rate(m=10)
    _hist_burst_taildist(d.max_rate[i]*1e-3, bins, pdf, weights=weights,
                         color=color, plot_style=plot_style)
    plt.xlabel('Peak rate (kcps)')


def hist_burst_delays(d, i=0, bins=(0, 10, 0.2), pdf=False, weights=None,
                      color=None, plot_style=None):
    """Histogram of waiting times between bursts.
    """
    weights = weights[i] if weights is not None else None
    bdelays = np.diff(d.mburst[i].start*d.clk_p)

    _hist_burst_taildist(bdelays, bins, pdf, weights=weights, color=color,
                         plot_style=plot_style)
    plt.xlabel('Delays between bursts (s)')

## Burst internal "symmetry"
def hist_asymmetry(d, i=0, bin_max=2, binwidth=0.1, stat_func=np.median):
    burst_asym = bext.asymmetry(d, ich=i, func=stat_func)
    bins_pos = np.arange(0, bin_max+binwidth, binwidth)
    bins = np.hstack([-bins_pos[1:][::-1], bins_pos])
    izero = (bins.size - 1)/2.
    assert izero == np.where(np.abs(bins) < 1e-8)[0]

    counts, _ = np.histogram(burst_asym, bins=bins)
    asym_counts_neg = counts[:izero] - counts[izero:][::-1]
    asym_counts_pos = counts[izero:] - counts[:izero][::-1]
    asym_counts = np.hstack([asym_counts_neg, asym_counts_pos])

    plt.bar(bins[:-1], width=binwidth, height=counts, fc=blue, alpha=0.5)
    plt.bar(bins[:-1], width=binwidth, height=asym_counts, fc=red,
            alpha=0.5)
    plt.grid(True)
    plt.xlabel('Time (ms)')
    plt.ylabel('# Bursts')
    plt.legend(['{func}$(t_D)$ - {func}$(t_A)$'.format(func=stat_func.__name__),
                'positive half - negative half'],
               frameon=False, loc='best')
    skew_abs = asym_counts_neg.sum()
    skew_rel = 100.*skew_abs/counts.sum()
    print('Skew: %d bursts, (%.1f %%)' % (skew_abs, skew_rel))

##
#  Scatter plots
#

def scatter_width_size(d, i=0):
    """Scatterplot of burst width versus size."""
    b = d.mburst[i]
    plot(b.width*d.clk_p*1e3, d.nt[i], 'o', mew=0, ms=3, alpha=0.7,
         color='blue')
    t_ms = arange(0, 50)
    plot(t_ms, ((d.m)/(d.T[i]))*t_ms*1e-3, '--', lw=2, color='k',
         label='Slope = m/T = min. rate = %1.0f cps' % (d.m/d.T[i]))
    plot(t_ms, d.bg_mean[Ph_sel('all')][i]*t_ms*1e-3, '--', lw=2, color=red,
         label='Noise rate: BG*t')
    xlabel('Burst width (ms)'); ylabel('Burst size (# ph.)')
    plt.xlim(0, 10); plt.ylim(0, 300)
    legend(frameon=False)

def scatter_rate_da(d, i=0):
    """Scatter of nd rate vs na rate (rates for each burst)."""
    b = d.mburst[i]
    Rate = lambda nX: nX[i]/b.width/d.clk_p*1e-3
    plot(Rate(d.nd), Rate(d.na), 'o', mew=0, ms=3, alpha=0.1, color='blue')
    xlabel('D burst rate (kcps)'); ylabel('A burst rate (kcps)')
    plt.xlim(-20, 100); plt.ylim(-20, 100)
    legend(frameon=False)

def scatter_fret_size(d, i=0, which='all', gamma=1, add_naa=False,
                      plot_style=None):
    """Scatterplot of FRET efficiency versus burst size.
    """
    if which == 'all':
        size = d.burst_sizes_ich(ich=i, gamma=gamma, add_naa=add_naa)
    else:
        assert which in d
        size = d[which][i]

    plot_style_ = dict(linestyle='', alpha=0.1, color=blue,
                       marker='o', markeredgewidth=0, markersize=3)
    plot_style_.update(_normalize_kwargs(plot_style, kind='line2d'))
    plot(d.E[i], size, **plot_style_)
    xlabel("FRET Efficiency (E)")
    ylabel("Corrected Burst size (#ph)")

def scatter_fret_nd_na(d, i=0, show_fit=False, no_text=False, gamma=1.,
                       **kwargs):
    """Scatterplot of FRET versus gamma-corrected burst size."""
    default_kwargs = dict(mew=0, ms=3, alpha=0.3, color=blue)
    default_kwargs.update(**kwargs)
    plot(d.E[i], gamma*d.nd[i]+d.na[i], 'o', **default_kwargs)
    xlabel("FRET Efficiency (E)")
    ylabel("Burst size (#ph)")
    if show_fit:
        _fitted_E_plot(d, i, F=1., no_E=no_text, ax=gca())
        if i == 0 and not no_text:
            plt.figtext(0.4, 0.01, _get_fit_E_text(d), fontsize=14)

def scatter_fret_width(d, i=0):
    """Scatterplot of FRET versus burst width."""
    b = d.mburst[i]
    plot(d.E[i], (b[:, 1]*d.clk_p)*1e3, 'o', mew=0, ms=3, alpha=0.1,
         color="blue")
    xlabel("FRET Efficiency (E)")
    ylabel("Burst width (ms)")

def scatter_da(d, i=0, alpha=0.3):
    """Scatterplot of donor vs acceptor photons (nd, vs na) in each burst."""
    plot(d.nd[i], d.na[i], 'o', mew=0, ms=3, alpha=alpha, color='blue')
    xlabel('# donor ph.'); ylabel('# acceptor ph.')
    plt.xlim(-5, 200); plt.ylim(-5, 120)

def scatter_naa_nt(d, i=0, alpha=0.5):
    """Scatterplot of nt versus naa."""
    plot(d.nt[i], d.naa[i], 'o', mew=0, ms=3, alpha=alpha, color='blue')
    plot(arange(200), color='k', lw=2)
    xlabel('Total burst size (nd+na+naa)'); ylabel('Accept em-ex BS (naa)')
    plt.xlim(-5, 200); plt.ylim(-5, 120)

def scatter_alex(d, i=0, **kwargs):
    """Scatterplot of E vs S. Keyword arguments passed to `plot`."""
    plot_style = dict(mew=1, ms=4, mec='black', color='purple',
                      alpha=0.1)
    plot_style = _normalize_kwargs(plot_style, 'line2d')
    plot_style.update(_normalize_kwargs(kwargs))
    plot(d.E[i], d.S[i], 'o', **plot_style)
    xlabel("E"); ylabel('S')
    plt.xlim(-0.2, 1.2); plt.ylim(-0.2, 1.2)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  High-level plot wrappers
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def _iter_plot(d, func, kwargs, iter_ch, nrows, ncols, figsize, AX,
               sharex, sharey, suptitle, grid, scale,
               top_adjust=0.95, hspace=0.15, left=0.08, right=0.96):
    if AX is None:
        fig, AX = plt.subplots(nrows, ncols, figsize=figsize, sharex=sharex,
                               sharey=sharey, squeeze=False)
        fig.subplots_adjust(left=0.08, right=0.96, top=top_adjust,
                            bottom=0.07, wspace=0.05)
        old_ax = False
    else:
        fig = AX[0, 0].figure
        old_ax = True

    for i, ich in enumerate(iter_ch):
        b = d.mburst[ich] if 'mburst' in d else None
        ax = AX.ravel()[i]
        if i == 0 and suptitle:
            fig.suptitle(d.status())
        s = '[%d]' % ich
        if 'bg_mean' in d:
            s += (' BG=%.1fk' % (d.bg_mean[Ph_sel('all')][ich] * 1e-3))
        if b is not None:
            s += (', #bu=%d' % b.num_bursts)
        ax.set_title(s)
        ax.grid(grid)
        plt.sca(ax)
        gui_status['first_plot_in_figure'] = (i == 0)
        func(d, ich, **kwargs)
        if ax.legend_ is not None:
            ax.legend_.remove()
    [a.set_xlabel('') for a in AX[:-1, :].ravel()]
    [a.set_ylabel('') for a in AX[:, 1:].ravel()]
    if sharex:
        plt.setp([a.get_xticklabels() for a in AX[:-1, :].ravel()],
                 visible=False)
        [a.set_xlabel('') for a in AX[:-1, :].ravel()]
        if not old_ax:
            fig.subplots_adjust(hspace=0.15)
    if sharey:
        if AX.shape[1] > 1:
            plt.setp([a.get_yticklabels() for a in AX[:, 1]], visible=False)
        fig.subplots_adjust(wspace=0.08)

        func_allows_autoscale = True
        if func.__name__ in _plot_status:
            func_allows_autoscale = _plot_status[func.__name__]['autoscale']
        if scale and func_allows_autoscale:
            ax.autoscale(enable=True, axis='y')
    return AX


def dplot_48ch(d, func, sharex=True, sharey=True, layout='horiz',
               grid=True, figsize=None, AX=None, suptitle=True,
               scale=True, top_adjust=0.95, **kwargs):
    """Plot wrapper for 48-spot measurements. Use `dplot` instead."""
    msg = "Wrong layout '%s'. Valid values: 'horiz', 'vert', '8x6'."
    assert (layout.startswith('vert') or layout.startswith('horiz') or
            layout == '8x6'), (msg % layout)
    global gui_status
    iter_ch = range(48)
    if layout == '8x6':
        nrows, ncols = 6, 8
        if figsize is None:
            figsize = (20, 16)
    else:
        nrows, ncols = 4, 12
        if layout.startswith('vert'):
            nrows, ncols = ncols, nrows
            iter_ch = np.arange(48).reshape(4, 12).T.ravel()
        if figsize is None:
            figsize = (20, 7)
            if layout.startswith('vert'):
                figsize = figsize[1], figsize[0]
    return _iter_plot(d, func, kwargs, iter_ch, nrows, ncols, figsize, AX,
                      sharex, sharey, suptitle, grid, scale,
                      top_adjust=top_adjust, hspace=0.15, left=0.08, right=0.96)


def dplot_16ch(d, func, sharex=True, sharey=True, ncols=8,
               pgrid=True, figsize=None, AX=None, suptitle=True,
               scale=True, top_adjust=0.93, **kwargs):
    """Plot wrapper for 16-spot measurements. Use `dplot` instead."""
    assert (ncols <= 16), '`ncols` needs to be <= 16.'
    global gui_status
    iter_ch = range(16)
    nrows = int(np.ceil(d.nch / ncols))
    if figsize is None:
        subplotsize = (3, 3)
        figsize = (subplotsize[0] * ncols, subplotsize[1] * nrows)
    return _iter_plot(d, func, kwargs, iter_ch, nrows, ncols, figsize, AX,
                      sharex, sharey, suptitle, grid, scale,
                      top_adjust=top_adjust, hspace=0.15, left=0.08, right=0.96)


def dplot_8ch(d, func, sharex=True, sharey=True,
              pgrid=True, figsize=(12, 9), nosuptitle=False, AX=None,
              scale=True, **kwargs):
    """Plot wrapper for 8-spot measurements. Use `dplot` instead."""
    global gui_status
    if AX is None:
        fig, AX = plt.subplots(4, 2, figsize=figsize, sharex=sharex,
                               sharey=sharey)
        fig.subplots_adjust(left=0.08, right=0.96, top=0.93, bottom=0.07,
                            wspace=0.05)
        old_ax = False
    else:
        fig = AX[0, 0].figure
        old_ax = True
    for i in range(d.nch):
        b = d.mburst[i] if 'mburst' in d else None
        if (func not in [timetrace, ratetrace, timetrace_single,
                         ratetrace_single, hist_bg_single, hist_bg,
                         timetrace_bg]) and np.size(b) == 0:
            continue
        ax = AX.ravel()[i]
        if i == 0 and not nosuptitle:
            fig.suptitle(d.status())
        s = u'[%d]' % (i+1)
        if 'bg_mean' in d:
            s += (' BG=%.1fk' % (d.bg_mean[Ph_sel('all')][i]*1e-3))
        if 'T' in d:
            s += (u', T=%dμs' % (d.T[i]*1e6))
        if b is not None: s += (', #bu=%d' %  b.num_bursts)
        ax.set_title(s, fontsize=12)
        ax.grid(pgrid)
        plt.sca(ax)
        gui_status['first_plot_in_figure'] = (i == 0)
        func(d, i, **kwargs)
        if i % 2 == 1: ax.yaxis.tick_right()
    [a.set_xlabel('') for a in AX[:-1, :].ravel()]
    [a.set_ylabel('') for a in AX[:, 1:].ravel()]
    if sharex:
        plt.setp([a.get_xticklabels() for a in AX[:-1, :].ravel()],
                 visible=False)
        [a.set_xlabel('') for a in AX[:-1, :].ravel()]
        if not old_ax: fig.subplots_adjust(hspace=0.15)
    if sharey:
        plt.setp([a.get_yticklabels() for a in AX[:, 1]], visible=False)
        fig.subplots_adjust(wspace=0.08)

        func_allows_autoscale = True
        if func.__name__ in _plot_status:
            func_allows_autoscale = _plot_status[func.__name__]['autoscale']
        if scale and func_allows_autoscale:
            ax.autoscale(enable=True, axis='y')
    return AX


def dplot_1ch(d, func, pgrid=True, ax=None,
              figsize=(9, 4.5), fignum=None, nosuptitle=False, **kwargs):
    """Plot wrapper for single-spot measurements. Use `dplot` instead."""
    global gui_status
    if ax is None:
        fig = plt.figure(num=fignum, figsize=figsize)
        ax = fig.add_subplot(111)
    else:
        fig = ax.figure
    s = d.name
    if 'bg_mean' in d:
        s += (' BG=%.1fk' % (d.bg_mean[Ph_sel('all')][0] * 1e-3))
    if 'T' in d:
        s += (u', T=%dμs' % (d.T[0] * 1e6))
    if 'mburst' in d:
        s += (', #bu=%d' % d.num_bursts[0])
    if not nosuptitle:
        ax.set_title(s, fontsize=12)
    ax.grid(pgrid)
    plt.sca(ax)
    gui_status['first_plot_in_figure'] = True
    func(d, **kwargs)
    return ax

def dplot(d, func, **kwargs):
    """Main plot wrapper for single and multi-spot measurements."""
    if d.nch == 1:
        return dplot_1ch(d=d, func=func, **kwargs)
    elif d.nch == 8:
        return dplot_8ch(d=d, func=func, **kwargs)
    elif d.nch == 16:
        return dplot_16ch(d=d, func=func, **kwargs)
    elif d.nch == 48:
        return dplot_48ch(d=d, func=func, **kwargs)

##
#  ALEX join-plot using seaborn
#
def _alex_plot_style(g, colorbar=True):
    """Set plot style and colorbar for an ALEX joint plot.
    """
    g.set_axis_labels(xlabel="E", ylabel="S")
    g.ax_marg_x.grid(True)
    g.ax_marg_y.grid(True)
    g.ax_marg_x.set_xlabel('')
    g.ax_marg_y.set_ylabel('')
    plt.setp(g.ax_marg_y.get_xticklabels(), visible=True)
    plt.setp(g.ax_marg_x.get_yticklabels(), visible=True)
    g.ax_marg_x.locator_params(axis='y', tight=True, nbins=3)
    g.ax_marg_y.locator_params(axis='x', tight=True, nbins=3)
    if colorbar:
        pos = g.ax_joint.get_position().get_points()
        X, Y = pos[:, 0], pos[:, 1]
        cax = plt.axes([1., Y[0], (X[1] - X[0]) * 0.045, Y[1] - Y[0]])
        plt.colorbar(cax=cax)

def _hist_bursts_marg(arr, dx, **kwargs):
    """Wrapper to call hist_burst_data() from seaborn plot_marginals().
    """
    vertical = kwargs.get('vertical', False)
    data_name = 'S' if vertical else 'E'
    hist_burst_data(dx, data_name=data_name, **kwargs)

def _alex_hexbin_vmax(patches, vmax_fret=True, Smax=0.8):
    """Return max the counts in the E-S hexbin histogram in `patches`.

    When `vmax_fret` is True, returns the max count for S < Smax.
    Otherwise returns the max count in all the histogram.
    """
    counts = patches.get_array()
    if vmax_fret:
        offset = patches.get_offsets()
        xoffset, yoffset = offset[:, 0], offset[:, 1]
        mask = yoffset < Smax
        vmax = counts[mask].max()
    else:
        vmax = counts.max()
    return vmax

def _calc_vmin(vmax, vmax_threshold, vmin_default):
    if vmax <= vmax_threshold:
        vmin = vmin_default - 0.5 * vmax
    elif vmax_threshold < vmax < 2 * vmax_threshold:
        vmin = vmin_default - 0.5 * vmax * ((2 * vmax_threshold - vmax) /
                                            (1 * vmax_threshold))
    else:
        vmin = vmin_default
    return vmin

def alex_jointplot(d, i=0, gridsize=50, cmap='Spectral_r', kind='hex',
                   vmax_fret=True, vmax_threshold=10,
                   vmin_default=0, vmin=None, cmap_compensate=False,
                   joint_kws=None, marginal_kws=None, histcolor_id=0,
                   rightside_text=False):
    """Plot an ALEX join plot: an E-S 2D histograms with marginal E and S.

    This function plots a jointplot: a main 2D histogram (hexbin plot)
    for E-S and the marginal histograms for E and S separately.
    The 2D histogram is an hexbin plot, i.e. the bin shape is hexagonal
    that has the advantage to reduce artifacts due to discretization.

    Arguments:
        d (Data object): the variable containing the bursts to plot
        i (int): the channel number. Default 0.
        gridsize (int): the grid size for the 2D histogram (hexbin)
        C (1D array or None): array of weights, it must have size equal to
            the number of bursts in channel `i` (d.num_bursts[i]).
            Passed to matplotlib hexbin().
        cmap (string): name of the colormap for the 2D histogram. In
            addition to matplotlib colormaps, FRETbursts defines
            these custom colormaps: 'alex_light', 'alex_dark' and 'alex_lv'.
            Default 'alex_light'.
        kind (string): kind of plot for the 2-D distribution. Valid values:
            'hex' for hexbin plots, 'kde' for kernel density estimation,
            'scatter' for scatter plot.
        vmax_fret (bool): if True, the colormap max value is equal to the
            max bin counts in the FRET region (S < 0.8). If False the
            colormap max is equal to the max bin counts.
        joint_kws (dict): keyword arguments passed to the function with plots
            the inner 2-D distribution (i.e matplotlib scatter or hexbin or
            seaborn kdeplot).
            and hence to matplolib hexbin to customize the plot style.
        marginal_kws (dict) : keyword arguments passed to :func:`hist_burst_`
            through seaborn plot_marginals() to customize the marginals plot
            style.
        histcolor_id (int): the colormap passes as `cmap` is divided in
            12 colors. `histcolor_id` is the index of the color to be used
            for the marginal 1D histogram. Default 1.
        rightside_text (bool): when True, print the measurement name on
            the right side of the figure. When False (default) no additional
            text is printed.

    .. seealso::
        The `Seaborn documentation <http://web.stanford.edu/~mwaskom/software/seaborn/index.html>`__
        has more info on plot customization:

        * http://web.stanford.edu/~mwaskom/software/seaborn/tutorial/axis_grids.html?highlight=jointgrid#plotting-bivariate-data-with-jointgrid-and-jointplot
        * http://web.stanford.edu/~mwaskom/software/seaborn/generated/seaborn.JointGrid.html

    """
    g = sns.JointGrid(x=d.E[i], y=d.S[i], ratio=3, space=0.2,
                      xlim=(-0.2, 1.2), ylim=(-0.2, 1.2))

    histcolor = sns.color_palette(cmap, 12)[histcolor_id]
    marginal_kws_ = dict(
        show_kde=True, bandwidth=0.03, binwidth=0.03,
        hist_bar_style={'facecolor': histcolor, 'edgecolor': 'k', 'linewidth': 0.2})
    if marginal_kws is not None:
        marginal_kws_.update(_normalize_kwargs(marginal_kws))

    if kind == "scatter":
        joint_kws_ = dict(s=40, color=histcolor, alpha=0.1, linewidths=0)
        if joint_kws is not None:
            joint_kws_.update(_normalize_kwargs(joint_kws))
        jplot = g.plot_joint(plt.scatter, **joint_kws_)
    elif kind.startswith('hex'):
        joint_kws_ = dict(edgecolor='none', linewidth=0.2, gridsize=gridsize,
                          cmap=cmap, extent=(-0.2, 1.2, -0.2, 1.2), mincnt=1)
        if joint_kws is not None:
            joint_kws_.update(_normalize_kwargs(joint_kws))
        jplot = g.plot_joint(plt.hexbin, **joint_kws_)

        # Set the vmin and vmax values for the colormap
        polyc = [c for c in jplot.ax_joint.get_children()
                 if isinstance(c, PolyCollection)][0]
        vmax = _alex_hexbin_vmax(polyc, vmax_fret=vmax_fret)
        if vmin is None:
            if cmap_compensate:
                # vmin = _calc_vmin(vmax, vmax_threshold, vmin_default)
                vmin = -vmax / 3
            else:
                vmin = vmin_default
        polyc.set_clim(vmin, vmax)
    elif kind.startswith("kde"):
        joint_kws_ = dict(shade=True, shade_lowest=False, n_levels=30,
                          cmap=cmap)
        if joint_kws is not None:
            joint_kws_.update(_normalize_kwargs(joint_kws))
        jplot = g.plot_joint(sns.kdeplot, **joint_kws_)

    g.plot_marginals(_hist_bursts_marg, dx=d, **marginal_kws_)
    g.annotate(lambda x, y: x.size, stat='# Bursts',
               template='{stat}: {val}')
    colorbar = kind.startswith('hex')
    _alex_plot_style(g, colorbar=colorbar)
    if rightside_text:
        plt.text(1.15, 0.6, d.name, transform=g.fig.transFigure, fontsize=14,
                 bbox=dict(edgecolor='r', facecolor='none', lw=1.3, alpha=0.5))
    return g

def _register_colormaps():
    import matplotlib as mpl
    import seaborn as sns

    c = sns.color_palette('nipy_spectral', 64)[2:43]
    cmap = mpl.colors.LinearSegmentedColormap.from_list('alex_lv', c)
    cmap.set_under(alpha=0)
    mpl.cm.register_cmap(name='alex_lv', cmap=cmap)

    c = sns.color_palette('YlGnBu', 64)[16:]
    cmap = mpl.colors.LinearSegmentedColormap.from_list('alex', c)
    cmap.set_under(alpha=0)
    mpl.cm.register_cmap(name='alex_light', cmap=cmap)
    mpl.cm.register_cmap(name='YlGnBu_crop', cmap=cmap)
    mpl.cm.register_cmap(name='alex_dark', cmap=mpl.cm.GnBu_r)

    # Temporary hack to workaround issue
    # https://github.com/mwaskom/seaborn/issues/855
    mpl.cm.alex_light = mpl.cm.get_cmap('alex_light')
    mpl.cm.alex_dark = mpl.cm.get_cmap('alex_dark')


# Register colormaps on import if not mocking
if not hasattr(sns, '_mock'):
    _register_colormaps()
