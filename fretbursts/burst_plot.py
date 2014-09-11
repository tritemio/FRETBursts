# encoding: utf-8
#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
This module defines all the plotting functions for the `Data` object.

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

"""

from __future__ import division
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
from matplotlib.collections import PatchCollection
#from matplotlib.collections import PathCollection
#from matplotlib.path import Path
#import matplotlib.patches as patches
#from matplotlib.lines import Line2D

# Local imports
from ph_sel import Ph_sel
import burstlib as bl
import burstlib_ext as bext
import background as bg
from utils.misc import clk_to_s, pprint
from scroll_gui import ScrollingToolQT
import gui_selection as gs


#ip = get_ipython()
#ip.magic("run -i scroll_gui.py")
#ip.magic("run -i gui_selection.py")
#ip.magic("run -i style.py")

params = {
        'font.size': 12,
        'legend.fontsize': 11,
        }
plt.rcParams.update(params)

##
#  Utility functions
#
def _normalize_kwargs(kwargs, kind='patch'):
    """Convert matplotlib keywords from short to long form."""
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
    plt.savefig(d.Name()+s)

##
#  Multi-channel plot functions
#

def mch_plot_bg(d, **kwargs):
    """Plot background vs channel for DA, D and A photons."""
    plot(r_[1:d.nch+1],[b.mean()*1e-3 for b in d.bg], lw=2, color='b',
            label=' T', **kwargs)
    plot(r_[1:d.nch+1],[b.mean()*1e-3 for b in d.bg_dd], color='g', lw=2,
            label=' D', **kwargs)
    plot(r_[1:d.nch+1],[b.mean()*1e-3 for b in d.bg_ad], color='r', lw=2,
            label=' A', **kwargs)
    xlabel("CH"); ylabel("kcps"); grid(True); legend(loc='best')
    title(d.name())

def mch_plot_bg_ratio(d):
    """Plot ratio of A over D background vs channel."""
    plot(r_[1:d.nch+1],[ba.mean()/bd.mean() for bd,ba in zip(d.bg_dd,d.bg_ad)],
            color='g', lw=2, label='A/D')
    xlabel("CH"); ylabel("BG Ratio A/D"); grid(True)
    title("BG Ratio A/D "+d.name())

def mch_plot_bsize(d):
    """Plot mean burst size vs channel."""
    CH = np.arange(1, d.nch+1)
    plot(CH, [b.mean() for b in d.nt], color='b', lw=2, label=' T')
    plot(CH, [b.mean() for b in d.nd], color='g', lw=2, label=' D')
    plot(CH, [b.mean() for b in d.na], color='r', lw=2, label=' A')
    xlabel("CH"); ylabel("Mean burst size")
    grid(True)
    legend(loc='best')
    title(d.name())


##
#  ALEX alternation period plots
#
def plot_alternation_hist(d, bins=100, **kwargs):
    plt.figure()
    ph_times_t, det_t, period = d.ph_times_t, d.det_t, d.alex_period
    d_ch, a_ch = d.det_donor_accept
    d_em_t = (det_t == d_ch)
    kwargs.update(bins=bins, alpha=0.2)
    hist(ph_times_t[d_em_t] % period, color='g', label='D', **kwargs)
    hist(ph_times_t[-d_em_t] % period, color='r', label='A', **kwargs)

    if d.D_ON[0] < d.D_ON[1]:
        plt.axvspan(d.D_ON[0], d.D_ON[1], color='g', alpha=0.1)
    else:
        plt.axvspan(0,d.D_ON[1], color='g', alpha=0.1)
        plt.axvspan(d.D_ON[0], period, color='g', alpha=0.1)

    if d.A_ON[0] < d.A_ON[1]:
        plt.axvspan(d.A_ON[0], d.A_ON[1], color='r', alpha=0.1)
    else:
        plt.axvspan(0, d.A_ON[1], color='r', alpha=0.1)
        plt.axvspan(d.A_ON[0], period, color='r', alpha=0.1)

    legend(loc='best')

def plot_alternation_hist_nsalex(d):
    nanotime_d = d.nanotime_t[d.det_t == d.det_donor_accept[0]]
    nanotime_a = d.nanotime_t[d.det_t == d.det_donor_accept[1]]
    hist(nanotime_d, bins=np.arange(4096), histtype='step', lw=1.2,
         alpha=0.5, color='g', label='Donor')
    hist(nanotime_a, bins=np.arange(4096), histtype='step', lw=1.2,
         alpha=0.5, color='r', label='Acceptor')
    plt.yscale('log')
    plt.axvspan(d.D_ON[0], d.D_ON[1], color='g', alpha=0.1)
    plt.axvspan(d.A_ON[0], d.A_ON[1], color='r', alpha=0.1)


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##  Multi-channel plots
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

##
#  Timetrace plots
#
def _plot_bursts(d, i, t_max_clk, pmax=1e3, pmin=0):
    """Highlights bursts in a timetrace plot."""
    b = d.mburst[i]
    if np.size(b) == 0: return
    pprint("CH %d burst..." % (i+1))
    bs = b[bl.b_start(b) < t_max_clk]
    start = bl.b_start(bs)*d.clk_p
    end = bl.b_end(bs)*d.clk_p
    R = []
    width = end-start
    ax = gca()
    for s,w in zip(start,width):
        r = Rectangle(xy=(s,pmin), height=pmax-pmin, width=w)
        r.set_clip_box(ax.bbox); r.set_zorder(0)
        R.append(r)
    ax.add_artist(PatchCollection(R, lw=0, color="#999999"))
    pprint("[DONE]\n")

def _gui_timetrace_burst_sel(d, i, func, fig, ax):
    """Add GUI burst selector via mouse click to the current plot."""
    if i == 0:
        func.burst_sel = gs.MultiAxPointSelection(fig, ax, d)
    else:
        func.burst_sel.ax_list.append(ax)

def _gui_timetrace_scroll(i, func, fig):
    """Add GUI to scroll a timetrace wi a slider."""
    if i == 0:
        func.scroll_gui = ScrollingToolQT(fig)


def _plot_rate_th(d, i, F, ph_sel, invert=False, bin_width=1,
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
    bg_dict = {Ph_sel('all'): d.bg[i],
               Ph_sel(Dex='Dem'): d.bg_dd[i],
               Ph_sel(Dex='Aem'): d.bg_ad[i],
               Ph_sel(Aex='Aem'): d.bg_aa[i],
               Ph_sel(Aex='Dem'): d.bg_da[i]}
    if ph_sel in bg_dict:
        rate_th_style_ = dict(plot_style_)
        rate_th_style_.update(linestyle='--', label='auto')
        rate_th_style_.update(_normalize_kwargs(rate_th_style, kind='line2d'))
        if rate_th_style_['label'] is 'auto':
            rate_th_style_['label'] = 'bg_rate*%d %s' % \
                                        (F, plot_style_['label'])
        x_rate = np.hstack(d.Ph_p[i])*d.clk_p
        y_rate = F*np.hstack([(rate, rate) for rate in bg_dict[ph_sel]])
        y_rate *= bin_width
        if invert:
            y_rate *= -1
        plot(x_rate, y_rate, **rate_th_style_)

def timetrace_single(d, i=0, bin_width=1e-3, bins=None, tmin=0, tmax=200,
                     ph_sel=Ph_sel('all'), invert=False, bursts=False,
                     burst_picker=True, scroll=False, cache_bins=True,
                     plot_style={}, show_rate_th=True, F=None,
                     rate_th_style={}):
    """Plot the timetrace (histogram) of timestamps for a photon selection.

    See :func:`timetrace` to plot multiple photon selections (i.e.
    Donor and Acceptor photons) in one step.
    """
    def _get_cache():
        return (timetrace_single.bins, timetrace_single.x,
                timetrace_single.bin_width,
                timetrace_single.tmin, timetrace_single.tmax)

    def _set_cache(bins, x, bin_width, tmin, tmax):
        cache = dict(bins=bins, x=x, bin_width=bin_width, tmin=tmin, tmax=tmax)
        for name, value in cache.items():
            setattr(timetrace_single, name, value)

    def _del_cache():
        names = ['bins', 'x', 'bin_width', 'tmin', 'tmax']
        for name in names:
            delattr(timetrace_single, name)

    def _has_cache():
        return hasattr(timetrace_single, 'bins')

    def _has_cache_for(bin_width, tmin, tmax):
        if _has_cache():
            return (bin_width, tmin, tmax) == _get_cache()[2:]
        return False

    # If cache_bins is False delete any previously saved attribute
    if not cache_bins and _has_cache:
        _del_cache()

    tmin_clk, tmax_clk = tmin/d.clk_p, tmax/d.clk_p
    bin_width_clk = bin_width/d.clk_p

    # If bins is not passed try to use the
    if bins is None:
        if cache_bins and _has_cache_for(bin_width, tmin, tmax):
            bins, x = timetrace_single.bins, timetrace_single.x
        else:
            bins = np.arange(tmin_clk, tmax_clk + 1, bin_width_clk)
            x = bins[:-1]*d.clk_p + 0.5*bin_width
            if cache_bins:
                _set_cache(bins, x, bin_width, tmin, tmax)

    # Compute histogram
    ph_times = d.get_ph_times(i, ph_sel=ph_sel)
    timetrace, _ = np.histogram(ph_times, bins=bins)
    if invert:
        timetrace *= -1

    # Plot bursts
    if bursts:
        t_max_clk = int(tmax/d.clk_p)
        _plot_bursts(d, i, t_max_clk, pmax=500, pmin=-500)

    # Plot timetrace
    color_dict = {Ph_sel('all'): 'k', Ph_sel(Dex='Dem'): 'g',
                  Ph_sel(Dex='Aem'): 'r', Ph_sel(Aex='Aem'): 'm',
                  Ph_sel(Aex='Dem'): 'c', }
    label_dict = {Ph_sel('all'): 'All-ph', Ph_sel(Dex='Dem'): 'DexDem',
                  Ph_sel(Dex='Aem'): 'DexAem', Ph_sel(Aex='Aem'): 'AexAem',
                  Ph_sel(Aex='Dem'): 'AexDem'}
    plot_style_ = dict(linestyle='-', linewidth=1.2, marker=None)
    if ph_sel in color_dict:
        plot_style_['color'] = color_dict[ph_sel]
        plot_style_['label'] = label_dict[ph_sel]
    plot_style_.update(_normalize_kwargs(plot_style, kind='line2d'))
    plot(x, timetrace, **plot_style_)

    # Plot burst-search rate-threshold
    if show_rate_th and 'bg' in d:
        _plot_rate_th(d, i, F=F, ph_sel=ph_sel, invert=invert,
                      bin_width=bin_width, plot_style_=plot_style_,
                      rate_th_style=rate_th_style)

    xlabel('Time (s)'); ylabel('# ph')#; plt.xlim(tmin, tmin + 1)
    if burst_picker:
        _gui_timetrace_burst_sel(d, i, timetrace_single, gcf(), gca())
    if scroll:
        _gui_timetrace_scroll(i, timetrace_single, gcf())

def timetrace(d, i=0, bin_width=1e-3, bins=None, tmin=0, tmax=200,
              bursts=False, burst_picker=True, scroll=False,
              show_rate_th=True, F=None, rate_th_style={'label': None},
              show_aa=True, legend=False,
              #dd_plot_style={}, ad_plot_style={}, aa_plot_style={}
              ):
    """Plot the timetraces (histogram) of photon timestamps.
    """
    # Plot bursts
    if bursts:
        t_max_clk = int(tmax/d.clk_p)
        _plot_bursts(d, i, t_max_clk, pmax=500, pmin=-500)

    # Plot multiple timetraces
    ph_sel_list = [Ph_sel(Dex='Dem'), Ph_sel(Dex='Aem')]
    invert_list = [False, True]
    if d.ALEX and show_aa:
         ph_sel_list.append(Ph_sel(Aex='Aem'))
         invert_list.append(True)

    for ph_sel, invert in zip(ph_sel_list, invert_list):
        if not bl.mask_empty(d.get_ph_mask(i, ph_sel=ph_sel)):
            timetrace_single(d, i, bin_width=bin_width, bins=bins, tmin=tmin,
                    tmax=tmax, ph_sel=ph_sel, invert=invert, bursts=False,
                    burst_picker=False, scroll=False, cache_bins=True,
                    show_rate_th=show_rate_th, F=F,
                    rate_th_style=rate_th_style)
    if legend:
        plt.legend(loc='best', fancybox=True)
    # Activate the burst picker
    if burst_picker:
        _gui_timetrace_burst_sel(d, i, timetrace, gcf(), gca())
    if scroll:
        _gui_timetrace_scroll(i, timetrace, gcf())

def ratetrace_single(d, i=0, m=None, max_num_ph=1e6, tmin=0, tmax=200,
                     ph_sel=Ph_sel('all'), invert=False, bursts=False,
                     burst_picker=True, scroll=False, plot_style={},
                     show_rate_th=True,  F=None, rate_th_style={}):
    """Plot the ratetrace of timestamps for a photon selection.

    See :func:`ratetrace` to plot multiple photon selections (i.e.
    Donor and Acceptor photons) in one step.
    """
    if m is None:
        m = d.m if m in d else 10

    # Compute ratetrace
    tmin_clk, tmax_clk = tmin/d.clk_p, tmax/d.clk_p
    ph_times = d.get_ph_times(i, ph_sel=ph_sel)
    iph1 = np.searchsorted(ph_times, tmin_clk)
    iph2 = np.searchsorted(ph_times, tmax_clk)
    if iph2 - iph1 > max_num_ph:
        iph2 = iph1 + max_num_ph
        tmax = ph_times[iph2]*d.clk_p
        warnings.warn(('Reached max number of photons, tmax reduced to %d s.',
                      '\nFor a wider time range increase `max_num_ph`'),
                      UserWarning)
    ph_times = ph_times[iph1:iph2]
    rates = bl.ph_rate(m, ph_times)/d.clk_p
    if invert:
        rates *= -1
    times = bl.ph_rate_t(m, ph_times)*d.clk_p

    # Plot ratetrace
    color_dict = {Ph_sel('all'): 'k', Ph_sel(Dex='Dem'): 'g',
                  Ph_sel(Dex='Aem'): 'r', Ph_sel(Aex='Aem'): 'm',
                  Ph_sel(Aex='Dem'): 'c', }
    label_dict = {Ph_sel('all'): 'All-ph', Ph_sel(Dex='Dem'): 'DexDem',
                  Ph_sel(Dex='Aem'): 'DexAem', Ph_sel(Aex='Aem'): 'AexAem',
                  Ph_sel(Aex='Dem'): 'AexDem'}
    plot_style_ = dict(linestyle='-', linewidth=1.2, marker=None)
    if ph_sel in color_dict:
        plot_style_['color'] = color_dict[ph_sel]
        plot_style_['label'] = label_dict[ph_sel]
    plot_style_.update(_normalize_kwargs(plot_style, kind='line2d'))
    plot(times, rates, **plot_style_)

    # Plot burst-search rate-threshold
    if show_rate_th and 'bg' in d:
        _plot_rate_th(d, i, F=F, ph_sel=ph_sel, invert=invert,
                      plot_style_=plot_style_, rate_th_style=rate_th_style)

    xlabel('Time (s)'); ylabel('# ph')#; plt.xlim(tmin, tmin + 1)
    if burst_picker:
        _gui_timetrace_burst_sel(d, i, ratetrace_single, gcf(), gca())
    if scroll:
        _gui_timetrace_scroll(i, ratetrace_single, gcf())

def ratetrace(d, i=0, m=None, max_num_ph=1e6, tmin=0, tmax=200,
              bursts=False, burst_picker=True, scroll=False,
              show_rate_th=True, F=None, rate_th_style={'label': None},
              show_aa=True, legend=False,
              #dd_plot_style={}, ad_plot_style={}, aa_plot_style={}
              ):
    """Plot the ratetraces of photon timestamps.
    """
    # Plot bursts
    if bursts:
        t_max_clk = int(tmax/d.clk_p)
        _plot_bursts(d, i, t_max_clk, pmax=1e6, pmin=-1e6)

    # Plot multiple timetraces
    ph_sel_list = [Ph_sel(Dex='Dem'), Ph_sel(Dex='Aem')]
    invert_list = [False, True]
    if d.ALEX and show_aa:
         ph_sel_list.append(Ph_sel(Aex='Aem'))
         invert_list.append(True)

    for ph_sel, invert in zip(ph_sel_list, invert_list):
        if not bl.mask_empty(d.get_ph_mask(i, ph_sel=ph_sel)):
            ratetrace_single(d, i, m=m, max_num_ph=max_num_ph, tmin=tmin,
                    tmax=tmax, ph_sel=ph_sel, invert=invert, bursts=False,
                    burst_picker=False, scroll=False,
                    show_rate_th=show_rate_th, F=F,
                    rate_th_style=rate_th_style)
    if legend:
        plt.legend(loc='best', fancybox=True)
    # Activate the burst picker
    if burst_picker:
        _gui_timetrace_burst_sel(d, i, ratetrace, gcf(), gca())
    if scroll:
        _gui_timetrace_scroll(i, ratetrace, gcf())

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

    style_kwargs = dict(marker='o', mew=0.5, color='b', mec='grey',
                        alpha=0.4, ls='')
    style_kwargs.update(**kwargs)

    t, E = bl.b_start(b)*d.clk_p, d.E[i]
    levels = sort_burst_sizes(bsizes)
    for ilev, level in enumerate(levels):
        plt.plot(t[level], E[level], ms=np.sqrt((ilev+1)*15),
                 **style_kwargs)
    plt.plot(bl.b_start(b)*d.clk_p, d.E[i], '-k', alpha=0.1, lw=1)
    xlabel('Time (s)'); ylabel('E')
    _gui_timetrace_burst_sel(d, i, timetrace_fret, gcf(), gca())

def timetrace_fret_scatter(d, i=0, gamma=1., **kwargs):
    """Timetrace of burst FRET vs time. Uses `scatter` (slow)."""
    b = d.mburst[i]
    bsizes = bl.select_bursts.get_burst_size(d, ich=i, gamma=gamma)

    style_kwargs = dict(s=bsizes, marker='o', alpha=0.5)
    style_kwargs.update(**kwargs)
    plt.scatter(bl.b_start(b)*d.clk_p, d.E[i], **style_kwargs)
    xlabel('Time (s)'); ylabel('E')

def timetrace_bg(d, i=0, nolegend=False, ncol=3):
    """Timetrace of background rates."""
    t = arange(d.bg[i].size)*d.bg_time_s
    plot(t, 1e-3*d.bg[i], 'k', lw=2, label="T: %d cps" % d.rate_m[i])
    plot(t, 1e-3*d.bg_dd[i], 'g', lw=2, label="DD: %d cps" % d.rate_dd[i])
    plot(t, 1e-3*d.bg_ad[i], 'r', lw=2, label="AD: %d cps" % d.rate_ad[i])
    if d.ALEX:
        plot(t, 1e-3*d.bg_aa[i], 'm', lw=2, label="AA: %d cps" % d.rate_aa[i])
    if not nolegend:
        legend(loc='best', frameon=False, ncol=ncol)
    xlabel("Time (s)"); ylabel("BG rate (kcps)"); grid(True)
    plt.ylim(ymin=0)

def timetrace_b_rate(d, i=0):
    """Timetrace of bursts-per-second in each period."""
    t = arange(d.bg[i].size)*d.bg_time_s
    b_rate = r_[[(d.bp[i] == p).sum() for p in xrange(d.bp[i].max()+1)]]
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
    ph_d = d.ph_times_m[i][SLICE][-d.A_em[i][SLICE]]
    ph_a = d.ph_times_m[i][SLICE][d.A_em[i][SLICE]]

    BSLICE = (bl.b_end(b) < ph_a[-1])
    start, end = bl.b_start(b[BSLICE]), bl.b_end(b[BSLICE])

    u = d.clk_p # time scale
    plt.vlines(ph_d*u, 0, 1, color='k', alpha=0.02)
    plt.vlines(ph_a*u, 0, 1, color='k', alpha=0.02)
    plt.vlines(start*u, -0.5, 1.5, lw=3, color='g', alpha=0.5)
    plt.vlines(end*u, -0.5, 1.5, lw=3, color='r', alpha=0.5)
    xlabel("Time (s)")


##
#  Histogram plots
#

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
                m1, s1, m2, s2, a1 =  d.fit_E_res[i, :]
                a2 = (1-a1)
            elif d.fit_E_res.shape[1] == 6:
                m1, s1, a1, m2, s2, a2 =  d.fit_E_res[i, :]
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
            print 'Fit Integral:', np.trapz(scale*y, x)

    ax2.axvline(d.E_fit[i], lw=3, color='r', ls='--', alpha=0.6)
    xtext = 0.6 if d.E_fit[i] < 0.6 else 0.2
    if d.nch > 1 and not no_E:
        ax2.text(xtext, 0.81,"CH%d: $E_{fit} = %.3f$" % \
                (i+1, d.E_fit[i]),
                transform = gca().transAxes, fontsize=16,
                bbox=dict(boxstyle='round', facecolor='#dedede', alpha=0.5))


def hist_width(d, i=0, bins=r_[0:10:0.025], yscale='log', density=True,
               **kwargs):
    b = d.mburst[i]
    histog, bins = np.histogram(bl.b_width(b)*d.clk_p*1e3, bins=bins,
                                density=density)
    #bins *= d.clk_p  # (ms)
    bins = bins[:-1]#+(bins[1]-bins[0])/2.
    plot_style = dict(color='red', lw=2)
    plot_style.update(**kwargs)
    plot(bins, histog, **plot_style)
    gca().set_yscale(yscale)
    #fill_between(bins,0,histog,color='red', alpha=0.5)
    xlabel('Burst width (ms)'); ylabel('# Burst')
    plt.xlim(xmin=0); plt.ylim(ymin=0)

def hist_size(d, i=0, vmax=600, binw=4, bins=None,
              which='all', gamma=1, add_naa=False,
              yscale='log', legend=True, plot_style={}):
    """Plot histogram of burst sizes.

    Parameters:
        d (Data): Data object
        i (int): channel index
        vmax (int/float): histogram max
        binw (int/float): histogram bin width
        bins (array or None): array of bin edges. If not NOne overrides `binw`
        which (string): which counts to consider. 'all' all-photon size
            computed with `d.burst_sizes()`; 'nd', 'na', 'naa' get counts from
            `d.nd`, `d.na`, `d.naa` (respectively Dex-Dem, Dex-Aem, Aex-Aem).
        yscale (string): 'log' or 'linear', sets the plot y scale.
        legend (bool): if True add legend to plot
        plot_style (dict): dict of matplotlib line style passed to `plot`.
    """
    valid_which = ["all", "nd", "na", "naa"]
    assert which in valid_which
    if which == 'all':
        size = d.burst_sizes_ich(ich=i, gamma=gamma, add_naa=add_naa)
        label = 'nd + na'
        if gamma != 1:
            label = "%.2f %s" % (gamma, label)
        if add_naa:
            label += " + naa"
    else:
        size = d[which][i]
        label = which

    colors = ['k', 'g', 'r', 'orange']
    colors_dict = {k: c for k, c in zip(valid_which, colors)}

    if bins is None:
        bins = np.arange(0, vmax+binw, binw)
    counts, bins = np.histogram(size, bins=bins)
    x = bins[:-1] + 0.5*(bins[1] + bins[0])
    plot_style_ = dict(linewidth=2, color=colors_dict[which])
    plot_style_.update(_normalize_kwargs(plot_style, kind='line2d'))
    plot(x, counts, label=label, **plot_style_)

    gca().set_yscale(yscale)
    xlabel('# Ph.'); ylabel('# Bursts')
    if legend: gca().legend(loc='best')

def hist_size_all(d, i=0, **kwargs):
    """Plot burst sizes for all the combinations of photons.

    Calls :func:`hist_size` multiple times with different `which` parameters.
    """
    fields = ['all', 'nd', 'na']
    if d.ALEX:
        fields.append('naa')
    for which in fields:
        hist_size(d, i, which=which, **kwargs)


def hist_fret(d, i=0, ax=None, binw=0.03, bins=None, pdf=True, hist_style='bar',
              weights=None, gamma=1., add_naa=False,            # weights args
              show_fit_stats=False, show_fit_value=False, fit_from='kde',
              show_kde=False, bandwidth=0.03, show_kde_peak=False,  # kde args
              show_model=False, show_model_peaks=True,
              hist_bar_style={}, hist_plot_style={}, model_plot_style={},
              kde_plot_style={}, verbose=False):

    """Plot FRET histogram and KDE.

    When `bins` is not None it overrides `binw (bin width).

    Histograms and KDE can be plotted on any Data variable after burst search.
    To show a model, a model must be fitted first by calling
    d.E_fitter.fit_histogram(). To show the KDE peaks position, they must be
    computed first with d.E_fitter.find_kde_max().
    """
    red = '#E41A1C'
    if ax is None:
        ax = gca()
    weights_tuple = (weights, float(gamma), add_naa)
    if not hasattr(d, 'E_fitter') or d.burst_weights != weights_tuple:
        if hasattr(d, 'E_fitter'):
            print ' - Overwriting the old E_fitter object with the new weights.'
            if verbose:
                print '   Old weights:', d.burst_weights
                print '   New weights:', weights_tuple
        bext.bursts_fitter(d, weights=weights, gamma=gamma, add_naa=add_naa)

    d.E_fitter.histogram(bin_width=binw, bins=bins, verbose=verbose)
    if pdf:
        ax.set_ylabel('PDF')
        hist_vals = d.E_fitter.hist_pdf[i]
    else:
        ax.set_ylabel('# Bursts')
        hist_vals = d.E_fitter.hist_counts[i]
    ax.set_xlabel('E')
    ax.set_xlim(-0.19, 1.19)

    hist_bar_style_ = dict(facecolor='#80b3ff', edgecolor='#5f8dd3',
                           linewidth=1.5, alpha=0.7, label='E Histogram')
    hist_bar_style_.update(**_normalize_kwargs(hist_bar_style))

    hist_plot_style_ = dict(linestyle='-', marker='o', markersize=6,
                            linewidth=2, alpha=0.6, label='E Histogram')
    hist_plot_style_.update(**_normalize_kwargs(hist_plot_style, kind='line2d'))
    if hist_style == 'bar':
        ax.bar(left = d.E_fitter.hist_bins[:-1], height=hist_vals,
                width = d.E_fitter.hist_bin_width, **hist_bar_style_)
    else:
        ax.plot(d.E_fitter.hist_axis, hist_vals, **hist_plot_style_)

    if show_model:
        model_plot_style_ = dict(color='k', alpha=0.8, label='Model')
        model_plot_style_.update(**_normalize_kwargs(model_plot_style,
                                                     kind='line2d'))
        fit_res = d.E_fitter.fit_res[i]
        x = d.E_fitter.x_axis
        ax.plot(x, fit_res.model.eval(x=x, **fit_res.values),
                 **model_plot_style_)
        if  fit_res.model.components is not None:
            for component in fit_res.model.components:
                model_plot_style_.update(ls = '--', label='Model component')
                ax.plot(x, component.eval(x=x, **fit_res.values),
                         **model_plot_style_)
        if show_model_peaks:
            for param in d.E_fitter.params:
                if param.endswith('center'):
                    ax.axvline(d.E_fitter.params[param][i], ls='--',
                                color=red)
    if show_kde:
        x = d.E_fitter.x_axis
        d.E_fitter.calc_kde(bandwidth=bandwidth)
        kde_plot_style_ = dict(linewidth=1.5, color='k', alpha=0.8, label='KDE')
        kde_plot_style_.update(**_normalize_kwargs(kde_plot_style,
                                                   kind='line2d'))
        ax.plot(x, d.E_fitter.kde[i](x), **kde_plot_style_)
    if show_kde_peak:
        ax.axvline(d.E_fitter.kde_max_pos[i], ls='--', color='orange')

    if show_fit_value or show_fit_stats:
        if fit_from == 'kde':
            fit_arr = d.E_fitter.kde_max_pos
        else:
            assert fit_from in d.E_fitter.params
            fit_arr = d.E_fitter.params[fit_from]

        if i == 0:
            if show_fit_stats:
                plt.figtext(0.4, 0.01, _get_fit_text_stats(fit_arr),
                            fontsize=16)
        if show_fit_value:
            _plot_fit_text_ch(fit_arr, i, ax=ax)

def _get_fit_text_stats(fit_arr, pylab=True):
    """Return a formatted string for mean E and max delta-E."""
    delta = (fit_arr.max() - fit_arr.min())*100
    fit_text = r'\langle{E}_{fit}\rangle = %.3f \qquad ' % fit_arr.mean()
    fit_text += r'\Delta E_{fit} = %.2f \%%' % delta
    if pylab: fit_text = r'$'+fit_text+r'$'
    return fit_text

def _plot_fit_text_ch(fit_arr, ich, fmt_str="CH%d: $E_{fit} = %.3f$", ax=None,
            bbox=dict(boxstyle='round', facecolor='#dedede', alpha=0.5),
            xtext_low=0.2, xtext_high=0.6, fontsize=16):
    """Plot a text box with ch and fit value."""
    if ax is None: ax = gca()
    xtext = xtext_high if fit_arr[ich] < xtext_high else xtext_low
    ax.text(xtext, 0.81, fmt_str % (ich+1, fit_arr[ich]),
            transform = ax.transAxes, fontsize=fontsize, bbox=bbox)

def hist_S(d, i=0, bins=None, binw=0.02, weights=None, gamma=1., normed=False,
           **kwargs):
    """Plot the Shoichiometry histogram and optionally the fitted model
    """
    if bins is None: bins = r_[-0.2:1.2:binw]
    style_kwargs = dict(bins=bins, normed=normed, histtype='stepfilled',
                        facecolor='#80b3ff', edgecolor='#5f8dd3',
                        linewidth=1.5, alpha=1)
    # kwargs overwrite style_kwargs
    style_kwargs.update(**_normalize_kwargs(kwargs))
    if weights is not None:
        w = bl.fret_fit.get_weights(d.nd[i], d.na[i], weights, gamma=gamma)
        w *= w.size/w.sum()
        style_kwargs.update(weights=w)
    hist(1.*d.S[i], **style_kwargs)
    xlabel('Stoichiometry'); ylabel('# Bursts')
    if normed: ylabel('PDF')
    plt.ylim(ymin=0); plt.xlim(-0.2,1.2)

def hist2d_alex(d, i=0, vmin=2, vmax=0, bin_step=None, S_max_norm=0.8,
                interp='bicubic', cmap='hot', under_color='white',
                over_color='white', scatter=True, scatter_ms=3,
                scatter_color='orange', scatter_alpha=0.2, gui_sel=False,
                ax=None, cbar_ax=None):
    """Plot 2-D E-S ALEX histogram with a scatterplot overlay.
    """
    if ax is None:
        ax = plt.gca()
    if bin_step is not None:
        d.calc_alex_hist(bin_step=bin_step)
    ES_hist, E_bins, S_bins, S_ax = d.ES_hist[i], d.E_bins, d.S_bins, d.S_ax

    colormap = plt.get_cmap(cmap)
    # Heuristic for colormap range
    if vmax <= vmin:
        S_range = (S_ax < S_max_norm)
        vmax = ES_hist[:, S_range].max()
        if vmax <= vmin: vmax = 10*vmin

    if scatter:
        ax.plot(d.E[i],d.S[i], 'o', mew=0, ms=scatter_ms, alpha=scatter_alpha,
                color=scatter_color)
    im = ax.imshow(ES_hist[:, ::-1].T, interpolation=interp,
            extent=(E_bins[0], E_bins[-1], S_bins[0], S_bins[-1]),
            vmin=vmin, vmax=vmax, cmap=colormap)
    im.cmap.set_under(under_color)
    im.cmap.set_over(over_color)
    if cbar_ax is None:
        gcf().colorbar(im)
    else:
        cbar_ax.colorbar(im)
    ax.set_xlim(-0.2, 1.2); ax.set_ylim(-0.2, 1.2)
    ax.set_xlabel('E');     ax.set_ylabel('S')
    ax.grid(color='gray')
    if gui_sel:
        # the selection object must be saved (otherwise will be destroyed)
        hist2d_alex.gui_sel = gs.rectSelection(gcf(), gca())


def plot_ES_selection(ax, E1, E2, S1, S2, rect=True, **kwargs):
    """Plot an overlay ROI on top of an E-S plot (i.e. ALEX histogram).

    This function plots a rectangle and inscribed ellipsis with x-axis limits
    (E1, E2) and y-axsis limits (S1, S2).

    Note that, a dict with keys (E1, E2, S1, S2, rect) can be both passed to
    :func:`fretbursts.select_bursts.ES` to apply a selection, and to
    `plot_ES_selection` to plot it.

    Parameters:
        ax (matplotlib axis): the axis where the rectangle is plotted.
            Typically you pass the axis of a previous E-S scatter plot
            or histogram.
        E1, E2, S1, S2 (floats): limits for E and S (X and Y axis respectively)
            used to plot the rectangle.
        rect (bool): if True, the rectangle is highlighted and the ellipsis is
            grey. The color are swapped otherwise.

    Any additional keyword argument specifies the matplotlib patch style
    for both the rectangle and the ellipsis.
    """
    if rect:
        rect_color, ellips_color = 'blue', 'gray'
    else:
        rect_color, ellips_color = 'gray', 'blue'
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
        print 'E1={E1:.3}, E2={E2:.3}, S1={S1:.3}, S2={S2:.3}'.format(**sel)
    return sel

def hist_sbr(d, ich=0, **hist_kwargs):
    """Histogram of per-burst Signal-to-Background Ratio (SBR).
    """
    if not 'sbr' in d:
        d.calc_sbr()
    style_kwargs = dict(bins=np.r_[0:20:0.5])  # default style
    style_kwargs.update(**hist_kwargs)
    hist(d.sbr[ich], **style_kwargs)
    xlabel('SBR'); ylabel('# Bursts')

def hist_bg_fit_single(d, i=0, bp=0, bg='bg_dd', bin_width_us=50, yscale='log',
        F=0.15, **kwargs):
    """Histog. of ph-delays compared with BG fitting in burst period 'bp'.
    """
    l1, l2 = d.Lim[i][bp][0], d.Lim[i][bp][1]
    ph = d.ph_times_m[i][l1:l2+1]*d.clk_p
    a_em = d.A_em[i][l1:l2+1]
    d_em = -a_em
    if bg=='bg': dph = np.diff(ph)
    if bg=='bg_dd': dph = np.diff(ph[d_em])
    if bg=='bg_ad': dph = np.diff(ph[a_em])
    hist_kwargs = dict(bins=r_[0:2000:bin_width_us], histtype='step',
                       color='k', lw=1)
    hist_kwargs.update(**kwargs)
    H = hist(dph*1e6, color='k', **hist_kwargs)
    gca().set_yscale('log')
    xlabel(u'Ph delay time (μs)')
    ylabel("# Ph")

    efun = lambda t, r: np.exp(-r*t)*r
    r = d[bg][i][bp]
    t = r_[0:2000]*1e-6
    #nF = 1 if 'normed' in hist_kwargs else H[0].sum()*(bin_width_us)

    bins = hist_kwargs['bins']
    ibin = bins.size*F
    t_min = bins[ibin]*1e-6
    C = H[0][ibin]/efun(t_min, r)
    plot(t*1e6, C*efun(t, r), lw=3, alpha=0.5, color='r',
        label="%s:  %d cps" % (bg, r))
    ym = 0.5
    if 'normed' in hist_kwargs and hist_kwargs['normed']: ym = 0.1/ph.size
    legend(loc='best', fancybox=True)
    plt.xlim(0,1500)
    plt.ylim(ymin=ym)

def hist_bg_fit(d, i=0, bp=0, bin_width_us=100, tmax=0.01, yscale='log',
                t_min_us=500, plot_style={}):
    """Histog. of ph-delays compared with BG fitting in burst period 'bp'.
    """
    # Compute histograms
    l1, l2 = d.Lim[i][bp][0], d.Lim[i][bp][1]
    ph = d.ph_times_m[i][l1:l2+1]*d.clk_p

    if d.ALEX:
        dd_mask = d.D_em[i][l1:l2+1]*d.D_ex[i][l1:l2+1]
        ad_mask = d.A_em[i][l1:l2+1]*d.D_ex[i][l1:l2+1]
        aa_mask = d.A_em[i][l1:l2+1]*d.A_ex[i][l1:l2+1]
    else:
        ad_mask = d.A_em[i][l1:l2+1]
        dd_mask = -ad_mask

    bins = np.arange(0, tmax*1e6, bin_width_us)
    t = bins[:-1] + 0.5*bin_width_us

    dph = np.diff(ph)
    dph_d, dph_a = np.diff(ph[dd_mask]), np.diff(ph[ad_mask])
    counts, _ = np.histogram(dph*1e6, bins=bins)
    counts_d, _  = np.histogram(dph_d*1e6, bins=bins)
    counts_a , _ = np.histogram(dph_a*1e6, bins=bins)
    if d.ALEX:
        counts_aa, _ = np.histogram(np.diff(ph[aa_mask])*1e6, bins=bins)

    # Plot histograms
    plot_style_ = dict(marker='o', markersize=5, linestyle='none', alpha=0.6)
    plot_style_.update(_normalize_kwargs(plot_style, kind='line2d'))
    plot(t, counts, color='k', label='T', **plot_style_)
    plot(t, counts_d, color='g', label='D', **plot_style_)
    plot(t, counts_a, color='r', label='A', **plot_style_)
    if d.ALEX:
        plot(t, counts_aa, color='m', label='AA', **plot_style_)

    # Plot fit
    r, rd, ra, = d.bg[i][bp], d.bg_dd[i][bp], d.bg_ad[i][bp]
    if d.ALEX:
        raa = d.bg_aa[i][bp]

    i_th = np.searchsorted(t, t_min_us)

    n_t = np.trim_zeros(counts).size
    n_td = np.trim_zeros(counts_d).size
    n_ta = np.trim_zeros(counts_a).size
    C = counts[i_th]/np.exp(-t[i_th]*1e-6*r)
    Cd = counts_d[i_th]/np.exp(-t[i_th]*1e-6*rd)
    Ca = counts_a[i_th]/np.exp(-t[i_th]*1e-6*ra)

    fit_style = dict(linewidth=4, alpha=0.5)
    plot(t[:n_t], C*np.exp(-t[:n_t]*1e-6*r), color='k',
         label="T:  %d cps" % r, **fit_style)
    plot(t[:n_td], Cd*np.exp(-t[:n_td]*1e-6*rd), color='g',
         label="DD:  %d cps" % rd, **fit_style)
    plot(t[:n_ta], Ca*np.exp(-t[:n_ta]*1e-6*ra), color='r',
         label="AD:  %d cps" % ra, **fit_style)
    if d.ALEX:
        n_taa = np.trim_zeros(counts_aa).size
        Caa = counts_aa[i_th]/np.exp(-t[i_th]*1e-6*raa)
        plot(t[:n_taa], Caa*np.exp(-t[:n_taa]*1e-6*raa), color='m',
             label="AA:  %d cps" % raa, **fit_style)
    if yscale == 'log':
        gca().set_yscale(yscale)
        plt.ylim(1)
    xlabel(u'Ph delay time (μs)'); ylabel("# Ph")
    plt.legend(loc='best', fancybox=True)

def hist_ph_delays(d, i=0, time_min_s=0, time_max_s=30, bin_width_us=10,
        mask=None, yscale='log', hfit_bin_ms=1, efit_tail_min_us=1000,
        **kwargs):
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
    plot(t*1e6, 0.65*F*efun(t,rc)*1e-6, lw=3, alpha=0.5, color='m',
            label="%d cps - Exp CDF (tail_min_p=%.2f)" % (rc, efit_tail_min_us))
    plot(t*1e6, 0.65*F*efun(t,re)*1e-6, lw=3, alpha=0.5, color='r',
            label="%d cps - Exp ML (tail_min_p=%.2f)" % (re, efit_tail_min_us))
    plot(t*1e6, 0.68*F*efun(t,rg)*1e-6, lw=3, alpha=0.5, color='g',
            label=u"%d cps - Hist (bin_ms=%d) [Δ=%d%%]" % (hfit_bin_ms, rg,
                                                           100*(rg-re)/re))
    plt.legend(loc='best', fancybox=True)

def hist_mdelays(d, i=0, m=10, bins_s=(0, 10, 0.02), bp=0, no_bg_fit=True,
                 hold=False, bg_ppf=0.01, ph_sel=Ph_sel('all'), spline=True,
                 s=1., bg_fit=True, bg_F=0.8):
    """Histogram of m-ph delays (all ph vs in-burst ph)."""
    ax = gca()
    if not hold:
        #ax.clear()
        for _ind in range(len(ax.lines)): ax.lines.pop()

    results = bext.calc_mdelays_hist(
                        d=d, ich=i, m=m, bp=bp, bins_s=bins_s, ph_sel=ph_sel,
                        bursts=True, bg_fit=bg_fit, bg_F=bg_F)
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
    #print com, com_b

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
    plot(bin_x, mdelays_pdf_y, lw=2, color='b', alpha=0.5,
         label="Delays dist.")
    plot(bin_x, mdelays_b_pdf_y, lw=2, color='r', alpha=0.5,
         label="Delays dist. (in burst)")
    plt.axvline(max_delay_th_P, color='k',
                label="BG ML dist. @ %.1f%%" % (bg_ppf*100))
    plt.axvline(max_delay_th_F, color='m',
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

def hist_mrates(d, i=0, m=10, bins=r_[0:20e3:20], yscale='log', normed=True,
        dense=True):
    """Histogram of m-photons rates."""
    ph = d.ph_times_m[i]
    if dense:
        ph_mrates = 1.*m/((ph[m-1:]-ph[:ph.size-m+1])*d.clk_p*1e3)
    else:
        ph_mrates = 1.*m/(np.diff(ph[::m])*d.clk_p*1e3)
    #gauss = lambda M,std: gaussian(M,std)/gaussian(M,std).sum()
    H = np.histogram(ph_mrates, bins=bins, normed=normed)
    epdf_x = H[1][:-1]
    epdf_y = H[0]
    plot(epdf_x, epdf_y, '.-')
    #plot(epdf_x, convolve(epdf_y, gauss(30,4),'same'), lw=2)
    gca().set_yscale(yscale)
    xlabel("Rates (kcps)")

## Bursts stats
def hist_rate_in_burst(d, i=0, bins=20):
    """Histogram of total photon rate in each burst."""
    b = d.mburst[i]
    rate = 1e-3*d.nt[i]/(bl.b_width(b)*d.clk_p)
    hist(rate, bins=bins, color="blue")
    xlabel('In-burst ph. rate (kcps)'); ylabel('# Bursts')
    #xlim(xmin=d.L/2); ylim(ymin=0)

def hist_burst_delays(d, i=0, tmax_seconds=0.5, bins=100, **kwargs):
    """Histogram of waiting times between bursts."""
    b = d.mburst[i]
    bd = clk_to_s(np.sort(np.diff(b[:,0].astype(float))))[:-20]
    hist(bd[bd<tmax_seconds], bins=bins, **kwargs)
    xlabel('Delays between bursts (s)'); ylabel('# bursts')

## Burst internal "symmetry"
def hist_asymmetry(d, i=0, bin_max=2, binw=0.1, func=np.median):
    burst_asym = bext.asymmetry(d, ich=i, func=func)
    bins_pos = np.arange(0, bin_max+binw, binw)
    bins = np.hstack([-bins_pos[1:][::-1], bins_pos])
    izero = (bins.size - 1)/2.
    assert izero == np.where(np.abs(bins) < 1e-8)[0]

    counts, _ = np.histogram(burst_asym, bins=bins)
    asym_counts_neg = counts[:izero] - counts[izero:][::-1]
    asym_counts_pos = counts[izero:] - counts[:izero][::-1]
    asym_counts = np.hstack([asym_counts_neg, asym_counts_pos])

    plt.bar(bins[:-1], width=binw, height=counts, fc='b', alpha=0.5)
    plt.bar(bins[:-1], width=binw, height=asym_counts, fc='r', alpha=0.5)
    plt.grid(True)
    plt.xlabel('Time (ms)')
    plt.ylabel('# Bursts')
    plt.legend(['{func}$(t_D)$ - {func}$(t_A)$'.format(func=func.__name__),
                'positive half - negative half'],
                frameon=False, loc='best')
    skew_abs = asym_counts_neg.sum()
    skew_rel = 100.*skew_abs/counts.sum()
    print 'Skew: %d bursts, (%.1f %%)' % (skew_abs, skew_rel)

##
#  Scatter plots
#

def scatter_width_size(d, i=0):
    """Scatterplot of burst width versus size."""
    b = d.mburst[i]
    plot(bl.b_width(b)*d.clk_p*1e3, d.nt[i], 'o', mew=0, ms=3, alpha=0.7,
         color='blue')
    t_ms = arange(0,50)
    plot(t_ms,((d.m)/(d.T[i]))*t_ms*1e-3,'--', lw=2, color='k',
            label='Slope = m/T = min. rate = %1.0f cps' % (d.m/d.T[i]))
    plot(t_ms,d.rate_m[i]*t_ms*1e-3,'--', lw=2, color='r',
            label='Noise rate: BG*t')
    xlabel('Burst width (ms)'); ylabel('Burst size (# ph.)')
    plt.xlim(0,10); plt.ylim(0,300)
    legend(frameon=False)

def scatter_rate_da(d, i=0):
    """Scatter of nd rate vs na rate (rates for each burst)."""
    b = d.mburst[i]
    Rate = lambda nX: nX[i]/bl.b_width(b)/d.clk_p*1e-3
    plot(Rate(d.nd), Rate(d.na), 'o', mew=0, ms=3, alpha=0.1, color='blue')
    xlabel('D burst rate (kcps)'); ylabel('A burst rate (kcps)')
    plt.xlim(-20,100); plt.ylim(-20,100)
    legend(frameon=False)

def scatter_fret_size(d, i=0, which='all', gamma=1, add_naa=False,
                      plot_style={}):
    """Scatterplot of FRET efficiency versus burst size.
    """
    if which == 'all':
        size = d.burst_sizes_ich(ich=i, gamma=gamma, add_naa=add_naa)
    else:
        assert which in d
        size = d[which][i]

    plot_style_ = dict(linestyle='', alpha=0.1, color='b',
                       marker='o', markeredgewidth=0, markersize=3)
    plot_style_.update(_normalize_kwargs(plot_style, kind='line2d'))
    plot(d.E[i], size, **plot_style_)
    xlabel("FRET Efficiency (E)")
    ylabel("Corrected Burst size (#ph)")

def scatter_fret_nd_na(d, i=0, show_fit=False, no_text=False, gamma=1.,
                       **kwargs):
    """Scatterplot of FRET versus gamma-corrected burst size."""
    default_kwargs = dict(mew=0, ms=3, alpha=0.3, color="blue")
    default_kwargs.update(**kwargs)
    plot(d.E[i], gamma*d.nd[i]+d.na[i], 'o', **default_kwargs)
    xlabel("FRET Efficiency (E)")
    ylabel("Burst size (#ph)")
    if show_fit:
        _fitted_E_plot(d, i, F=1., no_E=no_text, ax=gca())
        if i==0 and not no_text:
            plt.figtext(0.4,0.01, _get_fit_E_text(d),fontsize=14)

def scatter_fret_width(d, i=0):
    """Scatterplot of FRET versus burst width."""
    b = d.mburst[i]
    plot(d.E[i],(b[:,1]*d.clk_p)*1e3, 'o', mew=0, ms=3, alpha=0.1,
         color="blue")
    xlabel("FRET Efficiency (E)")
    ylabel("Burst width (ms)")

def scatter_da(d, i=0, alpha=0.3):
    """Scatterplot of donor vs acceptor photons (nd, vs na) in each burst."""
    plot(d.nd[i], d.na[i],'o', mew=0,ms=3, alpha=alpha, color='blue')
    xlabel('# donor ph.'); ylabel('# acceptor ph.')
    plt.xlim(-5,200); plt.ylim(-5,120)

def scatter_naa_nt(d, i=0, alpha=0.5):
    """Scatterplot of nt versus naa."""
    plot(d.nt[i], d.naa[i],'o', mew=0,ms=3, alpha=alpha, color='blue')
    plot(arange(200), color='k', lw=2)
    xlabel('Total burst size (nd+na+naa)'); ylabel('Accept em-ex BS (naa)')
    plt.xlim(-5,200); plt.ylim(-5,120)

def scatter_alex(d, i=0, **kwargs):
    """Scatterplot of E vs S. Keyword arguments passed to `plot`."""
    plot_style = dict(mew=1, ms=4, mec='black', color='purple',
                      alpha=0.1)
    plot_style = _normalize_kwargs(plot_style, 'line2d')
    plot_style.update(**_normalize_kwargs(kwargs))
    plot(d.E[i], d.S[i], 'o', **plot_style)
    xlabel("E"); ylabel('S')
    plt.xlim(-0.2,1.2); plt.ylim(-0.2,1.2)


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  High-level plot wrappers
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def dplot_48ch(d, fun=scatter_width_size, sharex=True, sharey=True,
        pgrid=True, figsize=None, AX=None, nosuptitle=False,
        scale=True, ich=None, **kwargs):
    """Plot wrapper for 48-spot measurements. Use `dplot` instead."""
    # Some function need an index of the number of calls so they are sure
    # that when idx == 0 it is thet first call of a series
    idx_funcs = [timetrace, timetrace_single, ratetrace, ratetrace_single]

    if ich is None:
        iter_ch = xrange(d.nch)
        if d.nch == 48:
            top_adjust = 0.95
            ax_ny, ax_nx = 6, 8
            if figsize is None:
                figsize = (20, 16)
        elif d.nch == 8:
            top_adjust = 0.93
            ax_ny, ax_nx = 4, 2
            if figsize is None:
                figsize = (12, 9)
    else:
        top_adjust = 0.9
        iter_ch = [ich]
        ax_ny, ax_nx = 1, 1
        if figsize is None:
            figsize = (8, 5)

    if AX is None:
        fig, AX = plt.subplots(ax_ny, ax_nx, figsize=figsize, sharex=sharex,
                               sharey=sharey, squeeze=False)
        fig.subplots_adjust(left=0.08, right=0.96, top=top_adjust,
                            bottom=0.07, wspace=0.05)
        old_ax = False
    else:
        fig = AX[0,0].figure
        old_ax = True

    for i, ich in enumerate(iter_ch):
        b = d.mburst[ich] if 'mburst' in d else None
        ax = AX.ravel()[i]
        if i == 0 and not nosuptitle:
            fig.suptitle(d.status())
        s = u'[%d]' % (ich+1)
        if 'rate_m' in d: s += (' BG=%.1fk' % (d.rate_m[ich]*1e-3))
        if b is not None: s += (', #bu=%d' %  b.shape[0])
        ax.set_title(s, fontsize=12)
        ax.grid(pgrid)
        plt.sca(ax)
        if fun in idx_funcs: kwargs.update(idx=i)
        fun(d, ich, **kwargs)
    [a.set_xlabel('') for a in AX[:-1,:].ravel()]
    [a.set_ylabel('') for a in AX[:,1:].ravel()]
    if sharex:
        plt.setp([a.get_xticklabels() for a in AX[:-1,:].ravel()], visible=False)
        [a.set_xlabel('') for a in AX[:-1,:].ravel()]
        if not old_ax: fig.subplots_adjust(hspace=0.15)
    if sharey:
        if AX.shape[1] > 1:
             plt.setp([a.get_yticklabels() for a in AX[:, 1]], visible=False)
        fig.subplots_adjust(wspace=0.08)
        if scale: ax.autoscale(enable=True, axis='y')
    return AX

def dplot_8ch(d, fun=scatter_width_size, sharex=True, sharey=True,
        pgrid=True, figsize=(12, 9), nosuptitle=False, AX=None,
        scale=True, **kwargs):
    """Plot wrapper for 8-spot measurements. Use `dplot` instead."""
    if AX is None:
        fig, AX = plt.subplots(4,2,figsize=figsize, sharex=sharex,
                               sharey=sharey)
        fig.subplots_adjust(left=0.08, right=0.96, top=0.93, bottom=0.07,
                wspace=0.05)
        old_ax = False
    else:
        fig = AX[0,0].figure
        old_ax = True
    for i in xrange(d.nch):
        b = d.mburst[i] if 'mburst' in d else None
        if (not fun in [timetrace, ratetrace, hist_bg_fit_single, hist_bg_fit,
            timetrace_bg]) and np.size(b) == 0:
            continue
        ax = AX.ravel()[i]
        if i == 0 and not nosuptitle:
            fig.suptitle(d.status())
        s = u'[%d]' % (i+1)
        if 'rate_m' in d: s += (' BG=%.1fk' % (d.rate_m[i]*1e-3))
        if 'T' in d: s += (u', T=%dμs' % (d.T[i]*1e6))
        if b is not None: s += (', #bu=%d' %  b.shape[0])
        ax.set_title(s, fontsize=12)
        ax.grid(pgrid)
        plt.sca(ax)
        fun(d, i, **kwargs)
        if i % 2 == 1: ax.yaxis.tick_right()
    [a.set_xlabel('') for a in AX[:-1,:].ravel()]
    [a.set_ylabel('') for a in AX[:,1:].ravel()]
    if sharex:
        plt.setp([a.get_xticklabels() for a in AX[:-1,:].ravel()], visible=False)
        [a.set_xlabel('') for a in AX[:-1,:].ravel()]
        if not old_ax: fig.subplots_adjust(hspace=0.15)
    if sharey:
        plt.setp([a.get_yticklabels() for a in AX[:,1]], visible=False)
        fig.subplots_adjust(wspace=0.08)
        if scale: ax.autoscale(enable=True, axis='y')
    return AX

def dplot_1ch(d, fun, pgrid=True, ax=None,
              figsize=(9, 4.5), fignum=None, nosuptitle=False, **kwargs):
    """Plot wrapper for single-spot measurements. Use `dplot` instead."""
    if ax is None:
        fig = plt.figure(num=fignum, figsize=figsize)
        ax = fig.add_subplot(111)
    else:
        fig = ax.figure
    s = d.name()
    if 'rate_m' in d: s += (' BG=%.1fk' % (d.rate_m[0]*1e-3))
    if 'T' in d: s += (u', T=%dμs' % (d.T[0]*1e6))
    if 'mburst' in d: s += (', #bu=%d' %  d.num_bursts()[0])
    if not nosuptitle: ax.set_title(s, fontsize=12)
    ax.grid(pgrid)
    plt.sca(ax)
    fun(d, **kwargs)
    return ax

def dplot(d, fun, **kwargs):
    """Main plot wrapper for single and multi-spot measurements."""
    if d.nch == 1:
        return dplot_1ch(d=d, fun=fun, **kwargs)
    elif d.nch == 8:
        return dplot_8ch(d=d, fun=fun, **kwargs)
    elif d.nch == 48:
        return dplot_48ch(d=d, fun=fun, **kwargs)

##
#  Other plot wrapper functions
#

def wplot(*args, **kwargs):
    AX, s = dplot_8ch(*args, **kwargs)
    kwargs.update(AX=AX)
    q = gs.mToolQT(gcf(), dplot_8ch, *args, **kwargs)
    return AX, q

def splot(d, fun=scatter_width_size,
        scroll=False, pgrid=True, figsize=(10, 8), nosuptitle=False, ax=None,
        scale=True, **kwargs):
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
    else:
        fig = ax.figure
    for i in xrange(d.nch):
        try:
            if i == 0 and not nosuptitle: fig.set_title(d.status())
        except:
            print "WARNING: No title in plots."
        ax.grid(pgrid)
        fun(d, i, **kwargs)
    s = None
    if scroll: s = ScrollingToolQT(fig)
    return ax, s

##
#  Other misc plot functions
#
def bplot(d, ich, b_index,  ph0=True, pad=0):
    """Plot photons in a burst as vertical bars. Burst: d.mburst[ich][b_index].
    """
    br = bl.b_irange(d.mburst[ich], b_index, pad=pad)
    accept = d.A_em[ich][br]
    donor = -accept
    ph = d.ph_times_m[ich][br]
    if ph0: ph -= ph[0]
    dtime = (ph[donor])*d.clk_p*1e6
    atime = (ph[accept])*d.clk_p*1e6
    plt.vlines(dtime,0,1, lw=2, color='g', alpha=0.8)
    plt.vlines(atime,0,1, lw=2, color='r', alpha=0.8)
    #plot(dtime, ones_like(ph[donor]), '^', color='g', alpha=0.5)
    #plot(atime, -ones_like(ph[accept]), 'o', color='r', alpha=0.5)
    xlabel("Time (us)")
    nd, na, nt = donor.sum(), accept.sum(), ph.size
    E = float(na)/(nt)
    title("#ph = %d, #D-ph = %d, #A-ph = %d, E = %.2f" % (nt,nd,na,E))
    plt.ylim(-10,10)


def bg_legend_8ch(d):
    ax = gca()
    L = ax.get_lines()[1::2]
    for i,l in enumerate(L):
        ich = i/3
        x = i%3
        s = ['Tot', 'D', 'A']
        r = [d.rate_m[ich], d.rate_dd[ich], d.rate_ad[ich]]
        l.set_label("CH%d, %s %d cps" % (ich+1, s[x], r[x]))
    ax.legend()
    plt.draw()

def bg_legend_alex(d):
    ax = gca()
    L = ax.get_lines()[1::2]
    for i,l in enumerate(L):
        ich = i/4
        x = i%4
        s = ['Tot', 'DD', 'AD', 'AA']
        r = [d.rate_m[ich], d.rate_dd[ich], d.rate_ad[ich], d.rate_aa[ich]]
        l.set_label("CH%d, %s %d cps" % (ich+1, s[x], r[x]))
    ax.legend()
    plt.draw()
