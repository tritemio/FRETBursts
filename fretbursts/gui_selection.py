#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
GUI selection of a range in a matplotlib figure.

Used by `hist2d_alex() and other functions in `burst_plot.py.
"""

from __future__ import absolute_import
from builtins import range, zip

import numpy as np
from scipy.stats import mode
from matplotlib.patches import Rectangle, Ellipse

from .utils.misc import pprint
from . import burstlib_ext as bext


class GuiSelection(object):
    """Abstract class for range selection in a matplotlib axis.

    Methods on_press_draw(), on_motion_draw() and on_release_print() must
    be overloaded by children classes.
    """
    def __init__(self, fig, ax, debug=False):
        self.ax = ax
        self.fig = fig
        self.pressed = False
        self.debug = debug
        self.id_press = fig.canvas.mpl_connect('button_press_event',
                                                self.on_press)
        if self.debug:
            pprint('Figure: ' + str(fig) + '\nAxis: ' + str(ax) + '\n')

    def on_press(self, event):
        if event.inaxes != self.ax: return
        self.pressed = True
        self.xs, self.ys = event.xdata, event.ydata
        if self.debug:
            pprint('PRESS button=%d, x=%d, y=%d, xdata=%f, ydata=%f\n' % (
                event.button, event.x, event.y, event.xdata, event.ydata))
        self.on_press_draw()
        self.fig.canvas.draw()
        self.id_motion = self.fig.canvas.mpl_connect('motion_notify_event',
                                                     self.on_motion)
        self.fig.canvas.mpl_connect('button_release_event',
                                             self.on_release)

    def on_motion(self, event):
        if event.inaxes != self.ax: return
        if self.debug:
            pprint('MOTION x=%d, y=%d, xdata=%f, ydata=%f\n' % (
                event.x, event.y, event.xdata, event.ydata))
        self.xe, self.ye = event.xdata, event.ydata
        self.on_motion_draw()
        self.fig.canvas.draw()

    def on_release(self, event):
        if not self.pressed: return
        self.pressed = False
        if self.debug:
            pprint('RELEASE button=%d, x=%d, y=%d, xdata=%f, ydata=%f\n' % (
                event.button, event.x, event.y, event.xdata, event.ydata))
        self.fig.canvas.mpl_disconnect(self.id_motion)
        self.on_release_print()

    def on_press_draw(self):
        pass

    def on_motion_draw(self):
        pass

    def on_release_print(self):
        pass


class xspanSelection(GuiSelection):
    """Interactive selection of an x range on the axis"""
    def on_press_draw(self):
        if 'r' in self.__dict__:
            self.r.set_width(0)
            self.r.set_xy((self.xs, self.ax.get_ylim()[0]))
            self.r.set_height(self.ax.get_ylim()[1] - self.ax.get_ylim()[0])
        else:
            self.r = Rectangle(xy=(self.xs, self.ax.get_ylim()[0]),
                               height=self.ax.get_ylim()[1] - \
                                      self.ax.get_ylim()[0],
                               width=0, fill=True, lw=2, alpha=0.5,
                               color='blue')
            self.ax.add_artist(self.r)
            self.r.set_clip_box(self.ax.bbox)
            self.r.set_zorder(10)

    def on_motion_draw(self):
        self.r.set_width(self.xe - self.xs)

    def on_release_print(self):
        pprint('X Span: (%d, %d)\n' % (self.xs, self.xe))


class rectSelection(GuiSelection):
    """Interactive selection of a rectangular region on the axis.

    Used by hist2d_alex().
    """
    def on_press_draw(self):
        if 'r' in self.__dict__:
            self.r.set_height(0)
            self.r.set_width(0)
            self.r.set_xy((self.xs, self.ys))
            self.e.height = 0
            self.e.width = 0
            self.e.center = (self.xs, self.ys)
        else:
            self.r = Rectangle(xy=(self.xs, self.ys), height=0, width=0,
                               fill=False, lw=2, alpha=0.5, color='blue')
            self.e = Ellipse(xy=(self.xs, self.ys), height=0, width=0,
                    fill=False, lw=2, alpha=0.6, color='blue')
            self.ax.add_artist(self.r)
            self.ax.add_artist(self.e)
            self.r.set_clip_box(self.ax.bbox)
            self.r.set_zorder(10)
            self.e.set_clip_box(self.ax.bbox)
            self.e.set_zorder(10)

    def on_motion_draw(self):
        self.r.set_height(self.ye - self.ys)
        self.r.set_width(self.xe - self.xs)
        self.e.height = (self.ye - self.ys)
        self.e.width = (self.xe - self.xs)
        self.e.center = (np.mean([self.xs, self.xe]),
                         np.mean([self.ys, self.ye]))
        self.fig.canvas.draw()

    def on_release_print(self):
        # This is the only custom method for hist2d_alex()
        E1, E2 = min((self.xs, self.xe)), max((self.xs, self.xe))
        S1, S2 = min((self.ys, self.ye)), max((self.ys, self.ye))
        self.selection = dict(E1=E1, E2=E2, S1=S1, S2=S2)
        pprint("Selection: \nE1=%.2f, E2=%.2f, S1=%.2f, S2=%.2f\n" %\
                (E1,E2,S1,S2))


class MultiAxPointSelection(object):
    """Class for point selection on a multi-axes plot.

    Used to select/print bursts by clicking in timetrace_ and ratetrace_ plots.
    """
    def __init__(self, fig, ax, d, debug=False):
        self.ax_list = [ax]
        self.fig = fig
        self.d = d
        self.debug = debug
        self._asymmetry_dict = {}
        self.id_press = fig.canvas.mpl_connect('button_press_event',
                                                self.on_press)

    def on_press(self, event):
        if self.debug:
            pprint('PRESS button=%d, x=%d, y=%d, xdata=%f, ydata=%f\n' % (
                event.button, event.x, event.y, event.xdata, event.ydata))

        iax = [i for i, ax in enumerate(self.ax_list) if ax == event.inaxes]
        if len(iax) == 0:
            if self.debug:
                pprint('NO axis found. event.inaxes "%s".\n' % event.inaxes)
                pprint('self.ax_list: ' + str(self.ax_list))
            return

        self.ich = iax[0]
        self.xp, self.yp = event.xdata, event.ydata
        self.on_press_print()

    def asymmetry(self, ich):
        if ich not in self._asymmetry_dict:
            self._asymmetry_dict[ich] = bext.asymmetry(self.d, ich)
        return self._asymmetry_dict[ich]

    def on_press_print(self):
        if self.debug:
            pprint("%s %s %s\n" % (self.xp, self.yp, self.ich))
        bursts = self.d.mburst[self.ich]
        t_clk = self.xp / self.d.clk_p
        mask = (bursts.start < t_clk) * (bursts.stop > t_clk)
        if mask.any():
            burst_index = np.where(mask)[0][0]
            burst = bursts[burst_index]
            params = dict(ich=self.ich, b_index=burst_index,
                start_ms = float(burst.start) * self.d.clk_p * 1e3,
                width_ms = float(burst.width) * self.d.clk_p * 1e3,
                asym = self.asymmetry(self.ich)[burst_index],
                nt = self.d.nt[self.ich][burst_index],
                nd = self.d.nd[self.ich][burst_index],
                na = self.d.na[self.ich][burst_index],
                E = self.d.E[self.ich][burst_index])
            msg = ("Burst [{b_index}-CH{ich}]: t = {start_ms:7.2f} ms"
                   "   width={width_ms:4.2f} ms"
                   "   size=(T{nt:3.0f}, D{nd:3.0f}, A{na:3.0f}")
            if self.d.ALEX:
                msg += ", AA{naa:3.0f}"
                params['naa'] = self.d.naa[self.ich][burst_index]
            msg += ")   E={E:4.2%}"
            if self.d.ALEX:
                msg += "   S={S:4.2%}"
                params['S'] = self.d.S[self.ich][burst_index]
            if 'particles' in self.d:
                msg = "P{par:2d} " + msg
                particles = self.d.particles[self.ich]
                params['par'] = mode(particles[burst.istart:burst.istop])[0][0]
            msg += "   Asym(D-A)={asym:5.2f} ms"
            pprint((msg + '\n').format(**params))
