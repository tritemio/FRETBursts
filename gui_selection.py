"""
GUI selection of a range in a matplotlib figure.

Used by hist2d_alex() and other functions in burst_plot.py.
"""

from matplotlib.patches import Rectangle, Ellipse
from utils.misc import pprint


class GuiSelection(object):
    """Abstract class for range selection.
    
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
    def on_press(self, event):
        if event.inaxes != self.ax: return
        self.pressed = True
        self.xs, self.ys = event.xdata, event.ydata
        if self.debug: 
            print 'PRESS button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % (
                event.button, event.x, event.y, event.xdata, event.ydata)
        self.on_press_draw()    
        self.fig.canvas.draw()     
        self.id_motion = self.fig.canvas.mpl_connect('motion_notify_event', 
                                                     self.on_motion)
        self.fig.canvas.mpl_connect('button_release_event', 
                                             self.on_release)

    def on_motion(self, event):
        if event.inaxes != self.ax: return
        if self.debug: 
            print 'MOTION x=%d, y=%d, xdata=%f, ydata=%f' % (
                event.x, event.y, event.xdata, event.ydata)
        self.xe, self.ye = event.xdata, event.ydata
        self.on_motion_draw()
        self.fig.canvas.draw()
     
    def on_release(self, event):
        if not self.pressed: return
        self.pressed = False
        if self.debug: 
            print 'RELEASE button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % (
                event.button, event.x, event.y, event.xdata, event.ydata)
        self.fig.canvas.mpl_disconnect(self.id_motion)
        self.on_release_print()
    
    def on_press_draw(self):
        pass

    def on_motion_draw(self):
        pass
    
    def on_release_print(self):
        pass
        

class xspanSelection(GuiSelection):
    """Select an x range on the figure"""
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
    """Select a rectangular region on the figure (for hist2d_alex())"""
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
        self.e.center = (mean([self.xs, self.xe]), mean([self.ys, self.ye]))
        self.fig.canvas.draw()

    def on_release_print(self):
        E1, E2 = min((self.xs, self.xe)), max((self.xs, self.xe))
        S1, S2 = min((self.ys, self.ye)), max((self.ys, self.ye))
        pprint("Selection: \nE1=%.2f, E2=%.2f, S1=%.2f, S2=%.2f\n" %\
                (E1,E2,S1,S2))
    

class Point_Selection:
    """Select a point in a plot"""
    xp = 0
    yp = 0
    AX = []
    d = None
    connected = False

PSel = Point_Selection()

def on_press_point(event):
    axb = r_[[ax == event.inaxes for ax in PSel.AX]]
    if not axb.any(): return
    ich = find(axb)
    ax = PSel.AX[ich]

    PSel.xp, PSel.yp, PSel.ich = event.xdata, event.ydata, ich
    #print 'PRESS button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % (
    #        event.button, event.x, event.y, event.xdata, event.ydata)
    pprint("%s %s %s\n" %(PSel.xp, PSel.yp, PSel.ich))
    mburst = PSel.d.mburst[ich]
    t_clk = PSel.xp/PSel.d.clk_p
    mask = (b_start(mburst) < t_clk)*(b_end(mburst) > t_clk)    
    if mask.any():
        ib = find(mask)[0]
        ts = b_start(mburst)[ib]*PSel.d.clk_p
        skew = bleaching1(PSel.d, ich, ib, use_median=True, normalize=True) 
        pprint("Burst [%d-CH%d]: t = %d us   nt = %5.1f   E = %4.2f   "
               "Skew = %.2f\n" % (ib, ich+1, ts*1e6, 
                   PSel.d.nt[ich][ib], PSel.d.E[ich][ib], skew))
        
