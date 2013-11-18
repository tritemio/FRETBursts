"""
Hack to allow GUI selection of bursts in a 2D ALEX histogram.

See hist2d_alex() in burst_plot.py.
"""


from matplotlib.patches import Rectangle, Ellipse

class GUI_Selection:
    id_motion = 0
    xs = 0
    ys = 0
    xe = 0
    ye = 0
    ax = 0
    pressed = False

GSel = GUI_Selection()

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
        
def _on_press(event):
    if event.inaxes != GSel.ax: return
    GSel.pressed = True
    GSel.xs, GSel.ys = event.xdata, event.ydata
    #print 'PRESS button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % (
    #        event.button, event.x, event.y, event.xdata, event.ydata)
    if hasattr(GSel, 'r'):
        GSel.r.set_height(0); GSel.r.set_width(0)
        GSel.r.set_xy((GSel.xs,GSel.ys))
        GSel.e.height = 0; GSel.e.width = 0
        GSel.e.center = (GSel.xs,GSel.ys)
    else:
        GSel.r = Rectangle(xy=(GSel.xs,GSel.ys), height=0, width=0, 
                fill=False, lw=2, alpha=0.6, color='blue')
        GSel.e = Ellipse(xy=(GSel.xs,GSel.ys), height=0, width=0, 
                fill=False, lw=2, alpha=0.6, color='blue')
        GSel.ax.add_artist(GSel.r)
        GSel.r.set_clip_box(GSel.ax.bbox); GSel.r.set_zorder(10)
        GSel.ax.add_artist(GSel.e)
        GSel.e.set_clip_box(GSel.ax.bbox); GSel.e.set_zorder(10)
    GSel.fig.canvas.draw()
    GSel.id_motion = GSel.fig.canvas.mpl_connect('motion_notify_event', 
            _on_motion)

def _on_release(event):
    if not GSel.pressed: return
    GSel.pressed = False
    #print 'RELEASE button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % (
    #        event.button, event.x, event.y, event.xdata, event.ydata)
    GSel.fig.canvas.mpl_disconnect(GSel.id_motion)
    #print "Selection: \nE1=%.3f; E2=%.3f; S1=%.3f; S2=%.3f" % \
    #        (GSel.xs, GSel.xe, GSel.ys, GSel.ye)
    E1, E2 = min((GSel.xs, GSel.xe)), max((GSel.xs, GSel.xe))
    S1, S2 = min((GSel.ys, GSel.ye)), max((GSel.ys, GSel.ye))
    print "Selection: \nE1=%.2f, E2=%.2f, S1=%.2f, S2=%.2f" % (E1,E2,S1,S2)

def _on_motion(event):
    if event.inaxes != GSel.ax: return
    #print 'MOTION x=%d, y=%d, xdata=%f, ydata=%f' % (
    #        event.x, event.y, event.xdata, event.ydata)
    GSel.xe, GSel.ye = event.xdata, event.ydata
    GSel.r.set_height(GSel.ye-GSel.ys); GSel.r.set_width(GSel.xe-GSel.xs)
    GSel.e.height=(GSel.ye-GSel.ys); GSel.e.width = (GSel.xe-GSel.xs)
    GSel.e.center = (mean([GSel.xs,GSel.xe]),mean([GSel.ys,GSel.ye]))
    GSel.fig.canvas.draw()



