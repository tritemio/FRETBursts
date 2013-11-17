from scipy import *
from matplotlib.patches import Rectangle, Ellipse

xs, ys, xe, ye = 0,0,0,0
id_motion = 0
r = 0

Shape = Rectangle
#Shape = Ellipse

def on_press(event):
    global id_motion, xs, ys, r
    print 'INAXES: ', event.inaxes
    if event.inaxes!=ax: return
    xs, ys = event.xdata, event.ydata
    print 'PRESS button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % (
            event.button, event.x, event.y, event.xdata, event.ydata)
    r = Rectangle(xy=(xs,ys), height=0, width=0, fill=False, lw=2, alpha=0.2)
    ax.add_artist(r)
    r.set_clip_box(ax.bbox)
    draw()
    id_motion = fig.canvas.mpl_connect('motion_notify_event', on_motion)

def on_release(event):
    #print 'RELEASE button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % (
    #        event.button, event.x, event.y, event.xdata, event.ydata)
    print "Selection: ", xs, ys, xe, ye
    fig.canvas.mpl_disconnect(id_motion)

def on_motion(event):
    global xe, ye
    #print 'MOTION x=%d, y=%d, xdata=%f, ydata=%f' % (
    #        event.x, event.y, event.xdata, event.ydata)
    xe, ye = event.xdata, event.ydata
    print xe,ye
    r.set_height(ye-ys)
    r.set_width(xe-xs)
    draw()

sx = normal(size=5000)
sy = normal(size=5000)

scatter(sx,sy, linewidth=0, alpha=0.2)

fig = gcf()
id_press = fig.canvas.mpl_connect('button_press_event', on_press)
id_rls = fig.canvas.mpl_connect('button_release_event', on_release)
#id_motion = fig.canvas.mpl_connect('motion_notify_event', on_release)

ax = gca()
print ax
#r = Rectangle(xy=(0,0), height=1, width=2, fill=False, lw=2, alpha=0.2)
#ax.add_artist(r)
#r.set_clip_box(ax.bbox)

show()
