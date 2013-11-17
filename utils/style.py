# USO TIPICO (copia-incolla)
# encoding: utf-8  μ αβγδ εζηθ κλμν ξοπρ ςστυ φχψω
# _ip.magic("run tesifig.py")
# axes(axis_pos)

from pylab import figure, rcParams, gca, draw

def pos_assi():
    """Posizione assi nella forma [left, bottom, width, height]"""
    p = gca().get_position().get_points().ravel()
    p[2] -= p[0]
    p[3] -= p[1]
    return p

def set_assi(pos):
    """Setta posizione assi con pos = [left, bottom, width, height]"""
    gca().set_position(pos)
    draw()

pylabfigure = figure

inch = 2.54
ratio = 3.15/2.5
x_size_in = 9/inch
y_size_in = 9/ratio/inch

grid_color='#b4b4b4' # grigio 30%

params = {
        'axes.axisbelow': True, # per griglia sotto linee plot
        'axes.grid': True,
        #'font.serif': 'FreeSerif',
        #'font.family': 'serif',
        'font.size': 6,
        'axes.titlesize': 6,
        'axes.labelsize': 6,
        'ytick.labelsize': 6,
        'xtick.labelsize': 6,
        'legend.fontsize': 6,
        'legend.fancybox': True,
        'legend.loc': 'best',

        'text.usetex': False,

        'figure.figsize': (x_size_in,y_size_in),
        #'figure.dpi': 200,
        'lines.markeredgewidth': 0.3, # regola anche lo spessore dei tick
        'lines.markersize': 4,
        #'lines.linewidth': 0.7,
        'axes.linewidth': 0.5,
        'grid.linewidth': 0.3,
        'grid.linestyle': '-',
        'grid.color': grid_color,
        'patch.linewidth': 0.5,    # spessore frame legenda
        #'mathtext.mit' : 'cmtt10.ttf',
    }
rcParams.update(params)

params = {
        'font.size': 6,
        'axes.titlesize': 6,
        #'axes.labelsize': 6,
        #'ytick.labelsize': 6,
        #'xtick.labelsize': 6,
        'legend.fontsize': 6,
        'legend.fancybox': True,
        'legend.loc': 'best',
    }
rcParams.update(params)


axis_pos = [0.16, 0.155, 0.78, 0.78]
#axis_pos = [0.15, 0.15, 0.80, 0.8]

###
#axes([0.15, 0.2, 0.80, 0.7])
###


#def tesifig(n=None, **kwargs):
#    
#    rcParams.update(params)
#    
#    if kwargs.has_key('figsize'):
#        if n is None:
#            f = pylabfigure(**kwargs)
#        else:
#            f = pylabfigure(n, **kwargs)
#    else:
#        if n is None:
#            f = pylabfigure(figsize=rcParams['figure.figsize'], **kwargs)
#        else:
#            f = pylabfigure(n, figsize=rcParams['figure.figsize'], **kwargs)
#    
#    return f

#figure = tesifig


