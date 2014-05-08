#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
A global plot style for matplotlib.
"""

from matplotlib.pyplot import rcParams


# This will fail on ReadTheDocs so I use a try-except
try:
    fontsize = fs = 12

    font1 = {'fontname':'Liberation Sans','fontsize':16}
    font2 = {'fontname':'Arial','fontsize':16}

    #rcParams["font.family"] = font2['fontname']
    rcParams["font.sans-serif"] = ['Arial', 'Liberation Sans']
    rcParams["font.size"] = fontsize

    rcParams['xtick.labelsize'] = fontsize
    rcParams['ytick.labelsize'] = fontsize
    rcParams['axes.labelsize'] = fontsize
    rcParams['legend.fontsize'] = fontsize - 1
    rcParams['lines.linewidth'] = 2

    # Define 8 different colors (1 per ch)
    rcParams['axes.color_cycle'] = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange']
except:
    pass
