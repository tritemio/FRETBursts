## GLOBAL Plot style
#import matplotlib as mpl

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
