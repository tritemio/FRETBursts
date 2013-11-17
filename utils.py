#
# Utility functions
#
import sys

def clk_to_s(t_ck, clk_p=12.5*1e-9): return t_ck*clk_p

def pprint(s):
    """Print immediately, even if inside a busy loop."""
    sys.stdout.write(s)
    sys.stdout.flush()

#
# Function to select a file graphically
#
from PySide import QtCore, QtGui

def gui_fname(dir=None):
    """Select a file via a dialog and returns the file name."""
    if dir is None: dir ='./'
    fname = QtGui.QFileDialog.getOpenFileName(None, "Select data file...", 
            dir, filter="All files (*);; SM Files (*.sm)")
    return fname[0]

