"""
GUI related helper functions.
"""

from PySide import QtCore, QtGui


def gui_fname(dir=None):
    """Select a file via a dialog and returns the file name.
    """
    if dir is None: dir ='./'
    fname = QtGui.QFileDialog.getOpenFileName(None, "Select data file...", 
            dir, filter="All files (*);; SM Files (*.sm)")
    return fname[0]

