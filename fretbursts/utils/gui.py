#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
GUI related helper functions.
"""

from __future__ import print_function

try:
    from PyQt4 import QtGui, QtCore
except ImportError:
    from PySide import QtGui, QtCore


def gui_fname(dir=None):
    """Select a file via a dialog and returns the file name.
    """
    if dir is None: dir ='./'
    fname = QtGui.QFileDialog.getOpenFileName(None, "Select data file...",
            dir, filter="All files (*);; SM Files (*.sm)")

    if type(fname) is tuple:
        fname = fname[0]
    elif type(fname) is QtCore.QString:
        fname = str(fname)

    print(fname)
    return fname

