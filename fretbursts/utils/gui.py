#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
GUI related helper functions.
"""

from sys import executable
from subprocess import check_output


def OpenFileDialog():
    file = check_output([executable, __file__])
    return file.strip()

def gui_fname(dir=None):
    """
    Select a file via a dialog and return the file name.
    """
    try:
        from PySide import QtGui
    except ImportError:
        from PyQt4 import QtGui

    if dir is None:
        dir ='./'

    app = QtGui.QApplication([dir])
    fname = QtGui.QFileDialog.getOpenFileName(None, "Select a file...",
            dir, filter="All files (*)")

    if isinstance(fname, tuple):
        return fname[0]
    else:
        return str(fname)

if __name__ == "__main__":
    print(gui_fname())
