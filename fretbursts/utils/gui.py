#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
# Credit for the out-of-process trick to Mikkel (github @mkrog):
# https://github.com/ipython/ipython/issues/5798#issuecomment-93851739
#
"""
GUI related helper functions.
"""

from sys import executable
from subprocess import check_output
import contextlib
import os


@contextlib.contextmanager
def chdir(dirname=None):
    """
    Source: http://www.astropython.org/snippet/2009/10/chdir-context-manager
    """
    curdir = os.getcwd()
    try:
        if dirname is not None:
            os.chdir(dirname)
        yield
    finally:
        os.chdir(curdir)


def OpenFileDialog(dirname=None):
    with chdir(dirname):
        file = check_output([executable, __file__])
    return file.strip()


def gui_fname(dir=None):
    """
    Select a file via a dialog and return the file name.
    """
    try:
        from PyQt5.QtWidgets import QApplication, QFileDialog
    except ImportError:
        try:
            from PyQt4.QtGui import QApplication, QFileDialog
        except ImportError:
            from PySide.QtGui import QApplication, QFileDialog

    if dir is None:
        dir = './'

    app = QApplication([dir])
    fname = QFileDialog.getOpenFileName(None, "Select a file...",
                                        dir, filter="All files (*)")

    if isinstance(fname, tuple):
        return fname[0]
    else:
        return str(fname)

if __name__ == "__main__":
    print(gui_fname())
