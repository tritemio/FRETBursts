"""
Example on how to add widgets the toolbar of a Matplotlib figure using the
QT backend.

No QT application is created, only the toolbar of the native MPL figure is
modified.
"""


from PySide import QtGui, QtCore
import matplotlib


def test():
    plot([1,2,3], lw=2)
    q = qt4_interface(gcf())
    return q   # WARNING: it's paramount to return the object otherwise, with 
               # no references, python deletes it and the GUI doesn't respond!
    
class qt4_interface:
    def __init__(self,fig):
        self.fig = fig
        toolbar = fig.canvas.toolbar
        
        self.line_edit = QtGui.QLineEdit()
        toolbar.addWidget(self.line_edit)
        self.line_edit.editingFinished.connect(self.do_something)

        self.spinbox = QtGui.QDoubleSpinBox()
        toolbar.addWidget(self.spinbox)
        self.spinbox.valueChanged.connect(self.do_something2)

    def do_something(self, *args):
        self.fig.axes[0].set_title(self.line_edit.text())
        self.fig.canvas.draw()
        #f = open('l','a'); f.write('yes\n'); f.flush(); f.close()    
    
    def do_something2(self, *args):
        self.fig.axes[0].set_xlim(0, self.spinbox.value())
        self.fig.canvas.draw()
        #f = open('l','a'); f.write('yes\n'); f.flush(); f.close()
