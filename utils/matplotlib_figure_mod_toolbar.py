from PySide import QtGui, QtCore
import matplotlib
#mpl_navtbar = matplotlib.backends.backend_qt4.NavigationToolbar2QT

def test():
    plot([1,2,3], lw=2)
    q = qt4_interface(gcf())
    return q   # WARNING: it's paramount to return the object otherwise, with 
               # no references, python deletes it and the GUI doesn't respond!
    
class qt4_interface:
    def __init__(self,fig):
        self.fig = fig

        #QMainWin = fig.canvas.parent()
        #toolbar = [w for w in QMainWin.children() if type(w) is mpl_navtbar][0]
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
        #self.fig.axes[0].set_xmin(self.line_edit.text())
        #self.fig.canvas.draw()
        print "ciao"
        klsjd
        f = open('l','a'); f.write('yes\n'); f.flush(); f.close()
