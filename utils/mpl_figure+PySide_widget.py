from PySide import QtGui, QtCore
from PySide import QtGui, QtCore

def test():
    plot([1,2,3], lw=2)
    q = qt4_interface(gcf())
    return q   # WARNING: it's paramount to return the object otherwise, with 
               # no references, python deletes it and the GUI doesn't respond!
    
class qt4_interface:
    def __init__(self,fig):
        self.fig = fig

        QMainWin = fig.canvas.parent()
        toolbar = QtGui.QToolBar(QMainWin)
        QMainWin.addToolBar(QtCore.Qt.BottomToolBarArea, toolbar)
    
        self.line_edit = QtGui.QLineEdit()#parent=toolbar)
        toolbar.addWidget(self.line_edit)
        self.line_edit.editingFinished.connect(self.do_something) 

    def do_something(self, *args):
        self.fig.axes[0].set_title(self.line_edit.text())
        self.fig.canvas.draw()
        #f = open('l','a'); f.write('yes\n'); f.flush(); f.close()


