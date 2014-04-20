#!/usr/bin/env python
#
# Demonstration file: how to embed a matplotlib figure in a QT4 window
#

import sys
from PyQt4 import QtGui, QtCore

from matplotlib.figure import Figure

# Import from matplotlib the FigureCanvas with QTAgg backend
from matplotlib.backends.backend_qt4agg \
        import FigureCanvasQTAgg as FigureCanvas

# Import the matplotlib Toolbar2
from matplotlib.backends.backend_qt4agg \
        import NavigationToolbar2QT as NavigationToolbar

class PlotWindow(QtGui.QWidget):
    def __init__(self, show=True, create_axis=False):
        super(PlotWindow, self).__init__()
        
        # Create the plot canvas ...
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)        
        if create_axis: self.axis = self.figure.add_subplot(111)
        #self.canvas.setParent(self)
        
        # ... and the toolbar
        self.toolbar = NavigationToolbar(self.canvas, self)

        # Pack all in QVBoxLayout
        self.vbox = QtGui.QVBoxLayout(self)
        self.vbox.addWidget(self.canvas)
        self.vbox.addWidget(self.toolbar)

        if show: self.show()

class ScrollPlotWindow(PlotWindow):
    def __init__(self):
        super(ScrollPlotWindow, self).__init__(show=False, create_axis=True)
        self.slider = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        self.slider.valueChanged.connect(self.changeValue)
        self.slider.setTickPosition(QtGui.QSlider.TicksBelow)
        self.slider.setTickInterval(10)
        self.vbox.insertWidget(1, self.slider)
        self.show()

    def changeValue(self, value):
        self.axis.set_xlim(value,value+1)
        self.canvas.draw()

def main():
    #win = PlotWindow()
    win = ScrollPlotWindow()
    #win = TestWindow()
    sys.exit(app.exec_())    

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    #pass
    main()    
