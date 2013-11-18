"""
A PyGraphQT timetrace plot figure with a slider to scroll the time axis 
back and forth.

Adapted from:
http://stackoverflow.com/questions/16824718/python-matplotlib-pyside-fast-timetrace-scrolling
"""


from PySide import QtGui, QtCore
import numpy as np
import pyqtgraph as pg

N_SAMPLES = 1e6

def test_plot():
    time = np.arange(N_SAMPLES)*1e-3
    sample = np.random.randn(N_SAMPLES)

    plt = pg.PlotWidget(title="Use the slider to scroll and the spin-box to set the width")
    plt.addLegend()
    plt.plot(time, sample, name="Gaussian noise")
    q = ScrollingToolQT(plt)
    return q   # WARNING: it's important to return this object otherwise
            # python will delete the reference and the GUI will not respond!


class ScrollingToolQT(object):
    def __init__(self, fig):
        # Setup data range variables for scrolling
        self.fig = fig
        self.xmin, self.xmax = fig.plotItem.vb.childrenBounds()[0]
        self.step = 1 # axis units

        self.scale = 1e3 # conversion betweeen scrolling units and axis units

        # Retrive the QMainWindow used by current figure and add a toolbar
        # to host the new widgets
        self.win = QtGui.QMainWindow()
        self.win.show()
        self.win.resize(800,600)
        self.win.setCentralWidget(fig)
        self.toolbar = QtGui.QToolBar()
        self.win.addToolBar(QtCore.Qt.BottomToolBarArea, self.toolbar)

        # Create the slider and spinbox for x-axis scrolling in toolbar
        self.set_slider(self.toolbar)
        self.set_spinbox(self.toolbar)

        # Set the initial xlimits coherently with values in slider and spinbox
        self.set_xlim = self.fig.setXRange
        self.set_xlim(0, self.step)

    def set_slider(self, parent):
        # Slider only support integer ranges so use ms as base unit
        smin, smax = self.xmin*self.scale, self.xmax*self.scale

        self.slider = QtGui.QSlider(QtCore.Qt.Horizontal, parent=parent)
        self.slider.setTickPosition(QtGui.QSlider.TicksAbove)
        self.slider.setTickInterval((smax-smin)/10.)
        self.slider.setMinimum(smin)
        self.slider.setMaximum(smax-self.step*self.scale)
        self.slider.setSingleStep(self.step*self.scale/5.)
        self.slider.setPageStep(self.step*self.scale)
        self.slider.setValue(0)  # set the initial position
        self.slider.valueChanged.connect(self.xpos_changed)
        parent.addWidget(self.slider)

    def set_spinbox(self, parent):
        self.spinb = QtGui.QDoubleSpinBox(parent=parent)
        self.spinb.setDecimals(3)
        self.spinb.setRange(0.001, 3600.)
        self.spinb.setSuffix(" s")
        self.spinb.setValue(self.step)   # set the initial width
        self.spinb.valueChanged.connect(self.xwidth_changed)
        parent.addWidget(self.spinb)

    def xpos_changed(self, pos):
        #pprint("Position (in scroll units) %f\n" %pos)
        #        self.pos = pos/self.scale
        pos /= self.scale
        self.set_xlim(pos, pos + self.step, padding=0)

    def xwidth_changed(self, xwidth):
        #pprint("Width (axis units) %f\n" % step)
        if xwidth <= 0: return
        self.step = xwidth
        self.slider.setSingleStep(self.step*self.scale/5.)
        self.slider.setPageStep(self.step*self.scale)
        old_xlim = self.fig.plotItem.vb.viewRange()[0]
        self.xpos_changed(old_xlim[0] * self.scale)

if __name__ == "__main__":
    app = pg.mkQApp()
    q = test_plot()
    #app.exec_()
