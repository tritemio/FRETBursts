"""
A Matplotlib timetrace plot figure with the addition of a QT slider to scroll
the time axis back and forth.
"""

from PySide import QtGui, QtCore
import pylab as plt
import numpy as np

N_SAMPLES = 1e6

def test_plot():
    time = np.arange(N_SAMPLES)*1e-3
    sample = np.random.randn(N_SAMPLES)
    plt.plot(time, sample, label="Gaussian noise")
    plt.title("1000s Timetrace \n (use the slider to scroll and the spin-box "
              "to set the width)")
    plt.xlabel('Time (s)')
    plt.legend(fancybox=True)
    q = ScrollingToolQT(plt.gcf(), scroll_step=1)
    return q   # WARNING: it's important to return this object otherwise
               # python will delete the reference and the GUI will not respond!


class ScrollingToolQT(object):
    def __init__(self, fig, scroll_step=10):
        # Setup data range variables for scrolling
        self.fig = fig
        self.scroll_step = scroll_step
        self.xmin, self.xmax = fig.axes[0].get_xlim()
        self.width = 1 # axis units
        self.pos = 0   # axis units
        self.scale = 1e3 # conversion betweeen scrolling units and axis units
        
        # Save some MPL shortcuts
        self.draw = self.fig.canvas.draw
        self.draw_idle = self.fig.canvas.draw_idle
        self.ax = self.fig.axes[0]

        # Retrive the QMainWindow used by current figure and add a toolbar
        # to host the new widgets
        QMainWin = fig.canvas.parent()
        toolbar = QtGui.QToolBar(QMainWin)
        QMainWin.addToolBar(QtCore.Qt.BottomToolBarArea, toolbar)

        # Create the slider and spinbox for x-axis scrolling in toolbar
        self.set_slider(toolbar)
        self.set_spinbox(toolbar)

        # Set the initial xlimits coherently with values in slider and spinbox
        self.ax.set_xlim(self.pos,self.pos+self.width)
        self.draw()

    def set_slider(self, parent):
        self.slider = QtGui.QSlider(QtCore.Qt.Horizontal, parent=parent)
        self.slider.setTickPosition(QtGui.QSlider.TicksAbove)
        self.slider.setTickInterval((self.xmax-self.xmin)/10.*self.scale)
        self.slider.setMinimum(self.xmin*self.scale)
        self.slider.setMaximum((self.xmax-self.width)*self.scale)
        self.slider.setSingleStep(self.width*self.scale/4.)
        self.slider.setPageStep(self.scroll_step*self.width*self.scale)
        self.slider.setValue(self.pos*self.scale) # set the initial position
        self.slider.valueChanged.connect(self.xpos_changed)
        parent.addWidget(self.slider)

    def set_spinbox(self, parent):
        self.spinb = QtGui.QDoubleSpinBox(parent=parent)
        self.spinb.setDecimals(3)
        self.spinb.setRange(0.001, 3600.)
        self.spinb.setSuffix(" s")
        self.spinb.setValue(self.width)   # set the initial width
        self.spinb.valueChanged.connect(self.xwidth_changed)
        parent.addWidget(self.spinb)

    def xpos_changed(self, pos):
        #pprint("Position (in scroll units) %f\n" %pos)
        pos /= self.scale
        self.ax.set_xlim(pos, pos + self.width)
        self.draw()

    def xwidth_changed(self, width):
        #pprint("Width (axis units) %f\n" % step)
        if width <= 0: return
        self.width = width
        self.slider.setSingleStep(self.width*self.scale/5.)
        self.slider.setPageStep(self.scroll_step*self.width*self.scale)
        old_xlim = self.ax.get_xlim()
        self.xpos_changed(old_xlim[0] * self.scale)


if __name__ == "__main__":
    q = test_plot()
    plt.show()
