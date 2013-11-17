from PySide import QtGui, QtCore
import pylab as plt
import numpy as np

N_SAMPLES = 1e6

def test_plot():
    dx = 1e-3 # this is the sampling period (or grid spacing) of `sample`
    sample = np.random.randn(N_SAMPLES)
    line, = plt.plot([], [], label="Gaussian noise")
    plt.legend(fancybox=True)
    plt.title("Use the slider to scroll and the spin-box to set the width")
    q = ScrollingPlotQT(fig=plt.gcf(), line=line, ydata=sample, dx=1e-3)
    return q   # WARNING: it's important to return this object otherwise
               # python will delete the reference and the GUI will not respond!


class ScrollingPlotQT(object):
    def __init__(self, fig, line, ydata, dx, scroll_step=10):
        # Setup data range variables for scrolling
        self.fig, self.line, self. ydata, self.dx = fig, line, ydata, dx
        self.scroll_step = scroll_step
        self.xmin, self.xmax = 0, dx*ydata.size
        self.plot_width = 1 # axis units
        self.scale = 1e3    # conversion betweeen scrolling units and axis units
        self.disp_points = self.plot_width/self.dx

        # Retrive the QMainWindow used by current figure and add a toolbar
        # to host the new widgets
        QMainWin = fig.canvas.parent()
        toolbar = QtGui.QToolBar(QMainWin)
        QMainWin.addToolBar(QtCore.Qt.BottomToolBarArea, toolbar)

        # Create the slider and spinbox for x-axis scrolling in toolbar
        self.set_slider(toolbar)
        self.set_spinbox(toolbar)

        # Set the initial xlimits coherently with values in slider and spinbox
        self.draw = self.fig.canvas.draw
        self.ax = self.fig.axes[0]
        self.ax.set_xlim(0, self.plot_width)
        self.ax.set_ylim(ydata.min(),ydata.max())
        
        # Setup the initial plot
        self.line.set_data(np.arange(self.disp_points)*dx,
                ydata[:self.disp_points])

        text0 = self.ax.text(0.01,0.02, "T = ", transform=fig.transFigure)
        self.text = self.ax.text(0.05,0.02, "0", transform=fig.transFigure)
        self.fig.canvas.draw()

    def set_slider(self, parent):
        self.slider = QtGui.QSlider(QtCore.Qt.Horizontal, parent=parent)
        self.slider.setTickPosition(QtGui.QSlider.TicksAbove)
        self.slider.setTickInterval((self.xmax-self.xmin)/10.*self.scale)
        self.slider.setMinimum(self.xmin*self.scale)
        self.slider.setMaximum((self.xmax-self.plot_width)*self.scale)
        self.slider.setSingleStep(self.plot_width*self.scale/4.)
        self.slider.setPageStep(self.scroll_step*self.plot_width*self.scale)
        self.slider.setValue(0)  # set the initial position
        self.slider.valueChanged.connect(self.xpos_changed)
        parent.addWidget(self.slider)

    def set_spinbox(self, parent):
        self.spinb = QtGui.QDoubleSpinBox(parent=parent)
        self.spinb.setDecimals(3)
        self.spinb.setRange(0.001, 3600.)
        self.spinb.setSuffix(" s")
        self.spinb.setValue(self.plot_width)   # set the initial width
        self.spinb.valueChanged.connect(self.xwidth_changed)
        parent.addWidget(self.spinb)

    def xpos_changed(self, pos):
        #print("Position (in scroll units) %f\n" %pos)
        pos /= (self.dx*self.scale) # pos converted in index units for ydata
        self.line.set_ydata(self.ydata[pos:pos+self.disp_points])
        self.text.set_text(pos*self.dx)
        self.draw()

    def xwidth_changed(self, xwidth):
        if xwidth <= 0: return
        self.plot_width = xwidth
        self.ax.set_xlim(0, self.plot_width)
        self.disp_points = self.plot_width/self.dx
        self.line.set_xdata(np.arange(self.disp_points)*self.dx)
        self.slider.setMaximum((self.xmax-self.plot_width)*self.scale)
        self.slider.setSingleStep(self.plot_width*self.scale/5.)
        self.slider.setPageStep(self.plot_width*self.scale)
        self.xpos_changed(self.slider.value()*self.scale)


if __name__ == "__main__":
    q = test_plot()
    plt.show()
