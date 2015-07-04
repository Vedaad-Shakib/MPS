import sys
import random
import types
import re
import sys
import os
import numpy
import matplotlib.pyplot as plt

from PyQt4                              import QtGui, QtCore
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib                         import pyplot as plt
from time                               import sleep

class Window(QtGui.QDialog):
    def __init__(self, increment=10, nSteps=500, timestep=0.0004, parent=None):
        super(Window, self).__init__(parent)

        # frame increment
        self.increment = increment
        self.nSteps = nSteps
        self.timestep = timestep # in seconds
        
        # a figure instance to plot on
        self.figure = plt.figure()
        self.ax1 = self.figure.add_subplot(1, 1, 1)

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        # Just some button connected to `plot` method
        self.button = QtGui.QPushButton('Plot')
        self.button.clicked.connect(self.start)

        # the current frame
        self.index = self.increment

        # cache the wallpoints and ghost points x and y
        # and get initial data set at t=0
        file = open("ghost_points.dat", "r")
        ghostX = []
        ghostY = []
        for line in file:
            tmp = [float(i) for i in line.split()]
            ghostX.append(tmp[0])
            ghostY.append(tmp[1])

        file = open("wall_points.dat", "r")
        wallX = []
        wallY = []
        for line in file:
            tmp = [float(i) for i in line.split()]
            wallX.append(tmp[0])
            wallY.append(tmp[1])

        file = open("mps.0.out", "r")
        x = []
        y = []
        for line in file:
            tmp = [float(i) for i in line.split()]
            x.append(tmp[0])
            y.append(tmp[1])

        self.line1, = self.ax1.plot(x, y, 'bo', alpha=0.8)
        self.line2, = self.ax1.plot(ghostX, ghostY, 'go', alpha=0.8)
        self.line3, = self.ax1.plot(wallX, wallY, 'ro', alpha=0.8)

        # idk why, but this works here
        self.ax1.cla()

        self.xmin, self.xmax, self.ymin, self.ymax = plt.axis()
        
        #self.figure.show()
        self.figure.canvas.draw()

        self.background = self.figure.canvas.copy_from_bbox(self.ax1.bbox)

        # set the layout
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        layout.addWidget(self.button)
        self.setLayout(layout)

    def start(self):
        # connects timer to dynamic plot
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.updatePlot)
        self.timer.start(1)

    def updatePlot(self):
        print self.index
        if self.index > self.increment*self.nSteps:
            print "asdfadf"
            self.timer.stop()
            return
        
        data = ("mps.%s.out") % self.index
        self.index += self.increment
        
        self.ax1.cla() # clear axes
        
        title = "Time = %.4g" % (self.timestep * self.index)
        self.ax1.set_title(title)

        file = open(data, "r")
        x = []
        y = []
        for line in file:
            tmp = [float(i) for i in line.split()]
            x.append(tmp[0])
            y.append(tmp[1])

        plt.axis([self.xmin, self.xmax, self.ymin, self.ymax])
        
        #self.ax1.plot(x, y, 'ro')
        self.line1.set_xdata(x)
        self.line1.set_ydata(y)
        self.figure.canvas.restore_region(self.background)
        self.ax1.draw_artist(self.line1)
        self.ax1.draw_artist(self.line2)
        self.ax1.draw_artist(self.line3)
        self.figure.canvas.blit(self.ax1.bbox)
        '''
        self.ax1.plot(numpy.arange(1422), numpy.arange(1422), 'ro')
        self.ax1.plot(numpy.arange(1007), numpy.arange(1007), 'ro')
        self.ax1.plot(numpy.arange(252), numpy.arange(252), 'ro')
        '''
        #self.ax1.plot(self.ghostX, self.ghostY, 'ro')
        #self.ax1.plot(self.wallX, self.wallY, 'ro')

        #self.canvas.draw()

def mpsPlot():
    app = QtGui.QApplication(sys.argv)

    main = Window()
    main.show()

    sys.exit(app.exec_())

if __name__ == "__main__":
    mpsPlot()
