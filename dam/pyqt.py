import sys
import random
import types
import re
import sys
import os
import matplotlib.pyplot as plt

from PyQt4                              import QtGui, QtCore
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib                         import pyplot as plt
from time                               import sleep


def findColor(i):
    colors = ["blue", "green", "red", "cyan", "magenta", "yellow", "black", "white"]
    for elem in i:
        elem = str(elem)
        if elem in colors or (len(elem) == 7 and elem[0] == "#"):
            return elem
    return "blue"

def findSize(i):
    for elem in i:
        if type(elem) == types.IntType:
            return elem
    return 20

def findMarker(i):
    for elem in i:
        if str(elem) in ". , o v ^ < > s p * + x D d h H | _ ":
            return elem.strip()
    return "o"

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
        self.index = 0

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
        if self.index > self.increment*self.nSteps:
            print "asdfadf"
            self.timer.stop()
            return
        
        data = [["ghost_points.dat", "green"], ["wall_points.dat", "red"], [("mps.%s.out") % self.index, "blue"]]
        self.index += self.increment
        
        self.ax1.cla() # clear axes
        
        for i in data:
            color = findColor(i)
            markersize = findSize(i)
            marker = findMarker(i)
            
            m = re.match("mps.([0-9]+).out", i[0])
            if m:
                step = int(m.group(1))
                title = "Time = %.4g" % (self.timestep * step)
                self.ax1.set_title(title)
            file = open(i[0], "r")
            x = []
            y = []
            for line in file:
                tmp = [float(j) for j in line.split()]
                x.append(tmp[0])
                y.append(tmp[1])
                
            self.ax1.scatter(x, y, c=color, s=markersize, marker=marker)

        '''self.ax1.cla()
        self.ax1.scatter(range(10), [random.randint(1, 10) for i in range(10)])'''
        self.canvas.draw()

def mpsPlot():
    app = QtGui.QApplication(sys.argv)

    main = Window()
    main.show()

    sys.exit(app.exec_())

if __name__ == "__main__":
    mpsPlot()
