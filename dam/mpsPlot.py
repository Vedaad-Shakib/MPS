#******************************************************************************
#** Copyright 2014-2014 Vedaad Shakib Inc.
#******************************************************************************

#******************************************************************************
#* 
#* "mpsPlot.py": Extension of matplotlib plotting software
#*
#******************************************************************************

from matplotlib import pyplot as plt
from time       import sleep

import types
import re
import sys
import os

_block = False

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

def mpsPlot(data):
    step = None # whether a mps.*.out is being plotted; if it is, then the plot is also saved into a pngx
    title = ""
    # data[i][0] is fileName
    # options include "line", "point", "circle", "x", "dot", "red", "blue", "green"
    plt.clf() # clear plt

    if type(data) == types.StringType:
	data = [data]

    if type(data) == types.ListType and len(data) > 0 and type(data[0]) == types.StringType:
        data = [data]

    for i in data:
        # error handling
	if type(i) == types.StringType:
	    i = [i, "point"]

        # parse options
          # color
        color = findColor(i)
          # markersize
        markersize = findSize(i)
          # marker type
        marker = findMarker(i)

        m = re.match("mps.([0-9]+).out", i[0])
        if m:
	    step = int(m.group(1))
	    title = "Time = %.4g" % (0.001 * step)
        file = open(i[0], "r")
        x = []
        y = []
        for line in file:
            tmp = [float(j) for j in line.split()]
            x.append(tmp[0])
            y.append(tmp[1])
        
        if "line" in i:
            plt.plot(x, y, c=color, linewidth=markersize, marker=marker)
        if "point" in i:
            plt.scatter(x, y, c=color, s=markersize, marker=marker)
        if not "line" in i and not "point" in i:
            plt.scatter(x, y, c=color, s=markersize, marker=marker)
    
    plt.suptitle(title)
    plt.axes().set_aspect('equal')
    plt.show(block = _block)
    if step: plt.savefig("mps.%d.png"%step)

# saves data into a png file
def mpsSave(data, outFile):
    step = None # whether a mps.*.out is being plotted; if it is, then the plot is also saved into a pngx
    title = ""
    # data[i][0] is fileName
    # options include "line", "point", "circle", "x", "dot", "red", "blue", "green"
    plt.clf() # clear plt

    if type(data) == types.StringType:
	data = [data]

    if type(data) == types.ListType and len(data) > 0 and type(data[0]) == types.StringType:
        data = [data]

    for i in data:
        # error handling
	if type(i) == types.StringType:
	    i = [i, "point"]

        # parse options
          # color
        color = findColor(i)
          # markersize
        markersize = findSize(i)
          # marker type
        marker = findMarker(i)

        m = re.match("mps.([0-9]+).out", i[0])
        if m:
	    step = int(m.group(1))
	    title = "Time = %.4g" % (0.001 * step)
        file = open(i[0], "r")
        x = []
        y = []
        for line in file:
            tmp = [float(j) for j in line.split()]
            x.append(tmp[0])
            y.append(tmp[1])
        
        if "line" in i:
            plt.plot(x, y, c=color, linewidth=markersize, marker=marker)
        if "point" in i:
            plt.scatter(x, y, c=color, s=markersize, marker=marker)
        if not "line" in i and not "point" in i:
            plt.scatter(x, y, c=color, s=markersize, marker=marker)

    
    plt.suptitle(title)
    plt.axes().set_aspect('equal')
    plt.axis((-0.5, 4.5, -0.5, 4.0))
    plt.savefig(outFile)

def mpsPlotDensity(data, density):
    plt.clf() # clear plt

    x = []
    y = []
    z = []
    for j in open(data):
        j = j.replace("\n", "")
        try: tmp = [float(k) for k in j.split(" ") if k != ""]
        except: print j.split(" ")
        x.append(tmp[0])
        y.append(tmp[1])
    for j in open(density):
        z.append(float(j))

    #fig, ax = plt.subplots()
    #ax.scatter(x, y, c=z, s=80, edgecolor='')
    plt.scatter(x, y, c=z, s=80, edgecolor='')
    plt.axes().set_aspect('equal')
    plt.show(block = _block)
    
#******************************************************************************
#* Test
#******************************************************************************

if __name__ == "__main__":
    inc         = int(sys.argv[1])
    nSteps      = int(sys.argv[2])

    _block	= True
    wall	= ["wall_points.dat", "green", "x"]
    ghost	= ["ghost_points.dat", "magenta", "+"]
    fluid	= ["fluid_points.dat", "blue"]

    for i in range(0, (nSteps+1)*inc, inc):
        print "Processing step %s" % i
        num     = (len(str(nSteps*inc))-len(str(i)))*"0"+str(i) # pad with 0's
        inFile  = "mps.%s.out" % i
	outFile	= "mps.%s.png" % num
	mpsSave([[inFile, "blue"], wall, ghost], outFile)
    print "Converting to gif"
    os.system("convert -delay 5 -loop 0 *png solution.gif")
    os.system("rm -rf *png")

