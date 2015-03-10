from mpsPlot import *
import time

_colors = ["blue", "green", "red", "cyan", "magenta", "yellow", "black", "white"]
_nColors = len(_colors)

def plot(*steps):
    pyplot.close()
    ps = []
    ps.append( ["wall_points.dat","cyan","."] )
    ps.append( ["ghost_points.dat","green","+"] )
    i = 0
    for step in steps:
	ps.append( ["mps.%d.out" % step, _colors[i] ] )
	i = (i+1) % _nColors
    mpsPlot( ps )

def anim( end, inc=20, beg=0, wait=0.2 ):
    for i in range(beg,end,inc):
	plot(i)
	time.sleep(wait)

