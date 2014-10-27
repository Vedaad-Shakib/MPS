from matplotlib import pyplot
import types

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
    # data[i][0] is fileName
    # options include "line", "point", "circle", "x", "dot", "red", "blue", "green"
    pyplot.clf() # clear pyplot

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
        
        file = open(i[0], "r")
        x = []
        y = []
        for line in file:
            tmp = [float(j) for j in line.split()]
            x.append(tmp[0])
            y.append(tmp[1])
        
        if "line" in i:
            pyplot.plot(x, y, c=color, linewidth=markersize, marker=marker)
        if "point" in i:
            pyplot.scatter(x, y, c=color, s=markersize, marker=marker)
        if not "line" in i and not "point" in i:
            pyplot.scatter(x, y, c=color, s=markersize, marker=marker)
                
    pyplot.show( block = False )
