#*******************************************************************************          
#** Copyright 2014-2014 Vedaad Shakib Inc.
#*******************************************************************************

#*******************************************************************************          
#** "makefile": makefile
#*******************************************************************************          

all:: install

include make/make.all

DIRS = src dam

install:: install_dirs

clean:: clean_dirs
	rm -rf *~
	rm -rf Log
	rm -rf \#*
	rm -rf *.dat
	rm -rf *.out
	rm -rf bin/*.dSYM

