#*******************************************************************************          
#** Copyright 2014-2014 Vedaad Shakib Inc.
#*******************************************************************************

#*******************************************************************************          
#** "makefile": makefile
#*******************************************************************************          

all:: install

include ../../make/make.all

INCS	= $(INCDIR)/sys.h $(INCDIR)/mem.h $(INCDIR)/que.h

install:: $(INCDIR)/mps.h $(BINDIR)/mps

$(INCDIR)/mps.h: $(INCDIR) mps.h
	$(INSTALL) mps.h $(INCDIR)/mps.h

mps.o: mps.h mps.c $(INCS)
	$(CC) -c $(CCFLAGS) -I$(INCDIR) mps.c

mpsCorners.o: mps.h mpsCorners.c $(INCS)
	$(CC) -c $(CCFLAGS) -I$(INCDIR) mpsCorners.c

mpsFluidPoints.o: mps.h mpsFluidPoints.c $(INCS)
	$(CC) -c $(CCFLAGS) -I$(INCDIR) mpsFluidPoints.c

mpsWallPoints.o: mps.h mpsWallPoints.c $(INCS)
	$(CC) -c $(CCFLAGS) -I$(INCDIR) mpsWallPoints.c

mpsGhostPoints.o: mps.h mpsGhostPoints.c $(INCS)
	$(CC) -c $(CCFLAGS) -I$(INCDIR) mpsGhostPoints.c

MPSOBJS	= mps.o mpsCorners.o mpsWallPoints.o mpsGhostPoints.o mpsFluidPoints.o

mps: $(MPSOBJS)
	$(LD) -o mps $(LDFLAGS) $(MPSOBJS) -lque $(LDLIBS)

$(BINDIR)/mps: $(BINDIR) mps
	$(INSTALL) mps $(BINDIR)/mps

clean::
	$(RM) *.o *.a mps

