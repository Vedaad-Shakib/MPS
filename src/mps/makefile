#*******************************************************************************          
#** Copyright 2014-2014 Vedaad Shakib Inc.
#*******************************************************************************

#*******************************************************************************          
#** "makefile": makefile
#*******************************************************************************          

all:: install

include ../../make/make.all

INCS	= $(INCDIR)/sys.h $(INCDIR)/mem.h $(INCDIR)/que.h $(INCDIR)/stn.h $(INCDIR)/slv.h

install:: $(INCDIR)/mps.h $(BINDIR)/mps

$(INCDIR)/mps.h: $(INCDIR) mps.h
	$(INSTALL) mps.h $(INCDIR)/mps.h

mps.o: mps.h mps.c $(INCS)
	$(CC) -c $(CCFLAGS) -I$(INCDIR) mps.c

mpsPoints.o: mps.h mpsPoints.c $(INCS)
	$(CC) -c $(CCFLAGS) -I$(INCDIR) mpsPoints.c

MPSOBJS	= mps.o mpsPoints.o

mps: $(MPSOBJS)
	$(LD) -o mps $(LDFLAGS) $(MPSOBJS) -lque -lstn -lslv -lcxs $(CXSPARSELIB) $(LDLIBS)

$(BINDIR)/mps: $(BINDIR) mps
	$(INSTALL) mps $(BINDIR)/mps

clean::
	$(RM) *.o *.a mps

