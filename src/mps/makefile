#*******************************************************************************          
#** Copyright 2014-2014 Vedaad Shakib Inc.
#*******************************************************************************

#*******************************************************************************          
#** "makefile": makefile
#*******************************************************************************          

all:: install

include ../../make/make.all

INCS	= $(INCDIR)/sys.h $(INCDIR)/mem.h $(INCDIR)/que.h $(INCDIR)/stn.h $(INCDIR)/slv.h
LIBS	= $(LIBDIR)/libque.a $(LIBDIR)/libstn.a $(LIBDIR)/libslv.a $(LIBDIR)/libcxs.a

install:: $(INCDIR)/mps.h $(BINDIR)/mps

$(INCDIR)/mps.h: $(INCDIR) mps.h
	$(INSTALL) mps.h $(INCDIR)/mps.h

mps.o: mps.h mps.c $(INCS)
	$(CC) -c $(CCFLAGS) -I$(INCDIR) mps.c

mpsPoints.o: mps.h mpsPoints.c $(INCS)
	$(CC) -c $(CCFLAGS) -I$(INCDIR) mpsPoints.c

mpsRead.o: mps.h mpsRead.c $(INCS)
	$(CC) -c $(CCFLAGS) -I$(INCDIR) mpsRead.c

mpsUtil.o: mps.h mpsUtil.c $(INCS)
	$(CC) -c $(CCFLAGS) -I$(INCDIR) mpsUtil.c

MPSOBJS	= mps.o mpsPoints.o mpsRead.o mpsUtil.o

mps: $(MPSOBJS) $(LIBS)
	$(LD) -o mps $(LDFLAGS) $(MPSOBJS) -lque -lstn -lslv -lcxs $(CXSPARSELIB) $(LDLIBS)

$(BINDIR)/mps: $(BINDIR) mps
	$(INSTALL) mps $(BINDIR)/mps

clean::
	$(RM) *.o *.a mps

