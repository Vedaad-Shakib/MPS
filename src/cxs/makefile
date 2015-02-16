#*******************************************************************************          
#** Copyright 2014-2014 Vedaad Shakib Inc.
#*******************************************************************************

#*******************************************************************************          
#** "makefile": makefile
#*******************************************************************************          

all:: install

include ../../make/make.all

install:: $(INCDIR)/cxs.h $(LIBDIR)/libcxs.a

$(INCDIR)/cxs.h: $(INCDIR) cxs.h
	$(INSTALL) cxs.h $(INCDIR)/cxs.h

$(LIBDIR)/libcxs.a: $(LIBDIR) libcxs.a
	$(INSTALL) libcxs.a $(LIBDIR)/libcxs.a

libcxs.a: cxs.o
	$(AR) libcxs.a cxs.o

cxs.o: cxs.h cxs.c 
	$(CC) -c $(CCFLAGS) -I$(INCDIR) cxs.c

clean::
	$(RM) cxs.o libcxs.a

