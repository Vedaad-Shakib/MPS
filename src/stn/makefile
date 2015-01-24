#*******************************************************************************          
#** Copyright 2014-2014 Vedaad Shakib Inc.
#*******************************************************************************

#*******************************************************************************          
#** "makefile": makefile
#*******************************************************************************          

all:: install

include ../../make/make.all

install:: $(INCDIR)/stn.h $(LIBDIR)/libstn.a

$(INCDIR)/stn.h: $(INCDIR) stn.h
	$(INSTALL) stn.h $(INCDIR)/stn.h

$(LIBDIR)/libstn.a: $(LIBDIR) libstn.a
	$(INSTALL) libstn.a $(LIBDIR)/libstn.a

libstn.a: stn.o
	$(AR) libstn.a stn.o

stn.o: stn.h stn.c $(INCDIR)/mem.h
	$(CC) -c $(CCFLAGS) -I$(INCDIR) stn.c

clean::
	$(RM) stn.o libstn.a

