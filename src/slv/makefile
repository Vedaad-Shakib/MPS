#*******************************************************************************          
#** Copyright 2014-2014 Vedaad Shakib Inc.
#*******************************************************************************

#*******************************************************************************          
#** "makefile": makefile
#*******************************************************************************          

all:: install

include ../../make/make.all

install:: $(INCDIR)/slv.h $(LIBDIR)/libslv.a

$(INCDIR)/slv.h: $(INCDIR) slv.h
	$(INSTALL) slv.h $(INCDIR)/slv.h

$(LIBDIR)/libslv.a: $(LIBDIR) libslv.a
	$(INSTALL) libslv.a $(LIBDIR)/libslv.a

libslv.a: slv.o
	$(AR) libslv.a slv.o

slv.o: slv.h slv.c $(INCDIR)/mem.h
	$(CC) -c $(CCFLAGS) -I$(INCDIR) slv.c

clean::
	$(RM) slv.o libslv.a

