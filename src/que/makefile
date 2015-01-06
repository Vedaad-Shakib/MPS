#*******************************************************************************          
#** Copyright 2014-2014 Vedaad Shakib Inc.
#*******************************************************************************

#*******************************************************************************          
#** "makefile": makefile
#*******************************************************************************          

all:: install

include ../../make/make.all

install:: $(INCDIR)/que.h $(LIBDIR)/libque.a

$(INCDIR)/que.h: $(INCDIR) que.h
	$(INSTALL) que.h $(INCDIR)/que.h

$(LIBDIR)/libque.a: $(LIBDIR) libque.a
	$(INSTALL) libque.a $(LIBDIR)/libque.a

libque.a: que.o
	$(AR) libque.a que.o

que.o: que.h que.c $(INCDIR)/mem.h
	$(CC) -c $(CCFLAGS) -I$(INCDIR) que.c

clean::
	$(RM) que.o libque.a

