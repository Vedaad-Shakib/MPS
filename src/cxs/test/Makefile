#*******************************************************************************          
#** Copyright 2014-2014 Vedaad Shakib Inc.
#*******************************************************************************

#*******************************************************************************          
#** "makefile": makefile
#*******************************************************************************          

all:: run

include ../../../make/make.all

main.o: main.c
	$(CC) -c $(CCFLAGS) -I$(INCDIR) main.c

test: main.o
	$(LD) -o test $(LDFLAGS) main.o -lcxs $(CXSPARSELIB) $(LDLIBS)

run: test
	./test

clean::
	$(RM) main.o test

