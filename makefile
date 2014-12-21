all:: mps run

CCFLAGS	= -std=c99 -I ./include/ -g -static -O0

STDINC	= sys.h mem.h

mpsCorners.o: mpsCorners.c mps.h $(STDINC)
	gcc $(CCFLAGS) -c mpsCorners.c

que.o: que.c que.h $(STDINC)
	gcc $(CCFLAGS) -c que.c

mpsWallPoints.o: mpsWallPoints.c mps.h $(STDINC)
	gcc $(CCFLAGS) -c mpsWallPoints.c

mpsGhostPoints.o: mpsGhostPoints.c mps.h $(STDINC)
	gcc $(CCFLAGS) -c mpsGhostPoints.c

mpsFluidPoints.o: mpsFluidPoints.c mps.h $(STDINC)
	gcc $(CCFLAGS) -c mpsFluidPoints.c

mps.o: mps.c mps.h $(STDINC)
	gcc $(CCFLAGS) -c mps.c

mps: mps.o mpsCorners.o mpsWallPoints.o mpsGhostPoints.o mpsFluidPoints.o que.o
	gcc -o mps que.o mps.o mpsCorners.o mpsWallPoints.o mpsGhostPoints.o mpsFluidPoints.o -lm

run:: mps
	mps

clean::
	rm -rf mps *.o *.dat mps.out *.pyc mps.dSYM Log *~
