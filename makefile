all:: mps run

CCFLAGS	= -std=c99 -I ./include/

mpsCorners.o: mpsCorners.c mps.h
	gcc $(CCFLAGS) -c mpsCorners.c

mpsWallPoints.o: mpsWallPoints.c mps.h
	gcc $(CCFLAGS) -c mpsWallPoints.c

mpsGhostPoints.o: mpsGhostPoints.c mps.h
	gcc $(CCFLAGS) -c mpsGhostPoints.c

mpsFluidPoints.o: mpsFluidPoints.c mps.h
	gcc $(CCFLAGS) -c mpsFluidPoints.c

mps.o: mps.c mps.h
	gcc $(CCFLAGS) -c mps.c

mps: mps.o mpsCorners.o mpsWallPoints.o mpsGhostPoints.o mpsFluidPoints.o
	gcc -o mps mps.o mpsCorners.o mpsWallPoints.o mpsGhostPoints.o mpsFluidPoints.o -lm

run:: mps
	mps

clean::
	rm -rf mps *.o *.dat mps.out *.pyc
