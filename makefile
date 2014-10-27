all:: mps run

mpsCorners.o: mpsCorners.c mps.h
	gcc -c mpsCorners.c

mpsWallPoints.o: mpsWallPoints.c mps.h
	gcc -c mpsWallPoints.c

mps.o: mps.c mps.h
	gcc -c mps.c

mps: mps.o mpsCorners.o mpsWallPoints.o
	gcc -o mps mps.o mpsCorners.o mpsWallPoints.o

run:: mps
	mps

clean::
	rm -rf mps *.o *.dat mps.out *.pyc
