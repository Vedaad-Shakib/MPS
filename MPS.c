#include <stdio.h>
#include <stdlib.h>

typedef struct {
  double x, y;
  int type;
} point;

int main() {

  double radius, dWallPoints, *wallSegments;
  int nWallSegments;
  FILE *fin;
  FILE *fout;
  fin = fopen("MPS.in", "r");
  fout = fopen("MPS.out", "w");

  // input
  fscanf(fin, "%lf", &radius);
  fscanf(fin, "%lf", &dWallPoints);
  fscanf(fin, "%d",  &nWallSegments);
  
  wallSegments = (double *) malloc(nWallSegments*sizeof(double)*4);

  for (int i = 0; i < nWallSegments; i++) {
    fscanf(fin, "%le %le %le %le", &wallSegments[4*i+0], &wallSegments[4*i+1], &wallSegments[4*i+2], &wallSegments[4*i+3]); 
  }

  return 0;
}
