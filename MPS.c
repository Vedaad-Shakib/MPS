#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>

typedef struct {
  double x, y;
  int type;
} point;

/*******************************************************************************
 * Memory allocation macros
 *******************************************************************************
 */
#define	memNew(T,N)	(T*) malloc(sizeof(T) * N)

#define	memResize(TYPE, ARRAY, OLDSIZE, NEWSIZE, MAXSIZE, NDIMS)	\
  if ((NEWSIZE) >= (MAXSIZE)) {						\
    TYPE *ary;								\
    (MAXSIZE) = 2 * (MAXSIZE) + 100;					\
    if ((MAXSIZE) < (NEWSIZE)) (MAXSIZE) = (NEWSIZE);			\
    ary = memNew(TYPE, (MAXSIZE)*(NDIMS));				\
    if ((OLDSIZE) != 0) {						\
      memcpy(ary, ARRAY, sizeof(TYPE)*(OLDSIZE)*(NDIMS));		\
      free(ARRAY);							\
    }									\
    ARRAY = ary;							\
  }

/*******************************************************************************
 * "contains": checks if array corners already contains point x,y
 *******************************************************************************
 */
int contains(double *corners, int nCorners, double x, double y) {
  for (int i = 0; i < nCorners; i++) {
    if (corners[2*i+0] == x && corners[2*i+1] == y) 
      return 1;
  }
  return 0;
}

/*******************************************************************************
 * "dist": returns the distance between two points
 *******************************************************************************
 */
double dist(int x1, int y1, int x2, int y2) {
  return sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
}

/*******************************************************************************
 * "main": main routine
 *******************************************************************************
 */
int main() {
  double	 radius;	/* the radius of influence of points		*/
  double	 wallSpacing;	/* the rounded spacing between wall points	*/
  double	*wallSegments;	/* an array of wall segments                    */
  double        *corners;	/* an array of corner points		        */
  double	*wallPoints;	/* an array of wall points		        */
  double	 tmp;		/* a temporary real		                */
  int		 nWallPoints;	/* number of current wall points		*/
  int		 maxWallPoints;	/* maximum wall points		                */
  double	 dx;		/* change in x			                */
  double	 dy;		/* change in y			                */
  int		 i;		/* index 		                        */
  int		 j;		/* index		                        */
  int		 nSegs;	        /* number of local wall intervals between points*/
  int		 nWallSegments; /* the number of glocal wall segmnets           */
  int            nCorners;      /* the number of global corners                 */
  FILE		*fin;
  FILE		*fout;
  fin  = fopen("MPS.in", "r");
  fout = fopen("MPS.out", "w");

  // input
  fscanf(fin, "%lf", &radius);
  fscanf(fin, "%lf", &wallSpacing);
  fscanf(fin, "%d",  &nWallSegments);
  
  // [x1, y1, x2, y2, x3, y3, x4, y4 ...]
  // initialize wall points
  nWallPoints	 = 0;
  maxWallPoints	 = 0;
  wallPoints	 = NULL;
  nCorners	 = 0;

  wallSegments	 = memNew(double, sizeof(double) * 4);

  corners	 = memNew(double, sizeof(double) * 4);	// max size of corners

  for (i = 0; i < nWallSegments; i++) {
    fscanf(fin, "%le %le %le %le", 
	    &wallSegments[4*i+0], &wallSegments[4*i+1], 
	    &wallSegments[4*i+2], &wallSegments[4*i+3]); 
  }

  // find number of points needed (to allocate)
  for (i = 0; i < nWallSegments; i++) {
      tmp	 = dist(wallSegments[4*i+0], wallSegments[4*i+1], 
	                wallSegments[4*i+2], wallSegments[4*i+3]);
      nSegs	 = ceil(tmp / wallSpacing);

      memResize(double, wallPoints, nWallPoints, nWallPoints+nSegs+1, maxWallPoints, 4);
      if (!contains(corners, nCorners, wallSegments[4*i+0], wallSegments[4*i+1])) {
          corners[2*nCorners+0] = wallSegments[4*i+0];
          corners[2*nCorners+1] = wallSegments[4*i+1];
          nCorners++;
          wallPoints[2*nWallPoints+0] = wallSegments[4*i+0];
          wallPoints[2*nWallPoints+1] = wallSegments[4*i+1];
          nWallPoints++;
      }
      dx	 = (wallSegments[4*i+2] - wallSegments[4*i+0]) / nSegs;
      dy	 = (wallSegments[4*i+3] - wallSegments[4*i+1]) / nSegs;
      for (j = 1; j < nSegs; j++) {
          wallPoints[2*nWallPoints+0] = wallSegments[4*i+0] + j * dx;
          wallPoints[2*nWallPoints+1] = wallSegments[4*i+1] + j * dy;
          nWallPoints++;
      }
      if (!contains(corners, nCorners, wallSegments[4*i+2], wallSegments[4*i+3])) {
          corners[2*nCorners+0] = wallSegments[4*i+2];
          corners[2*nCorners+1] = wallSegments[4*i+3];
          wallPoints[2*nWallPoints+0] = wallSegments[4*i+2];
          wallPoints[2*nWallPoints+1] = wallSegments[4*i+3];
          nCorners++;
          nWallPoints++;
      }
  }
  printf("wall points = %d\n", nWallPoints);

  for (int i = 0; i < nWallPoints; i++) {
    printf("%d %le, %le\n", i, wallPoints[2*i+0], wallPoints[2*i+1]);
  }

  return 0;
}




















