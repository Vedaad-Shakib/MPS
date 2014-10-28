#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>

#include "mps.h"

/*******************************************************************************
 * "outCrd": output the coordinates into a file
 *******************************************************************************
 */
void outCrd(char *fileName,double *crd, int nPoints) {
    FILE *fout;

    fout = fopen(fileName, "w") ;
    for (int i = 0; i < nPoints; i++) {
	fprintf(fout, "%.16e %.16e\n", crd[2*i+0], crd[2*i+1]);
    }
    fclose(fout) ;

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
  double	 radius;	/* the radius of influence of points		  */
  double	 wallSpacing;	/* the rounded spacing between wall points	  */
  double	*wallSegments;	/* an array of wall segments                      */
  double	*corners;	/* an array of corner points		          */
  double	*wallPoints;	/* an array of wall points		          */
  double	*ghostPoints;   /* an array of ghost points			  */
  double	*cornersNDirs;	/* normal directions at corner points		  */
  double	 nx;		/* normal in x					  */
  double	 ny;		/* normal in y					  */
  double	 tmp;		/* a temporary real		                  */
  int		 n1;		/* node of a segment				  */
  int		 n2;		/* node of a segment				  */
  int		 nWallPoints;	/* number of current wall points		  */
  int		 maxWallPoints;	/* maximum wall points		                  */
  int		 maxGhostPoints;/* maximum ghost points 		          */
  double	 dx;		/* change in x			                  */
  double	 dy;		/* change in y			                  */
  double         dxGhost;	/* change in x in ghost points                    */
  double	 dyGhost;	/* change in y in ghost points                    */
  int		 nSegs;	        /* number of local wall intervals between points  */
  int		 nGhostSegs;	/* number of local intervals between ghost points */
  int		 nWallSegments; /* the number of glocal wall segmnets             */
  int		 nCorners;      /* the number of global corners                   */
  int		 nGhostPoints;  /* the number of global ghost points              */ 
  MpsCornersHd	 cornersHd;	/* corners sturcture				  */
  MpsWallPointsHd wallPointsHd; /* wallPoints structure                           */
  FILE		*fin;           /* input file                                     */
  FILE		*fout;          /* output file                                    */

  fin  = fopen("mps.in", "r");
  fout = fopen("mps.out", "w");

  // input
  fscanf(fin, "%lf", &radius);
  fscanf(fin, "%lf", &wallSpacing);
  fscanf(fin, "%d",  &nWallSegments);
  
  // initialize wall and ghost points
  nWallPoints	 = 0;
  maxWallPoints	 = 0;
  nGhostPoints   = 0;
  maxGhostPoints = 0;
  wallPoints	 = NULL;
  ghostPoints    = NULL;
  nCorners	 = 0;

  wallSegments	 = memNew(double, sizeof(double) * 4);	// initially, four points; resize later

  cornersHd	 = mpsNewCorners();

  wallPointsHd   = mpsNewWallPoints();
  
  ghostPoints	 = memNew(double, sizeof(double) * 2);	// nDims is 2 for ghost poinnnnnts

  for (int i = 0; i < nWallSegments; i++) {
    fscanf(fin, "%le %le %le %le", 
	   &wallSegments[4*i+0], &wallSegments[4*i+1], 
	   &wallSegments[4*i+2], &wallSegments[4*i+3]); 
  }

/*------------------------------------------------------------------------------
 * Create the corner points and compute their normal directions
 *------------------------------------------------------------------------------
 */

  for (int i = 0; i < nWallSegments; i++) {
    mpsGetCornerId(cornersHd, wallSegments[4*i+0], wallSegments[4*i+1]);
    mpsGetCornerId(cornersHd, wallSegments[4*i+2], wallSegments[4*i+3]);
  }
  outCrd("corners.dat", cornersHd->cornerCrds, cornersHd->nCorners);

  nCorners	= mpsGetNCorners(cornersHd);
  cornersNDirs	= memNewZero(double, 2 * nCorners);

  for (int i = 0; i < nWallSegments; i++) {
    n1	= mpsGetCornerId(cornersHd, wallSegments[4*i+0], wallSegments[4*i+1]);
    n2	= mpsGetCornerId(cornersHd, wallSegments[4*i+2], wallSegments[4*i+3]);
    dx	= (wallSegments[4*i+2] - wallSegments[4*i+0]);
    dy	= (wallSegments[4*i+3] - wallSegments[4*i+1]);
    tmp	= sqrt(dx*dx + dy*dy);
    cornersNDirs[2*n1+0]	+= -dy / tmp;
    cornersNDirs[2*n1+1]	+= +dx / tmp;
    cornersNDirs[2*n2+0]	+= -dy / tmp;
    cornersNDirs[2*n2+1]	+= +dx / tmp;
  }

  // normalize the direction vector
  for (int i = 0; i < nCorners; i++) {
    dx	= cornersNDirs[2*i+0];
    dy	= cornersNDirs[2*i+1];
    tmp	= sqrt(dx*dx + dy*dy);
    cornersNDirs[2*i+0]	/= tmp;
    cornersNDirs[2*i+1]	/= tmp;
  }

  // the following is a hack for plotting purposes.  Must be removed
  for (int i = 0; i < nCorners; i++) {
    cornersNDirs[2*i+0]	= cornersHd->cornerCrds[2*i+0] + radius * cornersNDirs[2*i+0];
    cornersNDirs[2*i+1]	= cornersHd->cornerCrds[2*i+1] + radius * cornersNDirs[2*i+1];
  }
  outCrd("ndirs.dat", cornersNDirs, cornersHd->nCorners);

  // initiate the wall points with all corners
  memResize(double, wallPointsHd->wallPointCrds, 0, cornersHd->nCorners+1, maxWallPoints, 2);
  for (int i = 0; i < cornersHd->nCorners; i++) {
    wallPointsHd->wallPointCrds[2*i+0] = cornersHd->cornerCrds[2*i+0];
    wallPointsHd->wallPointCrds[2*i+1] = cornersHd->cornerCrds[2*i+1];
    wallPointsHd->nWallPoints++;
  }

  // for each wall segment, add the intermediate points
  for (int i = 0; i < nWallSegments; i++) {
    tmp   = dist(wallSegments[4*i+0], wallSegments[4*i+1], 
		 wallSegments[4*i+2], wallSegments[4*i+3]);
    nSegs = ceil(tmp / wallSpacing);
    dx	  = (wallSegments[4*i+2] - wallSegments[4*i+0]) / nSegs;
    dy	  = (wallSegments[4*i+3] - wallSegments[4*i+1]) / nSegs;

    memResize(double, wallPointsHd->wallPointCrds, wallPointsHd->nWallPoints, wallPointsHd->nWallPoints + nSegs+1, maxWallPoints, 2);

    for (int j = 1; j < nSegs; j++) {
      wallPointsHd->wallPointCrds[2*wallPointsHd->nWallPoints+0] = wallSegments[4*i+0] + j * dx;
      wallPointsHd->wallPointCrds[2*wallPointsHd->nWallPoints+1] = wallSegments[4*i+1] + j * dy;
      wallPointsHd->nWallPoints++;
    }
  }

  outCrd("wallPoints.dat", wallPointsHd->wallPointCrds, wallPointsHd->nWallPoints);

  return 0;
}




















