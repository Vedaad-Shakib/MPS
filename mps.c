#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>

#include "mps.h"

/*******************************************************************************
 * "outCrd": output the coordinates into a file
 *
 * Parameters:
 *            fileName: the fileName of the file to be written
 *            crd: the list of coordinates to be outputted
 *            nPoints: the length of crd
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
 * "containsPoint": checks if array corners already contains point x,y
 *
 * Parameters:
 *            corners: a pointer to a list of corners through which the function searches
 *            nCorners: the length of corners
 *            x: the x coordinate of the point to search for
 *            y: the y coordinate of the point to search for
 *******************************************************************************
 */
int containsPoint(double *corners, int nCorners, double x, double y) {
  for (int i = 0; i < nCorners; i+=2) {
    if (corners[2*i+0] == x && corners[2*i+1] == y) 
      return 1;
  }
  return 0;
}

/*******************************************************************************
 * "containsLine": checks if array corners already contains line x1 ,y1, x2, y2
 *
 * Parameters:
 *            corners: a pointer to a list of corners through which the function searches
 *            nCorners: the length of corners
 *            x1: the x coordinate of the first point
 *            y1: the y coordinate of the first point
 *            x2: the x coordinate of the second point
 *            y2: the y coordinate of the second point
 *******************************************************************************
 */
int containsLine(double *wallSegments, int nWallSegments, double x1, double y1, double x2, double y2) {
  for (int i = 0; i < nWallSegments; i++) {
    if (wallSegments[4*i+0] == x1 && wallSegments[4*i+1] == y1 && wallSegments[4*i+2] == x2 && wallSegments[4*i+3] == y2) 
      return 1;
  }
  return 0;
}

/*******************************************************************************
 * "dist": returns the distance between two points
 *
 * Parameters:
 *            x1: the x coordinate of the first point
 *            y1: the y coordinate of the first point
 *            x2: the x coordinate of the second point
 *            y2: the y coordinate of the second point
 *******************************************************************************
 */
double dist(int x1, int y1, int x2, int y2) {
  return sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
}

/*******************************************************************************
 * "mpsGhostCorners": Compute the ghost corners of a specific set of wall points
 *
 * Parameters:
 *            cornersHd: a pointer to a struct of corner points
 *            wallSegments: the specified list of wallsegments from which to calculate the corner points
 *            nWallSegments: the number of wall segments
 *            ghostPointsHd: a pointer to a struct of ghost point where the points will be stored
 *            nGhostPoints: the number of ghost points needed to be stored 
 *            radius: the radius of influence
 *******************************************************************************
 */
void mpsGhostCorners(MpsCornersHd cornersHd, double *wallSegments, int nWallSegments,
		     MpsGhostPointsHd ghostPointsHd, int nGhostPoints, double radius) {

    double	*cornerCrds;	/* corner coordinates		*/
    double	*cornerNDirs;	/* corner normal directions	*/
    double	 dx;		/* delta x			*/
    double	 dy;		/* delta y			*/
    double	 nx;		/* normal dir x			*/
    double	 ny;		/* normal dir y			*/
    double	 nx1;		/* normal dir x			*/
    double	 nx2;		/* normal dir x			*/
    double	 ny1;		/* normal dir y			*/
    double	 ny2;		/* normal dir y			*/
    double	 tmp;		/* a temporary var		*/
    int		 i1;		/* corner index 1		*/
    int		 i2;		/* corner index 2		*/

    // localize data
    cornerCrds	= cornersHd->cornerCrds;

    // allocate memory
    cornerNDirs	= memNewZero(double, 4 * nGhostPoints);

    // computer corner normal directions
    for (int i = 0; i < nWallSegments; i++) {
	i1  = mpsGetCornerId(cornersHd, wallSegments[4*i+0], wallSegments[4*i+1]);
	i2  = mpsGetCornerId(cornersHd, wallSegments[4*i+2], wallSegments[4*i+3]);
	dx  = (wallSegments[4*i+2] - wallSegments[4*i+0]);
	dy  = (wallSegments[4*i+3] - wallSegments[4*i+1]);
	tmp = sqrt(dx*dx + dy*dy);
	nx  = -dy / tmp;
	ny  = +dx / tmp;
	if (cornerNDirs[4*i1+0] == 0 && cornerNDirs[4*i1+1] == 0) {
	    cornerNDirs[4*i1+0]	= nx;
	    cornerNDirs[4*i1+1]	= ny;
	} else {
	    cornerNDirs[4*i1+2]	= nx;
	    cornerNDirs[4*i1+3]	= ny;
	}

	if (cornerNDirs[4*i2+0] == 0 && cornerNDirs[4*i2+1] == 0) {
	    cornerNDirs[4*i2+0]	= nx;
	    cornerNDirs[4*i2+1]	= ny;
	} else {
	    cornerNDirs[4*i2+2]	= nx;
	    cornerNDirs[4*i2+3]	= ny;
	}
    }

    // computer ghost corners
    for (int i = 0; i < nGhostPoints; i++) {
	nx1 = cornerNDirs[4*i+0];
	ny1 = cornerNDirs[4*i+1];
	nx2 = cornerNDirs[4*i+2];
	ny2 = cornerNDirs[4*i+3];

	if (nx2 == 0 && ny2 == 0) {
	    nx = nx1;
	    ny = ny1;
	} else {
	    tmp	= nx1 * ny2 - nx2 * ny1;
	    nx	= (+ny2 - ny1) / tmp;
	    ny	= (-nx2 + nx1) / tmp;
	}
	ghostPointsHd->ghostPointCrds[2*i+0] = cornerCrds[2*i+0] + radius * nx;
	ghostPointsHd->ghostPointCrds[2*i+1] = cornerCrds[2*i+1] + radius * ny;
	ghostPointsHd->nGhostPoints++;
    }

    // free memory
    memFree(cornerNDirs);
}

/*******************************************************************************
 * "constructIntermediatePoints": adds the intermediate points between corner points,
 *                                separated by wallSpacing units, to the outPoints array
 *
 * Parameters:
 *            cornerPoints: a pointer to an array of points that surround the intermediate points
 *            nCornerPoints: the number of points in the array of cornerPoints
 *            outPoints: a pointer to a pointer to an array that will store the outputted points
 *            nOutPoints: the number of current points in the outPoints array
 *            wallSpacing: the spacing between the points
 *            containsStart: if the first corner is included in the list of wall points - 0 or 1
 *            containsEnd: if the second corner is included in the list of wall points - 0 or 1
 *******************************************************************************
 */
void constructIntermediatePoints(double *cornerPoints, int nCornerPoints, double **outPoints,
				 int *nOutPoints, int *maxOutPoints, double wallSpacing,
				 int containsStart, int containsEnd) {
  double	tmp;		/* a temporary variable */
  int		nSegs;		/* the number of wall segments */
  double	dx;		/* the change in x */
  double	dy;		/* the change in y */
  int		maxPoints;      /* localized data */
  int           nPoints ;       /* number of points */
  double*	points ;	/* array of the points	*/

  maxPoints = *maxOutPoints;
  nPoints = *nOutPoints;
  points = *outPoints;
  
  for (int i = 0; i < nCornerPoints; i++) {
    tmp   = dist(cornerPoints[4*i+0], cornerPoints[4*i+1], 
		 cornerPoints[4*i+2], cornerPoints[4*i+3]);
    nSegs = ceil(tmp / wallSpacing);
    dx	  = (cornerPoints[4*i+2] - cornerPoints[4*i+0]) / nSegs;
    dy	  = (cornerPoints[4*i+3] - cornerPoints[4*i+1]) / nSegs;

    memResize(double, points, nPoints, nPoints + nSegs+1, maxPoints, 2);

    for (int j = 1-containsStart; j < nSegs+containsEnd; j++) {
      points[2*nPoints+0] = cornerPoints[4*i+0] + j * dx;
      points[2*nPoints+1] = cornerPoints[4*i+1] + j * dy;
      nPoints++;
    }
  }

  *maxOutPoints = maxPoints;
  *nOutPoints = nPoints;
  *outPoints = points;
}

/*******************************************************************************
 * "main": main function
 *******************************************************************************
 */
int main() {
  double		 r;	        /* the radius of influence of points		  */
  double		 wallSpacing;	/* the rounded spacing between wall points	  */
  double		*wallSegments;	/* an array of wall segments                      */
  double		 tmp;	        /* a temporary real		                  */
  double		 tmp1;	        /* a temporary real */
  double		 tmp2;	        /* a temporary real */
  double		 tmpArray[4];      /* a temporary array		                  */
  double		 dx;	        /* change in x			                  */
  double		 dx1;	        /* change in x			                  */
  double		 dx2;	        /* change in x			                  */
  double		 dy;	        /* change in y			                  */
  double		 dy1;	        /* change in x			                  */
  double		 dy2;	        /* change in x			                  */
  int			 nSegs;	        /* number of local wall intervals between points  */
  int			 nWallSegments;	/* the number of glocal wall segmnets             */
  int			 nCorners;	/* the number of global corners                   */
  MpsCornersHd		 cornersHd;	/* corners sturcture				  */
  MpsWallPointsHd	 wallPointsHd;	/* wallPoints structure                           */
  MpsGhostPointsHd       ghostPointsHd;	/* ghostPoints structure                          */
  FILE			*fin;           /* input file                                     */
  FILE			*fout;          /* output file                                    */

  fin  = fopen("mps.in", "r");
  fout = fopen("mps.out", "w");

  // input
  fscanf(fin, "%lf", &r);
  fscanf(fin, "%lf", &wallSpacing);
  fscanf(fin, "%d",  &nWallSegments);
  
  // initialize wall and ghost points
  wallSegments	 = memNew(double, sizeof(double) * 4);	// initially, four points; resize later

  cornersHd	 = mpsNewCorners();

  wallPointsHd   = mpsNewWallPoints();
  
  ghostPointsHd = mpsNewGhostPoints();
  
  for (int i = 0; i < nWallSegments; i++) {
    fscanf(fin, "%le %le %le %le", 
	   &wallSegments[4*i+0], &wallSegments[4*i+1], 
	   &wallSegments[4*i+2], &wallSegments[4*i+3]); 
  }
 
  outCrd("wall_segments.dat", wallSegments, nWallSegments*2);

  // create the corners from wall segments list
  for (int i = 0; i < nWallSegments; i++) {
    mpsGetCornerId(cornersHd, wallSegments[4*i+0], wallSegments[4*i+1]);
    mpsGetCornerId(cornersHd, wallSegments[4*i+2], wallSegments[4*i+3]);
  }

  // remove later
  outCrd("corners.dat", cornersHd->cornerCrds, cornersHd->nCorners);

  nCorners = mpsGetNCorners(cornersHd);

  // initialize the boundary corners of the ghost points
  mpsGhostCorners(cornersHd, wallSegments, nWallSegments, ghostPointsHd,
      nCorners, r);

  // remove later
  outCrd("ghost_corners.dat", ghostPointsHd->ghostPointCrds, nCorners);
  //printf("%d", ghostPointsHd->nGhostPoints);

  // initialize the wall points with all corners
  memResize(double, wallPointsHd->wallPointCrds, 0, cornersHd->nCorners+1, wallPointsHd->maxWallPoints, 2);

  for (int i = 0; i < cornersHd->nCorners; i++) {
    wallPointsHd->wallPointCrds[2*i+0] = cornersHd->cornerCrds[2*i+0];
    wallPointsHd->wallPointCrds[2*i+1] = cornersHd->cornerCrds[2*i+1];
    wallPointsHd->nWallPoints++;
  }

  // add the intermediate points between the wall corners
  constructIntermediatePoints(wallSegments, nWallSegments,
			      &(wallPointsHd->wallPointCrds),
			      &(wallPointsHd->nWallPoints),
			      &(wallPointsHd->maxWallPoints),
			      wallSpacing, 0, 0);

  // fill ghostPointsHd with the ghost points
  tmp = ghostPointsHd->nGhostPoints-1; // because it is changing in the for loop

  for (int i = 0; i < tmp; i++) {
    // check if the ghost segment is an actual ghost segment or just a segment with two points from different segments
    if (!containsLine(wallSegments, nWallSegments, cornersHd->cornerCrds[2*i+0], cornersHd->cornerCrds[2*i+1], cornersHd->cornerCrds[2*(i+1)+0], cornersHd->cornerCrds[2*(i+1)+1]))
	  continue;

    tmp1 = dist(cornersHd->cornerCrds[2*i+0], cornersHd->cornerCrds[2*i+1],
		ghostPointsHd->ghostPointCrds[2*i+0], ghostPointsHd->ghostPointCrds[2*i+1]);
    tmp2 = dist(cornersHd->cornerCrds[2*(i+1)+0], cornersHd->cornerCrds[2*(i+1)+1],
		ghostPointsHd->ghostPointCrds[2*(i+1)+0], ghostPointsHd->ghostPointCrds[2*(i+1)+1]);

    nSegs = fmax(ceil(tmp1 / wallSpacing), ceil(tmp2 / wallSpacing));

    dx1 = (ghostPointsHd->ghostPointCrds[2*i+0] - cornersHd->cornerCrds[2*i+0]) / nSegs;
    dy1 = (ghostPointsHd->ghostPointCrds[2*i+1] - cornersHd->cornerCrds[2*i+1]) / nSegs;
    dx2 = (ghostPointsHd->ghostPointCrds[2*(i+1)+0] - cornersHd->cornerCrds[2*(i+1)+0]) / nSegs;
    dy2 = (ghostPointsHd->ghostPointCrds[2*(i+1)+1] - cornersHd->cornerCrds[2*(i+1)+1]) / nSegs;

    // for each two corner and ghostPoint pairs, there are nSegs lines between them which comprise the ghost points
    for (int j = 1; j < nSegs+1; j++) {
      tmpArray[0] = cornersHd->cornerCrds[2*i+0] + j*dx1;
      tmpArray[1] = cornersHd->cornerCrds[2*i+1] + j*dy1;
      tmpArray[2] = cornersHd->cornerCrds[2*(i+1)+0] + j*dx2; // extend by dx of cornerCrds
      tmpArray[3] = cornersHd->cornerCrds[2*(i+1)+1] + j*dy2;

      constructIntermediatePoints(tmpArray, 4,
				  &(ghostPointsHd->ghostPointCrds),
				  &(ghostPointsHd->nGhostPoints),
				  &(ghostPointsHd->maxGhostPoints),
				  wallSpacing, 1, 1);
    }
  }

  outCrd("ghost_points.dat", ghostPointsHd->ghostPointCrds, ghostPointsHd->nGhostPoints);
 
  // remove later
  outCrd("wall_points.dat", wallPointsHd->wallPointCrds, wallPointsHd->nWallPoints);

  return 0;
}


