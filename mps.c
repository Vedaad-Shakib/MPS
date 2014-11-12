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
 * "contains": checks if array corners already contains point x,y
 *
 * Parameters:
 *            corners: a pointer to a list of corners through which the function searches
 *            nCorners: the length of corners
 *            x: the x coordinate of the point to search for
 *            y: the y coordinate of the point to search for
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
 *            ghostCorners: a pointer to an array where the determined ghost points will be stored
 *            nCorners: the number of corners 
 *            radius: the radius of influence
 *******************************************************************************
 */
void mpsGhostCorners(MpsCornersHd cornersHd, double *wallSegments, int nWallSegments,
		     double *ghostCorners, int nCorners, double	radius) {

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
    cornerNDirs	= memNewZero(double, 4 * nCorners);

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
    for (int i = 0; i < nCorners; i++) {
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
	ghostCorners[2*i+0] = cornerCrds[2*i+0] + radius * nx;
	ghostCorners[2*i+1] = cornerCrds[2*i+1] + radius * ny;
    }

    // free memory
    memFree(cornerNDirs);
}

/*******************************************************************************
 * "main": main function
 *******************************************************************************
 */
int main() {
  double		 r;	        /* the radius of influence of points		  */
  double		 wallSpacing;	/* the rounded spacing between wall points	  */
  double		*wallSegments;	/* an array of wall segments                      */
  double		*corners;	/* an array of corner points		          */
  double		*wallPoints;	/* an array of wall points		          */
  double		*ghostPoints;	/* an array of ghost points			  */
  double		*ghostCorners;	/* an array of ghost corners		*/
  double		*cornersNDirs;	/* normal directions at corner points		  */
  double		 nx;	        /* normal in x					  */
  double		 ny;      	/* normal in y					  */
  double		 tmp;	        /* a temporary real		                  */
  int			 n1;	        /* node of a segment				  */
  int			 n2;	        /* node of a segment				  */
  int			 i1;            /* the id of a corner */
  int			 i2;            /* the id of a corner */
  int			 nWallPoints;	/* number of current wall points		  */
  int			 maxWallPoints;	/* maximum wall points		                  */
  int			 maxGhostPoints;/* maximum ghost points 		          */
  double		 dx;	        /* change in x			                  */
  double		 dy;	        /* change in y			                  */
  double		 px;	        /* the x coordinate of a point                    */
  double		 py;	        /* the y coordinate of a point                    */
  int			 nSegs;	        /* number of local wall intervals between points  */
  int			 nGhostSegs;	/* number of local intervals between ghost points */
  int			 nWallSegments;	/* the number of glocal wall segmnets             */
  int			 nCorners;	/* the number of global corners                   */
  int			 nGhostPoints;	/* the number of global ghost points              */ 
  MpsCornersHd		 cornersHd;	/* corners sturcture				  */
  MpsWallPointsHd	 wallPointsHd;	/* wallPoints structure                           */
  FILE			*fin;           /* input file                                     */
  FILE			*fout;          /* output file                                    */

  fin  = fopen("mps.in", "r");
  fout = fopen("mps.out", "w");

  // input
  fscanf(fin, "%lf", &r);
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

  // remove
  outCrd("corners.dat", cornersHd->cornerCrds, cornersHd->nCorners);

  nCorners	= mpsGetNCorners(cornersHd);

  ghostCorners	= memNewZero(double, 2 * nCorners);
  mpsGhostCorners(cornersHd, wallSegments, nWallSegments, ghostCorners,
      nCorners, r);

  // remove
  outCrd("ghost_corners.dat", ghostCorners, nCorners);

  cornersNDirs	= memNewZero(double, 2 * nCorners);

  for (int i = 0; i < nWallSegments; i++) {
    i1	= mpsGetCornerId(cornersHd, wallSegments[4*i+0], wallSegments[4*i+1]);
    i2	= mpsGetCornerId(cornersHd, wallSegments[4*i+2], wallSegments[4*i+3]);
    dx	= (wallSegments[4*i+2] - wallSegments[4*i+0]);
    dy	= (wallSegments[4*i+3] - wallSegments[4*i+1]);
    tmp	= sqrt(dx*dx + dy*dy);

    // calculate normal direction or concatenate existing normal direction
    if (cornersNDirs[2*i1+0] == 0 && cornersNDirs[2*i1+1]) {
      cornersNDirs[2*i1+0] += -dy / tmp;
      cornersNDirs[2*i1+1] += +dx / tmp;
    } else {
      double nx1, ny1, nx2, ny2, d;
      nx1 = cornersNDirs[2*i1+0];
      ny1 = cornersNDirs[2*i1+1];
      nx2 = -dy / tmp;
      ny2 = +dx / tmp;
      d = nx1*ny2 - nx2*ny1;
      
      cornersNDirs[2*i1+0] = (+ny2*r - ny1*r) / d;
      cornersNDirs[2*i1+1] = (-nx2*r + nx1*r) / d;
    }

    if (cornersNDirs[2*i2+0] == 0 && cornersNDirs[2*i2+1] == 0) {
      cornersNDirs[2*i2+0] += -dy / tmp;
      cornersNDirs[2*i2+1] += +dx / tmp;
    } else {
      double nx1, ny1, nx2, ny2, d;
      nx1 = cornersNDirs[2*i2+0];
      ny1 = cornersNDirs[2*i2+1];
      nx2 = -dy / tmp;
      ny2 = +dx / tmp;
      d = nx1*ny2 - nx2*ny1;

      cornersNDirs[2*i1+0] = (+ny2*r - ny1*r) / d;
      cornersNDirs[2*i1+1] = (-nx2*r + nx1*r) / d;
    }
  }

  /*  // normalize the direction vector
  for (int i = 0; i < nCorners; i++) {
    dx	= cornersNDirs[2*i+0];
    dy	= cornersNDirs[2*i+1];
    tmp	= sqrt(dx*dx + dy*dy);
    cornersNDirs[2*i+0]	/= tmp;
    cornersNDirs[2*i+1]	/= tmp;
    }*/

  // the following is a hack for plotting purposes.  Must be removed
  for (int i = 0; i < nCorners; i++) {
    cornersNDirs[2*i+0]	= cornersHd->cornerCrds[2*i+0] + r * cornersNDirs[2*i+0];
    cornersNDirs[2*i+1]	= cornersHd->cornerCrds[2*i+1] + r * cornersNDirs[2*i+1];
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


