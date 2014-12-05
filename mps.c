/*******************************************************************************          
** Copyright 2014-2014 Vedaad Shakib Inc.
*******************************************************************************/

/*******************************************************************************
** 
** "mps.c": Moving Particle Semi-implicit
**
*******************************************************************************/

/*******************************************************************************
 * Standard includes
*******************************************************************************/

#include "sys.h"
#include "mps.h"
#include "que.h"

/*******************************************************************************
 * "mpsOutCrd": output the coordinates into a file
 *
 * Parameters:
 *            fileName: the fileName of the file to be written
 *            crd: the list of coordinates to be outputted
 *            nPoints: the length of crd
 *******************************************************************************
 */
void mpsOutCrd(char *fileName,double *crd, int nPoints) {
    FILE *fout;

    fout = fopen(fileName, "w");
    for (int i = 0; i < nPoints; i++) {
	fprintf(fout, "%.16e %.16e\n", crd[2*i+0], crd[2*i+1]);
    }
    fclose(fout);
} 

/*******************************************************************************
 * "mpsContainsLine": checks if array corners already contains line x1 ,y1, x2, y2
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
int mpsContainsLine(double *wallSegments, int nWallSegments, double x1, double y1, double x2, double y2) {
  for (int i = 0; i < nWallSegments; i++) {
    if (wallSegments[4*i+0] == x1 && wallSegments[4*i+1] == y1 && wallSegments[4*i+2] == x2 && wallSegments[4*i+3] == y2) 
      return 1;
  }
  return 0;
}

/*******************************************************************************
 * "mpsDist": returns the mpsDistance between two points
 *
 * Parameters:
 *            x1: the x coordinate of the first point
 *            y1: the y coordinate of the first point
 *            x2: the x coordinate of the second point
 *            y2: the y coordinate of the second point
 *******************************************************************************
 */
double mpsDist(int x1, int y1, int x2, int y2) {
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

    // compute ghost corners, combining the normal directions of two wall segments on one endpoint
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
 * "mpsConstructIntermediatePoints": adds the intermediate points between corner points,
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
void mpsConstructIntermediatePoints(double *cornerPoints, int nCornerPoints, double **outPoints,
				    int *nOutPoints, int *maxOutPoints, double wallSpacing,
				    bool containsStart, bool containsEnd) {
  double	tmp;		/* a temporary variable */
  int		nSegs;		/* the number of wall segments */
  double	dx;		/* the change in x */
  double	dy;		/* the change in y */
  int		maxPoints;      /* localized data */
  int           nPoints;        /* number of points */
  double*	points;  	/* array of the points	*/

  maxPoints = *maxOutPoints;
  nPoints = *nOutPoints;
  points = *outPoints;
  
  for (int i = 0; i < nCornerPoints; i++) {
    tmp   = mpsDist(cornerPoints[4*i+0], cornerPoints[4*i+1], 
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
 * "mpsCheckClosure": checks whether a list of boundaries are closed by ensuring
 *                    that every start connects to one and only one end and that
 *                    every end connects to one and only one start
 *
 * Parameters:
 *            fluidBoundaries: a list of boundaries
 *            nFluidBoundaries: the number of boundaries
 *******************************************************************************
 */
bool mpsCheckClosure(double *fluidBoundaries, int nFluidBoundaries) {
  double	x;		/* x coordinate             */
  double 	y;		/* y coordinate             */
  int		counter;	/* a counter                */

  // check that every start point connects with only one end point
  for (int i = 0; i < nFluidBoundaries; i++) {
    x = fluidBoundaries[4*i+0];
    y = fluidBoundaries[4*i+1];
    counter = 0;
    for (int j = 0; j < nFluidBoundaries; j++) {
      if (x == fluidBoundaries[4*j+2] && y == fluidBoundaries[4*j+3])
	counter++;
    } 
    if (counter != 1) return false;
  }

  // check that every end point connects with only one start point
  for (int i = 0; i < nFluidBoundaries; i++) {
    x = fluidBoundaries[4*i+2];
    y = fluidBoundaries[4*i+3];
    counter = 0;
    for (int j = 0; j < nFluidBoundaries; j++) {
      if (x == fluidBoundaries[4*j+0] && y == fluidBoundaries[4*j+1])
	counter++;
    } 
    if (counter != 1) return false;
  }

  return true;
}

/*******************************************************************************
 * "mpsIntegerize": converts a double coordinate point into an integer value
 *
 * Parameters:
 *            min: the minimum coordinate
 *            val: the value to mpsIntegerize
 *            wallSpacing: the wall spacing
 *******************************************************************************
 */
int mpsIntegerize(double min, double val, double wallSpacing) {
  return (int) round((val-min)/wallSpacing);
}

/*******************************************************************************
 * "mpsCrossesFluidBoundary": checks whether a line crosses  one fluid boundary
 *
 * Parameters:
 *            a1x: the x coordinate of the first point of the line 
 *            a2x: the y coordinate of the first point of the line 
 *            a1y: the x coordinate of the second point of the line 
 *            a2y: the y coordinate of the second point of the line 
 *            a1x: the x coordinate of the first point of the boundary
 *            a2x: the y coordinate of the first point of the boundary
 *            a1y: the x coordinate of the second point of the boundary
 *            a2y: the y coordinate of the second point of the boundary
 *******************************************************************************
 */
bool mpsCrossesFluidBoundary(double *fluidBoundaries, int nFluidBoundaries,
			  double x1, double y1, double x2, double y2) {
  for (int i = 0; i < nFluidBoundaries; i++) {
    if (mpsLinesIntersect(x1, y1, x2, y2,
		       fluidBoundaries[4*i+0], fluidBoundaries[4*i+1],
		       fluidBoundaries[4*i+2], fluidBoundaries[4*i+3])) {
      return true;
    }
  }
  return false;
}

/*******************************************************************************
 * "mpsLinesIntersect": checks whether a line crosses another line
 *
 * Parameters:
 *            a1x: the x coordinate of the first point of the line 
 *            a1y: the x coordinate of the second point of the line 
 *            a2x: the y coordinate of the first point of the line 
 *            a2y: the y coordinate of the second point of the line 
 *            a1x: the x coordinate of the first point of the second line
 *            a1y: the x coordinate of the second point of the second line
 *            a2x: the y coordinate of the first point of the second line
 *            a2y: the y coordinate of the second point of the second line
 *******************************************************************************
 */
bool mpsLinesIntersect(double a1x, double a1y, double a2x, double a2y,
		       double b1x, double b1y, double b2x, double b2y) {
  double dax; /* the ax delta */
  double dbx; /* the bx delta */
  double day; /* the ay delta */
  double dby; /* the by delta */
  double dbax; /* b1x-a1x */
  double dbay; /* b1y-a1y */
  double det; /* the determinant */
  double alpha; /* mpsDistance coefficient */
  double beta; /* mpsDistance coefficient */
  double tol; /* check tolerance */

  tol = 1.e-12;

  dax = a2x-a1x;
  day = a2y-a1y;
  dbx = b1x-b2x;
  dby = b1y-b2y;
  dbax = b1x-a1x;
  dbay = b1y-a1y;
  det = dax*dby-dbx*day;

  // if determinant is 0, they're parallel (and the rest will return error)
  if (fabs(det) <= tol * (fabs(dax)+fabs(day))*(fabs(dbx)+fabs(dby))) {
      return(false);
  }

  alpha = (+dby*dbax-dbx*dbay)/det;
  beta  = (-day*dbax+dax*dbay)/det;
  return (alpha > -tol && alpha < 1+tol && beta > -tol && beta < 1+tol);
}

/*******************************************************************************
 * "mpsFloodfill": iteratively fills a boundary with points
 *
 * Parameters:
 *            fluidPointsHd: an MpsFluidPointsHd struct which will contain the created points
 *            fluidBoundaries: a pointer to an array of boundaries containing the fluid
 *            nFluidBoundaries: the number of fluid boundaries
 *            initX: the x value from where floodfill starts
 *            initY: the y value from where floodfill starts
 *            wallSpacing: the mpsDistance between the points
 *******************************************************************************
 */
void mpsFloodfill(MpsFluidPointsHd fluidPointsHd, double *fluidBoundaries, int nFluidBoundaries,
		  double initX, double initY, double wallSpacing) {
  double	 x;		/* the x coordinate */
  double	 y;		/* the y coordinate */
  int		 xInt;          /* an mpsIntegerized x */
  int		 yInt;          /* an mpsIntegerized y */
  double	 xMin;		/* the minimum x */
  double	 xMax;		/* the maximum x */
  double	 yMin;		/* the minimum y */
  double	 yMax;		/* the maximum y */
  bool          *visited;       /* a 2D array of visited points */
  int            xDim;          /* the number of x dimensions in the visited array */
  int            yDim;          /* the number of y dimensions in the visited array */
  QueHd        	 queHd;         /* queue to store points needed to be visited */

  // get the minimum and maximum x and y coordinates
  xMin = fluidBoundaries[0];
  xMax = fluidBoundaries[0];
  yMin = fluidBoundaries[1];
  yMax = fluidBoundaries[1];

  for (int i = 0; i < nFluidBoundaries; i++) {
    if (fluidBoundaries[4*i+0] < xMin) xMin = fluidBoundaries[4*i+0];
    if (fluidBoundaries[4*i+0] > xMax) xMax = fluidBoundaries[4*i+0];
    if (fluidBoundaries[4*i+1] < yMin) yMin = fluidBoundaries[4*i+1];
    if (fluidBoundaries[4*i+1] > yMax) yMax = fluidBoundaries[4*i+1];

    if (fluidBoundaries[4*i+2] < xMin) xMin = fluidBoundaries[4*i+2];
    if (fluidBoundaries[4*i+2] > xMax) xMax = fluidBoundaries[4*i+2];
    if (fluidBoundaries[4*i+3] < yMin) yMin = fluidBoundaries[4*i+3];
    if (fluidBoundaries[4*i+3] > yMax) yMax = fluidBoundaries[4*i+3];
  }

  // make a 2D array to store whether points are visited or unvisited
  xDim = (int) round((xMax-xMin)/wallSpacing)+1;
  yDim = (int) round((yMax-yMin)/wallSpacing)+1;

  visited = memNewZero(bool, xDim * yDim);

  // initialize queue
  // using queue and while loop because recursion with a depth of nFluidPoints is way too much overhead
  queHd = queNew();
  quePush(queHd, initX);
  quePush(queHd, initY);

  while (!queIsEmpty(queHd)) {
    x = quePop(queHd);
    y = quePop(queHd);
    xInt = mpsIntegerize(xMin, x, wallSpacing);
    yInt = mpsIntegerize(yMin, y, wallSpacing);
    visited[xInt+yInt*xDim] = true;

    // add the point
    memResize(double, fluidPointsHd->fluidPointCrds, fluidPointsHd->nFluidPoints, fluidPointsHd->nFluidPoints+1, fluidPointsHd->maxFluidPoints, 2);
    fluidPointsHd->fluidPointCrds[2*fluidPointsHd->nFluidPoints+0] = x;
    fluidPointsHd->fluidPointCrds[2*fluidPointsHd->nFluidPoints+1] = y;
    fluidPointsHd->nFluidPoints++;

    // check and add the four adjacent points to the queue
    if (!mpsCrossesFluidBoundary(fluidBoundaries, nFluidBoundaries, x, y, x+wallSpacing, y) &&
	!visited[(xInt+1)+yInt*xDim]) {
      quePush(queHd, x+wallSpacing);
      quePush(queHd, y);
      visited[(xInt+1)+yInt*xDim] = true;
    }
    if (!mpsCrossesFluidBoundary(fluidBoundaries, nFluidBoundaries, x, y, x, y+wallSpacing) &&
	!visited[xInt+(yInt+1)*xDim]) {
      quePush(queHd, x);
      quePush(queHd, y+wallSpacing);
      visited[xInt+(yInt+1)*xDim] = true;
    }
    if (!mpsCrossesFluidBoundary(fluidBoundaries, nFluidBoundaries, x, y, x-wallSpacing, y) &&
        !visited[(xInt-1)+(yInt*xDim)]) {
      quePush(queHd, x-wallSpacing);
      quePush(queHd, y);
      visited[(xInt-1)+(yInt*xDim)] = true;
    }
    if (!mpsCrossesFluidBoundary(fluidBoundaries, nFluidBoundaries, x, y, x, y-wallSpacing) &&
	!visited[xInt+(yInt-1)*xDim]) {
      quePush(queHd, x);
      quePush(queHd, y-wallSpacing);
      visited[xInt+(yInt-1)*xDim] = true;
    }
  }
  
  queFree(queHd);
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
  double		 tmp1;	        /* a temporary real                               */
  double		 tmp2;	        /* a temporary real                               */
  double		 tmpArray[4];   /* a temporary array		                  */
  double		 dx;	        /* change in x			                  */
  double		 dx1;	        /* change in x			                  */
  double		 dx2;	        /* change in x			                  */
  double		 dy;	        /* change in y			                  */
  double		 dy1;	        /* change in x			                  */
  double		 dy2;	        /* change in x			                  */
  double                 x0;            /* the initial x                                  */
  double                 y0;            /* the initial y                                  */
  int			 nSegs;	        /* number of local wall intervals between points  */
  int			 nWallSegments;	/* the number of glocal wall segmnets             */
  int			 nCorners;	/* the number of global corners                   */
  int                    nFluidBoundaries; /* the number of fluid boundaries              */
  double                *fluidBoundaries;  /* an array of the boundaries enclosing the fluid */
  MpsCornersHd		 cornersHd;	/* corners sturcture				  */
  MpsWallPointsHd	 wallPointsHd;	/* wallPoints structure                           */
  MpsGhostPointsHd       ghostPointsHd;	/* ghostPoints structure                          */
  MpsFluidPointsHd       fluidPointsHd; /* fluidPoints structure                          */
  FILE			*fin;           /* input file                                     */
  FILE			*fout;          /* output file                                    */

  fin  = fopen("mps.in", "r");
  fout = fopen("mps.out", "w");

  /*=======================================================================================
   * Input and create wall points and ghost points
   *=======================================================================================
   */

  // input radius, wallSpacing, wallSegments
  fscanf(fin, "%lf", &r);
  fscanf(fin, "%lf", &wallSpacing);
  fscanf(fin, "%d",  &nWallSegments);
  
  // initialize wall and ghost points
  wallSegments	 = memNew(double, nWallSegments * 4); 
  cornersHd	 = mpsNewCorners();
  wallPointsHd   = mpsNewWallPoints();
  ghostPointsHd  = mpsNewGhostPoints();

  for (int i = 0; i < nWallSegments; i++) {
    fscanf(fin, "%le %le %le %le", 
	   &wallSegments[4*i+0], &wallSegments[4*i+1], 
	   &wallSegments[4*i+2], &wallSegments[4*i+3]); 
  }
 
  mpsOutCrd("wall_segments.dat", wallSegments, nWallSegments*2);



  // create the corners from wall segments list
  for (int i = 0; i < nWallSegments; i++) {
    mpsGetCornerId(cornersHd, wallSegments[4*i+0], wallSegments[4*i+1]);
    mpsGetCornerId(cornersHd, wallSegments[4*i+2], wallSegments[4*i+3]);
  }

  // remove later
  mpsOutCrd("corners.dat", cornersHd->cornerCrds, cornersHd->nCorners);

  nCorners = mpsGetNCorners(cornersHd);

  // initialize the boundary corners of the ghost points
  mpsGhostCorners(cornersHd, wallSegments, nWallSegments, ghostPointsHd,
      nCorners, r);

  // remove later
  mpsOutCrd("ghost_corners.dat", ghostPointsHd->ghostPointCrds, nCorners);

  // initialize the wall points with all corners
  memResize(double, wallPointsHd->wallPointCrds, 0, cornersHd->nCorners+1, wallPointsHd->maxWallPoints, 2);

  for (int i = 0; i < cornersHd->nCorners; i++) {
    wallPointsHd->wallPointCrds[2*i+0] = cornersHd->cornerCrds[2*i+0];
    wallPointsHd->wallPointCrds[2*i+1] = cornersHd->cornerCrds[2*i+1];
    wallPointsHd->nWallPoints++;
  }

  // add the intermediate points between the wall corners
  mpsConstructIntermediatePoints(wallSegments, nWallSegments,
				 &(wallPointsHd->wallPointCrds),
				 &(wallPointsHd->nWallPoints),
				 &(wallPointsHd->maxWallPoints),
				 wallSpacing, false, false);

  // fill ghostPointsHd with the ghost points
  tmp = ghostPointsHd->nGhostPoints-1; // because it is changing in the for loop

  for (int i = 0; i < tmp; i++) {
    // check if the ghost segment is an actual ghost segment or just a segment with two points from different segments
    if (!mpsContainsLine(wallSegments, nWallSegments, cornersHd->cornerCrds[2*i+0], cornersHd->cornerCrds[2*i+1], cornersHd->cornerCrds[2*(i+1)+0], cornersHd->cornerCrds[2*(i+1)+1]))
	  continue;

    tmp1 = mpsDist(cornersHd->cornerCrds[2*i+0], cornersHd->cornerCrds[2*i+1],
		ghostPointsHd->ghostPointCrds[2*i+0], ghostPointsHd->ghostPointCrds[2*i+1]);
    tmp2 = mpsDist(cornersHd->cornerCrds[2*(i+1)+0], cornersHd->cornerCrds[2*(i+1)+1],
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

      mpsConstructIntermediatePoints(tmpArray, 4,
				     &(ghostPointsHd->ghostPointCrds),
				     &(ghostPointsHd->nGhostPoints),
				     &(ghostPointsHd->maxGhostPoints),
				     wallSpacing, true, true);
    }
  }

  mpsOutCrd("ghost_points.dat", ghostPointsHd->ghostPointCrds, ghostPointsHd->nGhostPoints);
 
  // remove later
  mpsOutCrd("wall_points.dat", wallPointsHd->wallPointCrds, wallPointsHd->nWallPoints);

  /*=======================================================================================
   * Input and create fluid points
   *=======================================================================================
   */

  // input fluid boundaries
  fscanf(fin, "%d", &nFluidBoundaries);

  // initialize fluidBoundaries and fluidPoints
  fluidBoundaries = memNew(double, nFluidBoundaries*4);
  fluidPointsHd	  = mpsNewFluidPoints();

  for (int i = 0; i < nFluidBoundaries; i++) {
    fscanf(fin, "%le %le %le %le", 
	   &fluidBoundaries[4*i+0], &fluidBoundaries[4*i+1], 
	   &fluidBoundaries[4*i+2], &fluidBoundaries[4*i+3]); 
  }

  // check if the boundary is closed
  if (!mpsCheckClosure(fluidBoundaries, nFluidBoundaries)) {
    printf("The given fluid boundaries are not closed.");
    exit(0);
  } 
  
  // add starting point for floodfill
  dx  = (fluidBoundaries[2] - fluidBoundaries[0]);
  dy  = (fluidBoundaries[3] - fluidBoundaries[1]);
  tmp = sqrt(dx*dx + dy*dy);
  x0  = fluidBoundaries[0] + dx/2 + wallSpacing*(-dy/tmp);
  y0  = fluidBoundaries[1] + dy/2 + wallSpacing*(+dx/tmp);

  // fill the boundary through iterative floodfill
  mpsFloodfill(fluidPointsHd, fluidBoundaries, nFluidBoundaries,
	       x0, y0, wallSpacing);

  mpsOutCrd("fluid_points.dat", fluidPointsHd->fluidPointCrds, fluidPointsHd->nFluidPoints);

  return 0;
}
 
 
