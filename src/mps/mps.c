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
#include "stn.h"
#include "slv.h"

/*******************************************************************************
 * "mpsOutCrd": output the coordinates into a file
 *
 * Parameters:
 *            fileName: the fileName of the file to be written
 *            crd: the list of coordinates to be outputted
 *            nPoints: the length of crd
 *            nDims: the number of dimensions to the data
 *******************************************************************************
 */
void mpsOutCrd(char *fileName, double *crd, int nPoints,
               int   nDims) {
    FILE *fout;
    
    fout = fopen(fileName, "w");
    for (int i = 0; i < nPoints; i++) {
        for (int j = 0; j < nDims; j++) {
            fprintf(fout, "%.16e ", crd[nDims*i+j]);
        }
        fprintf(fout, "\n");
    }
    fclose(fout);
} 

/*******************************************************************************
 * "mpsOutCrdXY": output the coordinates into a file given two arrays of x, y
 *
 * Parameters:
 *            fileName: the fileName of the file to be written
 *            xCrd: the list of x coordinates to be outputted
 *            yCrd: the list of y coordinates to be outputted
 *            nPoints: the length of crd
 *******************************************************************************
 */
void mpsOutCrdXY(char *fileName, double *xCrd, double *yCrd,
                 int   nPoints) {
    FILE *fout;
    
    fout = fopen(fileName, "w");
    for (int i = 0; i < nPoints; i++) {
        fprintf(fout, "%.16e ", xCrd[i]);
        fprintf(fout, "%.16e ", yCrd[i]);
        fprintf(fout, "\n");
    }
    fclose(fout);
} 

/*******************************************************************************
 * "mpsVecL2": compute the L2 of a vector, for printing
 *******************************************************************************
 */
double mpsVecL2(double *vec, int nDims) {
    double        l2Norm;                        /* L2 norm */

    l2Norm = 0;
    for (int i = 0; i < nDims; i++) {
        l2Norm += vec[i] * vec[i];
    }
    return sqrt(l2Norm/nDims) ;
} 

/*******************************************************************************
 * "mps2VecL2": compute the L2 of a vector, for printing
 *******************************************************************************
 */
double mps2VecL2(double *vec1, double *vec2, int nDims) {
    double        l2Norm;                        /* L2 norm */

    l2Norm = 0;
    for (int i = 0; i < nDims; i++) {
        l2Norm += vec1[i] * vec1[i] + vec2[i] * vec2[i];
    }
    return sqrt(l2Norm/nDims) ;
} 

/*******************************************************************************
 * "mpsContainsLine": checks if array corners already contains line x1 ,y1, x2, y2
 *
 * Parameters:
 *            corners: a pointer to a list of corners through which the function searches
 *            nPoints: the length of corners
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
double mpsDist(double x1, double y1, double x2, double y2) {
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
void mpsGhostCorners(MpsPointsHd cornersHd, double *wallSegments, int nWallSegments,
                     MpsPointsHd ghostPointsHd, int nGhostPoints, double radius) {
    
    double      *pointCrds;    /* corner coordinates */
    double      *cornerNDirs;   /* corner normal directions */
    double       dx;            /* delta x */
    double       dy;            /* delta y */
    double       nx;            /* normal dir x */
    double       ny;            /* normal dir y */
    double       nx1;           /* normal dir x */
    double       nx2;           /* normal dir x */
    double       ny1;           /* normal dir y */
    double       ny2;           /* normal dir y */
    double       tmp;           /* a temporary var */
    int          i1;            /* corner index 1 */
    int          i2;            /* corner index 2 */
    
    // localize data
    pointCrds = cornersHd->pointCrds;
    
    // allocate memory
    cornerNDirs = memNewZero(double, 4 * nGhostPoints);
    
    // computer corner normal directions
    for (int i = 0; i < nWallSegments; i++) {
        i1  = mpsGetPointId(cornersHd, wallSegments[4*i+0], wallSegments[4*i+1]);
        i2  = mpsGetPointId(cornersHd, wallSegments[4*i+2], wallSegments[4*i+3]);
        dx  = (wallSegments[4*i+2] - wallSegments[4*i+0]);
        dy  = (wallSegments[4*i+3] - wallSegments[4*i+1]);
        tmp = sqrt(dx*dx + dy*dy);
        nx  = -dy / tmp;
        ny  = +dx / tmp;
        if (cornerNDirs[4*i1+0] == 0 && cornerNDirs[4*i1+1] == 0) {
            cornerNDirs[4*i1+0] = nx;
            cornerNDirs[4*i1+1] = ny;
        } else {
            cornerNDirs[4*i1+2] = nx;
            cornerNDirs[4*i1+3] = ny;
        }
        
        if (cornerNDirs[4*i2+0] == 0 && cornerNDirs[4*i2+1] == 0) {
            cornerNDirs[4*i2+0] = nx;
            cornerNDirs[4*i2+1] = ny;
        } else {
            cornerNDirs[4*i2+2] = nx;
            cornerNDirs[4*i2+3] = ny;
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
            tmp = nx1 * ny2 - nx2 * ny1;
            nx  = (+ny2 - ny1) / tmp;
            ny  = (-nx2 + nx1) / tmp;
        }
        ghostPointsHd->pointCrds[2*i+0] = pointCrds[2*i+0] + radius * nx;
        ghostPointsHd->pointCrds[2*i+1] = pointCrds[2*i+1] + radius * ny;
        ghostPointsHd->nPoints++;
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
 *            nSegs: the number of points this segment should have; -1 if calculated from cornerPoints
 *            outPoints: a pointer to a pointer to an array that will store the outputted points
 *            nOutPoints: the number of current points in the outPoints array
 *            wallSpacing: the spacing between the points
 *            containsStart: if the first corner is included in the list of wall points - 0 or 1
 *            containsEnd: if the second corner is included in the list of wall points - 0 or 1
 *******************************************************************************
 */
void mpsConstructIntermediatePoints(double *cornerPoints, int nCornerPoints, int nSegs,
                                    MpsPointsHd pointsHd, 
                                    double wallSpacing, bool containsStart, bool containsEnd) {
    double       tmp;              /* a temporary variable */
    double       dx;               /* the change in x */
    double       dy;               /* the change in y */
    bool         isNSegsProvided;  /* if the number of points is calculated using nSegs */
    
    isNSegsProvided = (nSegs == -1) ? false : true;
    
    for (int i = 0; i < nCornerPoints-1; i++) {
        if (!isNSegsProvided) {
            tmp   = mpsDist(cornerPoints[2*i+0], cornerPoints[2*i+1], 
                            cornerPoints[2*(i+1)+0], cornerPoints[2*(i+1)+1]);
            nSegs = ceil(tmp / wallSpacing);
        }

        dx = (cornerPoints[2*(i+1)+0] - cornerPoints[2*i+0]) / nSegs;
        dy = (cornerPoints[2*(i+1)+1] - cornerPoints[2*i+1]) / nSegs;
        
        for (int j = 1-containsStart; j < nSegs+containsEnd; j++) {
            mpsGetPointId(pointsHd, cornerPoints[2*i+0] + j * dx, cornerPoints[2*i+1] + j * dy);
        }
    }
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
    double      x;              /* x coordinate             */
    double      y;              /* y coordinate             */
    int         counter;        /* a counter                */
    
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
 *            val: the value to integerize
 *            wallSpacing: the wall spacing
 *******************************************************************************
 */
int mpsIntegerize(double min, double val, double wallSpacing) {
    return (int) round((val-min)/wallSpacing);
}

/*******************************************************************************
 * "mpsCrossesFluidBoundary": checks whether a line segment crosses any fluid boundaries or wall segments
 *
 * Parameters:
 *            fluidBoundaries: the list of fluid boundaries
 *            nFluidBoundaries: the number of fluid boundaries
 *            wallSegments: the list of wall segments
 *            nWallSegments: the number of wall segments
 *            x1: the x coordinate of the first point of the line
 *            y1: the y coordinate of the first point of the line
 *            x2: the x coordinate of the second point of the line
 *            y2: the y coordinate of the second point of the line
 *******************************************************************************
 */
bool mpsCrossesFluidBoundary(double *fluidBoundaries, int nFluidBoundaries,
                             double *wallSegments, int nWallSegments,
                             double x1, double y1, double x2, double y2) {
    for (int i = 0; i < nFluidBoundaries; i++) {
        if (mpsLinesIntersect(x1, y1, x2, y2,
                              fluidBoundaries[4*i+0], fluidBoundaries[4*i+1],
                              fluidBoundaries[4*i+2], fluidBoundaries[4*i+3])) {
            return true;
        }
    }

    for (int i = 0; i < nWallSegments; i++) {
        if (mpsLinesIntersect(x1, y1, x2, y2,
                              wallSegments[4*i+0], wallSegments[4*i+1],
                              wallSegments[4*i+2], wallSegments[4*i+3])) {
            return true;
        }
    }

    return false;
}

/*******************************************************************************
 * "mpsLinesIntersect": checks whether a line segment crosses another line segment
 *
 * Parameters:
 *            a1x: the x coordinate of the first point of the first line 
 *            a1y: the y coordinate of the first point of the first line 
 *            a2x: the x coordinate of the second point of the first line 
 *            a2y: the y coordinate of the second point of the first line 
 *            b1x: the x coordinate of the first point of the second line
 *            b1y: the y coordinate of the first point of the second line
 *            b2x: the x coordinate of the second point of the second line
 *            b2y: the y coordinate of the second point of the second line
 *******************************************************************************
 */
bool mpsLinesIntersect(double a1x, double a1y, double a2x, double a2y,
                       double b1x, double b1y, double b2x, double b2y) {
    double      dax;            /* the ax delta */
    double      dbx;            /* the bx delta */
    double      day;            /* the ay delta */
    double      dby;            /* the by delta */
    double      dbax;           /* b1x-a1x */
    double      dbay;           /* b1y-a1y */
    double      det;            /* the determinant */
    double      alpha;          /* mpsDistance coefficient */
    double      beta;           /* mpsDistance coefficient */
    double      tol;            /* tolerance for floating-point error */
    
    tol = 1.e-12;
    
    dax  = a2x-a1x;
    day  = a2y-a1y;
    dbx  = b1x-b2x;
    dby  = b1y-b2y;
    dbax = b1x-a1x;
    dbay = b1y-a1y;
    det  = dax*dby-dbx*day;
    
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
void mpsFloodfill(MpsPointsHd fluidPointsHd, double *fluidBoundaries,
                  int nFluidBoundaries, double* wallSegments, int nWallSegments,
                  double initX, double initY, double wallSpacing) {
    double       x;             /* the x coordinate */
    double       y;             /* the y coordinate */
    int          xInt;          /* an mpsIntegerized x */
    int          yInt;          /* an mpsIntegerized y */
    double       xMin;          /* the minimum x */
    double       xMax;          /* the maximum x */
    double       yMin;          /* the minimum y */
    double       yMax;          /* the maximum y */
    bool        *visited;       /* a 2D array of visited points */
    int          xDim;          /* the number of x dimensions in the visited array */
    int          yDim;          /* the number of y dimensions in the visited array */
    QueHd        queHd;         /* queue to store points needed to be visited */
    
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
    // using queue and while loop because recursion with a depth of nPoints is way too much overhead
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
        memResize(double, fluidPointsHd->pointCrds, fluidPointsHd->nPoints, fluidPointsHd->nPoints+1, fluidPointsHd->maxPoints, 2);
        fluidPointsHd->pointCrds[2*fluidPointsHd->nPoints+0] = x;
        fluidPointsHd->pointCrds[2*fluidPointsHd->nPoints+1] = y;
        fluidPointsHd->nPoints++;
        
        // check and add the four adjacent points to the queue
        if (!mpsCrossesFluidBoundary(fluidBoundaries, nFluidBoundaries, wallSegments,
                                     nWallSegments, x, y, x+wallSpacing, y) &&
            !visited[(xInt+1)+yInt*xDim]) {
            quePush(queHd, x+wallSpacing);
            quePush(queHd, y);
            visited[(xInt+1)+yInt*xDim] = true;
        }
        if (!mpsCrossesFluidBoundary(fluidBoundaries, nFluidBoundaries, wallSegments,
                                     nWallSegments, x, y, x, y+wallSpacing) &&
            !visited[xInt+(yInt+1)*xDim]) {
            quePush(queHd, x);
            quePush(queHd, y+wallSpacing);
            visited[xInt+(yInt+1)*xDim] = true;
        }
        if (!mpsCrossesFluidBoundary(fluidBoundaries, nFluidBoundaries, wallSegments,
                                     nWallSegments, x, y, x-wallSpacing, y) &&
            !visited[(xInt-1)+(yInt*xDim)]) {
            quePush(queHd, x-wallSpacing);
            quePush(queHd, y);
            visited[(xInt-1)+(yInt*xDim)] = true;
        }
        if (!mpsCrossesFluidBoundary(fluidBoundaries, nFluidBoundaries, wallSegments,
                                     nWallSegments, x, y, x, y-wallSpacing) &&
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
    double       r;             /* the radius of influence of points */
    double       dt;            /* the time difference between every calculation */
    double       density;       /* the density of the fluid */
    double       wallSpacing;   /* the rounded spacing between wall points */
    double       viscosity;     /* the viscosity */
    double       beta;          /* beta constant for free surface point calculation */
    double       nTimeSteps;    /* the number of timesteps that the simulation runs for */
    double      *wallSegments;  /* an array of wall segments */
    double       tmp;           /* a temporary real */
    double       tmp1;          /* a temporary real */
    double       tmp2;          /* a temporary real */
    double       tmpArray[4];   /* a temporary array */
    double       dx1;           /* change in x */
    double       dx2;           /* change in x */
    double       dy1;           /* change in y */
    double       dy2;           /* change in y */
    double       x0;            /* the initial x */
    double       y0;            /* the initial y */
    int          nSegs;         /* number of local wall intervals between points */
    int          nWallSegments; /* the number of glocal wall segmnets */
    int          nPoints;       /* the number of global corners */
    double      *fluidBoundaries;   /* an array of the boundaries enclosing the fluid */
    int          nFluidBoundaries;  /* the number of fluid boundaries */
    MpsPointsHd  cornersHd;     /* corners sturcture */
    MpsPointsHd  wallPointsHd;  /* wallPoints structure */
    MpsPointsHd  ghostPointsHd; /* ghostPoints structure */
    MpsPointsHd  fluidPointsHd; /* fluidPoints structure */
    FILE        *fin;           /* input file */

    fin  = fopen("mps.in", "r");

    /*=======================================================================================
     * Input and create wall points and ghost points
     *=======================================================================================
     */
    
    // input radius, wallSpacing, wallSegments
    fscanf(fin, "%lf", &r);
    fscanf(fin, "%lf", &wallSpacing);
    fscanf(fin, "%lf", &dt);
    fscanf(fin, "%lf", &density);
    fscanf(fin, "%lf", &viscosity);
    fscanf(fin, "%lf", &beta);
    fscanf(fin, "%lf", &nTimeSteps);
    fscanf(fin, "%d",  &nWallSegments);

    /*=======================================================================================
     * initialize wall and ghost points
     *=======================================================================================
     */
    wallSegments  = memNew(double, nWallSegments * 4); 
    cornersHd     = mpsNewPoints();
    wallPointsHd  = mpsNewPoints();
    ghostPointsHd = mpsNewPoints();


    for (int i = 0; i < nWallSegments; i++) {
        fscanf(fin, "%le %le %le %le", 
               &wallSegments[4*i+0], &wallSegments[4*i+1], 
               &wallSegments[4*i+2], &wallSegments[4*i+3]); 
    }

    mpsOutCrd("wall_segments.dat", wallSegments, nWallSegments*2, 2);

    // create the corners from wall segments list
    for (int i = 0; i < nWallSegments; i++) {
        mpsGetPointId(cornersHd, wallSegments[4*i+0], wallSegments[4*i+1]);
        mpsGetPointId(cornersHd, wallSegments[4*i+2], wallSegments[4*i+3]);
    }

    // remove later
    mpsOutCrd("corners.dat", cornersHd->pointCrds, cornersHd->nPoints, 2);
    
    nPoints = mpsGetNPoints(cornersHd);
    
    // initialize the boundary corners of the ghost points
    mpsGhostCorners(cornersHd, wallSegments, nWallSegments, ghostPointsHd,
                    nPoints, r);
    
    // remove later
    mpsOutCrd("ghost_corners.dat", ghostPointsHd->pointCrds, nPoints, 2);
    
    // initialize the wall points with all corners
    memResize(double, wallPointsHd->pointCrds, 0, cornersHd->nPoints+1, wallPointsHd->maxPoints, 2);
    
    for (int i = 0; i < cornersHd->nPoints; i++) {
        wallPointsHd->pointCrds[2*i+0] = cornersHd->pointCrds[2*i+0];
        wallPointsHd->pointCrds[2*i+1] = cornersHd->pointCrds[2*i+1];
        wallPointsHd->nPoints++;
    }
    
    // add the intermediate points between the wall corners
    mpsConstructIntermediatePoints(wallSegments, nWallSegments*2, // accepts number of points, not segments
                                   -1, // no reference points
                                   wallPointsHd, wallSpacing, false, false);
    
    // fill ghostPointsHd with the ghost points
    tmp = ghostPointsHd->nPoints-1; // because it is changing in the for loop

    for (int i = 0; i < (int) tmp; i++) {
        // check if the ghost segment is an actual ghost segment or just a segment with two points from different segments
        if (!mpsContainsLine(wallSegments, nWallSegments, cornersHd->pointCrds[2*i+0], cornersHd->pointCrds[2*i+1],
                             cornersHd->pointCrds[2*(i+1)+0], cornersHd->pointCrds[2*(i+1)+1]))
            continue;
        
        tmp1 = mpsDist(cornersHd->pointCrds[2*i+0], cornersHd->pointCrds[2*i+1],
                       ghostPointsHd->pointCrds[2*i+0], ghostPointsHd->pointCrds[2*i+1]);
        tmp2 = mpsDist(cornersHd->pointCrds[2*(i+1)+0], cornersHd->pointCrds[2*(i+1)+1],
                       ghostPointsHd->pointCrds[2*(i+1)+0], ghostPointsHd->pointCrds[2*(i+1)+1]);
        nSegs = ceil(r / wallSpacing);
        dx1 = (ghostPointsHd->pointCrds[2*i+0] - cornersHd->pointCrds[2*i+0]) / nSegs;
        dy1 = (ghostPointsHd->pointCrds[2*i+1] - cornersHd->pointCrds[2*i+1]) / nSegs;
        dx2 = (ghostPointsHd->pointCrds[2*(i+1)+0] - cornersHd->pointCrds[2*(i+1)+0]) / nSegs;
        dy2 = (ghostPointsHd->pointCrds[2*(i+1)+1] - cornersHd->pointCrds[2*(i+1)+1]) / nSegs;
        
        // for each two corner and ghostPoint pairs, there are nSegs lines between them which comprise the ghost points
        for (int j = 1; j < nSegs+1; j++) {
            tmpArray[0] = cornersHd->pointCrds[2*i+0] + j*dx1;
            tmpArray[1] = cornersHd->pointCrds[2*i+1] + j*dy1;
            tmpArray[2] = cornersHd->pointCrds[2*(i+1)+0] + j*dx2; // extend by dx of pointCrds
            tmpArray[3] = cornersHd->pointCrds[2*(i+1)+1] + j*dy2;
            
            mpsConstructIntermediatePoints(tmpArray, 2,
                                           mpsDist(cornersHd->pointCrds[2*i+0], cornersHd->pointCrds[2*i+1],
                                                   cornersHd->pointCrds[2*(i+1)+0], cornersHd->pointCrds[2*(i+1)+1]) / wallSpacing,
                                           ghostPointsHd, wallSpacing, true, true);
        }
    }
    
    mpsOutCrd("ghost_points.dat", ghostPointsHd->pointCrds, ghostPointsHd->nPoints, 2);
    
    // remove later
    mpsOutCrd("wall_points.dat", wallPointsHd->pointCrds, wallPointsHd->nPoints, 2);
    
    /*=======================================================================================
     * Input and create fluid points
     *=======================================================================================
     */
    
    // input fluid boundaries
    fscanf(fin, "%d", &nFluidBoundaries);
    
    // initialize fluidBoundaries and fluidPoints
    fluidBoundaries = memNew(double, nFluidBoundaries*4);
    fluidPointsHd   = mpsNewPoints();
    
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

    fscanf(fin, "%lf %lf", &x0, &y0);

    mpsFloodfill(fluidPointsHd, fluidBoundaries,
                 nFluidBoundaries, wallSegments, nWallSegments, 
                 x0, y0, wallSpacing);

    mpsOutCrd("fluid_points.dat", fluidPointsHd->pointCrds, fluidPointsHd->nPoints, 2);
    
    /*=======================================================================================
     * Solve the problem
     *=======================================================================================
     */

     mpsDriver(fluidPointsHd, wallPointsHd, ghostPointsHd, 
               r, wallSpacing, beta, density, viscosity, nTimeSteps, dt);

    /*=======================================================================================
     * End
     *=======================================================================================
     */

     return 0;

} 

/*******************************************************************************
 * LIM: limit update
 *******************************************************************************
 */
#define        LIM(X,T)        MIN(MAX((X),-(T)),(T))

/*******************************************************************************
 * "mpsDriver": MPS time stepping driver
 *******************************************************************************
 */
int mpsDriver(MpsPointsHd fluidPointsHd, MpsPointsHd wallPointsHd, MpsPointsHd ghostPointsHd,
              double      radius,        double      wallSpacing,  double      beta, 
              double      density,       double      viscosity,    int         nTimeSteps, 
              double      dt) {
    StnHd        stnHd;         /* an adjacency structure */
    char         buffer[1024];  /* a string buffer */
    double       dNum0;         /* wall spacing based density number */
    double       dr;            /* distance */
    double       dx;            /* change in x */
    double       dy;            /* change in y */
    double       l2Norm;        /* L2 norm of vector */
    double      *presCurr;      /* the pressure for the current time step */
    double      *presNext;      /* the pressure for the next time step */
    double      *xPosCurr;      /* the x position for the current time step */
    double      *xPosNext;      /* the x position for the next time step */
    double      *xPosStar;      /* the x position for the next time step */
    double      *xVelCurr;      /* the x velocity for the current time step */
    double      *xVelNext;      /* the x velocity for the next time step */
    double      *xVelStar;      /* the x velocity correction based on the pressure */
    double      *yPosCurr;      /* the y position for the current time step */
    double      *yPosNext;      /* the y position for the next time step */
    double      *yPosStar;      /* the y position for the next time step */
    double      *yVelCurr;      /* the y velocity for the current time step */
    double      *yVelNext;      /* the y velocity for the next time step */
    double      *yVelStar;      /* the y velocity correction based on the pressure */
    double       limX;          /* limit the update */
    int          n;             /* a count */
    int          nPoints;       /* total number of points */
    int          offSet;        /* point offset */
    int          nFluidPoints;  /* the total number of fluid points */
    int          nWallPoints;   /* the total number of wallPoints */

/*---------------------------------------------------------------------------------------
 * Get the dimensions
 *---------------------------------------------------------------------------------------
 */
    nFluidPoints = fluidPointsHd->nPoints;
    nWallPoints  = wallPointsHd->nPoints + ghostPointsHd->nPoints;
    nPoints      = nFluidPoints + nWallPoints;

    limX = 0.5 * wallSpacing;

/*---------------------------------------------------------------------------------------
 * Get a measure of density number
 *---------------------------------------------------------------------------------------
 */
    dNum0 = 0;
    n = (int) (radius / wallSpacing);
    for (int i = -n; i <= n; i++) {
        for (int j = -n; j <= n; j++) {
            dx     = i * wallSpacing;
            dy     = j * wallSpacing;
            dr     = sqrt(dx*dx + dy*dy);
            dNum0 += stnWeight(dr, radius);
        }
    }
    printf("Wall spacing num dens. = %g\n", dNum0);

/*---------------------------------------------------------------------------------------
 * Allocate arrays
 *---------------------------------------------------------------------------------------
 */
    xVelCurr = memNewZero(double, nPoints);
    yVelCurr = memNewZero(double, nPoints);
    xVelNext = memNewZero(double, nPoints);
    yVelNext = memNewZero(double, nPoints);
    xVelStar = memNewZero(double, nPoints);
    yVelStar = memNewZero(double, nPoints);

    xPosCurr = memNewZero(double, nPoints);
    yPosCurr = memNewZero(double, nPoints);
    xPosNext = memNewZero(double, nPoints);
    yPosNext = memNewZero(double, nPoints);
    xPosStar = memNewZero(double, nPoints);
    yPosStar = memNewZero(double, nPoints);

    presCurr = memNewZero(double, nPoints);
    presNext = memNewZero(double, nPoints);

/*---------------------------------------------------------------------------------------
 * Initialize the position
 *---------------------------------------------------------------------------------------
 */
    offSet = 0;
    for (int i = 0; i < fluidPointsHd->nPoints; i++) {
        xPosCurr[offSet+i] = fluidPointsHd->pointCrds[2*i+0];
        yPosCurr[offSet+i] = fluidPointsHd->pointCrds[2*i+1];
    }
    offSet += fluidPointsHd->nPoints;
    for (int i = 0; i < wallPointsHd->nPoints; i++) {
        xPosCurr[offSet+i] = wallPointsHd->pointCrds[2*i+0];
        yPosCurr[offSet+i] = wallPointsHd->pointCrds[2*i+1];
    }
    offSet += wallPointsHd->nPoints;
    for (int i = 0; i < ghostPointsHd->nPoints; i++) {
        xPosCurr[offSet+i] = ghostPointsHd->pointCrds[2*i+0];
        yPosCurr[offSet+i] = ghostPointsHd->pointCrds[2*i+1];
    }
    offSet += ghostPointsHd->nPoints;

    snprintf(buffer, sizeof(buffer), "mps.%d.out", 0);
    mpsOutCrdXY(buffer, xPosCurr, yPosCurr, nFluidPoints);

    snprintf(buffer, sizeof(buffer), "mps.vel.%d.out", 0);
    mpsOutCrdXY(buffer, xVelCurr, yVelCurr, nFluidPoints);
    
/*---------------------------------------------------------------------------------------
 * Initialize the values
 *---------------------------------------------------------------------------------------
 */
    for (int i = 0; i < nPoints; i++) {
        xVelCurr[i] = 0;
        yVelCurr[i] = 0;
        presCurr[i] = 0;
    }
/*---------------------------------------------------------------------------------------
 * Create and initialize the search data structure
 *---------------------------------------------------------------------------------------
 */
    stnHd = stnNew(nFluidPoints, nWallPoints, radius, beta);

    stnPopulate(stnHd, xPosCurr, yPosCurr);
    stnHd->n0 = (stnHd->n0+dNum0)/2;
    printf("Initial number density = %g\n", stnHd->n0);

/*---------------------------------------------------------------------------------------
 * Smooth the initial condition
 *---------------------------------------------------------------------------------------
 */
//mpsOutCrdXY("x0.dat", xPosCurr, yPosCurr, nPoints);
//    slvSmoothInit(stnHd, xPosCurr, yPosCurr, xPosStar, yPosStar);
//mpsOutCrdXY("x1.dat", xPosCurr, yPosCurr, nPoints);

/*---------------------------------------------------------------------------------------
 * Loop over the time steps
 *---------------------------------------------------------------------------------------
 */
    for (int stepId = 0; stepId < nTimeSteps; stepId++) {

        printf("===============> Processing time step %d ; time increment %g\n", 
               stepId+1, dt);

/*---------------------------------------------------------------------------------------
 * Find the neighboring points
 *---------------------------------------------------------------------------------------
 */
        stnPopulate(stnHd, xPosCurr, yPosCurr);

        if (stepId == 0) mpsOutCrd("density.dat", stnHd->dNum, nFluidPoints, 1);

        snprintf(buffer, sizeof(buffer), "mps.dens.%d.out", stepId+1);
        if(stepId%10==9)mpsOutCrd(buffer, stnHd->dNum, nFluidPoints, 1);

/*---------------------------------------------------------------------------------------
 * Advance the explicit part
 *---------------------------------------------------------------------------------------
 */
        slvCalcExplicitVelocity(stnHd, xVelCurr, xVelStar, viscosity, dt, 0);
        slvCalcExplicitVelocity(stnHd, yVelCurr, yVelStar, viscosity, dt, -9.8);

        l2Norm = mps2VecL2(xVelStar, yVelStar, nFluidPoints);
        printf("Explicit vel. L2 norm  = %g\n", l2Norm);

	if ( l2Norm > 100. ) exit(1);

        for (int i = 0; i < nPoints; i++) {
            xPosStar[i] = xPosCurr[i] + LIM(dt*xVelStar[i], limX);
            yPosStar[i] = yPosCurr[i] + LIM(dt*yVelStar[i], limX);
        }
/*---------------------------------------------------------------------------------------
 * Compute the pressure
 *---------------------------------------------------------------------------------------
 */
        //stnRecalc(stnHd, xPosStar, yPosStar);

        slvCalcPressure(stnHd,    xPosStar, yPosStar,
                        xVelStar, yVelStar, presNext,
                        dt,       density);

        snprintf(buffer, sizeof(buffer), "mps.pres.%d.out", stepId+1);
        if(stepId%10==9)mpsOutCrd(buffer, presNext, nFluidPoints, 1);

/*---------------------------------------------------------------------------------------
 * Correct the velocity
 *---------------------------------------------------------------------------------------
 */
        slvCalcCorrection(stnHd,       presNext, xVelNext,
                          xPosCurr,    density,  dt);
        slvCalcCorrection(stnHd,       presNext, yVelNext,
                          yPosCurr,    density,  dt);
        
/*---------------------------------------------------------------------------------------
 * Update the values
 *---------------------------------------------------------------------------------------
 */
        for (int i = 0; i < nPoints; i++) {
            xVelNext[i] = xVelNext[i] + xVelStar[i];
            yVelNext[i] = yVelNext[i] + yVelStar[i];
            
            xPosNext[i] = xPosCurr[i] + LIM(dt*xVelNext[i],limX);
            yPosNext[i] = yPosCurr[i] + LIM(dt*yVelNext[i],limX);
        }

        l2Norm = mps2VecL2(xVelNext, yVelNext, nFluidPoints);
        printf( "Implicit vel. L2 norm  = %g\n", l2Norm );

	if ( l2Norm > 100. ) exit(1);

/*---------------------------------------------------------------------------------------
 * Output the data
 *---------------------------------------------------------------------------------------
 */
        snprintf(buffer, sizeof(buffer), "mps.%d.out", stepId+1);
        if(stepId%10==9)mpsOutCrdXY(buffer, xPosNext, yPosNext, nFluidPoints);

        snprintf(buffer, sizeof(buffer), "mps.vel.%d.out", stepId+1);
        if(stepId%10==9)mpsOutCrdXY(buffer, xVelCurr, yVelCurr, nFluidPoints);

/*---------------------------------------------------------------------------------------
 * Reset the vectors
 *---------------------------------------------------------------------------------------
 */
        for (int i = 0; i < nPoints; i++) {
            xPosCurr[i] = xPosNext[i];
            yPosCurr[i] = yPosNext[i];
            xVelCurr[i] = xVelNext[i];
            yVelCurr[i] = yVelNext[i];
            presCurr[i] = presNext[i];
        }
    }
        
    return 0;

}


