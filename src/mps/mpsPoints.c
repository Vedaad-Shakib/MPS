/*******************************************************************************
** Copyright 2014-2014 Vedaad Shakib Inc.
*******************************************************************************/

/*******************************************************************************
** 
** "mpsPoints.c": Point storage and manipulation functions
**
*******************************************************************************/

#include "sys.h"
#include "mps.h"

/*******************************************************************************
 * "mpsNewPoints": return a new points data structure
 *******************************************************************************
 */
MpsPointsHd mpsNewPoints() {
    MpsPointsHd pointsHd; // points structure
    
    /* allocate structure */
    pointsHd = memNew(MpsPoints, 1);
    
    pointsHd->maxPoints = 100;
    pointsHd->nPoints   = 0;
    pointsHd->pointCrds = memNew(double, 2*pointsHd->maxPoints);
    
    /* return structure */
    return(pointsHd);
} 

/*******************************************************************************
 * "mpsFreePoints": free memory of the points
 *******************************************************************************
 */
void mpsFreePoints(MpsPointsHd pointsHd) {
    /* free memory */
    free(pointsHd->pointCrds);
    free(pointsHd);
} 

/*******************************************************************************
 * "mpsGetPointId": get index of a point; create one if not found
 *******************************************************************************
 */
int mpsGetPointId(MpsPointsHd pointsHd, double x, double y) {
    double      *pointCrds;	/* localized points  */
    int          nPoints;	/* number of points  */
    
    // localize data
    nPoints   = pointsHd->nPoints;
    pointCrds = pointsHd->pointCrds;
    
    // search for point in existing points
    for (int i = 0; i < nPoints; i++) {
        if (pointCrds[2*i+0] == x && pointCrds[2*i+1] == y) {
            return (i);
        }
    }
    
    // if not found, create a new point
    memResize(double,                pointsHd->pointCrds,        
              pointsHd->nPoints,     pointsHd->nPoints+1,        
              pointsHd->maxPoints,   2                             );
    
    pointsHd->pointCrds[2*pointsHd->nPoints+0] = x;
    pointsHd->pointCrds[2*pointsHd->nPoints+1] = y;

    pointsHd->nPoints++;

    return (pointsHd->nPoints-1);
} 

