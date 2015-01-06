/*******************************************************************************
** Copyright 2014-2014 Vedaad Shakib Inc.
*******************************************************************************/

/*******************************************************************************
** 
** "mpsWallPoints.c": Wall point manipulation functions
**
*******************************************************************************/

#include "sys.h"
#include "mps.h"

/*******************************************************************************
 * "mpsNewWallPoints": return a wall points data structure
 *******************************************************************************
 */
MpsWallPointsHd mpsNewWallPoints() {
    MpsWallPointsHd wallPointsHd; // wallPoints structure
    
    /* allocate structure */
    wallPointsHd = memNew(MpsWallPoints, 1);
    
    wallPointsHd->maxWallPoints = 100;
    wallPointsHd->nWallPoints   = 0;
    wallPointsHd->wallPointCrds = memNew(double, 2*wallPointsHd->maxWallPoints);
    
    /* return structure */
    return(wallPointsHd);
} 

/*******************************************************************************
 * "mpsFreeWallPoints": free memory of the wallPoints
 *******************************************************************************
 */
void mpsFreeWallPoints(MpsWallPointsHd wallPointsHd) {
    /* free memory */
    free(wallPointsHd->wallPointCrds);
    free(wallPointsHd);
} 

/*******************************************************************************
 * "mpsGetWallPointId": get index of a wallPoint point; create one if not found
 *******************************************************************************
 */
int mpsGetWallPointId(MpsWallPointsHd wallPointsHd, double x, double y) {
    double      *wallPointCrds; /* localized wall points  */
    int          nWallPoints;   /* number of wall points  */
    
    // localize data
    nWallPoints   = wallPointsHd->nWallPoints;
    wallPointCrds = wallPointsHd->wallPointCrds;
    
    // search for wall point in existing wall points
    for (int i = 0; i < nWallPoints; i++) {
        if (wallPointCrds[2*i+0] == x && wallPointCrds[2*i+1] == y) {
            return (i);
        }
    }
    
    // if not found, create a new wall point
    memResize(double,                        wallPointsHd->wallPointCrds,        
              wallPointsHd->nWallPoints,     wallPointsHd->nWallPoints+1,        
              wallPointsHd->maxWallPoints,   2                             );
    
    wallPointsHd->wallPointCrds[2*wallPointsHd->nWallPoints+0] = x;
    wallPointsHd->wallPointCrds[2*wallPointsHd->nWallPoints+1] = y;

    wallPointsHd->nWallPoints++;

    return (wallPointsHd->nWallPoints-1);
} 

