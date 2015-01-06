/*******************************************************************************
** Copyright 2014-2014 Vedaad Shakib Inc.
*******************************************************************************/

/*******************************************************************************
** 
** "mpsFluidPoints.c": FluidPoint manipulation functions
**
*******************************************************************************/

#include "sys.h"
#include "mps.h"

/*******************************************************************************
 * "mpsNewFluidPoints": return a new fluidPoints data structure
 *******************************************************************************
 */
MpsFluidPointsHd mpsNewFluidPoints() {
    MpsFluidPointsHd fluidPointsHd; // fluidPoints structure
    
    /* allocate structure */
    fluidPointsHd = memNew(MpsFluidPoints, 1);
    
    fluidPointsHd->maxFluidPoints = 100;
    fluidPointsHd->nFluidPoints   = 0;
    fluidPointsHd->fluidPointCrds = memNew(double, 2*fluidPointsHd->maxFluidPoints);
    
    /* return structure */
    return(fluidPointsHd);
} 

/*******************************************************************************
 * "mpsFreeFluidPoints": free memory of the fluidPoints
 *******************************************************************************
 */
void mpsFreeFluidPoints(MpsFluidPointsHd fluidPointsHd) {
    /* free memory */
    free(fluidPointsHd->fluidPointCrds);
    free(fluidPointsHd);
} 

/*******************************************************************************
 * "mpsGetFluidPointId": gets the ID a fluidPoint point and creates a new one if not found
 *******************************************************************************
 */
int mpsGetFluidPointId(MpsFluidPointsHd fluidPointsHd, double x, double y) {
    double      *fluidPointCrds;	/* localized fluidPoint points  */
    int          nFluidPoints;          /* number of fluidPoints        */
    
    // localize data
    nFluidPoints   = fluidPointsHd->nFluidPoints;
    fluidPointCrds = fluidPointsHd->fluidPointCrds;
    
    // search for fluidPoint in existing fluidPoints

    for (int i = 0; i < nFluidPoints; i++) {
        if (fluidPointCrds[2*i+0] == x && fluidPointCrds[2*i+1] == y) {
            return (i);
        }
    }
    
    // if not found, create a new fluidPoint point
    memResize(double,                        fluidPointsHd->fluidPointCrds,        
              fluidPointsHd->nFluidPoints,   fluidPointsHd->nFluidPoints+1,        
              fluidPointsHd->maxFluidPoints, 2                      );
    
    fluidPointsHd->fluidPointCrds[2*fluidPointsHd->nFluidPoints+0] = x;
    fluidPointsHd->fluidPointCrds[2*fluidPointsHd->nFluidPoints+1] = y;
    
    fluidPointsHd->nFluidPoints++;
    
    return (fluidPointsHd->nFluidPoints-1);
} 

