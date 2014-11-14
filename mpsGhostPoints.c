/*******************************************************************************
** Copyright 2014-2014 Vedaad Shakib Inc.
*******************************************************************************/

/*******************************************************************************
** 
** "mpsGhostPoints.c": Wall point manipulation functions
**
*******************************************************************************/

#include "mps.h"

/*******************************************************************************
 * "mpsNewWallPoints": return a wall points data structure
 *******************************************************************************
 */
MpsGhostPointsHd mpsNewGhostPoints() {
  MpsGhostPointsHd ghostPointsHd; // ghostPoints structure

  /* allocate structure */
  ghostPointsHd	= memNew(MpsGhostPoints, 1);
  
  ghostPointsHd->maxGhostPoints	= 100;
  ghostPointsHd->nGhostPoints	= 0;
  ghostPointsHd->ghostPointCrds	= memNew(double, 2*ghostPointsHd->maxGhostPoints);
  
  /* return structure */
  return(ghostPointsHd);
} 

/*******************************************************************************
 * "mpsFreeGhostPoints": free memory of the ghostPoints
 *******************************************************************************
 */
void mpsFreeGhostPoints(MpsGhostPointsHd ghostPointsHd) {
  /* free memory */
  free(ghostPointsHd->ghostPointCrds);
  free(ghostPointsHd);
} 

/*******************************************************************************
 * "mpsGetGhostPointId": get index of a ghostPoint point; create one if not found
 *******************************************************************************
 */
int mpsGetGhostPointId(MpsGhostPointsHd	ghostPointsHd, double x, double y) {
  double	*ghostPointCrds; /* localized ghost points  */
    int		 nGhostPoints;	 /* number of ghost points        */

    // localize data
    nGhostPoints	  = ghostPointsHd->nGhostPoints;
    ghostPointCrds = ghostPointsHd->ghostPointCrds;

    // search for ghost point in existing ghost points
    for (int i = 0; i < nGhostPoints; i++) {
	if (ghostPointCrds[2*i+0] == x && ghostPointCrds[2*i+1] == y) {
	    return (i);
	}
    }

    // if not found, create a new ghost point
    memResize(double,				 ghostPointsHd->ghostPointCrds,	
	      ghostPointsHd->nGhostPoints, 	 ghostPointsHd->nGhostPoints+1,	
	      ghostPointsHd->maxGhostPoints,	 2		             );

    ghostPointsHd->ghostPointCrds[2*ghostPointsHd->nGhostPoints+0] = x;
    ghostPointsHd->ghostPointCrds[2*ghostPointsHd->nGhostPoints+1] = y;

    ghostPointsHd->nGhostPoints++;

    return (ghostPointsHd->nGhostPoints-1);
} 

