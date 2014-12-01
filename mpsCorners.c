/*******************************************************************************
** Copyright 2014-2014 Vedaad Shakib Inc.
*******************************************************************************/

/*******************************************************************************
** 
** "mpsCorners.c": Corner manipulation functions
**
*******************************************************************************/

#include "sys.h"
#include "mps.h"

/*******************************************************************************
 * "mpsNewCorners": return a new corners data structure
 *******************************************************************************
 */
MpsCornersHd mpsNewCorners() {
  MpsCornersHd	cornersHd; // corners structure

  /* allocate structure */
  cornersHd		= memNew(MpsCorners, 1);
  
  cornersHd->maxCorners	= 100;
  cornersHd->nCorners	= 0;
  cornersHd->cornerCrds	= memNew(double, 2*cornersHd->maxCorners);
  
  /* return structure */
  return (cornersHd);
} 

/*******************************************************************************
 * "mpsFreeCorners": free memory of the corners
 *******************************************************************************
 */
void mpsFreeCorners(MpsCornersHd cornersHd) {
  /* free memory */
  free(cornersHd->cornerCrds);
  free(cornersHd);
} 

/*******************************************************************************
 * "mpsGetCornerId": gets the ID a corner point and creates a new one if not found
 *******************************************************************************
 */
int mpsGetCornerId(MpsCornersHd	cornersHd, double x, double y) {
    int		 i;		/* index		    */
    double	*cornerCrds;	/* localized corner points  */
    int		 nCorners;	/* number of corners        */

    // localize data
    nCorners	= cornersHd->nCorners;
    cornerCrds	= cornersHd->cornerCrds;

    // search for corner in existing corners

    for (i = 0; i < nCorners; i++) {
	if (cornerCrds[2*i+0] == x && cornerCrds[2*i+1] == y) {
	    return (i);
	}
    }

    // if not found, create a new corner point
    memResize(double,			cornersHd->cornerCrds,	
	      cornersHd->nCorners,	cornersHd->nCorners+1,	
	      cornersHd->maxCorners,	2		      );

    cornersHd->cornerCrds[2*cornersHd->nCorners+0] = x;
    cornersHd->cornerCrds[2*cornersHd->nCorners+1] = y;

    cornersHd->nCorners++;

    return (cornersHd->nCorners-1);

} 

