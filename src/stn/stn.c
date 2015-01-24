/*******************************************************************************
 ** Copyright 2014-2014 Vedaad Shakib Inc.
 *******************************************************************************/

/*******************************************************************************
 ** 
 ** "stn.c": A module for finding and storing points in a certain radius of other
 **          points
 **
 *******************************************************************************/

/*******************************************************************************
 * Standard includes
 *******************************************************************************/

#include "sys.h"
#include "stn.h"

/*******************************************************************************
 * "stnNew": creates a new stn structure
 *******************************************************************************/
StnHd stnNew() {
    StnHd stnHd;

    stnHd = memNew(Stn, 1);

    stnHd->maxPoints   = 100;
    stnHd->maxAdjacent = 100;
    stnHd->col	       = memNew(int, stnHd->maxPoints);
    stnHd->nPoints     = 0;
    stnHd->row	       = memNew(double, stnHd->maxAdjacent);

    return stnHd;
}

void stnPopulate(StnHd stnHd, double *points, int nPoints, double radius) {
    double	dx;		/* change in x */
    double	dy;		/* change in y */
    double	dist;		/* distance between points */
    int		nAdjacent;	/* a count of the number of points within the radius for each point */

    for (int i = 0; i < nPoints; i++) {
	nAdjacent = 0;
	for (int j = 0; j < nPoints; j++) {
	    if (i == j) continue;
	    dx = points[2*i+0] - points[2*j+0];
	    dy = points[2*i+1] - points[2*j+1];
	    dist = sqrt(dx*dx + dy*dy);
	    
	    if (dist <= radius) {
		if (stnHd->nPoints == 0) {
		    memResize(double, stnHd->row, nAdjacent, nAdjacent+1, stnHd->maxAdjacent, 2);
		    stnHd->row[2*nAdjacent+0] = points[2*j+0];
		    stnHd->row[2*nAdjacent+1] = points[2*j+1];
		} else {
		    memResize(double, stnHd->row, stnHd->col[stnHd->nPoints-1]+nAdjacent, stnHd->col[stnHd->nPoints-1]+nAdjacent+1, stnHd->maxAdjacent, 2);
		    stnHd->row[2*(stnHd->col[stnHd->nPoints-1]+nAdjacent)+0] = points[2*j+0];
		    stnHd->row[2*(stnHd->col[stnHd->nPoints-1]+nAdjacent)+1] = points[2*j+1];
		}
		nAdjacent++;
	    }
	}
	memResize(int, stnHd->col, nPoints, nPoints+1, stnHd->maxPoints, 2);
	if (stnHd->nPoints == 0) stnHd->col[stnHd->nPoints] = nAdjacent;
	else stnHd->col[stnHd->nPoints] = nAdjacent + stnHd->col[stnHd->nPoints-1];

	stnHd->nPoints++;
    }
}
