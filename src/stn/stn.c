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
StnHd stnNew(int nPoints) {
    StnHd stnHd;

    stnHd	       = memNew(Stn, 1);
    stnHd->nPoints     = nPoints;

    stnHd->maxAdjacent = 100;
    stnHd->row	       = memNew(int, stnHd->maxAdjacent);
    stnHd->col	       = memNew(int, nPoints+1);
    stnHd->weights     = memNew(double, stnHd->maxAdjacent);
    stnHd->dist        = memNew(double, stnHd->maxAdjacent);
    stnHd->dNum        = memNew(double, nPoints);
    for (int i = 0; i < nPoints+1; i++) stnHd->dNum[i] = 0;
    stnHd->diagIndex   = memNew(int, nPoints);

    return stnHd;
}

/********************************************************************************
 * "weight": calculates the kernel function
 ********************************************************************************/
double weight(double dist, double radius) {
    double tol;			/* floating-point error */

    tol = 1.e-12;

    if (dist < tol || dist > radius+tol) return 0;
    if (dist > tol && dist < radius+tol) return radius/dist - 1;
    else return 0;
}

void stnPopulate(StnHd stnHd, double *points, int nFluidPoints, int nWallPoints, double radius) {
    double	dx;		/* change in x */
    double	dy;		/* change in y */
    double	dist;		/* distance between points */
    int		nAdjacent;	/* a count of the number of points within the radius for each point */
    int		count;          /* the count of adjacent points so far */
    double      tol;            /* the floating-points error */
    double      totDNum;        /* the total density number */

    tol = 1.e-12;
    totDNum = 0;
    count = 0;

    // populate the row and col adjacency list
    for (int i = 0; i < stnHd->nPoints; i++) {
	stnHd->col[i] = count;
	for (int j = 0; j < stnHd->nPoints; j++) {
	    dx = points[2*i+0] - points[2*j+0];
	    dy = points[2*i+1] - points[2*j+1];
	    dist = sqrt(dx*dx + dy*dy);
	    if (dist < radius+tol) {
		memResize(int, stnHd->row, count, count+1, stnHd->maxAdjacent, 1);
		memResize(double, stnHd->weights, count, count+1, stnHd->maxAdjacent, 1);
		memResize(double, stnHd->dist, count, count+1, stnHd->maxAdjacent, 1);
		stnHd->row[count] = j;
		stnHd->weights[count] = weight(dist, radius);
		stnHd->dist[count] = dist;
		if (i < nFluidPoints) // only fluid points
		    stnHd->dNum[i] += stnHd->weights[count];
		count++;
	    }
	    if (i == j)
		stnHd->diagIndex[i] = count;
	}
	if (i < nFluidPoints)
	    totDNum += stnHd->dNum[i];
    }
    stnHd->col[stnHd->nPoints] = count;
    stnHd->n0 = totDNum/nFluidPoints;
}
