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
StnHd stnNew(int nFluidPoints, int nWallPoints) {
    StnHd stnHd;

    stnHd		= memNew(Stn, 1);
    stnHd->nFluidPoints = nFluidPoints;
    stnHd->nWallPoints	= nWallPoints;
    stnHd->nPoints	= nFluidPoints + nWallPoints;

    stnHd->maxAdjacent = 100;
    stnHd->row	       = memNew(int, stnHd->maxAdjacent);
    stnHd->col	       = memNew(int, stnHd->nPoints+1);
    stnHd->weights     = memNew(double, stnHd->maxAdjacent);
    stnHd->dist        = memNew(double, stnHd->maxAdjacent);
    stnHd->dNum        = memNew(double, stnHd->nPoints);
    stnHd->diagIndex   = memNew(int, stnHd->nPoints);

    stnHd->d = 2;

    for (int i = 0; i < stnHd->nPoints+1; i++) stnHd->dNum[i] = 0;

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

void stnPopulate(StnHd  stnHd,  double *xCrd, double *yCrd,
		 double radius, double  beta) {
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

    stnHd->radius = radius;
    stnHd->beta   = beta;

    // populate the row and col adjacency list
    for (int i = 0; i < stnHd->nPoints; i++) {
	stnHd->col[i] = count;
	for (int j = 0; j < stnHd->nPoints; j++) {
	    dx = xCrd[i+0] - xCrd[j+0];
	    dy = xCrd[i+1] - xCrd[j+1];
	    dist = sqrt(dx*dx + dy*dy);
	    if (dist < radius+tol) {
		memResize(int, stnHd->row, count, count+1, stnHd->maxAdjacent, 1);
		memResize(double, stnHd->weights, count, count+1, stnHd->maxAdjacent, 1);
		memResize(double, stnHd->dist, count, count+1, stnHd->maxAdjacent, 1);
		stnHd->row[count] = j;
		stnHd->weights[count] = weight(dist, radius);
		stnHd->dist[count] = dist;
		if (i < stnHd->nFluidPoints) // only fluid points
		    stnHd->dNum[i] += stnHd->weights[count];
		count++;
	    }
	    if (i == j)
		stnHd->diagIndex[i] = count;
	}
	if (i < stnHd->nFluidPoints)
	    totDNum += stnHd->dNum[i];
    }
    stnHd->col[stnHd->nPoints] = count;
    stnHd->n0 = totDNum/stnHd->nFluidPoints;
}

/*******************************************************************************
 * "stnRecalc": recalculates the weight and dist (assumes filled row and col)
 *******************************************************************************
 */
void stnRecalc(StnHd stnHd, double *xCrd, double *yCrd) {
    double	dx;		/* the x distance between two points */
    double	dy;		/* the y distance between two points */
    double	dist;		/* the distance between two points */

    // reset density num array
    for (int i = 0; i < stnHd->nFluidPoints; i++) stnHd->dNum[i] = 0;

    // recalculate the weight and dist using same col and row
    for (int i = 0; i < stnHd->nFluidPoints; i++) {
	for (int k = stnHd->col[i]; k < stnHd->col[i+1]; k++) {
	    int j = stnHd->row[k];
	    if (i == j) continue;

	    dx	 = xCrd[i] - xCrd[j];
	    dy	 = yCrd[i] - yCrd[j];
	    dist = sqrt(dx*dx + dy*dy);

	    stnHd->weights[k] = weight(dist, stnHd->radius);
	    stnHd->dist[k]    = dist;
	    stnHd->dNum[i]   += stnHd->weights[k];
	}
    }
}

void stnFree(StnHd stnHd) {
    free(stnHd->col);
    free(stnHd->row);
    free(stnHd->weights);
    free(stnHd->dist);
    free(stnHd->dNum);
    free(stnHd->diagIndex);
    
    free(stnHd);
}
