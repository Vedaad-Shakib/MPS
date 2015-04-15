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
StnHd stnNew(int nFluidPoints, int nWallPoints,
	     double radius, double  beta) {
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
    stnHd->rhs	       = memNew(double, stnHd->nPoints);
    stnHd->lhs         = memNew(double, stnHd->maxAdjacent);
    stnHd->freeSurf    = memNew(int, stnHd->nPoints);

    stnHd->radius = radius;
    stnHd->beta   = beta;

    stnHd->d  = 2;
    stnHd->n0 = 0;

    for (int i = 0; i < stnHd->nPoints; i++) {
	stnHd->dNum[i] = 0;
	stnHd->freeSurf[i] = 0;
    }

    return stnHd;
}

/********************************************************************************
 * "stnWeight": calculates the kernel function
 ********************************************************************************/
double stnWeight(double dist, double radius) {
    double tol;			/* floating-point error */

    tol = 1.e-12;

    if (dist < tol || dist > radius+tol) return 0;
    if (dist > tol && dist < radius+tol) return radius/dist - 1;
    else return 0;
}

/********************************************************************************
 * "stnPopulate": populate the stn adjacency lists by searching through the points
 ********************************************************************************/
void stnPopulate(StnHd  stnHd,  double *xCrd, double *yCrd) {
    double	beta;		/* free surface beta */
    double	dist;		/* distance between points */
    double	dx;		/* change in x */
    double	dy;		/* change in y */
    double	maxDNum;	/* max density number */
    double	radius;		/* search radius */
    double      tol;            /* the floating-points error */
    int		count;          /* the count of adjacent points so far */
    int		nFluidPoints;	/* no. fluid points */
    int		nPoints;	/* no. point */
    int         tmp;            /* a temporary variable */

    tol			= 1.e-12;
    count		= 0;

    radius		= stnHd->radius;
    beta		= stnHd->beta;
    nPoints		= stnHd->nPoints;
    nFluidPoints	= stnHd->nFluidPoints;

    // populate the row and col adjacency list
    for (int i = 0; i < nPoints; i++) {
	stnHd->col[i] = count;
	stnHd->dNum[i] = 0;
	for (int j = 0; j < nPoints; j++) {
	    dx = xCrd[i] - xCrd[j];
	    dy = yCrd[i] - yCrd[j];
	    dist = sqrt(dx*dx + dy*dy);
	    if (dist < radius+tol) {
		memResize(int, stnHd->row, count, count+1, stnHd->maxAdjacent, 1);

		stnHd->row[count] = j;
	    }
	}
    }
    stnHd->col[nPoints] = count;

    /*printf("No. fluid nodes        = %d\n", nFluidPoints);
    printf("No. nodes              = %d\n", nPoints);
    printf("No. nonzeros           = %d\n", stnHd->col[nPoints]);
    printf("Ave. degree            = %g\n", ((double)stnHd->col[nPoints])/nPoints);
    printf("No. free-surface nodes = %d\n", count);
    printf("Number density         = %g\n", totDNum/nFluidPoints);
    printf("Max number density     = %g\n", maxDNum);*/
}

/*******************************************************************************
 * "stnCalc": Calculates the weight and other statistics assuming that the adjacency list is full
 *******************************************************************************
 */
void stnCalc(StnHd stnHd, double *xCrd, double *yCrd) {
    double	dx;		/* the x distance between two points */
    double	dy;		/* the y distance between two points */
    double	dist;		/* the distance between two points */
    int		nFluidPoints;	/* no. fluid points */
    int		nPoints;	/* no. points */
    double      radius;         /* the radius of influence */
    int         count;          /* a counter */
    int         tmp;            /* a temporary variable */
    double      totDNum;        /* the total density number */
    double	n0Lim;		/* free surface density number */
    double      maxDNum;        /* the maximum density number */

    radius	 = stnHd->radius;
    nPoints	 = stnHd->nPoints;
    nFluidPoints = stnHd->nFluidPoints;
    count	 = 0;
    tmp		 = 0;
    totDNum	 = 0;

    // resize the arrays to the length of the row adjacency list
    memResize(double, stnHd->weights, 0, stnHd->maxAdjacent, tmp, 1);
    tmp = 0;
    memResize(double, stnHd->dist, 0, stnHd->maxAdjacent, tmp, 1);
    
    // reset density num array
    for (int i = 0; i < nPoints; i++) stnHd->dNum[i] = 0;

    // recalculate the weight and dist using same col and row
    for (int i = 0; i < nPoints; i++) {
	for (int k = stnHd->col[i]; k < stnHd->col[i+1]; k++) {
	    int j = stnHd->row[k];
	    if (i == j) {
		stnHd->diagIndex[i] = count;
		continue;
	    }

	    dx	 = xCrd[i] - xCrd[j];
	    dy	 = yCrd[i] - yCrd[j];
	    dist = sqrt(dx*dx + dy*dy);

	    stnHd->weights[k] = stnWeight(dist, radius);
	    stnHd->dist[k]    = dist;
	    stnHd->dNum[i]   += stnHd->weights[k];
	    count++;
	}
	if (i < nFluidPoints) totDNum += stnHd->dNum[i];
    }

    if (stnHd->n0 == 0) {
	stnHd->n0 = totDNum/nFluidPoints;
    }

    n0Lim = stnHd->n0 * stnHd->beta;

    count   = 0;
    maxDNum = 0;
    
    for (int i = 0; i < nFluidPoints; i++) {
	maxDNum	= MAX(maxDNum, stnHd->dNum[i]);
	stnHd->freeSurf[i] = 0;
	if (stnHd->dNum[i] < n0Lim) {
	    stnHd->freeSurf[i] = 1;
	    count++;
	}
    }
    
    for (int i = nFluidPoints; i < nPoints; i++) {
	stnHd->freeSurf[i] = 0;
    }
}

/*******************************************************************************
 * "stnFree": free the structure
 *******************************************************************************
 */
void stnFree(StnHd stnHd) {
    free(stnHd->col);
    free(stnHd->row);
    free(stnHd->weights);
    free(stnHd->dist);
    free(stnHd->dNum);
    free(stnHd->diagIndex);
    
    free(stnHd);
}
