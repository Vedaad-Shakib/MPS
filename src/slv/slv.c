/*******************************************************************************
 ** Copyright 2014-2014 Vedaad Shakib Inc.
 *******************************************************************************/

/*******************************************************************************
 ** 
 ** "slv.c": A module for solving differential equations
 **
 *******************************************************************************/

/*******************************************************************************
 * Standard includes
 *******************************************************************************/

#include "sys.h"
#include "stn.h"
#include "slv.h"

/*******************************************************************************
 * "slvCalcDivergence": calculates the divergence for a point
 *******************************************************************************
 */

double slvCalcDivergence(StnHd   stnHd, double *xCrd, double *yCrd,
			 double *xVel,  double *yVel, int     i) {
    double	div;		/* the divergence of the specified point */

    div = 0;
    
    for (int k = stnHd->col[i]; k < stnHd->col[i+1]; k++) {
	int j = stnHd->row[k];
	if (i == j) continue;
	
	div += ((xVel[j]-xVel[i])*(xCrd[j]-xCrd[i]) + (yVel[j]-yVel[i])*(yCrd[j]-yCrd[i]))
	       * stnHd->weights[k] / stnHd->dist[k] / stnHd->dist[k];
    }

    return 2 * div / stnHd->n0;
}

/*******************************************************************************
 * "slvCalcLaplacian": calculates the laplacian for a point
 *******************************************************************************
 */
double slvCalcLaplacian(StnHd stnHd, double *vel, int i) {
    double	lap;		/* the laplacian of the specified point */

    lap = 0;

    for (int k = stnHd->col[i]; k < stnHd->col[i+1]; k++) {
	int j = stnHd->row[k];
	if (i == j) continue;
	else lap += (vel[j] - vel[i]) * stnHd->weights[k] / (stnHd->dist[k]*stnHd->dist[k]);
    }

    return 4 * lap / stnHd->n0;
}

/*******************************************************************************
 * "slvCalcExplicitVelocity": calculates the explicit velocity of the next step, or u_i*
 *******************************************************************************
 */
void slvCalcExplicitVelocity(StnHd   stnHd,   double *vel, double viscosity,
	         	     double *velStep, double  dt,  double force) {
    double	lap;		/* the laplacian of a certain point */

    for (int i = 0; i < stnHd->nPoints; i++) {
	if (i >= stnHd->nFluidPoints)
	    velStep[i] = 0;
	else {
	    lap = slvCalcLaplacian(stnHd, vel, i);
	    velStep[i] = vel[i] + dt*(viscosity*lap + force);
	}
    }
}

/*******************************************************************************
 * "slvIsSurfacePoint": returns whether a specified point is a free surface point
 *******************************************************************************
 */

bool slvIsSurfacePoint(StnHd stnHd, int i) {
    return stnHd->dNum[i] <= stnHd->beta*stnHd->n0;
}

/*******************************************************************************
 * "slvCalcInitialPressure": calculates the initial pressure of the problem
 *******************************************************************************
 */
void slvCalcPressure(StnHd   stnHd, double *xCrd, double *yCrd,
		     double *xVel,  double *yVel, double *pressure,
		     double  dt,    double  density) {
    double	lhs[stnHd->col[stnHd->nPoints]];	/* the lhs of the initial pressure equation */
    double	rhs[stnHd->nPoints];	                /* the rhs of the initial pressure equation */
    double	alpha1;		                        /* α1 constant */
    double	alpha2;		                        /* α2 constant */

    alpha1 = 0.8;
    alpha2 = 0.2;

    // calculate lhs
    for (int i = 0; i < stnHd->col[stnHd->nPoints]; i++) lhs[i] = 0;
    
    for (int i = 0; i < stnHd->nPoints; i++) {
	int l = stnHd->diagIndex[i]; // index of diagonals in CCS matrix format
	for (int k = stnHd->col[i]; k < stnHd->col[i+1]; k++) {
	    int j = stnHd->row[k];
	    if (i == j) continue;

	    int t = 4 * stnHd->weights[k] / stnHd->dist[k] / stnHd->dist[k] / stnHd->n0;
	    if (slvIsSurfacePoint(stnHd, i) || slvIsSurfacePoint(stnHd, j)) {
		lhs[k] = 0;
	    } else {
		lhs[k]  = -t;
	    }

	    if (slvIsSurfacePoint(stnHd, i)) {
		lhs[l] = 1;
	    } else {
		lhs[l] += t;
	    }
	}
    }

    // calculate rhs
    for (int i = 0; i < stnHd->nPoints; i++) {
	if (slvIsFreeSurface(i))
	    rhs[i] = 0;
	else
	    rhs[i] = -(alpha1 * density/dt * slvCalcDivergence(stnHd, xCrd, yCrd, xVel, yVel, i)) +
		      (alpha2 * density/(dt*dt) * (stnHd->dNum[i]-stnHd->n0)/stnHd->n0);
    }

    cxsSolveSym(stnHd->col, stnHd->row, lhs,
		rhs,        pressure,   stnHd->nPoints,
		stnHd->col[stnHd->nPoints])
}
