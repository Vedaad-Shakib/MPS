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
 * "slvCalcLaplacian": calculates the laplacian operator for a point
 *******************************************************************************
 */
double slvCalcLaplacian(StnHd stnHd, double *vel, int i) {
    double lap; /* the laplacian of the specified point */
    int    j;   /* adjacent point */

    lap = 0;

    for (int k = stnHd->col[i]; k < stnHd->col[i+1]; k++) {
	j = stnHd->row[k];
	if (i == j) continue;
	else lap += (vel[j] - vel[i]) * stnHd->weights[k] / (stnHd->dist[k]*stnHd->dist[k]);
    }

    return 4 * lap / stnHd->n0;
}

/*******************************************************************************
 * "slvCalcExplicitVelocity": calculates the explicit velocity of the next step, or u_i*
 *******************************************************************************
 */
double* slvCalcExplicitVelocity(StnHd  stnHd, double *vel, double viscosity,
				double dt,    double  force) {
    double *velStep; /* the explicit velocity of that time step */
    double  lap;     /* the laplacian of a certain point */

    velStep = memNew(double, stnHd->nPoints);

    for (int i = 0; i < stnHd->nPoints; i++) {
	if (i >= stnHd->nFluidPoints)
	    velStep[i] = 0;
	else {
	    lap = slvCalcLaplacian(stnHd, vel, i);
	    velStep[i] = vel[i] + dt*(viscosity*lap + force);
	}
    }

    return velStep;
}

/*******************************************************************************
 * "slvCalcInitialPressure": calculates the initial pressure of the problem
 *******************************************************************************
 */
double* slvCalcInitialPressure(StnHd stnHd, double dt, double density) {
    double lhs[stnHd->col[stnHd->nPoints]]; /* the lhs of the initial pressure equation */
    double rhs[stnHd->nPoints];             /* the rhs of the initial pressure equation */

    // calculate lhs
    for (int i = 0; i < stnHd->col[stnHd->nPoints]; i++) lhs[i] = 0;
    
    for (int i = 0; i < stnHd->nPoints; i++) {
	int l = stnHd->diagIndex[i];
	for (int k = stnHd->col[i]; k < stnHd->col[i+1]; k++) {
	    int j = stnHd->row[k];
	    if (i == j) continue;
	    int t = stnHd->weights[k]/stnHd->dist[k];
	    lhs[k] -= t;
	    lhs[l] += t;
	}
    }

    // calculate rhs
    for (int i = 0; i < stnHd->nPoints; i++) {
	rhs[i] = -density/(dt*dt) * (stnHd->dNum[i]-stnHd->n0)/stnHd->n0;
    }

    // also multiply by -1
    // count influence of the wall points or not?
    // solve matrix Ax=B

    return NULL; // to avoid error
}
