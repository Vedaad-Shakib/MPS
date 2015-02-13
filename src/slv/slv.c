/*******************************************************************************
 ** Copyright 2014-2014 Vedaad Shakib Inc.
 *******************************************************************************/

/*******************************************************************************
 ** 
 ** "slv.c": A module for solving equations related to MPS
 **
 *******************************************************************************/

/*******************************************************************************
 * Standard includes
 *******************************************************************************/

#include "sys.h"
#include "slv.h"
#include "stn.h"

double* calcInitialPressure(StnHd stnHd, double dt, double density) {
    double lhs[stnHd->col[stnHd->nPoints]]; /* the lhs of the initial pressure equation */
    double rhs[stnHd->nPoints];             /* the rhs of the initial pressure equation */

    // calculate lhs
    for (int i = 0; i < count; i++) stnHd->lhs[i] = 0;
    
    for (int i = 0; i < nFluidPoints+nWallPoints; i++) {
	l = diagIndex[i];
	for (int k = col[i]; k < col[i+1]; k++) {
	    j = row[k];
	    if (i == j) continue;
	    t = weights[k]/dist[k];
	    lhs[k] -= t;
	    lhs[l] += t;
	}
    }

    // calculate rhs
    for (int i = 0; i < stnHd->nPoints; i++) {
	rhs[i] = -density/(dt*dt) * (stnHd->dNum[i]-n0)/n0;
    }

    // solve matrix Ax=B
}
