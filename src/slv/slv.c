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
#include "stn.h"

double* calcInitialPressure(StnHd stnHd, double dt, double density) {
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
