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
#include "cxs.h"

/*******************************************************************************
 * "slvCalcGradient": calculates the gradient of a vector
 *******************************************************************************
 */
double slvCalcGradient(StnHd stnHd, double *vec, double *pos, int i) {
    double	grad;		/* the gradient */
    double      minVec;		/* the minimum value */

    grad = 0;
    minVec = vec[i];
    if (true) {
	for (int k = stnHd->col[i]; k < stnHd->col[i+1]; k++) {
	    int j = stnHd->row[k];
	    if (vec[j] < minVec) {
		minVec = vec[j];
	    }
	}
    }
    minVec = (1*minVec + vec[i])/2;
    for (int k = stnHd->col[i]; k < stnHd->col[i+1]; k++) {
	int j = stnHd->row[k];
	if (i == j) continue;

	grad += (vec[j]-minVec)*(pos[j]-pos[i])
	       * stnHd->weights[k] / (stnHd->dist[k]*stnHd->dist[k]);
    }

    return stnHd->d * grad / stnHd->n0;
}

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
	       * stnHd->weights[k] / (stnHd->dist[k]*stnHd->dist[k]);
    }

    return stnHd->d * div / stnHd->n0;
}

/*******************************************************************************
 * "slvCalcLaplacian": calculates the laplacian for a point
 *******************************************************************************
 */
double slvCalcLaplacian(StnHd stnHd, double *vec, int i) {
    double	lap;		/* the laplacian of the specified point */

    lap = 0;

    for (int k = stnHd->col[i]; k < stnHd->col[i+1]; k++) {
	int j = stnHd->row[k];
	if (i == j) continue;
	else lap += (vec[j] - vec[i]) * stnHd->weights[k] / (stnHd->dist[k]*stnHd->dist[k]);
    }

    return 2 * stnHd->d * lap / stnHd->n0;
}

/*******************************************************************************
 * "slvCalcExplicitVelocity": calculates the explicit velocity of the next step, or u_i*
 *******************************************************************************
 */
void slvCalcExplicitVelocity(StnHd   stnHd,     double *vel, double *velStar,
	         	     double  viscosity, double  dt,  double force) {
    double	lap;		/* the laplacian of a certain point */
    int		nPoints;	/* no. points */
    int		nFluidPoints;	/* no. fluid points */

    nPoints		= stnHd->nPoints;
    nFluidPoints	= stnHd->nFluidPoints;

    for (int i = 0; i < nFluidPoints; i++) {
	lap = slvCalcLaplacian(stnHd, vel, i);
	velStar[i] = vel[i] + dt*(viscosity*lap + force);
    }
    for (int i = nFluidPoints; i < nPoints; i++ ) {
	velStar[i] = 0;
    }
}

/*******************************************************************************
 * "slvCalcPressure": calculates the pressure vector and puts it in the
 *                           pressure vector
 *******************************************************************************
 */
void slvCalcPressure(StnHd   stnHd, double *xCrd, double *yCrd,
		     double *xVel,  double *yVel, double *pres,
		     double  dt,    double  density) {
    double	*dNum;		/* density number of each point */
    double	*dist;		/* distance between two points */
    double	*lhs;		/* the lhs of the initial pressure equation */
    double	*rhs;	        /* the rhs of the initial pressure equation */
    double	*weights;	/* kernel weight between two points */
    double	 alpha1;	/* α1 constant */
    double	 alpha2;	/* α2 constant */
    double	 divVel;	/* divergence of velocity */
    double	 fct;		/* lhs factor */
    double	 fct1;		/* rhs factor 1 */
    double	 fct2;		/* rhs factor 2 */
    double	 n0;		/* average density number */
    double	 delDNumLim;	/* maximum change in dNum */
    double	 delDNum;	/* change in dNum */
    int		*col;		/* adjacency col */
    int		*diagIndex;	/* index of diagonal entry */
    int		*freeSurf;	/* free surface flag */
    int		*row;		/* adjacency row */
    int		 nNonzeros;	/* no. nonzeros */
    int		 nPoints;	/* no. points */
    int		 nFluidPoints;	/* no. fluid points */
    int		 errCode;	/* return code */

    nPoints	 = stnHd->nPoints;
    nFluidPoints = stnHd->nFluidPoints;
    col		 = stnHd->col;
    row		 = stnHd->row;
    lhs		 = stnHd->lhs;
    rhs		 = stnHd->rhs;
    weights	 = stnHd->weights;
    dist	 = stnHd->dist;
    diagIndex	 = stnHd->diagIndex;
    freeSurf	 = stnHd->freeSurf;
    nNonzeros	 = col[nPoints];
    n0		 = stnHd->n0;
    dNum	 = stnHd->dNum;

    alpha1     = 0.8;
    alpha2     = 0.2;
    delDNumLim = 0.02 * n0;

    // calculate lhs
    for (int i = 0; i < nNonzeros; i++) lhs[i] = 0;
    
    fct = 2 * stnHd->d / n0 ;
    for (int i = 0; i < nPoints; i++) {
	int l = diagIndex[i]; // index of diagonals in CCS matrix format
	if (freeSurf[i]) {
	    lhs[l] = 1;
	    continue;
	}
	for (int k = col[i]; k < col[i+1]; k++) {
	    int j = row[k];
	    if (i == j) continue;

	    double t = fct * weights[k] / (dist[k] * dist[k]);
	    lhs[l] += t;
	    if (freeSurf[j] == 0) {
		lhs[k]  = -t;
	    }
	}
    }

    double reg = 1.e-8;
    /* regularize the equation a bit */
    for (int i = 0; i < nPoints; i++) {
	int l = diagIndex[i]; // index of diagonals in CCS matrix format
	lhs[l] *= 1+reg;
    }

    // calculate rhs
    fct1 = -alpha1 * density / dt;
    fct2 = alpha2 * density / (dt*dt) / n0;
    for (int i = 0; i < nFluidPoints; i++) {
	if (freeSurf[i]) {
	    rhs[i] = 0;
	} else {
	    divVel = slvCalcDivergence(stnHd, xCrd, yCrd, xVel, yVel, i) ;
    	    delDNum = dNum[i] - n0;
	    delDNum = MIN(delDNum, +delDNumLim);
	    delDNum = MAX(delDNum, -delDNumLim);
	    rhs[i] = fct1 * divVel + fct2 * delDNum;
	}
    }
    for (int i = nFluidPoints; i < nPoints; i++) {
	rhs[i] = 0;
    }

    //mpsOutCrd("rhs.dat", rhs, nPoints, 1);

    errCode = cxsSolveSym(col, row, lhs, rhs, pres, nPoints, nNonzeros);

    if (errCode) {
	printf("Error %d from cxsSolveSym()\n",errCode);
	exit(1);
    }

    // check solution
    if (0) {
	double sum = 0;
	for (int i = 0; i < nPoints; i++) {
	    double t = 0;
	    for (int k = col[i]; k < col[i+1]; k++) {
		int j = row[k];
		t += lhs[k] * pres[j];
	    }
	    printf("Check %d %g %g %g\n", i, t, rhs[i], t - rhs[i]);
	    sum += (t - rhs[i]) * (t - rhs[i]);
	}
	sum = sqrt(sum/nPoints);
	printf("Total solution error = %g\n", sum);
    }
}

/*******************************************************************************
 * "slvCalcCorrection": transform u* into u'
 *******************************************************************************
 */
void slvCalcCorrection(StnHd   stnHd, double *pres, double *velPrime,
		       double *pos,   double  density,  double  dt) {
    for (int i = 0; i < stnHd->nFluidPoints; i++) {
	velPrime[i] = -dt/density * slvCalcGradient(stnHd, pres, pos, i);
    }
    for (int i = stnHd->nFluidPoints; i < stnHd->nPoints; i++) {
	velPrime[i] = 0;
    }
}

/*******************************************************************************
 * "slvSmoothInit": smooth initial fluid position
 *******************************************************************************
 */
void slvSmoothInit(StnHd   stnHd, double *xPos, double *yPos, 
		   double *xPos2, double *yPos2) {

    double	x;
    double	y;

    for (int i = 0; i < stnHd->nFluidPoints; i++) {
	if ( stnHd->freeSurf[i] ) {
	    xPos2[i] = xPos[i];
	    yPos2[i] = yPos[i];
	    continue;
	}
	x = 0;
	y = 0;
	for (int k = stnHd->col[i]; k < stnHd->col[i+1]; k++) {
	    int j = stnHd->row[k];
	    x += xPos[j];
	    y += yPos[j];
	}
	xPos2[i] = x / (stnHd->col[i+1] - stnHd->col[i]);
	yPos2[i] = y / (stnHd->col[i+1] - stnHd->col[i]);
    }
    for (int i = 0; i < stnHd->nFluidPoints; i++) {
	xPos[i]	= xPos2[i];
	yPos[i]	= yPos2[i];
    }
}
