/*******************************************************************************          
 ** Copyright 2014-2014 Vedaad Shakib Inc.
 *******************************************************************************/

/*******************************************************************************
 ** 
 ** "mpsDriver.c": Moving Particle Semi-implicit driver
 **
 *******************************************************************************/

/*******************************************************************************
 * Standard includes
 *******************************************************************************/

#include "sys.h"
#include "mps.h"
#include "que.h"
#include "stn.h"
#include "slv.h"

/*******************************************************************************
 * LIM: limit update
 *******************************************************************************
 */
#define        LIM(X,T)        MIN(MAX((X),-(T)),(T))

/*******************************************************************************
 * "mpsDriver": MPS time stepping driver
 *******************************************************************************
 */
void mpsDriver(MpsHd mpsHd) {
    StnHd        stnHd;         /* an adjacency structure */
    char         buffer[1024];  /* a string buffer */
    double       dNum0;         /* wall spacing based mpsHd->density number */
    double       dr;            /* distance */
    double       dx;            /* change in x */
    double       dy;            /* change in y */
    double       l2Norm;        /* L2 norm of vector */
    double      *presCurr;      /* the pressure for the current time step */
    double      *presNext;      /* the pressure for the next time step */
    double      *xPosCurr;      /* the x position for the current time step */
    double      *xPosNext;      /* the x position for the next time step */
    double      *xPosStar;      /* the x position for the next time step */
    double      *xVelCurr;      /* the x velocity for the current time step */
    double      *xVelNext;      /* the x velocity for the next time step */
    double      *xVelStar;      /* the x velocity correction based on the pressure */
    double      *yPosCurr;      /* the y position for the current time step */
    double      *yPosNext;      /* the y position for the next time step */
    double      *yPosStar;      /* the y position for the next time step */
    double      *yVelCurr;      /* the y velocity for the current time step */
    double      *yVelNext;      /* the y velocity for the next time step */
    double      *yVelStar;      /* the y velocity correction based on the pressure */
    double       limX;          /* limit the update */
    int          n;             /* a count */
    double       m;             /* a count */
    int          nPoints;       /* total number of points */
    int          offSet;        /* point offset */
    int          nFluidPoints;  /* the total number of fluid points */
    int          nWallPoints;   /* the total number of wallPoints */

/*---------------------------------------------------------------------------------------
 * Get the dimensions
 *---------------------------------------------------------------------------------------
 */
    nFluidPoints = mpsHd->fluidPointsHd->nPoints;
    nWallPoints  = mpsHd->wallPointsHd->nPoints + mpsHd->ghostPointsHd->nPoints;
    nPoints      = nFluidPoints + nWallPoints;

    limX = 0.5 * mpsHd->wallSpacing;

/*---------------------------------------------------------------------------------------
 * Get a measure of mpsHd->density number of an interior point
 *---------------------------------------------------------------------------------------
 */
    dNum0 = 0;
    n = (int) ((2/sqrt(3)) * mpsHd->radius / mpsHd->wallSpacing);
    for (int i = -n; i <= n; i++) {
	m = (abs(i)%2==1) ? 0.5*mpsHd->wallSpacing : 0;
	while (sqrt((m*m)+(i*2*mpsHd->wallSpacing/sqrt(3)*i*2*mpsHd->wallSpacing/sqrt(3))) < mpsHd->radius) {
	    dx = m;
	    dy = i*mpsHd->wallSpacing/(2/sqrt(3));
	    dr = sqrt(dx*dx + dy*dy);
	    dNum0 += stnWeight(dr, mpsHd->radius);
	    
	    m += mpsHd->wallSpacing;
	}
	m = (abs(i)%2==1) ? -0.5*mpsHd->wallSpacing : -mpsHd->wallSpacing;
	while (sqrt((m*m)+(i*2*mpsHd->wallSpacing/sqrt(3)*i*2*mpsHd->wallSpacing/sqrt(3))) < mpsHd->radius) {
	    dx = m;
	    dy = i*mpsHd->wallSpacing/(2/sqrt(3));
	    dr = sqrt(dx*dx + dy*dy);
	    dNum0 += stnWeight(dr, mpsHd->radius);
	    
	    m -= mpsHd->wallSpacing;
	}
    }
    printf("Interior point num dens. = %g\n", dNum0);
    
/*---------------------------------------------------------------------------------------
 * Allocate arrays
 *---------------------------------------------------------------------------------------
 */
    xVelCurr = memNewZero(double, nPoints);
    yVelCurr = memNewZero(double, nPoints);
    xVelNext = memNewZero(double, nPoints);
    yVelNext = memNewZero(double, nPoints);
    xVelStar = memNewZero(double, nPoints);
    yVelStar = memNewZero(double, nPoints);

    xPosCurr = memNewZero(double, nPoints);
    yPosCurr = memNewZero(double, nPoints);
    xPosNext = memNewZero(double, nPoints);
    yPosNext = memNewZero(double, nPoints);
    xPosStar = memNewZero(double, nPoints);
    yPosStar = memNewZero(double, nPoints);

    presCurr = memNewZero(double, nPoints);
    presNext = memNewZero(double, nPoints);

/*---------------------------------------------------------------------------------------
 * Initialize the position
 *---------------------------------------------------------------------------------------
 */
    offSet = 0;
    for (int i = 0; i < mpsHd->fluidPointsHd->nPoints; i++) {
        xPosCurr[offSet+i] = mpsHd->fluidPointsHd->pointCrds[2*i+0];
        yPosCurr[offSet+i] = mpsHd->fluidPointsHd->pointCrds[2*i+1];
    }
    offSet += mpsHd->fluidPointsHd->nPoints;
    for (int i = 0; i < mpsHd->wallPointsHd->nPoints; i++) {
        xPosCurr[offSet+i] = mpsHd->wallPointsHd->pointCrds[2*i+0];
        yPosCurr[offSet+i] = mpsHd->wallPointsHd->pointCrds[2*i+1];
    }
    offSet += mpsHd->wallPointsHd->nPoints;
    for (int i = 0; i < mpsHd->ghostPointsHd->nPoints; i++) {
        xPosCurr[offSet+i] = mpsHd->ghostPointsHd->pointCrds[2*i+0];
        yPosCurr[offSet+i] = mpsHd->ghostPointsHd->pointCrds[2*i+1];
    }
    offSet += mpsHd->ghostPointsHd->nPoints;

    snprintf(buffer, sizeof(buffer), "mps.%d.out", 0);
    mpsOutCrdXY(buffer, xPosCurr, yPosCurr, nFluidPoints);

    snprintf(buffer, sizeof(buffer), "mps.vel.%d.out", 0);
    mpsOutCrdXY(buffer, xVelCurr, yVelCurr, nFluidPoints);
    
/*---------------------------------------------------------------------------------------
 * Initialize the values
 *---------------------------------------------------------------------------------------
 */
    for (int i = 0; i < nPoints; i++) {
        xVelCurr[i] = 0;
        yVelCurr[i] = 0;
        presCurr[i] = 0;
    }
/*---------------------------------------------------------------------------------------
 * Create and initialize the search data structure
 *---------------------------------------------------------------------------------------
 */
    stnHd = stnNew(nFluidPoints, nWallPoints, mpsHd->radius, mpsHd->beta);

    stnPopulate(stnHd, xPosCurr, yPosCurr);
    stnCalc(stnHd, xPosCurr, yPosCurr);
    stnHd->n0 = (stnHd->n0+dNum0)/2;
    printf("Initial number mpsHd->density = %g\n", stnHd->n0);

    // print stnHd
    /*for (int i = 0; i < nFluidPoints; i++) {
	printf("for point %d; pos %lf, %lf; dnum %lf\n", i, xPosCurr[i], yPosCurr[i], stnHd->dNum[i]);
	for (int k = stnHd->col[i]; k < stnHd->col[i+1]; k++) {
	    int j = stnHd->row[k];
	    printf("\t%lf, %lf; dist: %lf \n", xPosCurr[j], yPosCurr[j], stnHd->dist[k]);
	}
    }*/

/*---------------------------------------------------------------------------------------
 * Smooth the initial condition
 *---------------------------------------------------------------------------------------
 */
    // mpsOutCrdXY("x0.dat", xPosCurr, yPosCurr, nPoints);
    // slvSmoothInit(stnHd, xPosCurr, yPosCurr, xPosStar, yPosStar);
    // mpsOutCrdXY("x1.dat", xPosCurr, yPosCurr, nPoints);

/*---------------------------------------------------------------------------------------
 * Loop over the time steps
 *---------------------------------------------------------------------------------------
 */
    printf("Calculating for %d time-steps\n", mpsHd->nTimeSteps);
    for (int stepId = 0; stepId < mpsHd->nTimeSteps; stepId++) {

        printf("===============> Processing time step %d ; time increment %g\n", 
               stepId+1, mpsHd->dt);

/*---------------------------------------------------------------------------------------
 * Find the neighboring points
 *---------------------------------------------------------------------------------------
 */
        stnPopulate(stnHd, xPosCurr, yPosCurr);
	stnCalc(stnHd, xPosCurr, yPosCurr);

        if (stepId == 0) mpsOutCrd("mpsHd->density.dat", stnHd->dNum, nFluidPoints, 1);

        snprintf(buffer, sizeof(buffer), "mps.dens.%d.out", stepId+1);
        if (stepId%10==9) mpsOutCrd(buffer, stnHd->dNum, nFluidPoints, 1);

/*---------------------------------------------------------------------------------------
 * Advance the explicit part
 *---------------------------------------------------------------------------------------
 */
        slvCalcExplicitVelocity(stnHd, xVelCurr, xVelStar, mpsHd->viscosity, mpsHd->dt, 0);
        slvCalcExplicitVelocity(stnHd, yVelCurr, yVelStar, mpsHd->viscosity, mpsHd->dt, -9.8);

        l2Norm = mps2VecL2(xVelStar, yVelStar, nFluidPoints);
        // printf("Explicit vel. L2 norm  = %g\n", l2Norm);

	// if (l2Norm > 100.) exit(1);

        for (int i = 0; i < nPoints; i++) {
            xPosStar[i] = xPosCurr[i] + LIM(mpsHd->dt*xVelStar[i], limX);
            yPosStar[i] = yPosCurr[i] + LIM(mpsHd->dt*yVelStar[i], limX);
        }
/*---------------------------------------------------------------------------------------
 * Compute the pressure
 *---------------------------------------------------------------------------------------
 */
        //stnCalc(stnHd, xPosStar, yPosStar);

        slvCalcPressure(stnHd,    xPosStar, yPosStar,
                        xVelStar, yVelStar, presNext,
                        mpsHd->dt,       mpsHd->density);

        snprintf(buffer, sizeof(buffer), "mps.pres.%d.out", stepId+1);
        if (stepId%10==9) mpsOutCrd(buffer, presNext, nFluidPoints, 1);

/*---------------------------------------------------------------------------------------
 * Correct the velocity
 *---------------------------------------------------------------------------------------
 */
        slvCalcCorrection(stnHd,       presNext, xVelNext,
                          xPosCurr,    mpsHd->density,  mpsHd->dt);
        slvCalcCorrection(stnHd,       presNext, yVelNext,
                          yPosCurr,    mpsHd->density,  mpsHd->dt);
        
/*---------------------------------------------------------------------------------------
 * Update the values
 *---------------------------------------------------------------------------------------
 */
        for (int i = 0; i < nPoints; i++) {
            xVelNext[i] = xVelNext[i] + xVelStar[i];
            yVelNext[i] = yVelNext[i] + yVelStar[i];
            
            xPosNext[i] = xPosCurr[i] + LIM(mpsHd->dt*xVelNext[i],limX);
            yPosNext[i] = yPosCurr[i] + LIM(mpsHd->dt*yVelNext[i],limX);
        }

        l2Norm = mps2VecL2(xVelNext, yVelNext, nFluidPoints);
        // printf("Implicit vel. L2 norm  = %g\n", l2Norm);

	// if (l2Norm > 100.) exit(1);

/*---------------------------------------------------------------------------------------
 * Output the data
 *---------------------------------------------------------------------------------------
 */
        snprintf(buffer, sizeof(buffer), "mps.%d.out", stepId+1);
        if (stepId%10==9) mpsOutCrdXY(buffer, xPosNext, yPosNext, nFluidPoints);

        snprintf(buffer, sizeof(buffer), "mps.vel.%d.out", stepId+1);
        if (stepId%10==9) mpsOutCrdXY(buffer, xVelCurr, yVelCurr, nFluidPoints);

/*---------------------------------------------------------------------------------------
 * Reset the vectors
 *---------------------------------------------------------------------------------------
 */
        for (int i = 0; i < nPoints; i++) {
            xPosCurr[i] = xPosNext[i];
            yPosCurr[i] = yPosNext[i];
            xVelCurr[i] = xVelNext[i];
            yVelCurr[i] = yVelNext[i];
            presCurr[i] = presNext[i];
        }
    }
}
