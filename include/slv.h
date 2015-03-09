/*******************************************************************************
 ** Copyright 2014-2014 Vedaad Shakib Inc.
 *******************************************************************************/

/*******************************************************************************
 ** 
 ** "slv.h": A module for solving differential equations
 **
 *******************************************************************************/

#include "cxs.h"

/*******************************************************************************
 * Abstract functions
 *******************************************************************************
 */
void   slvCalcExplicitVelocity(StnHd   stnHd,     double *vel, double *velStep,
			       double  viscosity, double  dt,  double force);
void   slvCalcPressure(StnHd   stnHd, double *xCrd, double *yCrd,
		       double *xVel,  double *yVel, double *pressure,
		       double  dt,    double  density);
double slvCalcDivergence(StnHd   stnHd, double *xCrd, double *yCrd,
			 double *xVel,  double *yVel, int     i);
double slvCalcLaplacian(StnHd stnHd, double *vel, int i);
double slvCalcGradient(StnHd stnHd, double *pressure, double *pos, int i);
void   slvCalcCorrection(StnHd   stnHd, double *pressure, double *velCorrect,
			 double *pos,   double  density,  double  dt);
void slvSmoothInit(StnHd stnHd, double *xPos, double *yPos, 
		   double *xPos2, double *yPos2);
