/*******************************************************************************
 ** Copyright 2014-2014 Vedaad Shakib Inc.
 *******************************************************************************/

/*******************************************************************************
 ** 
 ** "slv.h": A module for solving differential equations
 **
 *******************************************************************************/


/*******************************************************************************
 * Abstract functions
 *******************************************************************************
 */
double* slvCalcExplicitVelocity(StnHd   stnHd,   double *vel, double viscosity,
				double *velStep, double  dt,  double force);
double* slvCalcPressure(StnHd   stnHd, double *xCrd, double *yCrd,
			double *xVel,  double *yVel, double *pressure,
			double  dt,    double  density);
