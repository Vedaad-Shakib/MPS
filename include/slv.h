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

double  slvCalcLaplacian(StnHd stnHd, double *vel, int i);
double* slvCalcExplicitVelocity(StnHd  stnHd, double *vel, double viscosity,
				double dt,    double  force);
double* slvCalcInitialPressure(StnHd stnHd, double dt, double density);
