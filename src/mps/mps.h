/*******************************************************************************
 ** Copyright 2014-2014 Vedaad Shakib Inc.
 *******************************************************************************/

/*******************************************************************************
 ** 
 ** "mps.h": Moving Particle Semi-implicit
 **
 *******************************************************************************/

/*******************************************************************************
 * Point struct data
 *******************************************************************************
 */
typedef struct _MpsPoints {
    int          maxPoints;/* capacity of the fluidPointCrds array */
    int          nPoints;  /* number of fluid points                 */ 
    double      *pointCrds;/* the coordinates of the fluid points  */
} MpsPoints;

typedef MpsPoints* MpsPointsHd;

typedef struct _Mps {
    MpsPointsHd	fluidPointsHd;	/* the struct pointer of fluid points */
    MpsPointsHd	wallPointsHd;	/* the struct pointer of wall points */
    MpsPointsHd	ghostPointsHd;	/* the struct pointer of ghost points */
    double	radius;		/* the radius of influence */
    double	wallSpacing;	/* the spacing between points */
    double	beta;		/* the free-surface constant beta */
    double	density;	/* the density of the fluid */
    double	viscosity;	/* the viscosity of the fluid */
    int 	nTimeSteps;	/* the number of time steps for which the simulation runs */
    double	dt;		/* the time delta between each timestep */
} Mps;

typedef Mps*     MpsHd;

MpsPointsHd	mpsNewPoints();
int		mpsGetPointId(MpsPointsHd pointsHd, double x, double y);
void		mpsFreePoints(MpsPointsHd PointsHd);
#define         mpsGetNPoints(C)	((C)->nPoints)

/*******************************************************************************
 * Mps function definitions
 *******************************************************************************
 */

void   mpsOutCrd(char *fileName, double *crd, int nPoints, int nDims);
void   mpsOutCrdXY(char *fileName, double *xCrd, double *yCrd,
                 int   nPoints);
double mpsVecL2(double *vec, int nDims);
double mps2VecL2(double *vec1, double *vec2, int nDims);
int    mpsContainsLine(double *wallSegments, int nWallSegments, double x1, double y1, double x2, double y2);
double mpsDist(double x1, double y1, double x2, double y2);
void   mpsGhostCorners(MpsPointsHd cornersHd,     double *wallSegments, int nWallSegments,
                       MpsPointsHd ghostPointsHd, int nGhostPoints,     double radius,
                       double      wallSpacing);
void   mpsConstructIntermediatePoints(double     *cornerPoints, int    nCornerPoints, int nSegs, 
				      MpsPointsHd pointsHd,     double wallSpacing,
                                      bool containsStart, bool containsEnd);
bool   mpsCheckClosure(double *fluidBoundaries, int nFluidBoundaries);
int    mpsIntegerize(double min, double val, double wallSpacing);
bool   mpsCrossesFluidBoundary(double *fluidBoundaries, int    nFluidBoundaries, double x1,
			       double  y1,              double x2,               double y2);
bool   mpsCrossesWallBoundary(double *wallSegments, int    nWallSegments, double x1, 
			      double  y1,           double x2,            double y2);
bool   mpsCrossesBoundary(double *fluidBoundaries, int nFluidBoundaries,
			  double *wallSegments,    int nWallSegments,
			  double  x1,              double y1,
			  double  x2,              double y2);
bool   mpsLinesIntersect(double a1x, double a1y, double a2x, double a2y,
                         double b1x, double b1y, double b2x, double b2y);
void   mpsFloodfill(MpsPointsHd fluidPointsHd,    double *fluidBoundaries, int    nFluidBoundaries,
                    double     *wallSegments,     int     nWallSegments,   double initX,
		    double      initY,            double  wallSpacing);
MpsHd  mpsRead();
int    main();
void   mpsDriver(MpsHd mpsHd);

