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

MpsPointsHd	mpsNewPoints();
int		mpsGetPointId(MpsPointsHd pointsHd, double x, double y);
void		mpsFreePoints(MpsPointsHd PointsHd);
#define         mpsGetNPoints(C)	((C)->nPoints)

/*******************************************************************************
 * Mps function definitions
 *******************************************************************************
 */

void   mpsOutCrd(char *fileName, double *crd, int nPoints);
int    mpsContainsLine(double *wallSegments, int nWallSegments, double x1, double y1, double x2, double y2);
double mpsDist(double x1, double y1, double x2, double y2);
void   mpsGhostCorners(MpsPointsHd cornersHd, double *wallSegments, int nWallSegments,
                       MpsPointsHd ghostPointsHd, int nGhostPoints, double radius);
void   mpsConstructIntermediatePoints(double *cornerPoints, int nCornerPoints, int nSegs, double **outPoints, 
                                      int *nOutPoints, int *maxOutPoints, double wallSpacing,
                                      bool containsStart, bool containsEnd);
bool   mpsCheckClosure(double *fluidBoundaries, int nFluidBoundaries);
int    mpsIntegerize(double min, double val, double wallSpacing);
bool   mpsCrossesFluidBoundaries(double *fluidBoundaries, int nFluidBoundaries,
                                 double *wallSegments, int nWallSegments,
                                 double x1, double y1, double x2, double y2);
bool   mpsLinesIntersect(double a1x, double a1y, double a2x, double a2y,
                         double b1x, double b1y, double b2x, double b2y);
void   mpsFloodfill(MpsPointsHd fluidPointsHd, double *fluidBoundaries,
                    int nFluidBoundaries, double* wallSegments, int nWallSegments,
                    double initX, double initY, double wallSpacing);
int    main();

