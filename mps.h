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
typedef struct _MpsCorners {
    int                 maxCorners;        /* capacity of the cornerCrds array  */
    int                 nCorners;        /* number of corners                     */
    double        *cornerCrds;        /* the coordinates of the corners    */
} MpsCorners;

typedef struct _MpsWallPoints {
    int                 maxWallPoints;        /* capacity of the wallPointCrds array */
    int                 nWallPoints;        /* number of wall points                */
    double        *wallPointCrds;        /* the coordinates of the wall points  */
} MpsWallPoints;

typedef struct _MpsGhostPoints {
    int                 maxGhostPoints;/* capacity of the ghostPointCrds array */
    int                 nGhostPoints;        /* number of ghost points               */
    double        *ghostPointCrds;/* the coordinates of the ghost points  */
} MpsGhostPoints;

typedef struct _MpsFluidPoints {
    int                 maxFluidPoints;/* capacity of the fluidPointCrds array */
    int                 nFluidPoints;        /* number of fluid points                 */ 
    double        *fluidPointCrds;/* the coordinates of the fluid points  */
} MpsFluidPoints;

typedef MpsCorners* MpsCornersHd;
typedef MpsWallPoints* MpsWallPointsHd;
typedef MpsGhostPoints* MpsGhostPointsHd;
typedef MpsFluidPoints* MpsFluidPointsHd;

/*******************************************************************************
 * Corner function definitions
 *******************************************************************************
 */
MpsCornersHd	mpsNewCorners();
int             mpsGetCornerId(MpsCornersHd cornerHd, double x, double y);
void            mpsFreeCorners(MpsCornersHd cornerHd);
#define         mpsGetNCorners(C)   ((C)->nCorners)

/*******************************************************************************
 * Wall point function definitions
 *******************************************************************************
 */

MpsWallPointsHd mpsNewWallPoints();
int             mpsGetWallPointId(MpsWallPointsHd wallPointsHd, double x, double y);
void            mpsFreeWallPoints(MpsWallPointsHd wallPointsHd);
#define         mpsGetNWallPoints(C)     ((C)->nWallPoints)

/*******************************************************************************
 * Ghost point function definitions
 *******************************************************************************
 */

MpsGhostPointsHd mpsNewGhostPoints();
int		 mpsGetGhostPointId(MpsGhostPointsHd ghostPointsHd, double x, double y);
void		 mpsFreeGhostPoints(MpsGhostPointsHd ghostPointsHd);
#define          mpsGetNGhostPoints(C)   ((C)->nGhostPoints)

/*******************************************************************************
 * Fluid point function definitions
 *******************************************************************************
 */

MpsFluidPointsHd mpsNewFluidPoints();
int		 mpsGetFluidPointId(MpsFluidPointsHd fluidPointsHd, double x, double y);
void		 mpsFreeFluidPoints(MpsFluidPointsHd fluidPointsHd);
#define          mpsGetNFluidPoints(C)   ((C)->nFluidPoints)

/*******************************************************************************
 * Mps function definitions
 *******************************************************************************
 */

void   mpsOutCrd(char *fileName, double *crd, int nPoints);
int    mpsContainsLine(double *wallSegments, int nWallSegments, double x1, double y1, double x2, double y2);
double mpsDist(int x1, int y1, int x2, int y2);
void   mpsGhostCorners(MpsCornersHd cornersHd, double *wallSegments, int nWallSegments,
		       MpsGhostPointsHd ghostPointsHd, int nGhostPoints, double radius);
void   mpsConstructIntermediatePoints(double *cornerPoints, int nCornerPoints, double **outPoints,
				   int *nOutPoints, int *maxOutPoints, double wallSpacing,
				   bool containsStart, bool containsEnd);
bool   mpsCheckClosure(double *fluidBoundaries, int nFluidBoundaries);
int    mpsIntegerize(double min, double val, double wallSpacing);
bool   mpsCrossesFluidBoundaries(double *fluidBoundaries, int nFluidBoundaries,
				 double x1, double y1, double x2, double y2);
bool   mpsLinesIntersect(double a1x, double a1y, double a2x, double a2y,
			 double b1x, double b1y, double b2x, double b2y);
void   mpsFloodfill(MpsFluidPointsHd fluidPointsHd, double *fluidBoundaries, int nFluidBoundaries,
		    double initX, double initY, double wallSpacing);
int    main();

