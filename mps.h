/*******************************************************************************
** Copyright 2014-2014 Vedaad Shakib Inc.
*******************************************************************************/

/*******************************************************************************
** 
** "mps.h": Moving Particle Semi-implicit
**
*******************************************************************************/

/*******************************************************************************
 * Standard includes
 *******************************************************************************
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include "mem.h"

/*******************************************************************************
 * Corner data
 *******************************************************************************
 */
typedef struct _MpsCorners {
  int		 maxCorners;	/* capacity of the cornerCrds array  */
  int		 nCorners;	/* number of corners		     */
  double	*cornerCrds;	/* the coordinates of the corners    */
} MpsCorners;

typedef struct _MpsWallPoints {
  int		 maxWallPoints;	/* capacity of the wallPointCrds array */
  int		 nWallPoints;	/* number of wall points 	       */
  double	*wallPointCrds;	/* the coordinates of the wall points  */
} MpsWallPoints;

typedef struct _MpsGhostPoints {
  int		 maxGhostPoints;/* capacity of the ghostPointCrds array */
  int		 nGhostPoints;	/* number of ghost points               */
  double	*ghostPointCrds;/* the coordinates of the ghost points  */
} MpsGhostPoints;

typedef struct _MpsFluidPoints {
  int		 maxFluidPoints;/* capacity of the fluidPointCrds array */
  int		 nFluidPoints;	/* number of fluid points 	        */ 
  double	*fluidPointCrds;/* the coordinates of the fluid points  */
} MpsFluidPoints;

typedef MpsCorners* MpsCornersHd;
typedef MpsWallPoints* MpsWallPointsHd;
typedef MpsGhostPoints* MpsGhostPointsHd;
typedef MpsFluidPoints* MpsFluidPointsHd;

/*******************************************************************************
 * Corner function definitions
 *******************************************************************************
 */
MpsCornersHd mpsNewCorners();
int	     mpsGetCornerId(MpsCornersHd cornerHd, double x, double y);
void	     mpsFreeCorners(MpsCornersHd cornerHd);
#define	     mpsGetNCorners(C)	((C)->nCorners)

/*******************************************************************************
 * Wall point function definitions
 *******************************************************************************
 */

MpsWallPointsHd mpsNewWallPoints();
int		mpsGetWallPointId(MpsWallPointsHd wallPointsHd, double x, double y);
void		mpsFreeWallPoints(MpsWallPointsHd wallPointsHd);
#define	        mpsGetNWallPoints(C)	((C)->nWallPoints)

/*******************************************************************************
 * Ghost point function definitions
 *******************************************************************************
 */

MpsGhostPointsHd mpsNewGhostPoints();
int		 mpsGetGhostPointId(MpsGhostPointsHd ghostPointsHd, double x, double y);
void	  	 mpsFreeGhostPoints(MpsGhostPointsHd ghostPointsHd);
#define	         mpsGetNGhostPoints(C)	((C)->nGhostPoints)

/*******************************************************************************
 * Fluid point function definitions
 *******************************************************************************
 */

MpsFluidPointsHd mpsNewFluidPoints();
int		 mpsGetFluidPointId(MpsFluidPointsHd fluidPointsHd, double x, double y);
void		 mpsFreeFluidPoints(MpsFluidPointsHd fluidPointsHd);
#define	         mpsGetNFluidPoints(C)	((C)->nFluidPoints)
