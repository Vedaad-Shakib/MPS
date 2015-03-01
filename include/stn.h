/*******************************************************************************
 ** Copyright 2014-2014 Vedaad Shakib Inc.
 *******************************************************************************/

/*******************************************************************************
 ** 
 ** "stn.h": A module for finding and storing points in a certain radius of other
 **          points
 **
 *******************************************************************************/

/*******************************************************************************
 * Queue struct definition
 *******************************************************************************
 */
// Stored in a harwell-boeing sparse matrix format instead of an adjacency matrix
// which uses much more memory
typedef struct {
    int		*col;		/* the number of entries in each row of the adjacency matrix */
    int  	*row;		/* the points adjacent to each point in turn */
    int          maxAdjacent;   /* the maximum capacity of the row array */
    double      *weights;       /* the weights of each adjacent point */
    double      *dist;          /* the distance between each point */
    double      *dNum;          /* the density number of each point */
    int         *diagIndex;     /* the index of the diagonals in the LHS matrix */
    int          nPoints;       /* the number of points that stnHd contains */
    int          nFluidPoints;  /* the number of fluid points that stnHd contains */
    int          nWallPoints;   /* the number of wall points that stnHd contains */
    double       n0;            /* the average density number */
    double       radius;        /* the radius of influence */
    double       beta;          /* the free surface constant */
    int          d;             /* the number of dimensions of the problem */
} Stn;

typedef Stn* StnHd;

/*******************************************************************************
 * Abstract functions
 *******************************************************************************
 */

StnHd stnNew(int nFluidPoints, int nWallPoints);
void  stnPopulate(StnHd  stnHd,   double *xCrd, double *yCrd,
		  double radius,  double  beta);
void  stnRecalc(StnHd stnHd, double *xCrd, double *yCrd);
void  stnFree(StnHd stnHd);
