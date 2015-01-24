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
    double	*row;		/* the points adjacent to each point in turn */
    int		 nPoints;	/* the number of points the stn holds */
    int          maxPoints;     /* the maximum capacity of the col array */
    int          maxAdjacent;   /* the maximum capacity of the row array */
} Stn;

typedef Stn* StnHd;

/*******************************************************************************
 * Abstract functions
 *******************************************************************************
 */

StnHd	stnNew();
void	stnPopulate(StnHd stnHd, double *points, int nPoints, double radius);
