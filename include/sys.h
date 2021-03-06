/*******************************************************************************
 ** Copyright 2014-2014 Vedaad Shakib Inc.
 *******************************************************************************/

/*******************************************************************************
 ** 
 ** "sys.h": standard system level include
 **
 *******************************************************************************/

/*******************************************************************************
 * System includes
 *******************************************************************************
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <assert.h>

/*******************************************************************************
 * Standard includes
 *******************************************************************************
 */
#include "mem.h"

/*******************************************************************************
 * Definitions and typedefs
 *******************************************************************************
 */
#define SYS_H_INCLUDED
#ifndef true
#define true (1)
#endif

#ifndef false
#define false (0)
#endif

#define	MAX(A,B)	((A)>(B)?(A):(B))
#define	MIN(A,B)	((A)<(B)?(A):(B))

typedef int bool;

#define debugArr(NAME,SIZE,NDIMS) printf("Debug array %s: \n", #NAME);        \
                                  for (int u = 0; u < SIZE; u++) {            \
                                      for (int w = 0; w < NDIMS; w++)         \
                                          printf("%lf ", NAME[u*NDIMS+w]);    \
                                      printf("\n");                           \
                                  }                                  


