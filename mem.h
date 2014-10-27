/*******************************************************************************
** Copyright 2014-2014 Vedaad Shakib Inc.
*******************************************************************************/

/*******************************************************************************
** 
** "mem.h": memory management routines
**
*******************************************************************************/

/*******************************************************************************
 * Memory allocation macros
 *******************************************************************************
 */
#define	memNew(T,N)	(T*) malloc(sizeof(T) * N)

#define	memNewZero(T,N)	(T*) memset(malloc(sizeof(T)*N),0,sizeof(T)*N)

#define	memResize(TYPE, ARRAY, OLDSIZE, NEWSIZE, MAXSIZE, NDIMS)	\
  if ((NEWSIZE) >= (MAXSIZE)) {						\
    TYPE *ary;								\
    (MAXSIZE) = 2 * (MAXSIZE) + 100;					\
    if ((MAXSIZE) < (NEWSIZE)) (MAXSIZE) = (NEWSIZE);			\
    ary = memNew(TYPE, (MAXSIZE)*(NDIMS));				\
    if ((OLDSIZE) != 0) {						\
      memcpy(ary, ARRAY, sizeof(TYPE)*(OLDSIZE)*(NDIMS));		\
      free(ARRAY);							\
    }									\
    ARRAY = ary;							\
  }

