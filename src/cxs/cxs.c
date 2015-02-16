/*******************************************************************************
 ** Copyright 2014-2014 Vedaad Shakib Inc.
 *******************************************************************************/

/*******************************************************************************
 ** 
 ** "cxs.c": Interface for CXSparse
 **
 *******************************************************************************/

#include "cs.h"


/*******************************************************************************
 ** 
 ** "cxsSolve": solve a symmetric, positive-definite Ax=b using CXSparse routines
 **           
 ** Parameters:
 **            col: Adjacency col[nDims+1] (col array of CCS format matrix)
 **            row: Adjacency row[nNonzeros] (row array of CCS format matrix)
 **            mtx: Adjacency mtx[nNonzeros] (value array of CCS format matrix)
 **            b: RHS vector[nDims]
 **            x: solution vector[nDims]
 **            nDims: dimension of problem
 **            nNonzeros: the number of nonzeros
 **
 ******************************************************************************
 */
int cxsSolveSym(int *col, int *row, double *mtx, double *b,
		double *x, int nDims, int nNonzeros) {
    double     *y;              /* pivoted vector                */
    cs          A;              /* matrix in CS structure        */
    css        *S;              /* symbolic factorized matrix    */
    csn        *N;              /* numerical factorized matrix   */

/*------------------------------------------------------------------------------
 * Build the CS matrix from input values
 *------------------------------------------------------------------------------
 */
    A.nzmax = nNonzeros;        /* maximum number of entries        */
    A.n     = nDims;            /* number of columns                */
    A.m     = nDims;            /* number of rows                */
    A.p     = col;              /* column pointers (size n+1)        */
    A.i     = row;              /* row indices, size nzmax        */
    A.x     = mtx;              /* numerical values, size nzmax        */
    A.nz    = -1;               /* -1 for compressed-col        */

/*------------------------------------------------------------------------------
 * Perform symbolic factorization
 *------------------------------------------------------------------------------
 */
    S = cs_schol(1, &A);

    if (S == NULL) {
        return(0);
    }
/*------------------------------------------------------------------------------
 * Perform numerical (Cholesky) factorization
 *------------------------------------------------------------------------------
 */
    N = cs_chol(&A, S);                /* Cholesky factorization        */

    if (N == NULL) {
        cs_sfree(S);
        return(0);
    }
/*------------------------------------------------------------------------------
 * Allocate working memory
 *------------------------------------------------------------------------------
 */
    y = cs_malloc(nDims, sizeof(double));

    if (y == NULL) {
        cs_sfree(S);
        cs_nfree(N);
        return(0);
    }
/*------------------------------------------------------------------------------
 * Solve the matrix
 *------------------------------------------------------------------------------
 */
    cs_ipvec(S->pinv, b, y, nDims);           /* y = P*b                 */
    cs_lsolve(N->L, y);                       /* y = L\y                 */
    cs_ltsolve(N->L, y);                      /* y = L'\y                */
    cs_pvec(S->pinv, y, x, nDims);            /* b = P'*y                */

/*------------------------------------------------------------------------------
 * Free memory
 *------------------------------------------------------------------------------
 */
    cs_sfree(S);
    cs_nfree(N);
    cs_free(y);

/*------------------------------------------------------------------------------
 * Return
 *------------------------------------------------------------------------------
 */
    return 1;

} 
