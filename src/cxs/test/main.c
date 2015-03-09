/*******************************************************************************
 ** Copyright 2014-2014 Vedaad Shakib Inc.
 *******************************************************************************/

/*******************************************************************************
 ** 
 ** "main.c": test cxs
 **
 *******************************************************************************/

#include "cxs.h"
#include "sys.h"

/*******************************************************************************
 ** 
 ** "main": main function
 **
 ******************************************************************************
 */
int main(void) {
    int         nDims;          /* size of the matrix           */
    int         nNonzeros;      /* No. non-zeros                */
    int        *col;            /* Adjacency col                */
    int        *row;            /* Adjacency row                */
    double     *mtx;            /* Adjacency matrix             */
    double     *b;              /* problem RHS                  */
    double     *x;              /* problem solution             */
    int         i;              /* a running index              */
    int         j;              /* a running index              */
    int         k;              /* a running index              */
    int         errCode;        /* return code                  */
    double      tmp;            /* a temporary real             */
    double      aTmp;           /* a temporary real             */
    double      bTmp;           /* a temporary real             */

/*------------------------------------------------------------------------------
 * Initialize
 *------------------------------------------------------------------------------
 */
    nDims     = 10000;
    nNonzeros = 3 * nDims;

    col = memNew(int, nDims+1);
    row = memNew(int, nNonzeros);
    mtx = memNew(double, nNonzeros);
    b   = memNew(double, nDims);
    x   = memNew(double, nDims);

/*------------------------------------------------------------------------------
 * Build a periodic matrix adjacency
 *------------------------------------------------------------------------------
 */
    col[0] = 0;
    for (i = 0; i < nDims; i++) {
        col[i+1] = col[i] + 3;
    }

    j        = 0;
    row[j++] = 0;
    row[j++] = 1;
    row[j++] = nDims - 1;

    for (i = 1; i < nDims-1; i++) {
        row[j++] = i - 1;
        row[j++] = i;
        row[j++] = i + 1;
    }
    row[j++] = 0;
    row[j++] = nDims - 2;
    row[j++] = nDims - 1;

/*------------------------------------------------------------------------------
 * Set the matrix to 2 in diagonal and -1 for off diagonals 
 *------------------------------------------------------------------------------
 */
    for (i = 0; i < nDims; i++) {
        for (k = col[i]; k < col[i+1]; k++) {
            j = row[k];
            if (j == i) {
                mtx[k] = 2;
            } else {
                mtx[k] = -1;
            }
        }
    }
    mtx[0] = 3; // to avoid singularity

/*------------------------------------------------------------------------------
 * Build a vector
 *------------------------------------------------------------------------------
 */
    tmp = 0;
    for (i = 0; i < nDims; i++) {
        b[i]  = 1 + i;
        tmp  += b[i];
    }
    tmp /= nDims;
    for (i = 0; i < nDims; i++) {
        b[i] -= tmp;
    }
/*------------------------------------------------------------------------------
 * Solve the problem
 *------------------------------------------------------------------------------
 */
    errCode = cxsSolveSym(col, row, mtx, b, x, nDims, nNonzeros);

    if (errCode) {
        printf("Error %d from cxsSolveSym()\n",errCode);
    }
/*------------------------------------------------------------------------------
 * Check the solution
 *------------------------------------------------------------------------------
 */
    aTmp = 0;
    bTmp = 0;
    for (i = 0; i < nDims; i++) {
        aTmp += b[i] * b[i];
        tmp   = -b[i];
        for (k = col[i]; k < col[i+1]; k++) {
            j    = row[k];
            tmp += mtx[k] * x[j];
        }
        bTmp += tmp * tmp;
    }

    printf("Computational Error = %g\n", sqrt(bTmp / aTmp)); 

}
