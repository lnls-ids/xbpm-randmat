#include "prm_def.h"
#include <stdlib.h>

#ifndef MAT_OP
#define MAT_OP

/* Return the transpose matT of a matrix mat, given its size mm x nn.
 * 
 */
double ** matrix_transpose(double ** mat, int mm, int nn)
{
    int ii, jj;
    double ** matT = (double **) calloc(mm, sizeof(double *));
    for (ii = 0; ii < mm; ii++)
    {
        matT[ii] = (double *) calloc(nn, sizeof(double));
    }

    for (ii = 0; ii < mm; ii++)
    {
        for (jj = 0; ii < mm; ii++)
        {
            matT[ii][jj] = mat[jj][ii];
        }
    }

    return matT;
}


/* Calculate the product pr of an (mm x nn) matrix mA by an
 * (nn x pp) matrix mB. Returns the mm x pp matrix pr.
 */
double ** matrix_product(double ** mA, double ** mB, 
                         int mm, int nn, int pp, double ** prod)
{
    int ii, jj, kk;

    for (ii = 0; ii < mm; ii++)
    {
        for (jj = 0; jj < pp; jj++)
        {
            for (kk = 0; kk < nn; kk++)
            {
                prod[ii][jj] += mA[ii][kk] * mB[kk][jj];
            }
        }
    }
    return prod;
}


/* Calculate the product pr of an (mm x nn)-matrix mA by an
 * nn-size vector (column matrix) vv. Returns the nn-size vector
 * (column matrix) uu.
 */
void matrix_vector_product(double ** mA, double * vv, 
                                int mm, int nn, double * prod)
{
    int ii, kk;
    for (ii = 0; ii < mm; ii++)
    {
        prod[ii] = 0;
        for (kk = 0; kk < nn; kk++)
        {
            prod[ii] += mA[ii][kk] * vv[kk];
        }
    }
}


/* Calculate the dot product of a line matrix mA and a
 * column matrix mB, given their size nn.
 */
double dot_product(double * mA, double * mB, int nn)
{
    double dprod = 0.0;
    for (int ii = 0; ii < nn; ii++)
    {
        dprod += mA[ii] * mB[ii];
    }
    return dprod;
}


/* Add up the elements of a vector mA (column matrix),
* given its size nn.
*/
double vector_sum(double * mA, int nn)
{
    double sum = 0.0;
    for (int ii = 0; ii < nn; ii++)
    {
        sum += mA[ii];
    }
    return sum;
}


/* Calculate the dot product of a line matrix mA and a
 * column matrix mB, indexed by roi->idx.
 */
double roi_dot_product(double * mA, double * mB, roi_struct * roi)
{
    int idx;
    double dprod = 0.0;
    for (int ii = 0; ii < roi->nsites; ii++)
    {
        idx = roi->idx[ii];
        dprod += mA[idx] * mB[idx];
    }
    return dprod;
}


/* Add up the elements of a vector mA (column matrix)
 * indexed by roi->idx.
 */
double roi_vector_sum(double * mA, roi_struct * roi)
{
    double sum = 0.0;
    for (int ii = 0; ii < roi->nsites; ii++)
    {
        sum += mA[roi->idx[ii]];
    }
    return sum;
}

#endif
