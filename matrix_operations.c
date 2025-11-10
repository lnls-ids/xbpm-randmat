#include "prm_def.h"
#include <stdlib.h>

#ifndef MAT_OP
#define MAT_OP


/* Return the transpose matT of a mm x nn matrix mat.
 * The matrix is provided as a flat row-major array of size mm * nn:
 * element (i,j) is mat[i * nn + j]. The result is returned as a
 * newly allocated flat row-major array of size nn * mm such that
 * element (j,i) is matT[j * mm + i]. Caller must free the returned
 * pointer with free().
 */
double *matrix_transpose(const double *mat, size_t mm, size_t nn)
{
    if (mm <= 0 || nn <= 0 || mat == NULL) return NULL;
    size_t totsz = mm * nn;
    double *matT = malloc(totsz * sizeof(double));
    if (!matT) return NULL;
    for (size_t i = 0; i < mm; ++i) {
        for (size_t j = 0; j < nn; ++j) {
            /* mat[i,j] -> matT[j,i] */
            matT[j * mm + i] = mat[i * nn + j];
        }
    }
    return matT;
}


/* Calculate the product pr of an (mm x nn) matrix mA by an
 * (nn x pp) matrix mB. All matrices are flat row-major arrays:
 * mA: mm x nn (index i*nn + k)
 * mB: nn x pp (index k*pp + j)
 * prod: mm x pp (index i*pp + j) and must be allocated by caller.
 * prod can be pre-zeroed or the function will accumulate into it.
 */
double *matrix_product(const double *mA, const double *mB,
                       size_t mm, size_t nn, size_t pp, double *prod)
{
    if (!mA || !mB || !prod || mm == 0 || nn == 0 || pp == 0) return prod;

    for (size_t i = 0; i < mm; ++i) {
        for (size_t j = 0; j < pp; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < nn; ++k) {
                sum += mA[i * nn + k] * mB[k * pp + j];
            }
            prod[i * pp + j] = sum;
        }
    }
    return prod;
}


/* Multiply mm x nn matrix mA (flat row-major) by vector vv (size nn).
 * prod must be mm-sized and will be filled with results.
 */
void matrix_vector_product(const double *mA, const double *vv,
                           size_t mm, size_t nn, double *prod)
{
    if (!mA || !vv || !prod || mm == 0 || nn == 0) return;

    for (size_t ii = 0; ii < mm; ++ii) {
        double sum = 0.0;
        for (size_t kk = 0; kk < nn; ++kk) {
            sum += mA[ii * nn + kk] * vv[kk];
        }
        prod[ii] = sum;
    }
}


/* Calculate the dot product of nn-size line matrix mA and
 * column matrix mB.
 */
double dot_product(const double *mA, const double *mB, size_t nn)
{
    double dprod = 0.0;
    for (size_t ii = 0; ii < nn; ii++)
    {
        dprod += mA[ii] * mB[ii];
    }
    return dprod;
}


/* Add up the elements of an nn sized vector mA (column matrix).
*/
double vector_sum(const double *mA, size_t nn)
{
    double sum = 0.0;
    for (size_t ii = 0; ii < nn; ii++)
    {
        sum += mA[ii];
    }
    return sum;
}


/* Calculate the dot product of a line matrix mA and a
 * column matrix mB, indexed by roi->idx. The ROI skips
 * certain elements.
 */
double roi_dot_product(const double *mA, const double *mB,
                       const roi_struct * roi)
{
    double dprod = 0.0;
    for (size_t ii = 0; ii < roi->nsites; ii++)
    {
        size_t idx = roi->idx[ii];
        dprod += mA[idx] * mB[idx];
    }
    return dprod;
}


/* Add up the elements of a vector mA (column matrix)
 * indexed by roi->idx. The ROI skips
 * certain elements.
 */
double roi_vector_sum(const double *mA, const roi_struct * roi)
{
    double sum = 0.0;
    for (size_t ii = 0; ii < roi->nsites; ii++)
    {
        sum += mA[roi->idx[ii]];
    }
    return sum;
}

#endif
