/* Header for matrix operations used in xbpm-randmat.
 * Implementations live in matrix_operations.c
 */
#ifndef MAT_OP
#define MAT_OP

#include "prm_def.h"

/* Transpose a mm x nn matrix provided as a flat row-major array.
 * Returns a newly allocated flat array of size nn * mm (caller frees).
 */
double *matrix_transpose(const double *mat, size_t mm, size_t nn);

/* Matrix multiply: mA (mm x nn) * mB (nn x pp) -> prod (mm x pp).
 * All arrays are flat row-major. prod must be allocated by caller.
 */
double *matrix_product(const double *mA, const double *mB,
                       size_t mm, size_t nn, size_t pp, double *prod);

/* Matrix-vector product: mA (mm x nn) * vv (nn) -> prod (mm).
 */
void matrix_vector_product(const double *mA, const double *vv,
                           size_t mm, size_t nn, double *prod);

/* Basic vector helpers. */
double dot_product(const double *mA, const double *mB, size_t nn);
double vector_sum(const double *mA, size_t nn);

/* ROI-aware helpers (operate on flat vectors indexed by roi->idx). */
double roi_dot_product(const double *mA, const double *mB, roi_struct *roi);
double roi_vector_sum(const double *mA, roi_struct *roi);

#endif
