// #include "prm_def.h"
#include "matrix_operations.h"
#include <stdlib.h>
#include <math.h>

// DEBUG
// #include <stdio.h>
// DEBUG


/* Calculate linear coefficients k and delta for xp position scaling
 * by least squares method, given real (nominal) positions yp.
 * The fitting is calculated within the roi.
 * Returns a struct with k and delta.
 */
kdelta positions_scaling (const double * xp, const double * yp,
                          const roi_struct * roi)
{
    kdelta kd;

    double sxx = roi_dot_product(xp, xp, roi);
    double sx  = roi_vector_sum(xp, roi);
    double det = (double)roi->nsites * sxx - sx * sx;

    double sy  = roi_vector_sum(yp, roi);
    double sxy = roi_dot_product(xp, yp, roi);

    // DEBUG
    // printf(" (position scaling) : nsites = %zu\n", roi->nsites);
    // for (size_t ii = 0; ii < roi->nsites; ii++)
    //     printf(">>> [%zu : %zu] xp = %lf, yp = %lf\n", ii, roi->idx[ii],
    //             xp[roi->idx[ii]], yp[roi->idx[ii]]);

    // printf(" (position scaling) :\n");
    // printf("  yp[15] = %lf\n", yp[15]);
    // printf("      sx = %lf\n", sx);
    // printf("     sx2 = %lf\n", sx2);
    // printf("      sy = %lf\n", sy);
    // printf("     sxy = %lf\n", sxy);
    // printf("     det = %lf\n", det);
    // DEBUG

    kd.delta = (sxx * sy - sx * sxy) / det;
    kd.k = (sy - (double)roi->nsites * kd.delta) / sx;

    // DEBUG
    // printf("       k = %lf\n", kd.k);
    // printf("   delta = %lf\n", kd.delta);
    // DEBUG

    return kd;
}


/* Calculate positions multiplying by coefficients of suppression matrix.
 */
void raw_positions_calc (const dataset * ds, const double * supmat,
                         double * pos)
{
    double delta, sigma;
    for (size_t ii = 0; ii < ds->nsites; ii++)
    {
        delta = supmat[0] * ds->to[ii]
              + supmat[1] * ds->ti[ii] 
              + supmat[2] * ds->bi[ii]
              + supmat[3] * ds->bo[ii]; 

        sigma = supmat[4] * ds->to[ii]
              + supmat[5] * ds->ti[ii] 
              + supmat[6] * ds->bi[ii]
              + supmat[7] * ds->bo[ii]; 

        pos[ii] = (delta / sigma);
    }
}


/* Calculate positions multiplying by coefficients of suppression matrix.
 * Then scale calculated positions based on nominal positions.
 */
kdelta positions_calc (const dataset * ds, const double * supmat,
                       double * pos)
{
    /* Calculate positions according to suppression matrix. */
    raw_positions_calc(ds, supmat, pos);

    // DEBUG
    // printf(" (pos calc h) RAW: pos H[0], H[5], H[10] = %lf / %lf / %lf\n\n",
    //     pos[0], pos[5], pos[10]);
    // DEBUG

    /* Scale positions. */
    kdelta kd = positions_scaling(pos, ds->nom_h, &(ds->roi));

    // DEBUG
    // printf(" (pos calc h) scaling: k, delta = %lf , %lf \n\n",
    //     kd.k, kd.delta);
    // DEBUG

    /* If scaling is not successful. */
    if (isnan(kd.k) || isnan(kd.delta))
    {
        kd.k = 1.0;
        kd.delta = 0.0;
        return kd;
    }

    for (size_t ii = 0; ii < ds->nsites; ii++)
    {
        pos[ii] = pos[ii] * kd.k + kd.delta;
    }
    return kd;
}

