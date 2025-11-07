#include "prm_def.h"
#include "matrix_operations.h"
#include <stdlib.h>
#include <math.h>

// DEBUG
// #include <stdio.h>
// DEBUG


/* Calculate linear coefficients k and delta for xp position scaling,
 * given real positions yp. The fitting is calculated within the roi.
 * Returns a struct with k and delta.
 */
kdelta positions_scaling (double * xp, double * yp, roi_struct * roi)
{
    kdelta kd;

    double sx2 = roi_dot_product(xp, xp, roi);
    double sx  = roi_vector_sum(xp, roi);
    double det = (double)roi->nsites * sx2 - sx * sx;

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

    kd.delta = (sx2 * sy - sx * sxy) / det;
    kd.k = (sy - (double)roi->nsites * kd.delta) / sx;

    // DEBUG
    // printf("       k = %lf\n", kd.k);
    // printf("   delta = %lf\n", kd.delta);
    // DEBUG

    return kd;
}


/* Calculate vertical positions vv based on pairwise formula.
 * The gain matrix gain_v is applied to the blades' values
 * provided by dataset ds.
 */
void raw_positions_calc (dataset * ds, double * supmat,
                         double * pos_h, double * pos_v)
{
    /* Pointer to blades' values at each site.  */
    double pblade[4];
    /* Delta over sigma vector.  */
    double deltasigma[4] = {0.0, 0.0, 0.0, 0.0};

    for (size_t ii = 0; ii < ds->nsites; ii++)
    {
        pblade[0] = ds->to[ii];
        pblade[1] = ds->ti[ii];
        pblade[2] = ds->bi[ii];
        pblade[3] = ds->bo[ii];

        matrix_vector_product(supmat + 8, (double *)pblade, 4, 4, deltasigma);
        pos_h[ii] = (deltasigma[0] / deltasigma[1]);
        pos_v[ii] = (deltasigma[2] / deltasigma[3]);

        // p0 = gain[0] * ds->to[ii];
        // p1 = gain[1] * ds->ti[ii];
        // p2 = gain[2] * ds->bi[ii];
        // p3 = gain[3] * ds->bo[ii];
        // pos[ii] = ((p0 + p3 - p1 - p2) /
        //            (fabs(p0) + fabs(p1) + fabs(p2) + fabs(p3)));

    }
}


/* Calculate horizontal positions hh based on pairwise formula.
 * The gain matrix gain_h is applied to the blades' values
 * provided by dataset ds.
 */
void raw_positions_calc_h (dataset * ds, double * gain, double * pos)
{
    double p1, p2;
    for (size_t ii = 0; ii < ds->nsites; ii++)
    {
        p1 = fabs(gain[0]) * ds->to[ii] + fabs(gain[3]) * ds->bo[ii];
        p2 = fabs(gain[1]) * ds->ti[ii] + fabs(gain[2]) * ds->bi[ii];
        pos[ii] = (p1 - p2) / (p1 + p2);

        // p0 = gain[0] * ds->to[ii];
        // p1 = gain[1] * ds->ti[ii];
        // p2 = gain[2] * ds->bi[ii];
        // p3 = gain[3] * ds->bo[ii];
        // pos[ii] = ((p0 + p3 - p1 - p2) /
        //            (fabs(p0) + fabs(p1) + fabs(p2) + fabs(p3)));

    }
}


void positions_calc_h (dataset * ds, double * gain, double * pos)
{
    /* Calculate positions according to pairwise formula. */
    raw_positions_calc_h(ds, gain, pos);

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

    if (isnan(kd.k) || isnan(kd.delta))
        return;


    for (size_t ii = 0; ii < ds->nsites; ii++)
    {
        pos[ii] = pos[ii] * kd.k - kd.delta;
    }
}


void positions_calc_v (dataset * ds, double * gain, double * pos)
{
    /* Calculate positions according to pairwise formula. */
    raw_positions_calc_v(ds, gain, pos);

    /* Scale positions. */
    kdelta kd = positions_scaling(pos, ds->nom_v, &(ds->roi));
    for (size_t ii = 0; ii < ds->nsites; ii++)
    {
        pos[ii] = pos[ii] * kd.k - kd.delta;
    }
}
