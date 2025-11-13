// #include "prm_def.h"
#include "matrix_operations.h"
#include <stdlib.h>
#include <math.h>

/* Calculate positions (pos) by multiplying blades' measurements in dataset
 * ds (to, ti, bi. bo) by the elements of the suppression matrix. Horizontal
 * positions are calculated from the 8 first elements of supmat, vertical
 * positions from the last 8 elements; the pointer to the first of those must
 * be sent in each case.
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


/* Calculate linear coefficients k and delta for xp position scaling
 * by least squares method, given real (nominal) positions yp.
 * The fitting is calculated within the roi.
 * Returns a struct with k and delta.
 */ 
kdelta positions_scaling (const double * xp, const double * yp,
                          const roi_struct * roi)
{                      
    kdelta kd;
    double nsites = (double) roi->nsites;

    double sxx = roi_dot_product(xp, xp, roi);
    double sx  = roi_vector_sum(xp, roi);
    double det = nsites * sxx - sx * sx;

    double sy  = roi_vector_sum(yp, roi);
    double sxy = roi_dot_product(xp, yp, roi);

    kd.delta = (sxx * sy - sx * sxy)    / det;
    kd.k     = (sy - nsites * kd.delta) / sx;

    return kd;
}    


/* Calculate positions multiplying by coefficients of suppression matrix.
 * Then scale calculated positions based on nominal positions.
 */
kdelta positions_calc (const dataset * ds, const double * supmat,
                       const double * nompos, double * pos)
{

    
    /* Calculate positions (whole grid) according to suppression matrix. */
    raw_positions_calc(ds, supmat, pos);

    /* Scale positions. Only the ROI is relevant for scaling. */
    kdelta kd = positions_scaling(pos, nompos, &(ds->roi));

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

