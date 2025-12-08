#include "prm_def.h"
#include "pcg_random.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Temperature constant. Analogous to 1/kB. */
#define Bk  1.0e7

/* Prototype. Calculate positions from suppression matrix and
 * blades' measurements.
 */
kdelta positions_calc (const dataset * ds, const double * supmat,
                    double * nom_positions, double * positions);


/* Change the value of an element of the gain array by a 'step'.
 */
void mat_walk (double * supmat, double step)
{
    /* Pick an element.
    */
    int isite   = (int) (pcg_double() * 16);
    /* Choose sign +/- (increase/decrease step).
    */
    double sign = 1.0 ? -1.0 : (pcg_double() > 0.5);
    /* Add up in chosen matrix element value.
    */
    supmat[isite] += sign * step;
}


/* Chi-square of two vectors, v1 and v2, indexed by roi->idx.
 */
double chi2_calc(const double * v1, const double * v2,
                 const roi_struct * roi)
{
    double df, c2 = 0.0;
    size_t idx;
    for (size_t ii = 0; ii < roi->nsites; ii++)
    {
        idx = roi->idx[ii];
        df  = v1[idx] - v2[idx];
        c2 += df *  df;
    }
    if (roi->nsites <= 1) return 0.0;
    return c2 / ((double)(roi->nsites - 1));
}


/* Initialize seed with real random value from urandom device.
 */
uint64_t seed_get()
{
    /* Open urandom device.
    */
    FILE * sd = fopen("/dev/urandom", "rb");
    size_t nread;
    uint64_t buffer;

    /* Read 8 bytes from urandom into buffer. */
    nread = fread(&buffer, 8, 1, sd);
    fclose(sd);

    /* Warn if nothing read. */
    if (nread == 0)
    {
        printf("\n WARNING : nothing read from urandom.\n");
    }
    return buffer;
}


/* Perform random walk to optimize suppression matrix.
 */
rw_stats random_walk(dataset * ds, xbpm_prm * prm, double * supmat,
                   double * pos_h, double * pos_v)
{
    /* Counters. */
    size_t ii = 0;
    size_t accept = 0;

    /* Inverse of temperature. */
    double beta = prm->beta;
    
    /* Chi2 analysis.*/
    double oldval;
    double chi2_h, chi2_v, chi2_h_aft, chi2_v_aft;
    double chi2, chi2_aft, dchi2;

    /* Probability.*/
    double prob;

    /* Scaling parameters. */
    kdelta kd;
    
    /* Initialize random seed. */
    uint64_t seed = seed_get();
    pcg32_init(seed);
    
    /* Calculate initial positions and deviation from nominal
     * positions (chi2). */
    positions_calc(ds, supmat, ds->nom_h, pos_h);
    chi2_h = chi2_calc(ds->nom_h, pos_h, &ds->roi);
    
    positions_calc(ds, supmat + 8, ds->nom_v, pos_v);
    chi2_v = chi2_calc(ds->nom_v, pos_v, &ds->roi);
    
    chi2_h_aft = chi2_h;
    chi2_v_aft = chi2_v;
    chi2       = (chi2_h + chi2_v);
    chi2_aft   = chi2;

    /* Try to change the matrix nrand times. */
    size_t imat_h = 0;
    size_t imat_v = 0;
    size_t isite = 0;
    double sign = 1.0;
    for (ii = 0; ii < prm->nrand; ii++)
    {
        /* Pick an element of the suppression matrix. */
        isite = (size_t) (pcg_double() * 16);

        /* Choose sign (increase/decrease step).      */
        sign = (pcg_double() > 0.5) ? -1.0 : 1.0;
        /* Add up in chosen matrix element value.     */
        oldval = supmat[isite];
        supmat[isite] += sign * prm->step;

        /* Skip if new value would be zero. */
        if (supmat[isite] == 0.0) 
        {
            supmat[isite] = oldval;
            continue;
        }

        /* Recalculate positions and chi2. Check h and v separately.
         * Minimization step takes only ROI into account. */ 
        if (isite < 8)
        {
            /* Horizontal changes. */
            kd = positions_calc(ds, supmat, ds->nom_h, pos_h);
            chi2_h_aft = chi2_calc(ds->nom_h, pos_h, &ds->roi);            
            imat_h++;
        }
        else
        {
            /* Vertical changes. */
            kd = positions_calc(ds, supmat + 8, ds->nom_v, pos_v);
            chi2_v_aft = chi2_calc(ds->nom_v, pos_v, &ds->roi);
            imat_v++;
        }
        
        /* If the scaling failed, reject change. */
        if (kd.k == 1.0) 
        {
            printf("\n");
            continue;
        }
           
        /* Calculate the change in chi2. */
        chi2_aft = (chi2_h_aft + chi2_v_aft);
        dchi2    = chi2_aft - chi2;

        /* Probability of acceptance. */
        prob = exp(-dchi2 * beta * Bk);

        /* probability = min(1, prob). */
        if (prob > 1.0) prob = 1.0;

        /* Accept or reject change. */
        if (pcg_double() <= prob)
        {
            /* If change is accomplished, copy state values 
             * to former variables and reduce temperature. */
            chi2 = chi2_aft;
            chi2_h = chi2_h_aft;
            chi2_v = chi2_v_aft;
            accept++;
            // beta *= 0.999;
            beta *= 1.001;
        }
        else
        {
            /* If change is rejected, restore former value. */
            supmat[isite] = oldval;
            chi2_aft = chi2;
            chi2_h_aft = chi2_h;
            chi2_v_aft = chi2_v;
            beta *= 0.999;
        }
    }

    /* Final positions. */
    positions_calc(ds, supmat,     ds->nom_h, pos_h);
    positions_calc(ds, supmat + 8, ds->nom_v, pos_v);

    // DEBUG : print ROI
    /*
    printf("\n\n##### (RANDOM WALK) ROI: #####\n");
    for (size_t ii = 0; ii < ds->roi.nsites; ii++)
    {
        printf(" %10.4f %10.4f",
               ds->nom_h[ds->roi.idx[ii]], ds->nom_v[ds->roi.idx[ii]]);
        printf(" %10.4f %10.4f\n",
               pos_h[ds->roi.idx[ii]], pos_v[ds->roi.idx[ii]]);
    }
    printf("\n");
    */
    // DEBUG

    return (rw_stats) {imat_h, imat_v, accept, beta};
}