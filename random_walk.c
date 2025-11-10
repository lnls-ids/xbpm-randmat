#include "prm_def.h"
#include "pcg_random.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Prototype. Calculate positions from suppression matrix and
 * blades' measurements.
 */
int positions_calc (const dataset * ds, const double * supmat,
                    double * positions);


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
double chi2_calc(double * v1, double * v2, roi_struct * roi)
{
    double df, c2 = 0.0;
    for (size_t ii = 0; ii < roi->nsites; ii++)
    {
        size_t idx = roi->idx[ii];
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

    /* Read 8 bytes from urandom into buffer.
    */
    nread = fread(&buffer, 8, 1, sd);
    fclose(sd);

    /* Warn if nothing read.
    */
    if (nread == 0)
    {
        printf("\n WARNING : nothing read from urandom.\n");
    }
    return buffer;
}


/* Perform random walk to optimize suppression matrix.
 */
size_t random_walk(dataset * ds, xbpm_prm * prm, double * supmat,
                   double * pos_h, double * pos_v)
{
    /* Counters. */
    size_t ii = 0;
    size_t accept = 0;
    
    /* Inverse of temperature. */
    double beta = 1.0 / prm->temp;
    
    /* Chi2 analysis.*/
    double newval;
    double chi2_h, chi2_v, chi2_h_aft, chi2_v_aft;
    double chi2, chi2_aft, dchi2;

    /* Probability.*/
    double prob;
    
    /* Initialize random seed. */
    uint64_t seed = seed_get();
    pcg32_init(seed);
    
    /* Calculate initial positions and deviation from nominal
     * positions (chi2). */
    positions_calc(ds, supmat, pos_h);
    chi2_h     = chi2_calc(ds->nom_h, pos_h, &ds->roi);
    
    positions_calc(ds, supmat, pos_v);
    chi2_v     = chi2_calc(ds->nom_v, pos_v, &ds->roi);
    
    chi2       = (chi2_h + chi2_v) * 0.5;
    chi2_aft   = chi2;
    chi2_h_aft = chi2_h;
    chi2_v_aft = chi2_v;
    
    // DEBUG
    /*
    printf("\n##### (random walk) BEFORE #####\n");
    printf("\n\n >> gain << \n (random walk) H : \n");
    for (int ii = 0; ii < 4 ; ii++)
    printf(" %10.4f  ", gain_h[ii]);
    printf("\n\n (random walk) V : \n ");
    for (int ii = 0; ii < 4 ; ii++)
    printf(" %10.4f  ", gain_v[ii]);
    printf("\n\n");
    printf(" (random walk) H, V = %lf / %lf\n\n", pos_h[5], pos_v[5]);
    */
    // DEBUG
   
    /* Temporary copy of suppression matrix. */
    double * supmat_tmp = calloc(16, sizeof(double));
    memcpy(supmat_tmp, supmat, 16 * sizeof(double));

    /* Try to change the matrix nrand times. */
    for (ii = 0; ii < prm->nrand; ii++)
    {
        /* Pick an element of the suppression matrix. */
        int isite   = (int) (pcg_double() * 16);
        /* Choose sign (increase/decrease step).      */
        double sign = 1.0 ? -1.0 : (pcg_double() > 0.5);
        /* Add up in chosen matrix element value.     */
        newval = supmat_tmp[isite] + sign * prm->step;

        // Skip if new value would be zero.
        if (newval == 0.0) continue;

        supmat_tmp[isite] = newval;

        /* Recalculate positions and chi2. Check h and v separately.
         * Minimization step takes only ROI into account. */ 
        if (isite < 8)
        {
            positions_calc(ds, supmat_tmp, pos_h);
            chi2_h_aft = chi2_calc(ds->nom_h, pos_h, &ds->roi);
        }
        else
        {
            positions_calc(ds, supmat_tmp + 8, pos_v);
            chi2_v_aft = chi2_calc(ds->nom_v, pos_v, &ds->roi);
        }
        
        // DEBUG
        /*
        printf(" pos 0: H, V = %lf / %lf\n ", pos_h[0], pos_v[0]);
        for (int jj = 0; jj < 4; jj++)
        {
            printf("h = %lf ", gain_h_tmp[jj]);
            printf("v = %lf ", gain_v_tmp[jj]);
        }

        printf(" (random walk) chi2 = %lf , %lf\n", chi2_h_aft, chi2_v_aft);
        */
        // DEBUG
           
        /* Calculate the change in chi2. */
        chi2_aft = (chi2_h_aft + chi2_v_aft) * 0.5;
        dchi2    = chi2_aft - chi2;
        /* Probability of acceptance. */
        prob     = exp(-dchi2 * beta);
        /* probability = min(1, prob). */
        if (prob > 1.0) prob = 1.0;

        // DEBUG
        /*
        printf(" (random walk) AFT >> dchi2 = %lf - %lf = %lf ←→ %lf\n\n",
        chi2_aft, chi2_bef, dchi2, exp(-dchi2 / prm->temp));
        prob = (dchi2 < 0.0) ? 1.0 : exp(-dchi2 / prm->temp);
        printf(" (random walk) >> prob = %lf ←→ %Lg\n", prob, pcg_double());
        */
        // DEBUG

        /* Accept or reject change. */
        if (pcg_double() <= prob)
        {
            /* If change is accomplished, copy state values 
             * to former variables and reduce temperature.
             */
            supmat[isite] = newval;
            prm->temp *= 0.999;
            chi2 = chi2_aft;
            chi2_h = chi2_h_aft;
            chi2_v = chi2_v_aft;
            accept++;
        }
        else
        {
            /* If change is rejected, restore former value. */
            supmat_tmp[isite] = supmat[isite];
            chi2_aft = chi2;
            chi2_h_aft = chi2_h;
            chi2_v_aft = chi2_v;
        }
    }

    // DEBUG
    /*
    printf(" gain: \n");
    for (size_t jj = 0; jj < 4; jj++)
    {
        printf("H gain %zu → %lf\t", jj, gain_h[jj]);
        printf("V gain %zu → %lf\n", jj, gain_v[jj]);
    }
    */
    // DEBUG

    /* Final positions. */
    positions_calc(ds, supmat, pos_h);
    positions_calc(ds, supmat + 8, pos_v);

    // DEBUG
    /*
    printf("\n\n\n##### (random walk) AFTER #####\n");
    printf("\n\n >> gain << \n");
    for (int ii = 0; ii < 4 ; ii++)
    printf(" %10.4f  ", gain_h[ii]);
    printf("\n");
    for (int ii = 0; ii < 4 ; ii++)
    printf(" %10.4f  ", gain_v[ii]);
    printf("\n\n");
    printf("pos H, V = %lf / %lf\n\n", pos_h[5], pos_v[5]);
    */
    // DEBUG

    free(supmat_tmp);

    return accept;
}