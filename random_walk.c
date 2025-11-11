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
    double oldval;
    double chi2_h, chi2_v, chi2_h_aft, chi2_v_aft;
    double chi2, chi2_aft, dchi2;

    /* Probability.*/
    double prob;
    
    /* Initialize random seed. */
    uint64_t seed = seed_get();
    pcg32_init(seed);
    
    /* Calculate initial positions and deviation from nominal
     * positions (chi2). */
    positions_calc(ds, supmat, ds->nom_h, pos_h);
    chi2_h = chi2_calc(ds->nom_h, pos_h, &ds->roi);
    
    positions_calc(ds, supmat, ds->nom_v, pos_v);
    chi2_v = chi2_calc(ds->nom_v, pos_v, &ds->roi);
    
    chi2_h_aft = chi2_h;
    chi2_v_aft = chi2_v;
    chi2       = (chi2_h + chi2_v) * 0.5;
    chi2_aft   = chi2;
    
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
   
    /* Try to change the matrix nrand times. */
    size_t isite = 0;
    double sign = 1.0;
    for (ii = 0; ii < prm->nrand; ii++)
    {
        /* Pick an element of the suppression matrix. */
        isite   = (size_t) (pcg_double() * 16);
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
            positions_calc(ds, supmat, ds->nom_h, pos_h);
            chi2_h_aft = chi2_calc(ds->nom_h, pos_h, &ds->roi);
        }
        else
        {
            positions_calc(ds, supmat + 8, ds->nom_v, pos_v);
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
        prob = exp(-dchi2 * beta);
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
            prm->temp *= 0.999;
            chi2 = chi2_aft;
            chi2_h = chi2_h_aft;
            chi2_v = chi2_v_aft;
            accept++;
        }
        else
        {
            /* If change is rejected, restore former value. */
            supmat[isite] = oldval;
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
    positions_calc(ds, supmat, ds->nom_h, pos_h);
    positions_calc(ds, supmat + 8, ds->nom_v, pos_v);

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

    return accept;
}