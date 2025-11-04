#include "prm_def.h"
#include "pcg_random.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

int positions_calc_h (dataset * ds, double * gain_h, double * hh);
int positions_calc_v (dataset * ds, double * gain_v, double * vv);


/* Change the value of an element of the gain array by a 'step'.
 */
void mat_walk (double * gain, double step)
{
    int isite   = (int) (pcg_double() * 4);           /* Pick an element.  */
    double sign = 1.0 ? -1.0 : (pcg_double() > 0.5);  /* Choose sign +/-.  */
   gain[isite] += sign * step;          /* Increase step in element value. */
}


/* Chi-square of two vectors, v1 and v2, indexed by roi->idx.
 */
double chi2_calc(double * v1, double * v2, roi_struct * roi)
{
    int ii, idx;
    double df, c2 = 0.0;
    for (ii = 0; ii < roi->nsites; ii++)
    {
        idx = roi->idx[ii];
        df  = v1[idx] - v2[idx];
        c2 += df *  df;
    }
    return c2 / (roi->nsites - 1);
}


/* Initialize seed with real random value from urandom device.
 */
uint64_t seed_get()
{
    FILE * sd = fopen("/dev/urandom", "rb");
    size_t nread;
    uint64_t buffer;
    nread = fread(&buffer, 8, 1, sd);
    fclose(sd);
    if (nread == 0)
    {
        printf("\n WARNING : nothing read from urandom.\n");
    }
    return buffer;
}


int random_walk(dataset * ds, xbpm_prm * prm,
                 double * gain_h, double * gain_v,
                 double * pos_h, double * pos_v)
{
    double * gain_h_tmp = calloc(4, sizeof(double));
    double * gain_v_tmp = calloc(4, sizeof(double));

    int ii = 0;
    int accept = 0;

    double chi2_h_bef, chi2_v_bef, chi2_h_aft, chi2_v_aft;
    double chi2_bef, chi2_aft, dchi2;
    double prob;

    uint64_t seed = seed_get();
    pcg32_init(seed);

    /* Calculate positions. */
    positions_calc_h(ds, gain_h, pos_h);
    positions_calc_v(ds, gain_v, pos_v);

    /* Measrurement of deviation: Chi2. */
    chi2_h_bef = chi2_calc(ds->nom_h, pos_h, &ds->roi);
    chi2_v_bef = chi2_calc(ds->nom_v, pos_v, &ds->roi);
    chi2_bef   = (chi2_h_bef + chi2_v_bef) * 0.5;
    chi2_h_aft = chi2_h_bef;
    chi2_v_aft = chi2_v_bef;

    // DEBUG
    // printf("\n##### (random walk) BEFORE #####\n");
    // printf("\n\n >> gain << \n (random walk) H : \n");
    // for (int ii = 0; ii < 4 ; ii++)
    //     printf(" %10.4f  ", gain_h[ii]);
    // printf("\n\n (random walk) V : \n ");
    // for (int ii = 0; ii < 4 ; ii++)
    //     printf(" %10.4f  ", gain_v[ii]);
    // printf("\n\n");
    // printf(" (random walk) H, V = %lf / %lf\n\n", pos_h[5], pos_v[5]);
    // DEBUG

    /* Try to change the matrix nrand times. */
    for (ii = 0; ii < prm->nrand; ii++)
    {
        /* Choose between a horizontal or vertical change. */
        if (pcg_double() < 0.5)
        {
            memcpy(gain_h_tmp, gain_h, 4 * sizeof(double));
            mat_walk(gain_h_tmp, prm->step);
            positions_calc_h(ds, gain_h_tmp, pos_h);
            chi2_h_aft = chi2_calc(ds->nom_h, pos_h, &ds->roi);
        } 
        else 
        {
            memcpy(gain_v_tmp, gain_v, 4 * sizeof(double));
            mat_walk(gain_v_tmp, prm->step);
            positions_calc_v(ds, gain_v_tmp, pos_v);
            chi2_v_aft = chi2_calc(ds->nom_v, pos_v, &ds->roi);
        }

        // DEBUG
        // printf(" pos 0: H, V = %lf / %lf\n ", pos_h[0], pos_v[0]);

        // for (int jj = 0; jj < 4; jj++)
        // {
        //     printf("h = %lf ", gain_h_tmp[jj]);
        //     printf("v = %lf ", gain_v_tmp[jj]);
        // }

        // printf(" (random walk) chi2 = %lf , %lf\n", chi2_h_aft, chi2_v_aft);
        // DEBUG

        /* Calculate the change in chi2. */
        chi2_aft = (chi2_h_aft + chi2_v_aft) * 0.5;
        dchi2    = chi2_aft - chi2_bef;
        prob     = exp(-dchi2 / prm->temp);
        if (prob > 1.0)
            prob = 1.0;

        // DEBUG
        // printf(" (random walk) AFT >> dchi2 = %lf - %lf = %lf ←→ %lf\n\n",
        //         chi2_aft, chi2_bef, dchi2, exp(-dchi2 / prm->temp));
        // DEBUG

        // prob = (dchi2 < 0.0) ? 1.0 : exp(-dchi2 / prm->temp);
        // DEBUG
        // printf(" (random walk) >> prob = %lf ←→ %Lg\n", prob, pcg_double());
        // DEBUG

        /* Accept or reject change. */
        if (pcg_double() < prob)
        {
            /* If change is accomplished, copy state values 
             * to former variables and reduce temperature.
             */
            memcpy(gain_h, gain_h_tmp, 4 * sizeof(double));
            memcpy(gain_v, gain_v_tmp, 4 * sizeof(double));
            prm->temp *= 0.999;
            chi2_h_bef = chi2_h_aft;
            chi2_v_bef = chi2_v_aft;
            chi2_bef = chi2_aft;
            accept++;
        }
    }

    // DEBUG
    // printf(" gain: \n");
    // for (int jj = 0; jj < 4; jj++)
    // {
    //     printf("H gain %d → %lf\t", jj, gain_h[jj]);
    //     printf("V gain %d → %lf\n", jj, gain_v[jj]);
    // }
    // DEBUG

    /* Final posistions. */
    positions_calc_h(ds, gain_h, pos_h);
    positions_calc_v(ds, gain_v, pos_v);

    // DEBUG
    // printf("\n\n\n##### (random walk) AFTER #####\n");
    // printf("\n\n >> gain << \n");
    // for (int ii = 0; ii < 4 ; ii++)
    //     printf(" %10.4f  ", gain_h[ii]);
    // printf("\n");
    // for (int ii = 0; ii < 4 ; ii++)
    //     printf(" %10.4f  ", gain_v[ii]);
    // printf("\n\n");
    // printf("pos H, V = %lf / %lf\n\n", pos_h[5], pos_v[5]);
    // DEBUG

    free(gain_h_tmp);
    free(gain_v_tmp);

    return accept;
}