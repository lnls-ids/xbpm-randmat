#include "prm_def.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NLIN 4
#define NCOL 4

/* Prototypes.      */
/* Read parameters. */
xbpm_prm parameters_read(int argc, char **argv);

/* Read matrix from file. */
double ** matrix_read(char * matfile, double ** supmat);

/* Read data from file. */
dataset data_read(xbpm_prm prm);

/* Perform a random walk with the gain matrix. */
void random_walk(dataset * ds, xbpm_prm * prm,
                 double * gain_h, double * gain_v,
                 double * pos_h, double * pos_v);


/* Print coordinates of sites. */
void positions_print(dataset * ds, double * hh, double * vv);

/* Calculate the positions. */
int positions_calc_h (dataset * ds, double * gain_h, double * hh);
int positions_calc_v (dataset * ds, double * gain_v, double * vv);


/* Free up allocated memory. */
void dataset_free (dataset * ds)
{
    free(ds->nom_h);
    free(ds->nom_v);
    free(ds->to);
    free(ds->ti);
    free(ds->bi);
    free(ds->bo);
    free(ds->sto);
    free(ds->sti);
    free(ds->sbi);
    free(ds->sbo);
    free(ds->ord_sites);
    free(ds->roi.idx);
}


/* Print out matrix.
 */
void matrix_show (double ** mat, int nn, int mm)
{
    int ii, jj;
    for (ii = 0; ii < nn; ii++)
    {
        for (jj = 0; jj < mm; jj++)
        {
            printf(" %12.6lf  ", mat[ii][jj]);
        }
        printf("\n");
    }
}


/*
 */
int main(int argc, char **argv)
{
    /* Read parameters from command line. */
    xbpm_prm prm = parameters_read(argc, argv);

    /* Read XBPM data from file. */
    dataset ds = data_read(prm);

    int ii;
    double gain_h[4] = {1., 1., 1., 1.};
    double gain_v[4] = {1., 1., 1., 1.};
    double * hh  = calloc(prm.nsites, sizeof(double));
    double * vv  = calloc(prm.nsites, sizeof(double));

    /* The suppression matrix. */
    double ** supmat = calloc(4, sizeof(double *));
    for (ii = 0; ii < 4; ii++)
    {
        supmat[ii] = calloc(4, sizeof(double));
    }

    if (strlen(prm.matfile) != 0)
    {   
        matrix_read(prm.matfile, supmat);
        for (ii = 0; ii < 4; ii++)
        {
            gain_h[ii] = supmat[0][ii];
            gain_v[ii] = supmat[2][ii];
        }
        printf("\n##### Input matrix:\n");
        matrix_show(supmat, 4, 4);
        printf("#####\n");
    }

    /* Perform the random walk of gain matrix's elements. */
    random_walk(&ds, &prm, gain_h, gain_v, hh, vv);

    /* Show modified matrix. */
    for (ii = 0; ii < 4; ii++)
    {
        supmat[0][ii] = gain_h[ii];
        supmat[1][ii] = fabs(gain_h[ii]);
        supmat[2][ii] = gain_v[ii];
        supmat[3][ii] = fabs(gain_v[ii]);
    }
    printf("\n##### Modified matrix:\n");
    matrix_show(supmat, 4, 4);
    printf("#####\n\n");

    /* Rescale positions. */
    positions_calc_h(&ds, gain_h, hh);
    positions_calc_v(&ds, gain_v, vv);

    /* Print out final positions. */
    printf("# nom pos h, nom pos v, pos h, pos v\n");
    positions_print(&ds, hh, vv);

    /* Free up allocated memory. */
    if (supmat != NULL)
    {
        for (ii = 0; ii < 4; ii++)
        {
            free(supmat[ii]);
        }
        free(supmat);
    }
    dataset_free(&ds);
    free(hh);
    free(vv);
    return 0;
}
