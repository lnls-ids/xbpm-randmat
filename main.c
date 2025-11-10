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
double * matrix_read(char * matfile, double * supmat);

/* Read data from file. */
dataset data_read(xbpm_prm prm);

/* Perform a random walk with the gain matrix. */
size_t random_walk(dataset * ds, xbpm_prm * prm, double * supmat,
    double * pos_h, double * pos_v);

/* Print coordinates of sites. */
void positions_print(dataset * ds, double * pos_h, double * pos_v);

kdelta positions_calc(dataset * ds, double * supmat, double * pos);

/* Basic suppression matrix. */
double supmat_signs[16] = {
    1.0,   1.0,  -1.0,  -1.0,
    1.0,   1.0,   1.0,   1.0,
    1.0,  -1.0,  -1.0,   1.0,
    1.0,   1.0,   1.0,   1.0
};

/* Calculate the positions. */
/* Free up allocated memory. */
void dataset_free (dataset * ds, double * pos_h, double * pos_v,
                   double * supmat)
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
    free(supmat);
}


/* Print out matrix.
 */
void matrix_show (double * mat, size_t nn, size_t mm)
{
    size_t ii, jj;
    for (ii = 0; ii < nn; ii++)
    {
        for (jj = 0; jj < mm; jj++)
        {
            printf(" %12.6lf  ", mat[ii * mm + jj]);
        }
        printf("\n");
    }
    printf("\n");
}


/*
 */
int main(int argc, char **argv)
{
    size_t accept;

    /* Read parameters from command line. */
    xbpm_prm prm = parameters_read(argc, argv);

    /* Read XBPM data from file. */
    dataset ds = data_read(prm);


    /* Read initial suppression matrix from file if provided. */
    double * supmat = calloc(16, sizeof(double));
    if (supmat == NULL)
    {
        printf(" ERROR (main): could not allocate memory"
            " for suppression matrix. Aborting.\n");
        exit(-1);
    }


   if (strlen(prm.matfile) != 0)
   {   
        matrix_read(prm.matfile, supmat);
        printf("\n##### Input matrix:\n");
        matrix_show(supmat, 4, 4);
        printf("#####\n");
    } 
    else
    {
        /* Start with default suppression matrix. */
        memcpy(supmat, supmat_signs, 16 * sizeof(double));
    }

    /* Perform the random walk of suppression matrix's elements. */
    double * pos_h  = calloc(prm.nsites, sizeof(double));
    double * pos_v  = calloc(prm.nsites, sizeof(double));
    if (pos_h == NULL || pos_v == NULL)
    {
        printf(" ERROR (main): could not allocate memory"
            " for position arrays. Aborting.\n");
        exit(-1);
    }
    accept = random_walk(&ds, &prm, supmat, pos_h, pos_v);

    /* Show modified matrix. */
    printf("\n##### Modified matrix:\n");
    matrix_show(supmat, 4, 4);
    printf("#####\n\n");

    /* Rescale positions. */
    kdelta kdh = positions_calc(&ds, supmat, pos_h);
    kdelta kdv = positions_calc(&ds, supmat + 8, pos_v);

    /* Print out final positions. */
    printf("# nom pos h, nom pos v, pos h, pos v\n");
    positions_print(&ds, pos_h, pos_v);

    printf("\n Rescaling parameters:");
    printf("\n Horizontal: k = %lf, delta = %lf\n", kdh.k, kdh.delta);
    printf("\n Vertical:   k = %lf, delta = %lf\n", kdv.k, kdv.delta);
    printf("\n Acceptance rate: %lf %% \n",
           ((double) accept / (double) prm.nrand) * 100.0);

    /* Free up allocated memory. */
    dataset_free(&ds, pos_h, pos_v, supmat);
    return 0;
}
