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
dataset data_read(xbpm_prm * prm);

/* Perform a random walk with the gain matrix. */
rw_stats random_walk(dataset * ds, xbpm_prm * prm, double * supmat,
                     double * pos_h, double * pos_v);

/* Print coordinates of sites. */
void positions_print(const dataset * ds,
                     const double * pos_h, const double * pos_v,
                     const char * outfile);

kdelta positions_calc(const dataset * ds, const double * supmat,
                      const double * nompos, double * pos);

/* Basic suppression matrix.
 * It represents the usual delta/sigma calculation.
 */
double supmat_signs[16] = {
    1.0,   1.0,  -1.0,  -1.0,
    1.0,   1.0,   1.0,   1.0,
    1.0,  -1.0,  -1.0,   1.0,
    1.0,   1.0,   1.0,   1.0
};


/* Read suppresssion matrix.
 */
double * suppression_matrix_read(char * matfile)
{
    double * supmat = calloc(16, sizeof(double));
    if (supmat == NULL)
    {
        printf(" ERROR (main): could not allocate memory"
            " for suppression matrix. Aborting.\n");
        exit(-1);
    }

    if (strlen(matfile) != 0)
    {   
        matrix_read(matfile, supmat);
    } 
    else
    {
        /* Start with default suppression matrix. */
        memcpy(supmat, supmat_signs, 16 * sizeof(double));
    }
    return supmat;
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
            printf(" %10.6lf  ", mat[ii * mm + jj]);
        }
        printf("\n");
    }
    printf("\n");
}


/* Print scaling parameters.
 */
void scaling_params_print (kdelta kdh, kdelta kdv,
                           rw_stats rws, size_t nrand, double step)
{
    printf("\n##### Rescaling parameters:");
    printf("\n Horizontal:\n"
           "    k     = %12.6lf,\n    delta = %12.6lf", kdh.k, kdh.delta);
    printf("\n\n Vertical:\n"
           "    k     = %12.6lf,\n    delta = %12.6lf", kdv.k, kdv.delta);

    printf("\n\n##### Random walk statistics."
           "\n Total matrix changes H = %.2lf %%",
           (double) rws.imat_h / (double) nrand * 100.0);
    printf("\n Total matrix changes V = %.2lf %%",
           (double) rws.imat_v / (double) nrand * 100.0);
    printf("\n Acceptance rate        = %6.2lf %%",
           ((double) rws.accept / (double) nrand) * 100.0);
    printf("\n Final temperature      = %.4g \t(beta = %.4g)",
           1/rws.beta, rws.beta);
    printf("\n Final step size        = %.4g \n\n", step);

    // printf("\n Acceptance rate: %12.4lf %% \n\n",
    //        ((double) accept / (double) nrand) * 100.0);
    printf("\n");
}


/* Free up allocated memory. */
void dataset_free (dataset * ds, double * supmat,
                   double * pos_h, double * pos_v)
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
    free(pos_h);
    free(pos_v);
}


/* Main program.
 */
int main(int argc, char **argv)
{
    /* Read parameters from command line. */
    xbpm_prm prm = parameters_read(argc, argv);

    /* Read XBPM data from file. */
    dataset ds = data_read(&prm);

    /* Read initial suppression matrix from file if provided. */
    double * supmat = suppression_matrix_read(prm.matfile);
    printf("##### Input matrix:\n");
    matrix_show(supmat, 4, 4);

    /* Perform the random walk of suppression matrix's elements. */
    double * pos_h  = calloc(prm.nsites, sizeof(double));
    double * pos_v  = calloc(prm.nsites, sizeof(double));
    if (pos_h == NULL || pos_v == NULL)
    {
        printf(" ERROR (main): could not allocate memory"
               " for position arrays. Aborting.\n");
        exit(-1);
    }
    rw_stats rws = random_walk(&ds, &prm, supmat, pos_h, pos_v);

    /* Show modified matrix. */
    printf("##### Modified matrix:\n");
    matrix_show(supmat, 4, 4);

    /* Rescale positions. */
    kdelta kdh = positions_calc(&ds, supmat, ds.nom_h, pos_h);
    kdelta kdv = positions_calc(&ds, supmat + 8, ds.nom_v, pos_v);

    /* Print out final positions. */
    positions_print(&ds, pos_h, pos_v, prm.outfile);

    /* Print final scaling parameters. */
    scaling_params_print(kdh, kdv, rws, prm.nrand, prm.step);

    /* Free up allocated memory. */
    dataset_free(&ds, supmat, pos_h, pos_v);
    return 0;
}
