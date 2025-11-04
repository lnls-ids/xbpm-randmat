#ifndef PRM
#include "prm_def.h"
#define PRM
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pcg_random.h"

#define NLIN 4
#define NCOL 4

/* Prototypes.      */
/* Read parameters. */
xbpm_prm parameters_read(int argc, char **argv);

/* Read data from file. */
dataset data_read(char * datafile, int nsites);

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
}

double mat[4][4] = {
    {1.0,  1.0, -1.0, -1.0},
    {1.0,  1.0,  1.0,  1.0},
    {1.0, -1.0, -1.0,  1.0},
    {1.0,  1.0,  1.0,  1.0},
};

/* Calculate the product of an (mm x nn) matrix mA by an
 * (nn x pp) matrix mB. Returns a mm x pp matrix pr. */
double ** mat_prod(double ** mA, double ** mB,
                  size_t mm, size_t nn, size_t pp);

int main(int argc, char **argv)
{
    /* Read parameters from command line. */
    xbpm_prm prm = parameters_read(argc, argv);

    /* Read XBPM data from file. */
    dataset ds = data_read(prm.datafile, prm.nsites);

    // for (int ii=0; ii < prm.nsites; ii++)
    // {
    //     printf("to[%03d] = %lf\t %.4lf\n", ii, ds.nom_h[ii], ds.ti[ii]);
    // }

    double mA[1][4] = {{1.0, 2.0, 3.0, 4.0}};
    double mB[4][1] = {{1.0}, {1.0}, {1.0}, {1.0}};


    double ** pA = calloc(1, sizeof(double *));
    for (int ii=0; ii < 1; ii++)
    {
        pA[ii] = calloc(NCOL, sizeof(double));
        memcpy(pA[ii], mA[ii], NCOL * sizeof(double));
    }

    double ** pB = calloc(NLIN, sizeof(double *));
    for (int ii=0; ii < NLIN; ii++)
    {
        pB[ii] = calloc(1, sizeof(double));
        memcpy(pB[ii], mB[ii], sizeof(double));
    }

    printf("\n\n >>>>> (main) mA @ mB: res = \n ");
    double ** res = mat_prod(pA, pB, 1, 4, 1);
    for (int ii = 0; ii < 1; ii++)
    {
        for (int jj = 0; jj < 1; jj++)
        {            
            printf(">>>>> (main) [%d][%d] %.4lf  ", ii, jj, res[ii][jj]);
        }
        printf("\n");
    }
    // DEBUG
    printf("\n ##### ##### (main) END mA @ mB ");
    // DEBUG

    for (int ii = 0; ii < 1; ii++) free(res[ii]);
    free(res);
    // DEBUG
    printf("\n ##### ##### (main) FREED res ");
    // DEBUG


    printf("\n\n ##### NEXT ##### \n\n\n\n");
    fflush(stdout);

    double ** pmat = calloc(NLIN, sizeof(double *));
    printf(">>>>> (main) OK pmat allocated.\n");
    fflush(stdout);
    for (int ii = 0; ii < NLIN; ii++)
        pmat[ii] = mat[ii];
    
    // DEBUG
    printf("\n\n >>>>> (main) The 'mat' matrix: \n ");
    for (int ii = 0; ii < 4; ii++)
    {
        for (int jj = 0; jj < 4; jj++)
        {
            printf(" [%d][%d] = ", ii, jj);
            fflush(stdout);
            printf(" %.4lf  ", pmat[ii][jj]);
        }
        printf("\n");
    }
    // DEBUG
    

    printf("\n\n ##### NEXT ##### \n\n");
    fflush(stdout);

    printf("\n\n >>>>> (main) mat @ mA:  \n ");
    res = mat_prod(pA, pmat, 1, 4, 4);
    for (int ii = 0; ii < 1; ii++)
    {
        for (int jj = 0; jj < 4; jj++)
            printf(" (main) [%d][%d] %f ", ii, jj, res[ii][jj]);
        printf("\n");
    }

    // for (int ii = 0; ii < 4; ii++) 
    free(res[0]);
    free(res);

    printf("\n\n ##### NEXT ##### \n\n");




    printf("\n\n >>>>> (main) mat @ mB: res = \n ");
    res = mat_prod(pmat, pB, 4, 4, 1);
    for (int ii = 0; ii < 4; ii++)
        for (int jj = 0; jj < 1; jj++)
            printf(" (main) [%d][%d] %f\n ", ii, jj, res[ii][jj]);

    for (int ii = 0; ii < 4; ii++) free(res[ii]);
    free(res);


    printf("\n\n ##### END ##### \n\n");

    // for (int ii = 0; ii < 4; ii++)
    // free(pmat);
    free(pmat);
    free(pA[0]);
    free(pA);
    for (int ii=0; ii < NLIN; ii++) free(pB[ii]);
    free(pB);

    /* Free up allocated memory. */
    dataset_free(&ds);
    return 0;
}