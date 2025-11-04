#include "prm_def.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


/* Define a structure for minimum and maximum.
 */
typedef struct
{
    double min, max;
} minmax;


/* Return the minimum and maximum values of a vector vv of size nn.
 */
minmax min_and_max (double * vv, int nn)
{
    minmax mm;
    mm.min = vv[0];
    mm.max = vv[0];
    for (int ii = 1; ii < nn; ii++)
    {
        if (vv[ii] > mm.max)
            mm.max = vv[ii];
        if (vv[ii] < mm.min)
            mm.min = vv[ii];
    }
    return mm;
}


int site_visited (int * idx, int ic, int ivisit)
{
    for (int ii = 0; ii < ic; ii++)
    {
        if (ivisit == idx[ii])
        {
            printf(" True! ii = %d, ivisit = %d, idx = %d\n", 
                ii, ivisit, idx[ii]);
            return 1;
        }
    }
    printf(" False! ivisit = %d, idx = %d\n", ivisit, idx[ic]);
    return 0;
}


/* Exchange the values of two integer variables (a and b).
 */
void exchange_values(int * a, int * b)
{
    int tmp = *a;
    *a = *b;
    *b = tmp;
}


/* Given two arrays, hh and vv, with horizontal and vertical positions
 * of each site (out of nsites) _on a grid_, return an index array idx
 * of integers which maps the order through lines and columns
 * (from lower lines and columns to higher ones). It is supposed that 
 * positions given by the same index are correlated, namely, hh[i] and vv[i]
 * correspond to the same site.
 */
int * index_order_by_position (double * hh, double * vv, int nsites)
{
    int ic, jj, inext, icurr;
    int * idx = calloc(nsites, sizeof(int));
    double next_h, next_v;

    /* Initialize index array. */
    for (ic = 0; ic < nsites; ic++) { idx[ic] = ic; }

    /* Run through coordinate arrays. */
    idx[0] = 0;
    next_h = hh[0];
    next_v = vv[0];
    for (ic = 1; ic < nsites; ic++)
    {
        icurr = ic;
        next_h =  hh[ic];
        next_v =  vv[ic];
        inext  = idx[ic];
        for (jj = ic + 1; jj < nsites; jj++)
        {
            inext = idx[jj];
            /* Next coordinate x is lesser than previous one. */
            if (hh[inext] <= next_h)
            {
                /* Next coordinate y is lesser than previous one. */
                if (vv[inext] < next_v)  /* Guarantee no short circuit. */
                {
                    next_h = hh[inext];
                    next_v = vv[inext];
                    icurr = inext;
                }
            }
        }
        /* Rearrange order in index array. */
        exchange_values(&idx[ic], &idx[icurr]);
    }
    return idx;
}


/* Create an index structure for the ROI.
 */
roi_struct roi_indexation (dataset * ds, xbpm_prm prm)
{
    int   icount = 0;
    int * roi_sites = calloc(ds->nsites, sizeof(int));
    int   ord_idx;
    roi_struct roi;

    for (int ii = 0; ii < ds->nsites; ii++)
    {
        ord_idx = ds->ord_sites[ii];
        if (ds->nom_h[ord_idx] >= prm.roi_from && 
            ds->nom_h[ord_idx] <= prm.roi_to)
        {
            if (ds->nom_v[ord_idx] >= prm.roi_from &&
                ds->nom_v[ord_idx] <= prm.roi_to)
            {
                roi_sites[icount++] = ord_idx;
            }
        }

    }

    roi.nsites = icount;
    roi.idx    = calloc(roi.nsites, sizeof(int));
    memcpy(roi.idx, roi_sites, roi.nsites * sizeof(int));
    free(roi_sites);
    return roi;
}

/* Read matrix from file.
 */
void matrix_read(char * matfile, double ** mat)
{
    int ii;
    FILE * df = fopen(matfile, "r");

    for(ii = 0; ii < 4; ii++)
    {
        if(! fscanf(df, "%lf %lf %lf %lf",
            &mat[ii][0], &mat[ii][1], &mat[ii][2], &mat[ii][3]))
        {
            printf(" WARNING (at matrix read) : could not read"
                " from file %s at line (%d)", matfile, ii);
        }
    }
    fclose(df);
}


/* Read data from file.
 */
dataset data_read(xbpm_prm prm)
{
    FILE * df = fopen(prm.datafile, "r");
    char line[MAX_LINE];
    char * pd, * parse;
    dataset ds;
    int ii = 0;
    
    /* Check when opening data file. */
    if (df == NULL)
    {
        perror(prm.datafile);
        printf("##### (data_read) file: '%s'\n"
            "ERROR: Aborting.\n\n", prm.datafile);
            exit(-1);
    }

    /* Allocate space for data. */
    ds.nsites = prm.nsites;

    ds.nom_h = calloc(prm.nsites, sizeof(double));
    ds.nom_v = calloc(prm.nsites, sizeof(double));

    ds.to    = calloc(prm.nsites, sizeof(double));
    ds.ti    = calloc(prm.nsites, sizeof(double));
    ds.bi    = calloc(prm.nsites, sizeof(double));
    ds.bo    = calloc(prm.nsites, sizeof(double));

    ds.sto   = calloc(prm.nsites, sizeof(double));
    ds.sti   = calloc(prm.nsites, sizeof(double));
    ds.sbi   = calloc(prm.nsites, sizeof(double));
    ds.sbo   = calloc(prm.nsites, sizeof(double));


    // DEBUG
    // printf("\n>> (DATA READ) OK 1.\n");
    // DEBUG

    while (fgets(line, sizeof(line), df) != NULL)
    {
        parse = line;

        /* Skip empty lines. */
        if (line[0] == '\n')
            continue;
        
        ds.nom_h[ii] = atof(strtok_r(parse, " ", &pd));
        ds.nom_v[ii] = atof(strtok_r(NULL, " ", &pd));

        ds.to[ii]    = atof(strtok_r(NULL, " ", &pd));
        ds.sto[ii]   = atof(strtok_r(NULL, " ", &pd));

        ds.ti[ii]    = atof(strtok_r(NULL, " ", &pd));
        ds.sti[ii]   = atof(strtok_r(NULL, " ", &pd));

        ds.bi[ii]    = atof(strtok_r(NULL, " ", &pd));
        ds.sbi[ii]   = atof(strtok_r(NULL, " ", &pd));

        ds.bo[ii]    = atof(strtok_r(NULL, " ", &pd));
        ds.sbo[ii]   = atof(strtok_r(NULL, " ", &pd));


        // DEBUG
        // printf("\n>> (DATA READ) ii = %d\n", ii);
        // DEBUG

        ii++;
    }

    // DEBUG
    // printf("\n>> (DATA READ) AFTER WHILE OK.\n");
    // DEBUG

    ds.ord_sites = index_order_by_position(ds.nom_h, ds.nom_v, ds.nsites);
    ds.roi       = roi_indexation(&ds, prm);

    fclose(df);
    return ds;
}