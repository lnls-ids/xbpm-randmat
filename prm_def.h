#ifndef PRM
#define PRM
#define MAX_LINE 1024
#include <stddef.h>

/* Struct for parameters.
 */
typedef struct
{
    char datafile[256];         /* Data file name.                */
    char matfile[256];          /* Suppression matrix file name.  */
    int nrand;                  /* Number of random trials.       */
    size_t nsites;              /* Number of sites in the grid.   */
    double roi_from, roi_to;    /* Interval that defines the ROI. */
    double temp;                /* Temperature.                   */
    double step;                /* Random step size.              */
} xbpm_prm;


/* Struct for ROI parameters and sites' indices.
 */
typedef struct
{
    size_t * idx;               /* Pointer for ROI sites's indices.   */
    size_t   nsites;           /* Actual number of sites in the ROI. */
} roi_struct;


/* Struct for data.
 */
typedef struct
{
    size_t nsites;           /* Number of sites.                            */
    size_t * ord_sites;      /* Indices for the correct order of positions. */
    double * nom_h, * nom_v; /* Nominal values of positions (sites).        */

    /* Pointers to blades' currents and std devs.   */
    double * to, * sto;         /* Top, out.        */
    double * ti, * sti;         /* Top, in.         */
    double * bi, * sbi;         /* Bottom, in.      */
    double * bo, * sbo;         /* Bottom, out.     */

    roi_struct roi;             /* Index for the sites within the ROI. */
} dataset;


/* Scaling parameters.
 */
typedef struct 
{
    double k, delta;
} kdelta;


/* The basic suppression matrix, defining the signals for horizontal 
 * and vertical calculations of beam position. 
 */
double supmat_signals[16] = {
    1.0,   1.0,  -1.0,  -1.0,
    1.0,   1.0,   1.0,   1.0,
    1.0,  -1.0,  -1.0,   1.0,
    1.0,   1.0,   1.0,   1.0
};


/* The basic pairwise blades calculation matrix.
 */
extern double paircalc_mat[4][4];

#endif
