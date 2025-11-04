#include "prm_def.h"
#include <stdio.h>


/* Print the cartesian positions (hh, vv) of n sites (nsites) in the grid,
 * given by the idx index order array.
 */
void positions_print(dataset * ds, double * hh, double * vv)
{
    int ia;
    for (int ii = 0; ii < ds->nsites; ii++)
    {
        ia = ds->ord_sites[ii];
        printf("%12.4f %12.4f %12.4f %12.4f \n",
               ds->nom_h[ia], ds->nom_v[ia], hh[ia], vv[ia]);
    }
    // printf("\n");
}
