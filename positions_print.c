#include "prm_def.h"
#include <stdio.h>
#include <string.h>

/* Print the cartesian positions (hh, vv) of n sites (nsites) in the grid,
 * given by the idx index order array.
 */
void positions_print(const dataset * ds,
                    const double * hh, const double * vv,
                    const char * outfile)
{
    FILE * fout = NULL;
    if (outfile != NULL && strlen(outfile) > 0)
    {
        fout = fopen(outfile, "w");
        if (fout == NULL)
        {
            printf(" ERROR (positions_print): could not open"
                " output file '%s'.\n", outfile);
                return;
            }
        }
        else
    {
        fout = stdout;
    }

    fprintf(fout, "#    nom pos h,   nom pos v,"
                  "         pos h,       pos v\n");
    size_t ia;
    for (size_t ii = 0; ii < ds->nsites; ii++)
    {
        ia = ds->ord_sites[ii];
        fprintf(fout, "%12.4f %12.4f %12.4f %12.4f \n",
               ds->nom_h[ia], ds->nom_v[ia], hh[ia], vv[ia]);
    }

    if (fout != stdout)
    {
        fclose(fout);
    }    
}
