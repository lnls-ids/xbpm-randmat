#include "prm_def.h"
// #include "pcg_random.h"

#include <errno.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* prototype: help text. */
void help(void);

/* Initialize parameters.
 * Define default values for the general parameters.
 */
void parameters_initialize (xbpm_prm * prm)
{
    prm->nrand    = 1000;
    prm->temp     =  1.0;
    prm->step     =  0.1;
    prm->roi_from = -8.0;
    prm->roi_to   =  8.0;
    prm->nsites   =  0;
    strcpy(prm->datafile, "");
    strcpy(prm->matfile, "");
}


/* Get parameters form command line. Return a parameters struct.
 */
xbpm_prm parameters_read (int argc, char **argv)
{
    int opt;
    xbpm_prm prm;
  
    parameters_initialize(&prm);

    static struct option long_options[] =
    {
        /* "name", no_argument / required_argument, 0 /
            &verbose_flag, 'symbol' */
        {"help",  no_argument, 0, 'H'},
        //{"split",  no_argument, 0, 'S'},
        {0, 0, 0, 0}
    };
    /* getopt_long stores the option index here. */
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, "hHb:d:f:m:n:r:s:u:",
                            long_options, &option_index)) != -1)
    {
        switch (opt)
        {
        case 'b':                    /* Inverse of temperature. */
            prm.temp = optarg[0];
            break;
        
        case 'd':                   /* Input data file. */
            strcpy(prm.datafile, optarg);
            break;
        
        case 'f':                    /* ROI initial index. */
            prm.roi_from = atof(optarg);
            break;
            
        case 'h':  /* Help. */
            help();
            break;
            
        case 'H':  /* Help. */
            help();
            break;

        case 'm':                   /* Initial matrix file. */
            strcpy(prm.matfile, optarg);
            break;
        
        case 'n':                   /* Total number of sites. */
            prm.nsites = (size_t) strtoul(optarg, NULL, 10);
            break;
        
        case 'r':                    /* Number of random changes. */
            prm.nrand = (int) atoi(optarg);
            break;

        case 's':                    /* Step size. */
            prm.step = atof(optarg);
            break;

        case 'u':                  /* ROI last index. */
            prm.roi_to = atof(optarg);
            break;
        
            
        default:
            perror("No options given, aborting.\n");
            exit(-1);
    }
    }

    if (strlen(prm.datafile) == 0)
    {
        printf(" ERROR: no data file given."
               " Use '-h' to see options. Aborting.\n");
        exit(-1);
    }

    if (prm.nsites == 0)
    {
        printf(" ERROR: number of sites not defined."
            " Use '-h' to see options. Aborting.\n");
        exit(-1);
    }

    return prm;
}
