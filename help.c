#include <stdio.h>
#include <stdlib.h>

void help (void)
{
  printf(" mc_search - random walks in parameter's space to optimize \n"
         "   fitting for suppression matrix.\n");
  printf(
    "\n Usage:"
    "\n    ./mc_search [options] [args]"
    "\n where  "
    "\n  -h                : this help"
    "\n  -b <inv temp>     : the inverse of the temperature, beta = 1/T "
    "\n  -d <data file>    : data file name "
    "\n  -f <init index>   : ROI initial index "
    "\n  -u <last index>   : ROI last index "
    "\n  -n <nsites>       : total number of sites in the grid "
    "\n  -r <rand change>  : number of random changes "
    "\n  -m <matrix file>  : initial matrix to be update by annealing"
    "\n  -s <changes size> : step size (for random changes in the parameters)"
    "\n"
  );
 
  //Hb:d:f:u:n:s:t]
  exit(0);
}
