#include <stdio.h>
#include <stdlib.h>

void help (void)
{
  printf(" mc_search - random walks in parameter's space to optimize \n"
         "             fitting for suppression matrix.\n");
  printf(
    "\n Usage:"
    "\n    ./mc_search -d <data file> -n <sites> [options] [args]"
    "\n\n where  "
    "\n  -d <data file>    : data file name (obligatory)"
    "\n  -n <# sites>      : total number of sites in the grid (obligatory)"
    "\n\n with optional arguments"
    "\n  -h                : this help"
    "\n  -b <inv. temp>    : the inverse of the temperature, beta = 1/T"
    "\n  -f <init. index>  : ROI initial index (from)"
    "\n  -u <last index>   : ROI last index (up to)"
    "\n  -r <# rand.>      : number of random changes"
    "\n  -m <matrix file>  : initial matrix to be update by annealing"
    "\n  -s <changes size> : step size of random changes in the"
    "\n                      suppression matrix elements (default = 1e-5)"
    "\n"
    "\n The data is a 10-column text; the two first columns are the nominal"
    "\n positions. The other eight are the values measured by each blade"
    "\n with respective errors. The sequence of blades is:"
    "\n top out, top in, bottom in, bottom out."
    "\n"
    "\n The initial and last indices of the ROI are the positions of the"
    "\n ROI's boundaries. It defines the optimal adjustment domain."
    "\n"
    "\n If no initial matrix is provided, the program starts with a standard"
    "\n matrix, whose gains/suppressions are equal to 1."
    "\n"
    "\n Obs.: (Important) Check the signs of the initial matrix elements,"
    "\n       they must correspond to the usual delta/sigma calculation"
    "\n       (two first rows correspond to horizontal position calculation,"
    "\n       the last two to vertical)."
    "\n"
  );
 
  exit(0);
}
