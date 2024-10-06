#include <mpi-ext.h>

void getrank26d_(int* ix, int* iy, int* iz, int* ia, int* ib, int* ic, int* irank)
{
     int ierr;
     ierr = FJMPI_Topology_rel_rank2xyzabc(*irank, ix, iy, iz, ia, ib, ic);
     return;
}
