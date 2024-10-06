#include <mpi-ext.h>

void get6d2rank_(int* ix, int* iy, int* iz, int* ia, int* ib, int* ic, int* irank)
{
     int ierr;
     ierr = FJMPI_Topology_rel_xyzabc2rank(*ix, *iy, *iz, *ia, *ib, *ic, irank);
     return;
}

void get3d2rank_(int* ix, int* iy, int* iz, int* irank)
{
     int ierr;
     ierr = FJMPI_Topology_xyz2rank(*ix, *iy, *iz, irank);
     return;
}
