#include "cxxheaders.h"
#include "mpi.h"
#include "debugfunc.h"

void showCheckPoint(const char *msg)
{
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0) printf("%s\n", msg);
    MPI_Barrier(MPI_COMM_WORLD);
}
