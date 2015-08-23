#include "cxxheaders.h"
#include "mympi.h"
#include "mathfunc.h"

/* mympi::all_reduce */


/* Reduce an array
 * Arguments:
 *  nvar -- number of variables
 *  u -- reduced array
 *  uloc -- locoal array
 *  l2g -- local to global index
 *  comm -- the communicator
 * Note:
 *   Must reserve enough number of elements in u[] */
void mympi::all_reduce(int nvar, double *u, int nloc, double *uloc, int *l2g, MPI_Comm comm)
{
    int mpi_size;
    MPI_Comm_size(comm, &mpi_size);

    int *nloc_all, ntot;
    int *recvcounts, *displs;
    int *l2g_all;
    double *u_all;

    // Allocate working arrays
    nloc_all = new int[mpi_size];
    recvcounts = new int[mpi_size];
    displs = new int[mpi_size];

    // gather nloc
    MPI_Allgather(&nloc, 1, MPI_INT, nloc_all, 1, MPI_INT, comm);

    // Find the total number of points
    ntot = 0;
    for (int i = 0; i < mpi_size; ++i) ntot += nloc_all[i];

    // Allocate working arrays
    l2g_all = new int[ntot];
    u_all = new double[ntot*nvar];

    // Gather l2g_all
    for (int i = 0; i < mpi_size; ++i) 
        recvcounts[i] = nloc_all[i];

    displs[0] = 0;
    for (int i = 1; i < mpi_size; ++i) 
        displs[i] = displs[i-1] + recvcounts[i-1];

    MPI_Allgatherv(l2g, nloc, MPI_INT, l2g_all, recvcounts, displs, MPI_INT, comm);

    // Gather u_all
    for (int i = 0; i < mpi_size; ++i) 
        recvcounts[i] = nvar*nloc_all[i];

    displs[0] = 0;
    for (int i = 1; i < mpi_size; ++i)  
        displs[i] = displs[i-1] + recvcounts[i-1];

    MPI_Allgatherv(uloc, nloc*nvar, MPI_DOUBLE,
            u_all, recvcounts, displs, MPI_DOUBLE, comm);

    // Reduction
    for (int i = 0; i < ntot; i++) {
        int iglb = l2g_all[i];
        m_dadd(nvar, u_all + i*nvar, u + iglb*nvar);
    }

    // Delete working arrays
    delete [] nloc_all;
    delete [] recvcounts;
    delete [] displs;
    delete [] u_all;
    delete [] l2g_all;
}
