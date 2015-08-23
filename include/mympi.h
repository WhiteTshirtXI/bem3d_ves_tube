#ifndef MYMPI_H
#define MYMPI_H

#include "mpi.h"

namespace mympi {
    // MPI shortcuts
    inline int comm_size(MPI_Comm comm = MPI_COMM_WORLD)
    {
	int mpi_inited;
	MPI_Initialized(&mpi_inited);

	int mpi_size = 1;
	if (mpi_inited) MPI_Comm_size(comm, &mpi_size);
	return mpi_size;
    }


    inline int comm_rank(MPI_Comm comm = MPI_COMM_WORLD)
    {
	int mpi_inited;
	MPI_Initialized(&mpi_inited);

	int mpi_rank = 0;
	if (mpi_inited) MPI_Comm_rank(comm, &mpi_rank);
	return mpi_rank;
    }

    void all_reduce(int nvar, double *u, int nloc, double *uloc, int *l2g, 
    			MPI_Comm comm = MPI_COMM_WORLD);
};

#endif
