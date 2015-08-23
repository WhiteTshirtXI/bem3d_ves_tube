#include "cxxheaders.h"
#include "mypetsc.h"
#include "mathfunc.h"

/* myAllReduce
 * myVecSetValues
 * myVecScatterCreateToAll 
 * myVecScatter 
 * myMatMult
 * myMatMultAdd  */


/* Reduce an array
 * Arguments:
 *  nvar -- number of variables
 *  u -- reduced array
 *  uloc -- locoal array
 *  l2g -- local to global index
 *  comm -- the communicator
 * Note:
 *   Must reserve enough number of elements in u[] */
void myAllReduce(int nvar, double *u, int nloc, double *uloc, int *l2g, MPI_Comm comm)
{
    int numprocs;
    int *nloc_all, ntot;
    int *recvcounts, *displs;
    int *l2g_all;
    double *u_all;

    MPI_Comm_size(comm, &numprocs);

    // allocate working arrays
    nloc_all = new int[numprocs];
    recvcounts = new int[numprocs];
    displs = new int[numprocs];

    // gather nloc
    MPI_Allgather(&nloc, 1, MPI_INT, nloc_all, 1, MPI_INT, comm);

    // find the total number of points
    ntot = 0;
    for (int i = 0; i < numprocs; ++i) ntot += nloc_all[i];

    // allocate working arrays
    l2g_all = new int[ntot];
    u_all = new double[ntot*nvar];

    // gather l2g_all
    for (int i = 0; i < numprocs; ++i) 
        recvcounts[i] = nloc_all[i];

    displs[0] = 0;
    for (int i = 1; i < numprocs; ++i) 
        displs[i] = displs[i-1] + recvcounts[i-1];

    MPI_Allgatherv(l2g, nloc, MPI_INT, l2g_all, recvcounts, displs, MPI_INT, comm);

    // gather u_all
    for (int i = 0; i < numprocs; ++i) 
        recvcounts[i] = nvar*nloc_all[i];

    displs[0] = 0;
    for (int i = 1; i < numprocs; ++i)  
        displs[i] = displs[i-1] + recvcounts[i-1];

    MPI_Allgatherv(uloc, nloc*nvar, MPI_DOUBLE,
            u_all, recvcounts, displs, MPI_DOUBLE, comm);

    // reduction
    for (int i = 0; i < ntot; i++) {
        int iglb = l2g_all[i];
        m_dadd(nvar, u_all + i*nvar, u + iglb*nvar);
    }

    // delete working arrays
    delete [] nloc_all;
    delete [] recvcounts;
    delete [] displs;
    delete [] u_all;
    delete [] l2g_all;
}


void myVecSetValues(Vec vec, const double *val, InsertMode mode) 
{
    int rowmin, rowmax;
    VecGetOwnershipRange(vec, &rowmin, &rowmax);

    int n = rowmax - rowmin;
    int *idx = new int[n];
    for (int p = 0; p < n; p++) idx[p] = p + rowmin;

    VecSetValues(vec, n, idx, val+rowmin, mode);

    delete [] idx;
}


/* PETSc VecScatterCreateToAll */
void myVecScatterCreateToAll(Vec vec, VecScatter &scatter)
{
    Vec vecTmp;
    VecScatterCreateToAll(vec, &scatter, &vecTmp);
    VecDestroy(&vecTmp);
}


/* PETSc vec scatter
 * Arguments:
 *   -- Same as VecScatterBegin() */
void myVecScatter(VecScatter scatter, Vec x, Vec y, InsertMode addv, ScatterMode mode)
{
    VecScatterBegin(scatter, x, y, addv, mode);
    VecScatterEnd(scatter, x, y, addv, mode);
}


/* PETSc vec scatter
 * Arguments:
 *   -- Same as VecScatterBegin() except the target vector is of pointer type */
void myVecScatter(VecScatter scatter, Vec x, double *py, InsertMode addv, ScatterMode mode)
{
    int nrow;
    Vec y;

    VecGetSize(x, &nrow);
    VecCreateSeqWithArray(PETSC_COMM_SELF, PETSC_DECIDE, nrow, py, &y);

    VecScatterBegin(scatter, x, y, addv, mode);
    VecScatterEnd(scatter, x, y, addv, mode);

    VecDestroy(&y);
}


/* Another overloaded VecScatter 
 * Note:
 *   1. It is useful for reverse scattering, i.e. y->x, where x is an parallel vector */
void myVecScatter(VecScatter scatter, double *py, Vec x, InsertMode addv, ScatterMode mode)
{
    int nrow;
    Vec y;

    VecGetSize(x, &nrow);
    VecCreateSeqWithArray(PETSC_COMM_SELF, PETSC_DECIDE, nrow, py, &y);

    VecScatterBegin(scatter, y, x, addv, mode);
    VecScatterEnd(scatter, y, x, addv, mode);

    VecDestroy(&y);
}


/* Short hand for matmult with sequential matrix 
 * Arguments:
 *   y = A x */
void myMatMult(Mat A, double *x, double *y)
{
    int nrow, ncol;
    MatGetSize(A, &nrow, &ncol);
    if (nrow == 0 || ncol == 0) return;

    Vec x_vec, y_vec;
    VecCreateSeqWithArray(PETSC_COMM_SELF, PETSC_DECIDE, nrow, y, &y_vec);
    VecCreateSeqWithArray(PETSC_COMM_SELF, PETSC_DECIDE, ncol, x, &x_vec);

    MatMult(A, x_vec, y_vec);

    VecDestroy(&x_vec);
    VecDestroy(&y_vec);
}


/* Arguments:
 *   y = A x + y */
void myMatMultAdd(Mat A, double *x, double *y)
{
    int nrow, ncol;
    MatGetSize(A, &nrow, &ncol);
    if (nrow == 0 || ncol == 0) return;

    Vec x_vec, y_vec;
    VecCreateSeqWithArray(PETSC_COMM_SELF, PETSC_DECIDE, nrow, y, &y_vec);
    VecCreateSeqWithArray(PETSC_COMM_SELF, PETSC_DECIDE, ncol, x, &x_vec);

    MatMultAdd(A, x_vec, y_vec, y_vec);

    VecDestroy(&x_vec);
    VecDestroy(&y_vec);
}


/* Short hand for matmult with sequential matrix 
 * Arguments:
 *   y = alpha*A*x + y */
void myMatMultAdd(double alpha, Mat A, double *x, double *y)
{
    if (fabs(alpha) < 1.E-10) return;

    int nrow, ncol;
    MatGetSize(A, &nrow, &ncol);
    if (nrow == 0 || ncol == 0) return;

    double *tmp = new double[nrow];

    Vec x_vec, tmp_vec;
    VecCreateSeqWithArray(PETSC_COMM_SELF, PETSC_DECIDE, nrow, tmp, &tmp_vec);
    VecCreateSeqWithArray(PETSC_COMM_SELF, PETSC_DECIDE, ncol, x, &x_vec);

    MatMult(A, x_vec, tmp_vec);
    m_daxpy(nrow, alpha, tmp, y);

    VecDestroy(&x_vec);
    VecDestroy(&tmp_vec);
    delete [] tmp;
}
