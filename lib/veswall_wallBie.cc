#include "cxxheaders.h"
#include "veswall.h"
#include "ewald.h"
#include "mathfunc.h"
#include "mypetsc.h"
#include "debugfunc.h"
#include "param.h"

/* VesWall::solveWallBie
 * VesWall::initWallBieSolver
 * VesWall::getWallBieRhs
 * VesWall::getWallBieSol
 * VesWall::getWallBieScatter
 * solver::matmult */

/* Context for the matrix shell
 * Note: PETSc context only supports integer type, but we 
 *       need pointers, so use a local namespace  */
namespace {
namespace solver
{
    bool INITED = false;
    bool FIRST_SOLVE = true;

    VesWall *stks;

    Vec rhs_vec, sol_vec;
    Mat lhs;
    KSP ksp;
    PC pc;
    VecScatter scatter;

    void init();
    PetscErrorCode matmult(Mat, Vec, Vec);
}};


/* Solve the wall friction force 
 * Note:
 *  -- The rhs is stored in rhs_vec, and is calculated by calling calcRhs_wallBie
 *     before this function is called
 *  -- The result is stored in the f(:,:) array of each wall */
void VesWall::solveWallBie()
{
    assert(solver::INITED);

    using namespace solver;

    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    KSPSetTolerances(ksp, 1.E-5, PETSC_DEFAULT, PETSC_DEFAULT, 20);

    if (solver::FIRST_SOLVE) {
	KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
	solver::FIRST_SOLVE = false;
    } 
    else {
	KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    }

    ewald::phys::resetTiming();
    ewald::four::resetTiming();

    KSPSolve(ksp, rhs_vec, sol_vec);

    double wtime_phys, wtime_four;

    if (ewald::timing) {
        double wtime_phys, wtime_four;
	ewald::getWallTime(wtime_phys, wtime_four);

	if (mpi_rank == 0) {
	    printf("    Wall solver: wtime = %7.2f %7.2f s\n", wtime_phys, wtime_four);
	}
    }

    if (mpi_rank == 0) {
	int niter;
	KSPGetIterationNumber(ksp, &niter);

	double rnorm;
	KSPGetResidualNorm(ksp, &rnorm);

        printf("    Wall solver: niter = %d  rnorm = %.2E\n", niter, rnorm);
    }

    // Copy solution
    int nrow;
    VecGetSize(sol_vec, &nrow);
    MArray<double,1> sol(nrow);
    myVecScatter(solver::scatter, sol_vec, sol.data(), INSERT_VALUES, SCATTER_FORWARD);

    for (int iwall = 0; iwall < stks->numWalls(); iwall++) {
        Wall &wall = stks->walls[iwall];
	for (int ivert = 0; ivert < wall.numVerts(); ivert++) {
	    Point &vert = wall.verts[ivert];
	    m_dcopy(3, sol.data()+3*vert.Gindx, &wall.f(ivert,0));
	}
    }
}


void VesWall::initWallBieSolver()
{
    using solver::rhs_vec;
    using solver::sol_vec;
    using solver::scatter;
    using solver::lhs;
    using solver::matmult;
    using solver::ksp;
    using solver::pc;

    solver::stks = this;
    solver::INITED = true;

    // rhs and solution vectors
    int nrow = 3*vertList[WALL].size();
    VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, nrow, &rhs_vec);
    VecDuplicate(rhs_vec, &sol_vec);

    // Vector scatter
    myVecScatterCreateToAll(rhs_vec, scatter);

    // lhs matrix 
    int m;
    VecGetLocalSize(rhs_vec, &m);
    MatCreateShell(MPI_COMM_WORLD, m, m, nrow, nrow, NULL, &lhs);
    MatShellSetOperation(lhs, MATOP_MULT, (void(*)(void))matmult);

    // ksp solver
    KSPCreate(MPI_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, lhs, lhs, DIFFERENT_NONZERO_PATTERN);
    KSPSetType(ksp, KSPGMRES);
    KSPSetFromOptions(ksp);

    // No preconditioner
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCNONE);

}


Vec VesWall::getWallBieRhs()
{
    assert(solver::INITED);
    return solver::rhs_vec;
}


Vec VesWall::getWallBieSol()
{
    assert(solver::INITED);
    return solver::sol_vec;
}


VecScatter VesWall::getWallBieScatter()
{
    assert(solver::INITED);
    return solver::scatter;
}



/* Matrix-vector multiplication
 * Arguments:
 *  v_vec = lhs * src_vec */
PetscErrorCode solver::matmult(Mat lhs, Vec src_vec, Vec v_vec)
{
    // Each processor gets a copy of the global src array
    int nrow;
    VecGetSize(src_vec, &nrow);

    // src holds the local copy of the global src_vec vector
    MArray<double,1> src(nrow);
    myVecScatter(scatter, src_vec, src.data(), INSERT_VALUES, SCATTER_FORWARD);

    for (int iwall = 0; iwall < stks->numWalls(); iwall++) {
        Wall &wall = stks->walls[iwall];

	for (int ivert = 0; ivert < wall.numVerts(); ivert++) {
	    Point &vert = wall.verts[ivert];
	    m_dcopy(3, src.data()+3*vert.Gindx, vert.f);
	}
    }


    // Init
    VecZeroEntries(v_vec);

    double cf[] = { 0.0, 1.0/(8*M_PI) };
    double cg[] = { 0.0, 0.0 };
    stks->addBoundaryIntegrals(WALL, v_vec, cf, cg);

    return 0;
}
