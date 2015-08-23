// Projection step of the flow solver
#include "cxxheaders.h"
#include "stokesys.h"
#include "ewald.h"
#include "mathfunc.h"
#include "mypetsc.h"
#include "debugfunc.h"
#include "param.h"

/* VeSys::projectVel
 * Stokesys::initProjectSolver
 * Stokesys::finalizeProjectSolver 
 * ctx::matmult
 * ctx::lbd2vel */


/* Context for the matrix shell
 * Note: PETSc context only supports integer type, but we 
 *       need pointers, so use a local namespace  */
namespace {
namespace ctx
{
	double dt;
    StokeSys *stks;
    bool first_solve = true;

    const int NBLK = 3;
    int blk_size[NBLK], blk_offset[NBLK];

    double coeffSL[NBLK], coeffDL[NBLK];

    // rigid_xref -- coordiantes of the 1st triangle of every rigid particle
    Mat matSL_RigidSelf, matDL_RigidSelf;
    double (*rigid_xref)[3][3];

    Vec rhs_vec, sol_vec;
    Mat lhs, matSL[NBLK][NBLK], matDL[NBLK][NBLK];
    KSP ksp;
    VecScatter src_scatter;

    PetscErrorCode matmult(Mat, Vec, Vec);

    void calcLhs();  // Might be useful to put calcLhs and calcRhs in StokeSys
    void calcRhs();

    void addRigidSelfTerm(double, double, double *);

    void setSourceFromVec(Vec, double *, double *);
	void lbd2vel(const MArray<double1> &, MArray<double,2> &);  //Need to fill in this function
}};


/* VeSys::projectVel
 * VeSys::initProjectSolver
 * VeSys::finalizeProjectSolver 
 * ctx::matmult
 * ctx::lbd2vel */


/* Project the surface velocity to make it divergence free
 * Arguments:
 *  vold -- surface velocity before projection
 *  divTar -- target surface velocity divergence
 *  vnew -- surface velocity after projection
 *  lbd -- the Lagrange multiplier 
 *  initial_guess -- whether use lbd as initial guess */
void StokeSys::projectSol(const DoubArray2d &divTar, DoubArray1d &lbd,
		DoubArray2d &cell_v, DoubArray2d &rigid_v, DoubArray2d &wall_f,
		bool initial_guess)
{
    using namespace ctx;
    
	
	//Below here is copy/pasted and should be adjusted.  LHS can probably stay the same, but need to declare new RHS vector and update accordingly.
	
	
	//cout << "a\n";
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    
    //======================================================================
    // Single- and double-layer coefficients in the lhs
    ctx::coeffSL[CELL] = 1.0/(8*M_PI);
    ctx::coeffDL[CELL] = (stks->viscRat - 1.0)/(8*M_PI);

    ctx::coeffSL[RIGID] = 0.0;
    ctx::coeffDL[RIGID] = 1.0/(8*M_PI);

    ctx::coeffSL[WALL] = 1.0/(8*M_PI);
    ctx::coeffDL[WALL] = 0.0;
    //cout << "b\n";
    //======================================================================
    // Calculate lhs sub-matrices except wall-wall matrix, which is precomputed
    // 
    // THE LHS OPERATOR IS UNCHANGED, BUT THE ARGUMENT OF THE OPERATOR IS NOW
    //    u^(n+1) - u*,   where u* is the intermediate velocity
    ewald::phys::resetTiming();
    //cout << "c\n";
    ctx::calcLhs();
    //cout << "d\n";
    if (ewald::timing) {
        double wtime_phys, wtime_four;
	ewald::getWallTime(wtime_phys, wtime_four);

	if (mpi_rank == 0) {
	    printf("    calc lhs matrix, wtime = %7.2f s\n", wtime_phys);
	}
    }

    //======================================================================
    // THIS IS WHAT'S REALLY GOING TO CHANGE - THE RIGHT-HAND SIDE
    ewald::phys::resetTiming();
    ewald::four::resetTiming();
    //cout << "e\n";
    ctx::calcRhs();
	//======================================================================

    if (ewald::timing) {
        double wtime_phys, wtime_four;
	ewald::getWallTime(wtime_phys, wtime_four);

	if (mpi_rank == 0) {
	    printf("    calc rhs, wtime = %7.2f %7.2f s\n", wtime_phys, wtime_four);
	}
    }


    //======================================================================
    // Solve the equation
    KSPSetTolerances(ksp, 5.E-5, PETSC_DEFAULT, PETSC_DEFAULT, 25);

    if (first_solve) {
	KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
	first_solve = false;
    } 
    else {
        // Use last time step's solution as initial guess
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
	    printf("    matmult, wtime = %7.2f %7.2f s\n", wtime_phys, wtime_four);
	}
    }

    if (mpi_rank == 0) {
	int niter;
	KSPGetIterationNumber(ksp, &niter);

	double rnorm;
	KSPGetResidualNorm(ksp, &rnorm);

        printf("    GMRES: niter = %d  rnorm = %.2E\n", niter, rnorm);
    }

    // Copy solution
    int nrow;
    VecGetSize(sol_vec, &nrow);
    MArray<double,1> sol(nrow);
    myVecScatter(ctx::src_scatter, sol_vec, sol.data(), INSERT_VALUES, SCATTER_FORWARD);

    int offset = blk_offset[CELL];
    for (int icell = 0; icell < stks->numCells(); icell++) {
        Vesicle &cell = stks->cells[icell];
	for (int ivert = 0; ivert < cell.numVerts(); ivert++) {
	    Point &vert = cell.verts[ivert];
	    m_dcopy(3, sol.data()+offset+3*vert.Gindx, &cell.v(ivert,0));
	}
    }

    offset = blk_offset[RIGID];
    for (int irigid = 0; irigid < stks->numRigids(); irigid++) {
        Rigid &rigid = stks->rigids[irigid];
	for (int ivert = 0; ivert < rigid.numVerts(); ivert++) {
	    Point &vert = rigid.verts[ivert];
	    m_dcopy(3, sol.data()+offset+3*vert.Gindx, &rigid.v(ivert,0));
	}
    }

    offset = blk_offset[WALL];
    for (int iwall = 0; iwall < stks->numWalls(); iwall++) {
        Wall &wall = stks->walls[iwall];
	for (int ivert = 0; ivert < wall.numVerts(); ivert++) {
	    Point &vert = wall.verts[ivert];
	    m_dcopy(3, sol.data()+offset+3*vert.Gindx, &wall.f(ivert,0));
	}
    }
    
    
    
    
}





void StokeSys::initProjectSolver()
{
	//Most of this structure is correct, but references to rhs_vec should be updated to account for the fact that we are now solving for a lbd instead of velocity
	
	
    using namespace ctx;

    int nrow, m;
    PC pc;

    stks = this;

    // Submatrix size and index offset
    for (int i = 0; i < NBLK; i++) {
        blk_size[i] = 3*stks->vertList[i].size();
		blk_offset[i] = (i == 0)? 0 : blk_offset[i-1] + blk_size[i-1];
    }

    // Create rhs and solution vectors
    nrow = blk_size[0] + blk_size[1] + blk_size[2];
    VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, nrow, &rhs_vec);
    VecDuplicate(rhs_vec, &sol_vec);
    //cout << "if0\n";
    // Create lhs matrix and ksp solver
    VecGetLocalSize(rhs_vec, &m);
    MatCreateShell(MPI_COMM_WORLD, m, m, nrow, nrow, NULL, &lhs);
    MatShellSetOperation(lhs, MATOP_MULT, (void(*)(void))matmult);

    KSPCreate(MPI_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, lhs, lhs, DIFFERENT_NONZERO_PATTERN);
    KSPSetType(ksp, KSPGMRES);
    KSPSetFromOptions(ksp);

    // No preconditioner
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCNONE);

    // Create vector scatter
    myVecScatterCreateToAll(sol_vec, src_scatter);
    //cout << "if1\n";
    // Pre-compute rigid self-interaction and wall-wall interaction
    stks->calcFixedSTNList();

    if (ewald::phys::active) {
        for (int irigid = 0; irigid < numRigids(); irigid++) {
	    Rigid &rigid = stks->rigids[irigid];
	    rigid.updateGeometry();
	}

	SrcList &slist = stks->slist_RigidSelf;
	NbrList &nlist = stks->nlist_RigidSelf;

	// Compute the matrix
	int nrow = 3*nlist.verts.size();
	int ncol = 3*stks->vertList[RIGID].size();

	Mat &SL = matSL_RigidSelf;
	Mat &DL = matDL_RigidSelf;
	//cout << "if2\n";
	MatCreate(PETSC_COMM_SELF, &SL);
	MatSetType(SL, MATSEQAIJ);
	MatSetSizes(SL, PETSC_DECIDE, PETSC_DECIDE, nrow, ncol);
	//cout << "if3\n";
	MatCreate(PETSC_COMM_SELF, &DL);
	MatSetType(DL, MATSEQAIJ);
	MatSetSizes(DL, PETSC_DECIDE, PETSC_DECIDE, nrow, ncol);
	//cout << "if4\n";
	ewald::phys::calcSurfIntMat(nlist, SL, DL);

	// Save the coordinates of the 1st triangle of each particle for future
	// coordinate transform computation
	int nrigid = stks->numRigids();
	rigid_xref = new double[nrigid][3][3];
	//cout << "if5\n";
	for (int irigid = 0; irigid < nrigid; irigid++) {
	    Rigid &rigid = stks->rigids[irigid];
	    Tri &face = rigid.faces[0];

	    FOR_I3 FOR_J3 rigid_xref[irigid][i][j] = face.vert[i]->x[j];
	}
    }

    // Pre-compute wall-wall interaction
    if (ewald::phys::active) {
	SrcList &slist = stks->slist_phys[WALL];
	NbrList &nlist = stks->nlist_phys[WALL][WALL];
       
	int nrow = 3*nlist.verts.size();
	int ncol = 3*stks->vertList[WALL].size();
	//cout << "if6\n";
	Mat &SL = matSL[WALL][WALL];
	MatCreate(PETSC_COMM_SELF, &SL);
	MatSetType(SL, MATSEQAIJ);
	MatSetSizes(SL, PETSC_DECIDE, PETSC_DECIDE, nrow, ncol);
	//cout << "if7\n";
	ewald::phys::calcSurfIntMat(nlist, SL, NULL);
    }
}


void StokeSys::finalizeProjectSolver()  //Should be good
{
}


/* Calculate velocity induced by surface tension
 * Arguments:
 *  lbd -- the surface tension
 *  v -- the induced velocity */
void ctx::lbd2vel(const MArray<double,1> &lbd, MArray<double,2> &v)
{
	
	
	
    // Compute source density
    for (int ives= 0; ives < stks->numCells(); ives++) {
        Vesicle &vesicle = stks->cells[ives];
	if ( ! vesicle.isActive ) continue;  //Check this domain decomposition flag

	int nvert = vesicle.numVerts();
	MArray<double,1> lbdtmp(nvert);
	MArray<double,2> ftmp(nvert,3);

	int p = vesicle.verts[0].Gindx;
	m_dcopy(nvert, &lbd(p), lbdtmp.data());
	vesicle.tensionForce(lbdtmp, ftmp);

	for (int ivert = 0; ivert < nvert; ivert++) {
	    Point &vert = vesicle.verts[ivert];

	    m_dcopy(3, &ftmp(ivert,0), vert.f);
	    m_dclear(3, vert.g);
	}
    } 



    //This is copy pasted from vesys_project for vesys structures and needs to be checked.  Check that all these steps have analogous counterparts in stks

    // Compute -1/(8*pi*mu) N(f)
    VecZeroEntries(stks->rhs_vel);
    stks->addBoundaryIntegral(-1.0/(8*M_PI), 0.0, stks->rhs_vel); //Find analogous function
    VecAssemblyBegin(stks->rhs_vel);
    VecAssemblyEnd(stks->rhs_vel);

    // Calculate the induced velocity by solving the BIE
    stks->solveBie(0.0);
    myVecScatter(stks->scatter_vel, stks->sol_vel, v.data(), 
    			INSERT_VALUES, SCATTER_FORWARD);
}



