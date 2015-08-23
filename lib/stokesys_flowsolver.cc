#include "cxxheaders.h"
#include "stokesys.h"
#include "ewald.h"
#include "mathfunc.h"
#include "mypetsc.h"
#include "debugfunc.h"
#include "param.h"

/* StokeSys::solveFlow
 * StokeSys::initFlowSolver
 * StokeSys::finalizeFlowSolver
 * StokeSys::verifySolution
 * ctx::calcLhs
 * ctx::calcRhs 
 * ctx::matmult
 * ctx::addRigidSelfTerm
 * ctx::setSourceFromVec */


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

    void calcLhs();
    void calcRhs();

    void addRigidSelfTerm(double, double, double *);

    void setSourceFromVec(Vec, double *, double *);
}};


void StokeSys::solveFlow()
{
    using namespace ctx;
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
    ewald::phys::resetTiming();
    ewald::four::resetTiming();
    //cout << "e\n";
    ctx::calcRhs();

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
    
    
    //======================================
    // PROJECTION
    //======================================
    
    const double tau = param:exist("AREA_RELAX_TIME");
    param::getDoubleValue("AREA_RELAX_TME") : 1.E99;
    
    MArray<double,1> divTar(vertList[CELL].size());
    MArray<double,1> dlbd(vertList[CELL].size());
    
    
    for (int p =0, int icell=0; icell < stks->numCells(); icell++) 
	{	
		Vesicle &cell = stks->cells[icell];
		int nvert = cell.numVerts();
		
		double s = (cell.areaTar/cell.area - 1.0)/tau;
	
		for (int ivert = 0; ivert < cell.numVerts(); ivert++)
		{
				divTar(p) = s*cell.vertArea(ivert);
				p++; //Gindx
		}
		
	}
	
	//Project and update the solution
	dlbd = 0.0;
	projectSol(divTar,dlbd,cell.v,rigid.v,wall.f,cell.sigma,true);
}


void StokeSys::initFlowSolver()
{
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


void StokeSys::finalizeFlowSolver()
{
}


/*
void StokeSys::verifySolution()
{
    using namespace ctx;

    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    double rhs_L1, rhs_L2, rhs_Linf;
    double err_L1, err_L2, err_Linf;

    //======================================================================
    // RHS norm
    int nrow;
    VecGetSize(rhs_vec, &nrow);

    VecNorm(rhs_vec, NORM_1, &rhs_L1);		//rhs_L1 /= nrow;
    VecNorm(rhs_vec, NORM_2, &rhs_L2);		//rhs_L2 /= sqrt(1.0*nrow);
    VecNorm(rhs_vec, NORM_INFINITY, &rhs_Linf);


    //======================================================================
    // Set source
    double *src = new double[nrow];
    ctx::setSourceFromVec(sol_vec, src);
    delete [] src;

    // source for RHS
    for (int icell = 0; icell < numCells(); icell++){
	Cell &cell = cells[icell];
	for (int ivert = 0; ivert < cell.numVerts(); ivert++){
	    Point &vert = cell.verts[ivert];
	    FOR_J3 vert.f[j] = cell.f(ivert,j);
	}
    }

    // No need to set up rigid surface force because it is not modified during the
    // solution procedure

    for (int iwall = 0; iwall < numWalls(); iwall++){
	Wall &wall = walls[iwall];
	for (int ivert = 0; ivert < wall.numVerts(); ivert++){
	    Point &vert = wall.verts[ivert];
	    FOR_J3 vert.g[j] = 0.0;
	}
    }

    //======================================================================
    // Calculate residual
    // Init
    VecZeroEntries(rhs_vec);

    // Add physical sum
    if (ewald::phys::active){
	for (int iblk = 0; iblk < NBLK; iblk++)
	for (int jblk = 0; jblk < NBLK; jblk++){
	    // Build neighbor list
	    NbrList nlist;
	    nlist.build(vertList[iblk], slist_phys[jblk]);

	    vector<Point*> &tlist = nlist.verts;
	    int nloc = tlist.size();
	    double *vloc = new double[3*nloc];
	    int *idx = new int[3*nloc];

	    double c1, c2;

	    switch (jblk){
                case CELL:
                    c1 = 1.0/(8*M_PI);
                    c2 = -(1.0 - viscRat)/(8*M_PI);
                    break;
                case RIGID:
                    c1 = 1.0/(8*M_PI);
                    c2 = 1.0/(8*M_PI);
                    break;
                case WALL:
                    c1 = 1.0/(8*M_PI);
                    c2 = 0.0;
                    break;
	    }

	    m_dclear(3*nloc, vloc);
            ewald::phys::addSurfInt(nlist, c1, c2, vloc);

	    setRowIndices(tlist, blk_offset[iblk], idx);
	    VecSetValues(rhs_vec, 3*nloc, idx, vloc, ADD_VALUES);

	    delete [] vloc;
	    delete [] idx;
        }
    } // ewald::physical

    if (ewald::four::active) {
	ewald::four::clear_source();

	if (numCells() > 0) {
	    double c1 = 1.0/(8*M_PI);
	    double c2 =  -(1.0 - viscRat)/(8*M_PI);
	    ewald::four::add_source(c1, c2, slist_four[CELL]);
	}

	if (numRigids() > 0) {
	    double c1 = 1.0/(8*M_PI);
	    double c2 = 1.0/(8*M_PI);
	    ewald::four::add_source(c1, c2, slist_four[RIGID]);
	}

	if (numWalls() > 0) {
	    double c1 = 1.0/(8*M_PI);
	    double c2 = 0.0;
	    ewald::four::add_source(c1, c2, slist_four[WALL]);
	}

	ewald::four::transform();

	for (int iblk = 0; iblk < NBLK; iblk++) {
	    vector<Point*> &tlist = tlist_four[iblk];
	    if (tlist.size()==0) continue;

	    int nloc = tlist.size();
	    double *vloc = new double[3*nloc];
	    int *idx = new int[3*nloc];

	    m_dclear(3*nloc, vloc);
	    ewald::four::add_interp_vel(tlist, vloc);

	    setRowIndices(tlist, blk_offset[iblk], idx);
	    VecSetValues(rhs_vec, 3*nloc, idx, vloc, ADD_VALUES);

	    delete [] vloc;
	    delete [] idx;
	}
    } // ewald::four

    // Other terms
    if (ewald::phys::active && ewald::phys::comm_rank == 0) {
	double vlin[3];

	// Cell
	{
	    const int iblk = CELL;

	    vector<Point*> &tlist = vertList[iblk];
	    int nloc = tlist.size();
	    int *idx = new int[3*nloc];
	    double *vloc = new double[3*nloc];

	    m_dclear(3*nloc, vloc);

	    for (int iloc = 0; iloc < nloc; iloc++){
		Point &vert = *tlist[iloc];
		FOR_J3 {
		    int idx = 3*iloc+j;
		    vloc[idx] += vert.g[j];
		    FOR_K3 vloc[idx] += -0.5*(1-viscRat)*vert.C[j][k]*vert.g[k];
		    vloc[idx] += vlin[j] - vbkg[j];
	        }
	    }

	    setRowIndices(tlist, blk_offset[iblk], idx);
	    VecSetValues(rhs_vec, 3*nloc, idx, vloc, ADD_VALUES);

	    delete [] vloc;
	    delete [] idx;
	} // cells

	// Rigid
	{
	    int iblk = RIGID;

	    vector<Point*> &tlist = vertList[iblk];
	    int nloc = tlist.size();
	    int *idx = new int[3*nloc];
	    double *vloc = new double[3*nloc];
	    m_dclear(3*nloc, vloc);


	    Rigid *last_rigid = NULL;
	    double UT[3], OMG[3];

	    for (int iloc = 0; iloc < nloc; iloc++) {
		Point &vert = *tlist[iloc];
		Rigid *rigid = (Rigid*)(vert.mesh);

		// Re-compute the rigid body motion when reaching a new rigid body
		if (last_rigid != rigid) {
		    last_rigid = rigid;
		    int nvert = rigid->numVerts();
		    double (*u)[3] = new double[nvert][3];

		    for (int ivert = 0; ivert < nvert; ivert++) {
			FOR_J3 u[ivert][j] = rigid->verts[ivert].g[j];
		    }

		    rigid->calcTransRotatVel(u, UT, OMG);
		    delete [] u;
		}

		double xx[3], omgxx[3];
		FOR_I3 xx[i] = vert.x[i] - rigid->center[i];
		cross_product(OMG, xx, omgxx);

		FOR_J3 {
		    int idx = iloc*3 + j;
		    vloc[idx] += UT[j] + omgxx[j];
		    FOR_K3 vloc[idx] += 0.5*vert.C[j][k]*vert.g[k];
		    vloc[idx] += vlin[j] - vbkg[j];
	        }
	    } // iloc

	    setRowIndices(tlist, blk_offset[iblk], idx);
	    VecSetValues(rhs_vec, 3*nloc, idx, vloc, ADD_VALUES);

	    delete [] vloc;
	    delete [] idx;
	} // rigid

	{ // walls
	    int iblk = WALL;

	    vector<Point*> &tlist = vertList[iblk];
	    int nloc = tlist.size();
	    int *idx = new int[3*nloc];
	    double *vloc = new double[3*nloc];

	    m_dclear(3*nloc, vloc);

	    for (int iloc = 0; iloc < nloc; iloc++){
		FOR_J3 {
		    int idx = iloc*3 + j;
		    vloc[idx] += vlin[j] - vbkg[j];
	        }
	    }

	    setRowIndices(tlist, blk_offset[iblk], idx);
	    VecSetValues(rhs_vec, 3*nloc, idx, vloc, ADD_VALUES);

	    delete [] vloc;
	    delete [] idx;
	} // wall
    }  // extra terms

    VecAssemblyBegin(rhs_vec);
    VecAssemblyEnd(rhs_vec);

    // Error norm
    VecNorm(rhs_vec, NORM_1, &err_L1);		//err_L1 /= nrow;
    VecNorm(rhs_vec, NORM_2, &err_L2);		//err_L2 /= sqrt(1.0*nrow);
    VecNorm(rhs_vec, NORM_INFINITY, &err_Linf);

    if (mpi_rank == 0){
	printf("    Verify solution\n");
	printf("\tRHS norm = %9.2E %9.2E %9.2E\n", rhs_L1, rhs_L2, rhs_Linf);
	printf("\tERR norm = %9.2E %9.2E %9.2E\n", err_L1, err_L2, err_Linf);
	printf("\tRelative error = %9.2E %9.2E %9.2E\n",
			err_L1/rhs_L1, err_L2/rhs_L2, err_Linf/rhs_Linf);
    }
} */


/* Calculate lhs sub-block
 * Arguments:
 *   (iblk, jblk) -- block identifier
 * Note:
 *  1. The neighbor list must be constructed. */
void ctx::calcLhs()
{
    if (! ewald::phys::active) return;	// can't form explict matrix for Fourier sum

    ewald::phys::startTiming();

    for (int iblk = 0; iblk < NBLK; iblk++)
    for (int jblk = 0; jblk < NBLK; jblk++) {
	// Wall-wall interactions are pre-calcualted
//	if (iblk == WALL && jblk == WALL) continue;

	NbrList &nlist = stks->nlist_phys[iblk][jblk];

        int nrow, ncol;
        nrow = 3*nlist.verts.size();
        ncol = 3*stks->vertList[jblk].size();

        Mat &SL = matSL[iblk][jblk];
        Mat &DL = matDL[iblk][jblk];
        if (SL) MatDestroy(&SL);
        if (DL) MatDestroy(&DL);
	//cout << "lh1\n";
	MatCreate(PETSC_COMM_SELF, &SL);
	MatSetType(SL, MATSEQAIJ);
	MatSetSizes(SL, PETSC_DECIDE, PETSC_DECIDE, nrow, ncol);
	//cout << "lh2\n";
	MatCreate(PETSC_COMM_SELF, &DL);
	MatSetType(DL, MATSEQAIJ);
	MatSetSizes(DL, PETSC_DECIDE, PETSC_DECIDE, nrow, ncol);
	//cout << "lh3\n";
	// Trivial case for zero size matrix

	if (nrow == 0 || ncol == 0) {
	    //This causes a segfault in PETSC 3.4.  This may be correct for PETSC 3.2 and below
/*  	MatZeroEntries(SL);
	    MatAssemblyBegin(SL, MAT_FINAL_ASSEMBLY);
	    MatAssemblyEnd(SL, MAT_FINAL_ASSEMBLY);
	
	    MatZeroEntries(DL);
	    MatAssemblyBegin(DL, MAT_FINAL_ASSEMBLY);
	    MatAssemblyEnd(DL, MAT_FINAL_ASSEMBLY);
	    */
	} else {
	    switch (jblk) {
		case CELL:
		    ewald::phys::calcSurfIntMat(nlist, SL, DL);
		    break;

		case RIGID:
		    ewald::phys::calcSurfIntMat(nlist, SL, DL);
		    break;
		    
		case WALL:
		    ewald::phys::calcSurfIntMat(nlist, SL, NULL);
		    break;
	    }
        }
    }

    ewald::phys::addWallTime();
}


// Calculate rhs
void ctx::calcRhs()
{
    int nvert_tot = stks->vertList[0].size() 
                  + stks->vertList[1].size()
		  + stks->vertList[2].size();
    MArray<double,2> f(nvert_tot,3);
    f = 0.0;

    // Cell surface force
    int offset = blk_offset[CELL];

    for (int icell = 0; icell < stks->numCells(); icell++) {
        Vesicle &cell = stks->cells[icell];
        cell.vesaddtension();
	if ( ! cell.isActive) continue;

	for (int ivert = 0; ivert < cell.numVerts(); ivert++) {
	    Point &vert = cell.verts[ivert];
	    m_dcopy(3, &cell.f(ivert,0), vert.f);
	    m_dcopy(3, &cell.f(ivert,0), f.data()+offset+3*vert.Gindx);
	}
    }

    // Rigid repulsion force
    offset = blk_offset[RIGID];

    for (int irigid = 0; irigid < stks->numRigids(); irigid++) {
        Rigid &rigid = stks->rigids[irigid];
	if ( ! rigid.isActive) continue;

	for (int ivert = 0; ivert < rigid.numVerts(); ivert++) {
	    Point &vert = rigid.verts[ivert];
	    m_dclear(3, vert.f);
	}
    }

	//Note to Joe: Uncomment line below if you start mixing rigids into the code.  With no rigids, there's no need to consider their collisions
    //stks->addCollisionForce(RIGID, RIGID, 0.2);

    for (int irigid = 0; irigid < stks->numRigids(); irigid++) {
        Rigid &rigid = stks->rigids[irigid];
	if ( ! rigid.isActive) continue;

	for (int ivert = 0; ivert < rigid.numVerts(); ivert++) {
	    Point &vert = rigid.verts[ivert];
	    m_dcopy(3, vert.f, f.data()+offset+3*vert.Gindx);
	}
    }

    //======================================================================
    // Init
    VecZeroEntries(rhs_vec);

    // Physical sum
    if (ewald::phys::active) {
        using ewald::phys::comm;
        using ewald::phys::comm_size;
        using ewald::phys::comm_rank;

        // Add interactions
        for (int iblk = 0; iblk < NBLK; iblk++)
        for (int jblk = 0; jblk < NBLK; jblk++) {
	    if (jblk == WALL) continue;

            NbrList &nlist = stks->nlist_phys[iblk][jblk];
	    vector<Point*> &tlist = nlist.verts;
	    int n = tlist.size();

	    MArray<double,2> v(n,3);
	    v = 0.0;
	    myMatMultAdd(-1.0/(8*M_PI), matSL[iblk][jblk], f.data(), v.data());

	    MArray<int,2> idx(n,3);
	    for (int i = 0; i < n; i++) {
	        int p = blk_offset[iblk] + 3*tlist[i]->Gindx;
	        FOR_J3 idx(i,j) = p++;
	    }

	    VecSetValues(rhs_vec, 3*n, idx.data(), v.data(), ADD_VALUES);
        }

        // Add rigid self interaction
	{ 
            NbrList &nlist = stks->nlist_RigidSelf;
            vector<Point*> &tlist = nlist.verts;
	    int n = tlist.size();

	    MArray<double,2> v(n,3);
	    v = 0.0;
	    addRigidSelfTerm(-1.0/(8*M_PI), 0.0, v.data());

	    MArray<int,2> idx(n,3);
	    for (int i = 0; i < n; i++) {
	        int p = blk_offset[RIGID] + 3*tlist[i]->Gindx;
	        FOR_J3 idx(i,j) = p++;
	    }

	    VecSetValues(rhs_vec, 3*n, idx.data(), v.data(), ADD_VALUES);
        }
    }  // phys

    // Fourier sum
    if (ewald::four::active) {
	ewald::four::clear_source();
	ewald::four::add_source(-1.0/(8*M_PI), 0.0, stks->slist_four[CELL]);
	ewald::four::add_source(-1.0/(8*M_PI), 0.0, stks->slist_four[RIGID]);

	ewald::four::transform();

        for (int iblk = 0; iblk < NBLK; iblk++) {
	    vector<Point*> &tlist = stks->tlist_four[iblk];
	    int n = tlist.size();

	    MArray<double,2> v(n,3);
	    v = 0.0;
	    ewald::four::add_interp_vel(tlist, v.data());

	    MArray<int,2> idx(n,3);
	    for (int i = 0; i < n; i++) {
	        int p = blk_offset[iblk] + 3*tlist[i]->Gindx;
	        FOR_J3 idx(i,j) = p++;
	    }

	    VecSetValues(rhs_vec, 3*n, idx.data(), v.data(), ADD_VALUES);
	}
    } // ewald::four

    // Background velocity
    for (int blk = 0; blk < NBLK; blk++)
    for (int imesh = 0; imesh < stks->meshList[blk].size(); imesh++) {
	Mesh &mesh = *stks->meshList[blk][imesh];
	if ( ! mesh.isPrivate) continue;

	int n = mesh.numVerts();

	MArray<double,2> v(n,3);
	v = 0.0;
	for (int i = 0; i < n; i++) {
	    Point &vert = mesh.verts[i];
	    if (vert.pp) continue;	// Avoid dupicate points at periodic boundaries
	    FOR_J3 v(i,j) += stks->vbkg[j];
	}

	MArray<int,2> idx(n,3);
	for (int i = 0; i < n; i++) {
	    int p = blk_offset[blk] + 3*mesh.verts[i].Gindx;
	    FOR_J3 idx(i,j) = p++;
	}

	VecSetValues(rhs_vec, 3*n, idx.data(), v.data(), ADD_VALUES);
    }

    VecAssemblyBegin(rhs_vec);
    VecAssemblyEnd(rhs_vec);
}


/* Matrix vector product for the wall solver
 * Arguments:
 *  v_vec = lhs * src_vec */
PetscErrorCode ctx::matmult(Mat lhs, Vec src_vec, Vec v_vec)
{
    // Each processor gets a copy of the global src array
    int nrow;
    VecGetSize(src_vec, &nrow);

    int nvert_tot = stks->vertList[0].size() 
                  + stks->vertList[1].size()
                  + stks->vertList[2].size();
    MArray<double,2> f(nvert_tot,3), g(nvert_tot,3);
    setSourceFromVec(src_vec, f.data(), g.data());

    // Init
    VecZeroEntries(v_vec);

    if (ewald::phys::active) {
        using ewald::phys::comm;
        using ewald::phys::comm_size;
        using ewald::phys::comm_rank;

	ewald::phys::startTiming();

	for (int iblk = 0; iblk < NBLK; iblk++)
	for (int jblk = 0; jblk < NBLK; jblk++) {
	    NbrList &nlist = stks->nlist_phys[iblk][jblk];
	    vector<Point*> &tlist = nlist.verts;
	    int n = tlist.size();

	    MArray<double,2> v(n,3);

	    v = 0.0;

	    myMatMultAdd(coeffSL[jblk], matSL[iblk][jblk], f.data()+blk_offset[jblk], v.data());
	    myMatMultAdd(coeffDL[jblk], matDL[iblk][jblk], g.data()+blk_offset[jblk], v.data());

	    MArray<int,2> idx(n,3);
	    for (int i = 0; i < n; i++) {
	        int p = blk_offset[iblk] + 3*tlist[i]->Gindx;
	        FOR_J3 idx(i,j) = p++;
	    }

	    VecSetValues(v_vec, 3*n, idx.data(), v.data(), ADD_VALUES);
	} // iblk, jblk

	// Add rigid particle self interaction
	{ 
	    NbrList &nlist = stks->nlist_RigidSelf;
	    vector<Point*> &tlist = nlist.verts;

	    int n = tlist.size();
	    MArray<double,2> v(n,3);

	    v = 0.0;
	    addRigidSelfTerm(coeffSL[RIGID], coeffDL[RIGID], v.data());

	    MArray<int,2> idx(n,3);
	    for (int i = 0; i < n; i++) {
	        int p = blk_offset[RIGID] + 3*tlist[i]->Gindx;
	        FOR_J3 idx(i,j) = p++;
	    }

	    VecSetValues(v_vec, 3*n, idx.data(), v.data(), ADD_VALUES);
	}

	ewald::phys::addWallTime();
    } // phys

    if (ewald::four::active) {
        ewald::four::startTiming();

        ewald::four::clear_source();
	for (int blk = 0; blk < NBLK; blk++) {
	    ewald::four::add_source(coeffSL[blk], coeffDL[blk], stks->slist_four[blk]);
	}

	ewald::four::transform();

	for (int blk = 0; blk < NBLK; blk++) {
	    vector<Point*> &tlist = stks->tlist_four[blk];
	    int n = tlist.size();

	    MArray<double,2> v(n,3);
	    v = 0.0;
	    ewald::four::add_interp_vel(tlist, v.data());

	    MArray<int,2> idx(n,3);
	    for (int i = 0; i < n; i++) {
		int p = blk_offset[blk] + 3*tlist[i]->Gindx;
	        FOR_J3 idx(i,j) = p++;
	    }

	    VecSetValues(v_vec, 3*n, idx.data(), v.data(), ADD_VALUES);
	} // blk

        ewald::four::addWallTime();
    }


    // Diagonal term
    // Double-layer jump 
    // Double-layer linear term
    // Deflation term (for RIGID only)
    {
        // Compute the linear term in double layer potential
	double vlin[3];
	double my_vlin[3] = { 0.0 };

	for (int blk = 0; blk < NBLK; blk++)  {
	    if ( fabs(coeffDL[blk]) < 1.E-10 ) continue;

	    for (int imesh = 0; imesh < stks->meshList[blk].size(); imesh++) {
		Mesh *mesh = stks->meshList[blk][imesh];
		if ( ! mesh->isPrivate ) continue;

		double tmp[3];
		mesh->calcDoubleLayerLinearTerm(tmp);
		FOR_I3 my_vlin[i] += coeffDL[blk]*tmp[i];
	    }
	}

	MPI_Allreduce(my_vlin, vlin, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	// Cells
	for (int icell = 0; icell < stks->numCells(); icell++) {
	    Vesicle &cell = stks->cells[icell];
	    if ( ! cell.isPrivate ) continue;

	    int nvert = cell.numVerts();

	    MArray<double,2> v(nvert,3);
	    v = 0.0;

	    for (int ivert = 0; ivert < nvert; ivert++) {
	        Point &vert = cell.verts[ivert];
		if (vert.pp) continue;	// no duplicity

		double Cg[3] = {0.0};
		FOR_J3 FOR_K3 Cg[j] += vert.C[j][k]*vert.g[k];

		double s = 4*M_PI*coeffDL[CELL];
	        FOR_J3 v(ivert,j) += vert.g[j] + s*Cg[j] + vlin[j];
	    }

	    MArray<int,2> idx(nvert,3);
	    for (int ivert = 0; ivert < nvert; ivert++) {
	        int p = blk_offset[CELL] + 3*cell.verts[ivert].Gindx;
	        FOR_J3 idx(ivert,j) = p++;
	    }

	    VecSetValues(v_vec, 3*nvert, idx.data(), v.data(), ADD_VALUES);
	}


	// Rigids
	for (int irigid = 0; irigid < stks->numRigids(); irigid++) {
	    Rigid &rigid = stks->rigids[irigid];
	    if ( ! rigid.isPrivate ) continue;

	    int n = rigid.numVerts();

	    double UT[3], OMG[3];
	    MArray<double,2> gtmp(n,3);
	    for (int i = 0; i < n; i++) {
	        Point &vert = rigid.verts[i];
		FOR_J3 gtmp(i,j) = vert.g[j];
	    }
	    rigid.calcTransRotatVel(gtmp, UT, OMG);

	    MArray<double,2> v(n,3);
	    v = 0.0;

	    for (int i = 0; i < n; i++) {
	        Point &vert = rigid.verts[i];
		if (vert.pp) continue;	// no duplicity

		double xx[3], omgxx[3];
		FOR_J3 xx[j] = vert.x[j] - rigid.center[j];
	        cross_product(OMG, xx, omgxx);

		double s = 4*M_PI*coeffDL[RIGID];
		double Cg[3] = {0.0};
		FOR_J3 FOR_K3 Cg[j] += vert.C[j][k]*vert.g[k];

	        FOR_J3 v(i,j) += s*Cg[j] + vlin[j] + UT[j] + omgxx[j];
	    }

	    MArray<int,2> idx(n,3);
	    for (int i = 0; i < n; i++) {
	        int p = blk_offset[RIGID] + 3*rigid.verts[i].Gindx;
	        FOR_J3 idx(i,j) = p++;
            }

	    VecSetValues(v_vec, 3*n, idx.data(), v.data(), ADD_VALUES);
	}


	// Wall
	for (int iwall = 0; iwall < stks->numWalls(); iwall++) {
	    Wall &wall = stks->walls[iwall];
	    if ( ! wall.isPrivate ) continue;
	    int n = wall.numVerts();

	    MArray<double,2> v(n,3);
	    v = 0.0;
	    for (int i = 0; i < n; i++) {
	        Point &vert = wall.verts[i];
		if (vert.pp) continue;	// no duplicity

	        FOR_J3 v(i,j) += vlin[j];
	    }

	    MArray<int,2> idx(n,3);
	    for (int i = 0; i < n; i++) {
	        int p = blk_offset[WALL] + 3*wall.verts[i].Gindx;
	        FOR_J3 idx(i,j) = p++;
            }

	    VecSetValues(v_vec, 3*n, idx.data(), v.data(), ADD_VALUES);
        }
    } 

    VecAssemblyBegin(v_vec); 
    VecAssemblyEnd(v_vec); 

    return 0;
}


/* Calculate the self-interaction term of rigid body, assuming that
 * the sources are stored by vert.f[] and vert.g[] arrays
 * Argument:
 *   c1, c2 -- single- and double-layer coefficients
 *   v -- the integral */
void ctx::addRigidSelfTerm(double c1, double c2, double *v)
{
    bool c1_flag = (fabs(c1) > 1.E-10);
    bool c2_flag = (fabs(c2) > 1.E-10);

    NbrList &nlist = stks->nlist_RigidSelf;
    int nrigid = stks->numRigids();
    int nrow = 3*nlist.verts.size();
    int ncol = 3*stks->vertList[RIGID].size();

    // Trivial case
    if (nrow == 0 || ncol == 0) return;

    double (*T)[3][3] = new double[nrigid][3][3];
    double *f = new double[ncol];
    double *g = new double[ncol];
    double *vref = new double[nrow];

    // Compute the transform matrix of every particle
    for (int irigid = 0; irigid < nrigid; irigid++) {
        Rigid &rigid = stks->rigids[irigid];
	Tri &face = rigid.faces[0];

	double x[3][3];
	FOR_I3 FOR_J3 x[i][j] = face.vert[i]->x[j];
	calcRotateMatrix(rigid_xref[irigid], x, T[irigid]);
    }

    // Set up source term
    for (int irigid = 0; irigid < stks->numRigids(); irigid++) {
        Rigid &rigid = stks->rigids[irigid];

	// Convert the density to reference frame, i.e. multiply by transpose of T
        for (int ivert = 0; ivert < rigid.numVerts(); ivert++) {
	    Point &vert = rigid.verts[ivert];
	    int p = vert.Gindx;

	    if (c1_flag) {
	        FOR_I3 f[3*p+i] = T[irigid][0][i]*vert.f[0]
		                + T[irigid][1][i]*vert.f[1]
				+ T[irigid][2][i]*vert.f[2];
            }

	    if (c2_flag) {
	        FOR_I3 g[3*p+i] = T[irigid][0][i]*vert.g[0]
		                + T[irigid][1][i]*vert.g[1]
				+ T[irigid][2][i]*vert.g[2];
            }
	}
    } // irigid


    // Calculate the matmult in reference frame
    m_dclear(nrow, vref);
    if (c1_flag) myMatMultAdd(c1, matSL_RigidSelf, f, vref);
    if (c2_flag) myMatMultAdd(c2, matDL_RigidSelf, g, vref);

    // Trasform the matmult result to the current frame
    for (int p = 0; p < nlist.verts.size(); p++) {
        Point *vert = nlist.verts[p];
        int irigid = vert->mesh->Gindx;
	FOR_I3 v[3*p+i] += T[irigid][i][0]*vref[3*p+0]
	                 + T[irigid][i][1]*vref[3*p+1]
	                 + T[irigid][i][2]*vref[3*p+2];
    }

    delete [] f;
    delete [] g;
    delete [] vref;
    delete [] T;
}



/* Set up source from a distributed vector
 * Argument:
 *   src_vec -- PETSc distributed vector 
 *   f -- single-layer density
 *   g -- double-layer density */
void ctx::setSourceFromVec(Vec src_vec, double *f, double *g)
{
    int nrow;
    VecGetSize(src_vec, &nrow);

    // Init
    m_dclear(nrow, f);
    m_dclear(nrow, g);

    // src holds the local copy of the global src_vec vector
    double *src = new double[nrow];
    myVecScatter(src_scatter, src_vec, src, INSERT_VALUES, SCATTER_FORWARD);

    // Cell surface source
    int offset = blk_offset[CELL];
    
    for (int icell= 0; icell < stks->numCells(); icell++) {
        Vesicle &cell = stks->cells[icell];
	int nvert = cell.numVerts();

	if (cell.isActive) {
	    MArray<double,2> dx(nvert,3), df(nvert,3);

	    for (int ivert = 0; ivert < nvert; ivert++) {
		Point &vert = cell.verts[ivert];

		m_dcopy(3, src+offset+3*vert.Gindx, g+offset+3*vert.Gindx);
		m_dcopy(3, src+offset+3*vert.Gindx, vert.g);
		FOR_D3 dx(ivert,d) = stks->Ts*vert.g[d];
	    }

	    cell.diffBendForce(dx, df);

	    for (int ivert = 0; ivert < nvert; ivert++) {
		Point &vert = cell.verts[ivert];
		FOR_D3 vert.f[d] = df(ivert,d);
		m_dcopy(3, vert.f, f+offset+3*vert.Gindx);
	    } 
        } else {
	    for (int ivert = 0; ivert < nvert; ivert++) {
		Point &vert = cell.verts[ivert];
		FOR_J3 vert.f[j] = 0.0;
		FOR_J3 vert.g[j] = 0.0;
	    }
	}
    }

    // Rigid double layer
    offset = blk_offset[RIGID];

    for (int irigid = 0; irigid < stks->numRigids(); irigid++) {
        Rigid &rigid = stks->rigids[irigid];
	
	for (int ivert = 0; ivert < rigid.numVerts(); ivert++) {
	    Point &vert = rigid.verts[ivert];

	    m_dcopy(3, src+offset+3*vert.Gindx, g+offset+3*vert.Gindx);
	    m_dcopy(3, src+offset+3*vert.Gindx, vert.g);
	}
    }

    // Wall single layer
    offset = blk_offset[WALL];

    for (int iwall = 0; iwall < stks->numWalls(); iwall++) {
        Wall &wall = stks->walls[iwall];

	for (int ivert = 0; ivert < wall.numVerts(); ivert++) {
	    Point &vert = wall.verts[ivert];

	    m_dcopy(3, src+offset+3*vert.Gindx, f+offset+3*vert.Gindx);
	    m_dcopy(3, src+offset+3*vert.Gindx, vert.f);
	}
    }

    delete [] src;
}
