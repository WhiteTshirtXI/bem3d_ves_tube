#include "cxxheaders.h"
#include "shearsys.h"
#include "ewald.h"
#include "mathfunc.h"
#include "mypetsc.h"
#include "debugfunc.h"
#include "param.h"

/* ShearSys::solveFlow
 * ShearSys::initFlowSolver
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
    ShearSys *stks;
    bool first_solve = true;

    const int NBLK = 2;
    int blk_size[NBLK], blk_offset[NBLK];

    Vec rhs_vec, sol_vec;
    Mat lhs, matSL[NBLK][NBLK], matDL[NBLK][NBLK];
    double coeffSL[NBLK], coeffDL[NBLK];

    // rigid_xref -- coordiantes of the 1st triangle of every rigid particle
    Mat matSL_RigidSelf, matDL_RigidSelf;
    double (*rigid_xref)[3][3];

    KSP ksp;
    VecScatter src_scatter;

    PetscErrorCode matmult(Mat, Vec, Vec);

    void calcLhs();
    void calcLhs_RigidSelf();
    void calcRhs();

    void addExtraTerm_Cell(vector<Point*> &, double *);
    void addExtraTerm_Rigid(vector<Point*> &, double *);

    void addRigidSelfTerm(double, double, double *);

    void setSourceFromVec(Vec, double *f, double *g);
}};


void ShearSys::solveFlow()
{
    using namespace ctx;

    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);


    //======================================================================
    // Single- and double-layer coefficients in the lhs
    ctx::coeffSL[CELL] = 1.0/(8*M_PI);
    ctx::coeffDL[CELL] = (stks->viscRat - 1.0)/(8*M_PI);

    ctx::coeffSL[RIGID] = 0.0;
    ctx::coeffDL[RIGID] = 1.0/(8*M_PI);

    //======================================================================
    ewald::phys::resetTiming();

    ctx::calcLhs();

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

    calcRhs();

    if (ewald::timing) {
        double wtime_phys, wtime_four;
	ewald::getWallTime(wtime_phys, wtime_four);

	if (mpi_rank == 0) {
	    printf("    calc rhs, wtime = %7.2f %7.2f s\n", wtime_phys, wtime_four);
	}
    }


    //======================================================================
    // Solve the equation
    KSPSetTolerances(ksp, 1.E-3, PETSC_DEFAULT, PETSC_DEFAULT, 100);

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
        Cell &cell = stks->cells[icell];
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
}


void ShearSys::initFlowSolver()
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
    nrow = blk_size[0] + blk_size[1];
    VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, nrow, &rhs_vec);
    VecDuplicate(rhs_vec, &sol_vec);

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

    // Rigid self interaction matrix
    if (ewald::phys::active) {
	// Build the source list and neighbor list for rigid self-interaction
	SrcList &slist = stks->slist_RigidSelf;
	NbrList &nlist = stks->nlist_RigidSelf;

	slist.build(stks->faceList[RIGID]);
	for (int i = 0; i < slist.numFaces(); i++) {
	    slist.faces[i]->updateGeometry();
	}

	// NO_SELF = false, ONLY_SELF = true
	nlist.build(stks->vertList[RIGID], slist, false, true);

	// Compute the matrix
	int nrow = 3*nlist.verts.size();
	int ncol = 3*stks->vertList[RIGID].size();

	MatCreate(PETSC_COMM_SELF, &matSL_RigidSelf);
	MatSetType(matSL_RigidSelf, MATSEQAIJ);
	MatSetSizes(matSL_RigidSelf, PETSC_DECIDE, PETSC_DECIDE, nrow, ncol);

	MatCreate(PETSC_COMM_SELF, &matDL_RigidSelf);
	MatSetType(matDL_RigidSelf, MATSEQAIJ);
	MatSetSizes(matDL_RigidSelf, PETSC_DECIDE, PETSC_DECIDE, nrow, ncol);

	ewald::phys::calcSurfIntMat(nlist, matSL_RigidSelf, matDL_RigidSelf);

	// Save the coordinates of the 1st triangle of every particle to perform
	// coordinate transform in the future
	int nrigid = stks->numRigids();
	rigid_xref = new double[nrigid][3][3];

	for (int irigid = 0; irigid < nrigid; irigid++) {
	    Rigid &rigid = stks->rigids[irigid];
	    Tri &face = rigid.faces[0];

	    FOR_I3 FOR_J3 rigid_xref[irigid][i][j] = face.vert[i]->x[j];
	}
    }
}


/* Calculate Stokes interaction matrices */
void ctx::calcLhs()
{
    if (! ewald::phys::active) return;

    ewald::phys::startTiming();

    for (int iblk = 0; iblk < NBLK; iblk++)
    for (int jblk = 0; jblk < NBLK; jblk++) {
	NbrList &nlist = stks->nlist_phys[iblk][jblk];

	int nrow, ncol;
	nrow = 3*nlist.verts.size();
	ncol = 3*stks->vertList[jblk].size();

	Mat &SL = matSL[iblk][jblk];
	Mat &DL = matDL[iblk][jblk];
	if (SL) MatDestroy(&SL);
	if (DL) MatDestroy(&DL);

	MatCreate(PETSC_COMM_SELF, &SL);
	MatSetType(SL, MATSEQAIJ);
	MatSetSizes(SL, PETSC_DECIDE, PETSC_DECIDE, nrow, ncol);

	MatCreate(PETSC_COMM_SELF, &DL);
	MatSetType(DL, MATSEQAIJ);
	MatSetSizes(DL, PETSC_DECIDE, PETSC_DECIDE, nrow, ncol);

	if (nrow == 0 || ncol == 0) {
	    MatZeroEntries(SL);
	    MatAssemblyBegin(SL, MAT_FINAL_ASSEMBLY);
	    MatAssemblyEnd(SL, MAT_FINAL_ASSEMBLY);
	
	    MatZeroEntries(DL);
	    MatAssemblyBegin(DL, MAT_FINAL_ASSEMBLY);
	    MatAssemblyEnd(DL, MAT_FINAL_ASSEMBLY);
	}
	else {
	    ewald::phys::calcSurfIntMat(nlist, SL, DL);
	}

    }
    ewald::phys::addWallTime();
}


// Calculate rhs
void ctx::calcRhs()
{
    int nvert_tot = stks->vertList[0].size() + stks->vertList[1].size();
    MArray<double,2> f(nvert_tot,3);
    f = 0.0;

    // Cell surface force
    int offset = blk_offset[CELL];

    for (int icell = 0; icell < stks->numCells(); icell++) {
        Cell &cell = stks->cells[icell];
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

    stks->addCollisionForce(RIGID, RIGID, 0.2);

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

    // Back ground velocity
    for (int blk = 0; blk < NBLK; blk++) 
    for (int imesh = 0; imesh < stks->meshList[blk].size(); imesh++) {
	Mesh &mesh = *stks->meshList[blk][imesh];
	if ( ! mesh.isPrivate) continue;

	int n = mesh.numVerts();

	MArray<double,2> v(n,3);
	v = 0.0;
	for (int i = 0; i < n; i++) {
	    Point &vert = mesh.verts[i];

	    double vbkg[3];
	    stks->calcBkgVel(vert.x, vbkg);
	    FOR_J3 v(i,j) += vbkg[j];
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
/*
    VecZeroEntries(rhs_vec);

    // Back ground velocity
    for (int blk = 0; blk < NBLK; blk++) 
    for (int imesh = 0; imesh < stks->meshList[blk].size(); imesh++) {
	Mesh &mesh = *stks->meshList[blk][imesh];
	if ( ! mesh.isPrivate) continue;

	int n = mesh.numVerts();

	MArray<double,2> v(n,3);
	v = 0.0;
	for (int i = 0; i < n; i++) {
	    Point &vert = mesh.verts[i];

	    double vbkg[3];
	    stks->calcBkgVel(vert.x, vbkg);
	    FOR_J3 v(i,j) += vbkg[j];
	}

	MArray<int,2> idx(n,3);
	for (int i = 0; i < n; i++) {
	    double p = blk_offset[blk] + 3*mesh.verts[i].Gindx;
	    FOR_J3 idx(i,j) = p++;
	}

	VecSetValues(rhs_vec, 3*n, idx.data(), v.data(), ADD_VALUES);
    }

    VecAssemblyBegin(rhs_vec);
    VecAssemblyEnd(rhs_vec);
*/
}


/* Matrix vector product for the solver
 * Arguments:
 *  v_vec = lhs * src_vec */
PetscErrorCode ctx::matmult(Mat lhs, Vec src_vec, Vec v_vec)
{
    // Each processor gets a copy of global source arrays
    int nrow;
    VecGetSize(src_vec, &nrow);

    int nvert_tot = stks->vertList[0].size() + stks->vertList[1].size();
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
    // Double-layer liner term
    // Deflation term (for RIGID only)
    {
        // Compute the linear term in double layer potential
	double vlin[3];
	double my_vlin[3] = { 0.0 };

	for (int blk = 0; blk < NBLK; blk++) 
	for (int imesh = 0; imesh < stks->meshList[blk].size(); imesh++) {
	    Mesh *mesh = stks->meshList[blk][imesh];
	    if ( ! mesh->isPrivate ) continue;

	    double tmp[3];
	    mesh->calcDoubleLayerLinearTerm(tmp);
	    FOR_I3 my_vlin[i] += coeffDL[blk]*tmp[i];
	}
	MPI_Allreduce(my_vlin, vlin, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	// Cells
	for (int icell = 0; icell < stks->numCells(); icell++) {
	    Cell &cell = stks->cells[icell];
	    if ( ! cell.isPrivate ) continue;

	    int nvert = cell.numVerts();

	    MArray<double,2> v(nvert,3);
	    v = 0.0;
	    for (int ivert = 0; ivert < nvert; ivert++) {
	        Point &vert = cell.verts[ivert];

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
    } 

    VecAssemblyBegin(v_vec); 
    VecAssemblyEnd(v_vec); 

    return 0;
}


/* Add the Ewald physical sum due to rigid self interaction
 *   v --> v + c1 * N(f) + c2 * K(g)
 *   where N(f) is single layer and K(g) is double layer */
void ctx::addRigidSelfTerm(double c1, double c2, double *v)
{
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

    // Compute the transform matrix of every particle and set up source term
    for (int irigid = 0; irigid < nrigid; irigid++) {
        Rigid &rigid = stks->rigids[irigid];
	if ( ! rigid.isActive) continue;

	Tri &face = rigid.faces[0];
	double x[3][3];
	FOR_I3 FOR_J3 x[i][j] = face.vert[i]->x[j];
	calcRotateMatrix(rigid_xref[irigid], x, T[irigid]);

	// Convert the density to reference frame, i.e. multiply by the transpose of T
        for (int ivert = 0; ivert < rigid.numVerts(); ivert++) {
	    Point &vert = rigid.verts[ivert];
	    int p = vert.Gindx;

	    FOR_I3 f[3*p+i] = T[irigid][0][i]*vert.f[0]
			    + T[irigid][1][i]*vert.f[1]
			    + T[irigid][2][i]*vert.f[2];

	    FOR_I3 g[3*p+i] = T[irigid][0][i]*vert.g[0]
			    + T[irigid][1][i]*vert.g[1]
			    + T[irigid][2][i]*vert.g[2];
	}
    } // irigid


    // Calculate the matmult in reference frame
    m_dclear(nrow, vref);
    myMatMultAdd(c1, matSL_RigidSelf, f, vref);
    myMatMultAdd(c2, matDL_RigidSelf, g, vref);

    // Trasform the matmult result to the current frame
    for (int p = 0; p < nlist.verts.size(); p++) {
        Point *vert = nlist.verts[p];
	if ( ! vert->mesh->isActive) continue;	// Very important

        int irigid = vert->mesh->Gindx;
	FOR_I3 v[3*p+i] += T[irigid][i][0]*vref[3*p+0]
	                 + T[irigid][i][1]*vref[3*p+1]
	                 + T[irigid][i][2]*vref[3*p+2];
    }

    delete [] T;
    delete [] f;
    delete [] g;
    delete [] vref;
}



/* Set up source from a distributed vector
 * Argument:
 *   src_vec -- PETSc vector 
 *   f -- single-layer density
 *   g -- double-layer density */
void ctx::setSourceFromVec(Vec src_vec, double *f, double *g)
{
    int nrow;
    VecGetSize(src_vec, &nrow);

    // Init
    m_dclear(nrow, f);
    m_dclear(nrow, g);

    // Src holds the local copy of the global src_vec vector
    myVecScatter(src_scatter, src_vec, g, INSERT_VALUES, SCATTER_FORWARD);

    // Cell surface source
    int offset = blk_offset[CELL];

    for (int icell= 0; icell < stks->numCells(); icell++) {
        Cell &cell = stks->cells[icell];
	int nvert = cell.numVerts();

	if (cell.isActive) {
	    MArray<double,2> dx(nvert,3), df(nvert,3);

	    for (int ivert = 0; ivert < nvert; ivert++) {
		Point &vert = cell.verts[ivert];

		m_dcopy(3, g+offset+3*vert.Gindx, vert.g);
		FOR_D3 dx(ivert,d) = stks->Ts*vert.g[d];
	    }

	    cell.diffForce(dx, df);

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

    // Rigid surace source
    offset = blk_offset[RIGID];

    for (int irigid = 0; irigid < stks->numRigids(); irigid++) {
        Rigid &rigid = stks->rigids[irigid];

	if (rigid.isActive) {
	    for (int ivert = 0; ivert < rigid.numVerts(); ivert++) {
		Point &vert = rigid.verts[ivert];
		m_dcopy(3, g+offset+3*vert.Gindx, vert.g);
	    }
        } else {
	    for (int ivert = 0; ivert < rigid.numVerts(); ivert++) {
		Point &vert = rigid.verts[ivert];
		FOR_J3 vert.f[j] = 0.0;
		FOR_J3 vert.g[j] = 0.0;
	    }
	}
    }
}
