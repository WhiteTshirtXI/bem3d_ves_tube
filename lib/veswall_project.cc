// Projection step of the flow solver
#include "cxxheaders.h"
#include "veswall.h"

/* Context for the matrix shell
 * Note: PETSc context only supports integer type, but we 
 *       need pointers, so use a local namespace  */
namespace { 
namespace solver {
    bool INITED = false;

    VesWall *stks;

    Vec rhs_vec, sol_vec;
    Mat lhs;
    KSP ksp;
    PC pc;
    VecScatter scatter;

    PetscErrorCode matmult(Mat, Vec, Vec);
    void lbd2vel(const MArray<double,1> &, MArray<double,2> &);
}};


/* VesWall::projectVel
 * VesWall::initProjectSolver
 * solver::matmult
 * solver::lbd2vel */

/* Project the surface velocity to make it divergence free
 * Arguments:
 *  vold -- surface velocity before projection
 *  divTar -- target surface velocity divergence
 *  vnew -- surface velocity after projection
 *  lbd -- the Lagrange multiplier 
 *  initial_guess -- whether use lbd as initial guess */
void VesWall::projectVel(const DoubArray2d &vold, const DoubArray1d &divTar,
		DoubArray2d &vnew, DoubArray1d &lbd,
		bool initial_guess)
{
    using namespace solver;
    int nvert_tot = vertList[VESICLE].size();

    //======================================================================
    // rhs = divTar - div
    // Init
    VecZeroEntries(rhs_vec); 

    for (int ives = 0; ives < numVesicles(); ives++) {
        Vesicle &vesicle = vesicles[ives];
	if ( ! vesicle.isPrivate ) continue;

	int nvert = vesicle.numVerts();

	MArray<double,2> v_tmp(nvert,3);
	MArray<double,1> div_tmp(nvert), divTar_tmp(nvert), rhs_tmp(nvert);

	int p = vesicle.verts[0].Gindx;
	m_dcopy(3*nvert, &vold(p,0), v_tmp.data());
	m_dcopy(nvert, &divTar(p), divTar_tmp.data());

	vesicle.velDiv(v_tmp, div_tmp);

	rhs_tmp = divTar_tmp;
	rhs_tmp -= div_tmp;

	MArray<int,1> idx(nvert);
	for (int ivert = 0; ivert < nvert; ivert++) {
	   idx(ivert) = vesicle.verts[ivert].Gindx;
	}

	VecSetValues(rhs_vec, nvert, idx.data(), rhs_tmp.data(), ADD_VALUES);
    }


    VecAssemblyBegin(rhs_vec);
    VecAssemblyEnd(rhs_vec);

    //======================================================================
    // Solve the equation
    if (initial_guess) {
        myVecSetValues(sol_vec, lbd.data(), INSERT_VALUES);
	VecAssemblyBegin(sol_vec);
	VecAssemblyEnd(sol_vec);

        KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    }
    else {
        KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
    }

    double atol = 0.0;
    for (int ives = 0; ives < numVesicles(); ives++) {
        Vesicle &vesicle = vesicles[ives];
	for (int ivert = 0; ivert < vesicle.numVerts(); ivert++) {
	    atol += square(1.E-3*vesicle.vertArea(ivert));
	}
    }
    atol = sqrt(atol);

    KSPSetTolerances(ksp, PETSC_DEFAULT, atol, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSolve(ksp, rhs_vec, sol_vec);

    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank == 0) {
	int niter;
	double rnorm;
	KSPGetIterationNumber(ksp, &niter);
	KSPGetResidualNorm(ksp, &rnorm);
	rnorm /= sqrt(1.0*nvert_tot);
        printf("    GMRES: niter = %d  rnorm = %.2E\n", niter, rnorm);
    }

    //======================================================================
    // Calculate the projected velocity
    myVecScatter(scatter, sol_vec, lbd.data(), INSERT_VALUES, SCATTER_FORWARD);
    MArray<double,2> dv(nvert_tot,3);
    solver::lbd2vel(lbd, dv);

    vnew = vold;
    vnew += dv;
}


void VesWall::initProjectSolver()
{
    using solver::rhs_vec;
    using solver::sol_vec;
    using solver::lhs;
    using solver::matmult;
    using solver::ksp;
    using solver::pc;
    using solver::scatter;

    solver::stks = this;
    solver::INITED = true;

    // rhs and solution vectors
    int nrow = vertList[VESICLE].size();
    VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, nrow, &rhs_vec);
    VecDuplicate(rhs_vec, &sol_vec);

    // scatter
    myVecScatterCreateToAll(sol_vec, scatter);

    // lhs matrix 
    int m;
    VecGetLocalSize(rhs_vec, &m);
    MatCreateShell(MPI_COMM_WORLD, m, m, nrow, nrow, NULL, &lhs);
    MatShellSetOperation(lhs, MATOP_MULT, (void(*)(void))matmult);

    // ksp solver
    KSPCreate(MPI_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, lhs, lhs, DIFFERENT_NONZERO_PATTERN);
    KSPSetType(ksp, KSPGMRES);

    // preconditioner 
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCNONE);
}


/* Matrix vector product for the solver
 * Arguments:
 *  div_vec = lhs * src_vec */
PetscErrorCode solver::matmult(Mat lhs, Vec src_vec, Vec div_vec)
{
    // Init
    VecZeroEntries(div_vec);

    int nvert_tot = stks->vertList[VESICLE].size();
    MArray<double,1> lbd(nvert_tot);
    myVecScatter(scatter, src_vec, lbd.data(), INSERT_VALUES, SCATTER_FORWARD);

    MArray<double,2> v(nvert_tot,3);
    lbd2vel(lbd, v);

    for (int ives = 0; ives < stks->numVesicles(); ives++) {
        Vesicle &vesicle = stks->vesicles[ives];
	if ( ! vesicle.isPrivate ) continue;

	int nvert = vesicle.numVerts();

	MArray<double,2> v_tmp(nvert,3);
	MArray<double,1> div_tmp(nvert);

	int p = vesicle.verts[0].Gindx;
	m_dcopy(3*nvert, &v(p,0), v_tmp.data());
	vesicle.velDiv(v_tmp, div_tmp);

	MArray<int,1> idx(nvert);
	for (int ivert = 0; ivert < nvert; ivert++) {
	    idx(ivert) = vesicle.verts[ivert].Gindx;
	}

	VecSetValues(div_vec, nvert, idx.data(), div_tmp.data(), ADD_VALUES);
    }

    VecAssemblyBegin(div_vec);
    VecAssemblyEnd(div_vec);

    return 0;
}


/* Calculate velocity induced by surface tension
 * Arguments:
 *  lbd -- the surface tension
 *  v -- the induced velocity */
void solver::lbd2vel(const MArray<double,1> &lbd, MArray<double,2> &v)
{
    int nvert_tot = stks->vertList[VESICLE].size();
    MArray<double,2> f(nvert_tot,3);
    f = 0.0;

    // Compute source density
    for (int ives= 0; ives < stks->numVesicles(); ives++) {
        Vesicle &vesicle = stks->vesicles[ives];
	if ( ! vesicle.isActive ) continue;

	int nvert = vesicle.numVerts();
	MArray<double,1> lbdtmp(nvert);
	MArray<double,2> ftmp(nvert,3);

	int p = vesicle.verts[0].Gindx;
	m_dcopy(nvert, &lbd(p), lbdtmp.data());
	vesicle.tensionForce(lbdtmp, ftmp);

	for (int ivert = 0; ivert < nvert; ivert++) {
	    Point &vert = vesicle.verts[ivert];

	    m_dcopy(3, &ftmp(ivert,0), vert.f);
	    m_dcopy(3, &ftmp(ivert,0), f.data()+3*vert.Gindx);
	    m_dclear(3, vert.g);
	}
    } 

    Vec rhs_vel = stks->getVesicleBieRhs();
    Vec sol_vel = stks->getVesicleBieSol();
    VecScatter scatter_vel = stks->getVesicleBieScatter();

    // Compute -1/(8*pi*mu) N(f)
    VecZeroEntries(rhs_vel);

    double cf = -1/(8*M_PI);

    if (ewald::phys::active) {
	NbrList &nlist = stks->nlist_phys[VESICLE][VESICLE];
	vector<Point*> &tlist = nlist.verts;
	int n = tlist.size();

	MArray<double,2> v(n,3);
	v = 0.0;
	myMatMultAdd(cf, stks->matSL[VESICLE][VESICLE], f.data(), v.data());

	MArray<int,2> idx(n,3);
	for (int i = 0; i < n; i++) {
	    int offset = 3*tlist[i]->Gindx;
	    FOR_J3 idx(i,j) = offset++;
	}
	VecSetValues(rhs_vel, 3*n, idx.data(), v.data(), ADD_VALUES);
    }

    if (ewald::four::active) {
        vector<Tri*> &slist = stks->slist_four[VESICLE];
	vector<Point*> &tlist = stks->tlist_four[VESICLE];

        ewald::four::clear_source();
	ewald::four::add_source(cf, 0.0, slist);
	ewald::four::transform();

	int n = tlist.size();
	MArray<double,2> v(n,3);
	v = 0.0;
	ewald::four::add_interp_vel(tlist, v.data());

	MArray<int,2> idx(n,3);
	for (int i = 0; i < n; i++) {
	    int offset = 3*tlist[i]->Gindx;
	    FOR_J3 idx(i,j) = offset++;
	}
	VecSetValues(rhs_vel, 3*n, idx.data(), v.data(), ADD_VALUES);
    }

    VecAssemblyBegin(rhs_vel);
    VecAssemblyEnd(rhs_vel);

    // Calculate the induced velocity by solving the BIE
    stks->solveVesicleBie(0.0);

    myVecScatter(scatter_vel, sol_vel, v.data(), INSERT_VALUES, SCATTER_FORWARD);
}
