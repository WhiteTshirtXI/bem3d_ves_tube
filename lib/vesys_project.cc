// Projection step of the flow solver
#include "cxxheaders.h"
#include "vesys.h"
#include "mympi.h"

/* Context for the matrix shell
 * Note: PETSc context only supports integer type, but we 
 *       need pointers, so use a local namespace  */
namespace { 
namespace ctx {
    VeSys *vesys;

    Mat lhs;
    KSP ksp;

    PetscErrorCode matmult(Mat, Vec, Vec);
    void lbd2vel(const MArray<double,1> &, MArray<double,2> &);
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
void VeSys::projectVel(const DoubArray2d &vold, const DoubArray1d &divTar,
		DoubArray2d &vnew, DoubArray1d &lbd,
		bool initial_guess)
{
    using namespace ctx;
    int nvert_tot = vertList.size();

    //======================================================================
    // rhs = divTar - div
    // Init
    VecZeroEntries(rhs_lbd); 

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

	VecSetValues(rhs_lbd, nvert, idx.data(), rhs_tmp.data(), ADD_VALUES);
    }


    VecAssemblyBegin(rhs_lbd);
    VecAssemblyEnd(rhs_lbd);

    //======================================================================
    // Solve the equation
    if (initial_guess) {
        myVecSetValues(sol_lbd, lbd.data(), INSERT_VALUES);
	VecAssemblyBegin(sol_lbd);
	VecAssemblyEnd(sol_lbd);

        KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    }
    else {
        KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
    }

    double rtol = 0.0;
    double atol = 0.0;

    for (int ives = 0; ives < numVesicles(); ives++) {
        Vesicle &vesicle = vesicles[ives];
	for (int ivert = 0; ivert < vesicle.numVerts(); ivert++) {
	    atol += square(1.E-3*vesicle.vertArea(ivert));
	}
    }
    atol = sqrt(atol);

    KSPSetTolerances(ksp, rtol, atol, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSolve(ksp, rhs_lbd, sol_lbd);

    // Convergence info
    int niter;
    KSPGetIterationNumber(ksp, &niter);

    double rhs_norm;
    VecNorm(rhs_lbd, NORM_2, &rhs_norm);
    rhs_norm /= sqrt(vertList.size());

    double res_norm;
    KSPGetResidualNorm(ksp, &res_norm);
    res_norm /= sqrt(vertList.size());

    if (mympi::comm_rank() == 0) {
	printf("    Project: niter = %d  rhs_norm = %9.2E, res_norm = %9.2E\n", 
			niter, rhs_norm, res_norm);
    }

    //======================================================================
    // Calculate the projected velocity
    myVecScatter(scatter_lbd, sol_lbd, lbd.data(), INSERT_VALUES, SCATTER_FORWARD);
    MArray<double,2> dv(nvert_tot,3);
    ctx::lbd2vel(lbd, dv);

    vnew = vold;
    vnew += dv;
}


void VeSys::initProjectSolver()
{
    using namespace ctx;
    ctx::vesys = this;

    int nrow, m;

    // vectors
    nrow = vertList.size();
    VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, nrow, &rhs_lbd);
    VecDuplicate(rhs_lbd, &sol_lbd);

    // scatter
    myVecScatterCreateToAll(sol_lbd, scatter_lbd);

    // lhs matrix 
    VecGetLocalSize(rhs_lbd, &m);
    MatCreateShell(MPI_COMM_WORLD, m, m, nrow, nrow, NULL, &lhs);
    MatShellSetOperation(lhs, MATOP_MULT, (void(*)(void))matmult);

    // ksp solver
    KSPCreate(MPI_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, lhs, lhs, DIFFERENT_NONZERO_PATTERN);
    KSPSetType(ksp, KSPGMRES);

    // preconditioner (none)
    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCNONE);
}


void VeSys::finalizeProjectSolver()
{
}

/* Matrix vector product for the solver
 * Arguments:
 *  div_vec = lhs * src_vec */
PetscErrorCode ctx::matmult(Mat lhs, Vec src_vec, Vec div_vec)
{
    // Init
    VecZeroEntries(div_vec);

    int nvert_tot = vesys->vertList.size();
    MArray<double,1> lbd(nvert_tot);
    myVecScatter(vesys->scatter_lbd, src_vec, lbd.data(), INSERT_VALUES, SCATTER_FORWARD);

    MArray<double,2> v(nvert_tot,3);
    lbd2vel(lbd, v);

    for (int ives = 0; ives < vesys->numVesicles(); ives++) {
        Vesicle &vesicle = vesys->vesicles[ives];
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
void ctx::lbd2vel(const MArray<double,1> &lbd, MArray<double,2> &v)
{
    // Compute source density
    for (int ives= 0; ives < vesys->numVesicles(); ives++) {
        Vesicle &vesicle = vesys->vesicles[ives];
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
	    m_dclear(3, vert.g);
	}
    } 


    // Compute -1/(8*pi*mu) N(f)
    VecZeroEntries(vesys->rhs_vel);
    vesys->addBoundaryIntegral(-1.0/(8*M_PI), 0.0, vesys->rhs_vel);
    VecAssemblyBegin(vesys->rhs_vel);
    VecAssemblyEnd(vesys->rhs_vel);

    // Calculate the induced velocity by solving the BIE
    vesys->solveBie(0.0);
    myVecScatter(vesys->scatter_vel, vesys->sol_vel, v.data(), 
    			INSERT_VALUES, SCATTER_FORWARD);
}
