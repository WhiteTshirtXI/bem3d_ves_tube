// BIE solver for the vesicle system
#include "cxxheaders.h"
#include "vesys.h"
#include "mympi.h"

/* Context for the matrix shell
 * Note: PETSc context only supports integer type, but we 
 *       need pointers, so use a local namespace  */
namespace { namespace ctx {
    VeSys *vesys;
    double dt;	// time step length, for implicit time marching
    
    Mat lhs;
    KSP ksp;

    PetscErrorCode matmult(Mat, Vec, Vec);
}};


/* VeSys::solveBie
 * VeSys::initBieSolver
 * VeSys::finalizeBieSolver 
 * VeSys::addBoundaryIntegral
 * ctx::matmult */


/* Solve an implicit BIE
 * Arguments:
 *   myDt -- time step length (0 for explicit stepping) 
 * Note
 *   rhs_vel -- rhs
 *   sol_vel -- solution :*/
void VeSys::solveBie(double myDt)
{
    using namespace ctx;

    // Set implicit time stepping
    ctx::dt = myDt;	

    double rtol = 0.0;
    double atol = 0.0;
    double dtol = 100.0;
    int maxits = 100;

    if (fabs(myDt) < 1.E-10) {
        rtol = 1.E-2;
    } else {
	atol = myDt*sqrt(vertList.size());
    }
    KSPSetTolerances(ksp, rtol, atol, dtol, maxits);
    KSPSolve(ksp, rhs_vel, sol_vel);

    // Convergence info
    int niter;
    KSPGetIterationNumber(ksp, &niter);

    double rhs_norm;
    VecNorm(rhs_vel, NORM_2, &rhs_norm);
    rhs_norm /= sqrt(vertList.size());

    double res_norm;
    KSPGetResidualNorm(ksp, &res_norm);
    res_norm /= sqrt(vertList.size());

    if (mympi::comm_rank() == 0 && myDt > 1.E-10) {
        printf("    BIE: niter = %d  rhs_norm = %9.2E res_norm = %9.2E\n", 
			niter, rhs_norm, res_norm);
    }
}


void VeSys::initBieSolver()
{
    using namespace ctx;
    ctx::vesys = this;

    int nrow, m;

    // vector
    nrow = 3*vertList.size();
    VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, nrow, &rhs_vel);
    VecDuplicate(rhs_vel, &sol_vel);

    // scatter
    myVecScatterCreateToAll(sol_vel, scatter_vel);

    // lhs matrix 
    VecGetLocalSize(rhs_vel, &m);
    MatCreateShell(MPI_COMM_WORLD, m, m, nrow, nrow, NULL, &lhs);
    MatShellSetOperation(lhs, MATOP_MULT, (void(*)(void))matmult);

    // ksp solver
    KSPCreate(MPI_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, lhs, lhs, DIFFERENT_NONZERO_PATTERN);
    KSPSetType(ksp, KSPGMRES);
    KSPSetFromOptions(ksp);

    // preconditioner (none)
    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCNONE);
}


void VeSys::finalizeBieSolver()
{
}


/* Add boundary integral
 * Argumets:
 *   cf, cg -- single- and double-layer integral coefficients
 *   v -- the induced velocity 
 * Note: 
 *   -- v = v + cf*Nf + cg*Kg 
 *   -- The densities are stored in vert.f and vert.g arrays */
void VeSys::addBoundaryIntegral(double cf, double cg, Vec vec)
{
    int nvert_tot = vertList.size();
    const bool fflag = fabs(cf) > 1.E-10;
    const bool gflag = fabs(cg) > 1.E-10;

    // Set source
    MArray<double,2> f(nvert_tot,3), g(nvert_tot,3);
    f = 0.0;
    g = 0.0;

    for (int ives = 0; ives < numVesicles(); ives++) {
        Vesicle &vesicle = vesicles[ives];
	if ( ! vesicle.isActive ) continue;

	for (int ivert = 0; ivert < vesicle.numVerts(); ivert++) {
	    Point &vert = vesicle.verts[ivert];

	    int p = vert.Gindx;
	    m_dcopy(3, vert.f, &f(p,0));
	    m_dcopy(3, vert.g, &g(p,0));
	}
    }

    // Physical sum
    if (ewald::phys::active) {
	NbrList &nlist = nlist_phys;

	vector<Point*> &tlist = nlist.verts;
	int n = tlist.size();

	MArray<double,2> v(n,3);
	v = 0.0;
	if (fflag) myMatMultAdd(cf, matSL, f.data(), v.data());
	if (gflag) myMatMultAdd(cg, matDL, g.data(), v.data());

	MArray<int,2> idx(n,3);
	for (int i = 0; i < n; i++) {
	    int offset = 3*tlist[i]->Gindx;
	    FOR_J3 idx(i,j) = offset++;
	}
	VecSetValues(vec, 3*n, idx.data(), v.data(), ADD_VALUES);
    }


    // Ewald Fourier sum
    if (ewald::four::active) {
        ewald::four::clear_source();
	ewald::four::add_source(cf, cg, slist_four);
	ewald::four::transform();

	vector<Point*> &tlist = tlist_four;
	int n = tlist.size();

	MArray<double,2> v(n,3);
	v = 0.0;
	ewald::four::add_interp_vel(tlist, v.data());

	MArray<int,2> idx(n,3);
	for (int i = 0; i < n; i++) {
	    int offset = 3*tlist[i]->Gindx;
	    FOR_J3 idx(i,j) = offset++;
	}
	VecSetValues(vec, 3*n, idx.data(), v.data(), ADD_VALUES);
    }

    // Linear term in the double layer potential
    if (gflag) {
	double vlin[3];
	double my_vlin[3] = { 0.0 };

	for (int ives = 0; ives < numVesicles(); ives++) {
	    Vesicle &vesicle = vesicles[ives];
	    if ( ! vesicle.isPrivate ) continue;

	    double tmp[3];
	    vesicle.calcDoubleLayerLinearTerm(tmp);
	    FOR_I3 my_vlin[i] += cg*tmp[i];
	}
	MPI_Allreduce(my_vlin, vlin, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	for (int ives = 0; ives < numVesicles(); ives++) {
	    Vesicle &vesicle = vesicles[ives];
	    if ( ! vesicle.isPrivate ) continue;

	    int nvert = vesicle.numVerts();
	    MArray<double,2> v(nvert,3);
	    for (int ivert = 0; ivert < nvert; ivert++) {
	        FOR_J3 v(ivert,j) = vlin[j];
	    }

	    MArray<int,2> idx(nvert,3);
	    for (int ivert = 0; ivert < nvert; ivert++) {
	        int offset = 3*vesicle.verts[ivert].Gindx;
		FOR_J3 idx(ivert,j) = offset++;
	    }

	    VecSetValues(vec, 3*nvert, idx.data(), v.data(), ADD_VALUES);
	}
    }
}


/* Matrix vector product for the solver
 * Arguments:
 *   v_vec = lhs * src_vec */
PetscErrorCode ctx::matmult(Mat lhs, Vec src_vec, Vec v_vec)
{
    //======================================================================
    // Set source term
    int nvert_tot = vesys->vertList.size();
    MArray<double,2> g(nvert_tot,3);
    myVecScatter(vesys->scatter_vel, src_vec, g.data(), INSERT_VALUES, SCATTER_FORWARD);

    for (int ives= 0; ives < vesys->numVesicles(); ives++) {
        Vesicle &vesicle = vesys->vesicles[ives];
	int nvert = vesicle.numVerts();

	if (vesicle.isActive) {
	    MArray<double,2> dx(nvert,3);
	    for (int ivert = 0; ivert < nvert; ivert++) {
		Point &vert = vesicle.verts[ivert];

		m_dcopy(3, &g(vert.Gindx,0), vert.g);
		m_dcopy(3, &g(vert.Gindx,0), &dx(ivert,0));
	    }
	    dx *= dt;
	    
	    MArray<double,2> df(nvert,3);
	    if (fabs(dt) < 1.E-10) 
		df = 0.0;
	    else 
		vesicle.diffBendForce(dx, df);

	    for (int ivert = 0; ivert < nvert; ivert++) {
		Point &vert = vesicle.verts[ivert];
		m_dcopy(3, &df(ivert,0), vert.f);
	    }
	}
    }

    //======================================================================
    // Init
    VecZeroEntries(v_vec);

    double cf = 1.0/(8*M_PI);
    double cg = -(1.0 - vesys->viscRat)/(8*M_PI);
    vesys->addBoundaryIntegral(cf, cg, v_vec);

    // Diagonal term
    {
        for (int ives = 0; ives < vesys->numVesicles(); ives++) {
	    Vesicle &vesicle = vesys->vesicles[ives];
	    if ( ! vesicle.isPrivate ) continue;

	    int nvert = vesicle.numVerts();

	    MArray<double,2> v(nvert,3);
	    v = 0.0;
	    for (int ivert = 0; ivert < nvert; ivert++) {
	        Point &vert = vesicle.verts[ivert];

		double Cg[3] = {0.0};
		FOR_J3 FOR_K3 Cg[j] += vert.C[j][k]*vert.g[k];

		double s = 4*M_PI*cg;
	        FOR_J3 v(ivert,j) += vert.g[j] + s*Cg[j];
	    }

	    MArray<int,2> idx(nvert,3);
	    for (int ivert = 0; ivert < nvert; ivert++) {
	        int offset = 3*vesicle.verts[ivert].Gindx;
	        FOR_J3 idx(ivert,j) = offset++;
	    }

	    VecSetValues(v_vec, 3*nvert, idx.data(), v.data(), ADD_VALUES);
	}
    }

    VecAssemblyBegin(v_vec); 
    VecAssemblyEnd(v_vec); 

    return 0;
}
