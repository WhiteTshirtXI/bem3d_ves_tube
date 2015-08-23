// BIE solver for the vesicle system
#include "cxxheaders.h"
#include "veswall.h"

/* VesWall::solveVesicleBie
 * VesWall::initVesicleBieSolver
 * VesWall::getVesicleBieRhs
 * VesWall::getVesicleBieSol
 * VesWall::getVesicleBieScatter
 * solver::matmult */


/* Context for the matrix shell
 * Note: PETSc context only supports integer type, but we 
 *       need pointers, so use a local namespace  */
namespace { namespace solver {
    bool INITED = false;

    VesWall *stks;
    double dt;	// time step length
    
    Vec rhs_vec, sol_vec;
    Mat lhs;
    KSP ksp;
    PC pc;
    VecScatter scatter;

    PetscErrorCode matmult(Mat, Vec, Vec);
}};


/* Solve an implicit BIE
 * Arguments:
 *   myDt -- time step length (0 for explicit stepping) */
void VesWall::solveVesicleBie(double myDt)
{
    assert(solver::INITED);

    using namespace solver;

    // Implicit time stepping
    solver::dt = myDt;	

    double rtol = 1.E-5;
//    double atol = 1.E-3*sqrt(vertList.size());
    double atol = 0.0;
    double dtol = 100.0;
    KSPSetTolerances(ksp, rtol, atol, dtol, 100);
    KSPSolve(ksp, rhs_vec, sol_vec);
}


void VesWall::initVesicleBieSolver()
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
    int nrow = 3*vertList[VESICLE].size();
    VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, nrow, &rhs_vec);
    VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, nrow, &sol_vec);

    // vector scatter
    myVecScatterCreateToAll(rhs_vec, scatter);

    // lhs matrix 
    int m;
    VecGetLocalSize(rhs_vec, &m);
    MatCreateShell(MPI_COMM_WORLD, m, m, nrow, nrow, NULL, &lhs);
    MatShellSetOperation(lhs, MATOP_MULT, (void(*)(void))matmult);

    KSPCreate(MPI_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, lhs, lhs, DIFFERENT_NONZERO_PATTERN);
    KSPSetType(ksp, KSPGMRES);
    KSPSetFromOptions(ksp);

    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCNONE);
}


Vec VesWall::getVesicleBieRhs()
{
    assert(solver::INITED);
    return solver::rhs_vec;
}


Vec VesWall::getVesicleBieSol()
{
    assert(solver::INITED);
    return solver::sol_vec;
}

VecScatter VesWall::getVesicleBieScatter()
{
    assert(solver::INITED);
    return solver::scatter;
}


/* Matrix vector product for the solver
 * Arguments:
 *   v_vec = lhs * src_vec */
PetscErrorCode solver::matmult(Mat lhs, Vec src_vec, Vec v_vec)
{
    // Set source term
    int nvert_tot = stks->vertList[VESICLE].size();
    MArray<double,2> f(nvert_tot,3), g(nvert_tot,3);

    f = 0.0;
    g = 0.0;

    myVecScatter(scatter, src_vec, g.data(), INSERT_VALUES, SCATTER_FORWARD);

    for (int ives= 0; ives < stks->numVesicles(); ives++) {
        Vesicle &vesicle = stks->vesicles[ives];

	if (!vesicle.isActive) continue;

	int nvert = vesicle.numVerts();

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
	    m_dcopy(3, &df(ivert,0), &f(vert.Gindx,0));
	}
    }

    // Calculate boundary integrals
    // Init
    VecZeroEntries(v_vec);

    double cf = 1.0/(8*M_PI);
    double cg = -(1.0 - stks->viscRat)/(8*M_PI);

    bool fflag = fabs(cf) > 1.E-10;
    bool gflag = fabs(cg) > 1.E-10;

    // Physical sum
    if (ewald::phys::active) {
	NbrList &nlist = stks->nlist_phys[VESICLE][VESICLE];

	vector<Point*> &tlist = nlist.verts;
	int n = tlist.size();

	MArray<double,2> v(n,3);
	v = 0.0;
	if (fflag) myMatMultAdd(cf, stks->matSL[VESICLE][VESICLE], f.data(), v.data());
	if (gflag) myMatMultAdd(cg, stks->matDL[VESICLE][VESICLE], g.data(), v.data());

	MArray<int,2> idx(n,3);
	for (int i = 0; i < n; i++) {
	    int offset = 3*tlist[i]->Gindx;
	    FOR_J3 idx(i,j) = offset++;
	}
	VecSetValues(v_vec, 3*n, idx.data(), v.data(), ADD_VALUES);
    }


    // Ewald Fourier sum
    if (ewald::four::active) {
        ewald::four::clear_source();
	ewald::four::add_source(cf, cg, stks->slist_four[VESICLE]);
	ewald::four::transform();

	vector<Point*> &tlist = stks->tlist_four[VESICLE];
	int n = tlist.size();

	MArray<double,2> v(n,3);
	v = 0.0;
	ewald::four::add_interp_vel(tlist, v.data());

	MArray<int,2> idx(n,3);
	for (int i = 0; i < n; i++) {
	    int offset = 3*tlist[i]->Gindx;
	    FOR_J3 idx(i,j) = offset++;
	}
	VecSetValues(v_vec, 3*n, idx.data(), v.data(), ADD_VALUES);
    }

    // Linear term in the double layer potential
    if (gflag) {
	double vlin[3];
	double my_vlin[3] = { 0.0 };

	for (int ives = 0; ives < stks->numVesicles(); ives++) {
	    Vesicle &vesicle = stks->vesicles[ives];
	    if ( ! vesicle.isPrivate ) continue;

	    double tmp[3];
	    vesicle.calcDoubleLayerLinearTerm(tmp);
	    FOR_I3 my_vlin[i] += cg*tmp[i];
	}
	MPI_Allreduce(my_vlin, vlin, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	for (int ives = 0; ives < stks->numVesicles(); ives++) {
	    Vesicle &vesicle = stks->vesicles[ives];
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

	    VecSetValues(v_vec, 3*nvert, idx.data(), v.data(), ADD_VALUES);
	}
    }

    // Diagonal term
    {
        for (int ives = 0; ives < stks->numVesicles(); ives++) {
	    Vesicle &vesicle = stks->vesicles[ives];
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
