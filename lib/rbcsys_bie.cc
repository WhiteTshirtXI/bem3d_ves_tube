// BIE solver for the red cell system
#include "cxxheaders.h"
#include "rbcsys.h"
#include "mathfunc.h"

/* Context for the matrix shell
 * Note: PETSc context only supports integer type, but we 
 *       need pointers, so use a local namespace  */
namespace { namespace ctx {
    RbcSys *rbcsys;
    double dt;	// time step length, for implicit time marching
    
    Mat lhs;
    KSP ksp;

    PetscErrorCode matmult(Mat, Vec, Vec);
    void setRowIndices(vector<Point*> &, int *);
}};


/* RbcSys::solveBie
 * RbcSys::initBieSolver
 * RbcSys::finalizeBieSolver 
 * RbcSys::addBoundaryIntegral
 * RbcSys::calcResidual_Bie
 * ctx::matmult
 * ctx::setRowIndices */


/* Solve an implicit BIE
 * Arguments:
 *   myDt -- time step length (0 for explicit stepping) 
 * Note
 *   rhs_vel -- rhs
 *   sol_vel -- solution :*/
void RbcSys::solveBie(double myDt)
{
    using namespace ctx;

    // Set implicit time stepping
    ctx::dt = myDt;	

    double rtol = 1.E-3;
    double atol = 0.0;
    double dtol = 100.0;
    KSPSetTolerances(ksp, rtol, atol, dtol, 100);
    KSPSolve(ksp, rhs_vel, sol_vel);

    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank == 0) {
	int niter;
	double rnorm;
	KSPGetIterationNumber(ksp, &niter);
	KSPGetResidualNorm(ksp, &rnorm);
	rnorm /= sqrt(1.0*vertList.size());
        printf("    GMRES: niter = %d  rnorm = %.2E\n", niter, rnorm);
    }
}


void RbcSys::initBieSolver()
{
    using namespace ctx;
    ctx::rbcsys = this;

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


void RbcSys::finalizeBieSolver()
{
}


/* Add boundary integral
 *   v <- v + cf*Nf + cg*Kg */
void RbcSys::addBoundaryIntegral(double cf, MArray<double,2> &f,
		double cg, MArray<double,2> &g, Vec v)
{
    int nvert_tot = vertList.size();
    const bool fflag = fabs(cf) > 1.E-10;
    const bool gflag = fabs(cg) > 1.E-10;

    // Set source
    for (int i = 0; i < nvert_tot; i++) {
        Point *vert = vertList[i];
	m_dcopy(3, &f(i,0), vert->f);
	m_dcopy(3, &g(i,0), vert->g);
    }

    // Ewald physical sum
    if (ewald::phys::active) {
	NbrList &nlist = nlist_phys;

	vector<Point*> &tlist = nlist.verts;
	int nloc = tlist.size();
	MArray<double,2> vloc(nloc,3);
	MArray<int,2> idx(nloc,3);

	vloc = 0.0;
	if (fflag) myMatMultAdd(cf, matSL, f.data(), vloc.data());
	if (gflag) myMatMultAdd(cg, matDL, g.data(), vloc.data());

	ctx::setRowIndices(tlist, idx.data());
	VecSetValues(v, 3*nloc, idx.data(), vloc.data(), ADD_VALUES);
    }


    // Ewald Fourier sum
    if (ewald::four::active) {
        ewald::four::clear_source();
	ewald::four::add_source(cf, cg, slist_four);
	ewald::four::transform();

	vector<Point*> &tlist = tlist_four;
	int nloc = tlist.size();
	MArray<double,2> vloc(nloc,3);
	MArray<int,2> idx(nloc,3);

	vloc = 0.0;
	ewald::four::add_interp_vel(tlist, vloc.data());

	ctx::setRowIndices(tlist, idx.data());
	VecSetValues(v, 3*nloc, idx.data(), vloc.data(), ADD_VALUES);
    }

    // Linear term in the double layer potential
    if (gflag) {
	double vlin[3] = {0.0, 0.0, 0.0};
	for (int ives = 0; ives < numCells(); ives++) {
	    RedCell &cell = cells[ives];

	    double tmp[3];
	    cell.calcDoubleLayerLinearTerm(tmp);
	    m_daxpy(3, cg, tmp, vlin);
	}

	int mpi_size, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	int first, last;
	divideInterval(vertList.size(), mpi_size, mpi_rank, first, last);

	vector<Point*> tlist(vertList.begin() + first,
			    vertList.begin() + last);
	int nloc = last - first;
	MArray<double,2> vloc(nloc,3);
	MArray<int,2> idx(nloc,3);

	for (int i = 0; i < nloc; i++) {
	    m_dcopy(3, vlin, &vloc(i,0));
	}
	ctx::setRowIndices(tlist, idx.data());
	VecSetValues(v, 3*nloc, idx.data(), vloc.data(), ADD_VALUES);
    }
}


/* Calculate the residual in the BIE 
 * Arguments:
 *   vec -- rhs vector in the BIE
 * Note:
 *   rhs = b - A*x */
void RbcSys::calcResidual_Bie(Vec vec)
{
    int nvert_tot = vertList.size();
    MArray<double,2> f(nvert_tot,3), g(nvert_tot,3);

    f = 0.0;
    g = 0.0;

    for (int p = 0, icell = 0; icell < numCells(); icell++) {
	RedCell &cell = cells[icell];
	int nvert = cell.numVerts();

	// single-layer
	m_dadd(3*nvert, cell.f.data(), &f(p,0));

	p += nvert;
    }

    MArray<double,2> diag(nvert_tot,3);
    for (int i = 0; i < nvert_tot; i++) {
	Point *vert = vertList[i];

	double vbkg[3];
	calcBkgVel(vert->x, vbkg);

	FOR_J3 diag(i,j) = vbkg[j];
    }

    VecZeroEntries(vec);
    myVecSetValues(vec, diag.data(), ADD_VALUES);
    addBoundaryIntegral(-1.0/(8*M_PI), f, 0, g, vec);
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);
}


/* Matrix vector product for the solver
 * Arguments:
 *   v_vec = lhs * src_vec */
PetscErrorCode ctx::matmult(Mat lhs, Vec src_vec, Vec v_vec)
{
    //======================================================================
    // Set source term
    int nvert_tot = rbcsys->vertList.size();
    MArray<double,2> f(nvert_tot,3), g(nvert_tot,3);
    myVecScatter(rbcsys->scatter_vel, src_vec, g.data(), INSERT_VALUES, SCATTER_FORWARD);

    for (int p = 0, icell= 0; icell < rbcsys->numCells(); icell++) {
        RedCell &cell = rbcsys->cells[icell];
//	if (!cell.active) continue;

	int nvert = cell.numVerts();
	MArray<double,2> dx(nvert,3), ftmp(nvert,3);
	
	m_dcopy(3*nvert, &g(p,0), dx.data());
	dx *= dt;

	if (fabs(dt) < 1.E-10) 
	    ftmp = 0.0;
	else 
	    cell.diffElasBendForce(dx, ftmp);

        m_dcopy(3*nvert, ftmp.data(), &f(p,0));

	p += nvert;
    }

    //======================================================================
    // Init
    VecZeroEntries(v_vec);

    double cf = 1.0/(8*M_PI);
    double cg = -(1.0 - rbcsys->viscRat)/(8*M_PI);

    rbcsys->addBoundaryIntegral(cf, f, cg, g, v_vec);

    // Diagonal term
    {
	int mpi_size, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	vector<Point*> &vertList = rbcsys->vertList;

	int first, last;
	divideInterval(vertList.size(), mpi_size, mpi_rank, first, last);

	vector<Point*> tlist(vertList.begin() + first,
			    vertList.begin() + last);
	int nloc = last - first;
	MArray<double,2> vloc(nloc,3);
	MArray<int,2> idx(nloc,3);

        const double cdiag = 4*M_PI*cg;
	vloc = 0.0;
	for (int i = 0; i < nloc; i++) {
	    Point *vert = tlist[i];
	    FOR_D3 {
	        vloc(i,d) = vert->g[d];
	        double tmp = 0.0;
	        FOR_K3 tmp += vert->C[d][k]*vert->g[k];
	        vloc(i,d) += cdiag*tmp;
	    }
	}

	// Add to global vector
	setRowIndices(tlist, idx.data());
	VecSetValues(v_vec, 3*nloc, idx.data(), vloc.data(), ADD_VALUES);
    }

    VecAssemblyBegin(v_vec); 
    VecAssemblyEnd(v_vec); 

    return 0;
}


/* Create row index for a target list
 * Arguments:
 *   tlist -- list of target point
 *   idx -- the row indices */
void ctx::setRowIndices(vector<Point*> &tlist, int *idx)
{
    for (int i = 0; i < tlist.size(); i++) {
	int p = 3*i;
	int irow = 3*tlist[i]->Gindx;

	idx[p++] = irow++;
	idx[p++] = irow++;
	idx[p++] = irow++;
    }
}
