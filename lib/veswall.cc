#include "cxxheaders.h"
#include "veswall.h"
#include "ewald.h"
#include "mathfunc.h"
#include "quadrature.h"
#include "param.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "debugfunc.h"
#include "collision.h"

/* VesWall::VesWall 
 * VesWall::~VesWall 
 * VesWall::initialHook 
 * VesWall::domainDecomp
 * VesWall::updateNeighborList
 * VesWall::calcBieMat
 * VesWall::addBoundaryIntegrals
 * VesWall::calcRhs_vesicleBie
 * VesWall::calcRhs_WallBie
 * VesWall::timeInt
 * VesWall::rebox
 * VesWall::syncCoord
 * VesWall::calcStress
 * VesWall::calcProbeVel 
 * VesWall::writeAll 
 * VesWall::writeVesicles 
 * VesWall::writeWalls 
 * VesWall::writeRestart 
 * VesWall::readRestart 
 * VesWall::writeProbe
 * VesWall::writeVesicleCenter 
 * VesWall::volFraction */

// Subsystem identifies
namespace {
    const int NBLK = 2;
};


VesWall::VesWall()
{
    for (int iblk = 0; iblk < NBLK; iblk++)
    for (int jblk = 0; jblk < NBLK; jblk++) {
        matSL[iblk][jblk] = NULL;
        matDL[iblk][jblk] = NULL;
    }
}


VesWall::~VesWall()
{ 
}

/* -- Index meshes, vertices, and faces
 * -- Build vertex lists and face lists 
 * -- Allocate additional arrays */
void VesWall::initialHook()
{
    // Build mesh list and index meshes
    for (int i = 0; i < numVesicles(); i++) {
        Vesicle &vesicle = vesicles[i];

        meshList[VESICLE].push_back(&vesicle);

	vesicle.Gindx = i;
	vesicle.setInternalPointers();
	vesicle.connectPeriodicBoundaries();
	vesicle.buildControlPoints();

	int nvert = vesicle.numVerts();

	vesicle.vertArea.resize(nvert);
	vesicle.vertNrml.resize(nvert,3);
	vesicle.vertH.resize(nvert);

	vesicle.fbend.resize(nvert,3);
	vesicle.ftens.resize(nvert,3);
	vesicle.v.resize(nvert,3);
	vesicle.sigma.resize(nvert);
    }


    for (int i = 0; i < numWalls(); i++) {
	Wall &wall = walls[i];

	meshList[WALL].push_back(&wall);

	wall.Gindx = i;
	wall.setInternalPointers();
	wall.connectPeriodicBoundaries(ewald::L);

	int nvert = wall.numVerts();
	wall.f.resize(nvert,3);
    }

    // Build vertex list, vertex indices, and face list
    for (int blk = 0; blk < NBLK; blk++) {
        vector<Mesh*> &meshes = meshList[blk];
	vector<Point*> &vlist = vertList[blk];
	vector<Tri*> &flist = faceList[blk];
        int vert_cnt = 0;
	int face_cnt = 0;

	for (int imesh = 0; imesh < meshes.size(); imesh++) {
	    Mesh *mesh = meshes[imesh];

	    // Vertices
	    for (int ivert = 0; ivert < mesh->numVerts(); ivert++) {
		Point &vert = mesh->verts[ivert];
		if (vert.pp == NULL) {
		    vlist.push_back(&vert);
		    vert.Gindx = vert_cnt++;
		}
	    }

	    // Periodic boundaries
	    for (int ivert = 0; ivert < mesh->numVerts(); ivert++) {
		Point &vert = mesh->verts[ivert];
		if (vert.pp != NULL) vert.Gindx = vert.pp->Gindx;
	    }

	    // Faces
	    for (int iface = 0; iface < mesh->numFaces(); iface++) {
		Tri &face = mesh->faces[iface];
		flist.push_back(&face);
		face.Gindx = face_cnt++;
	    }
        } // imesh
    } // blk
}


/* Build source and target lists
 * Note:
 *   - The fixed lists are taken care by calcFixedSTNList */
void VesWall::domainDecomp()
{
    static bool INITED = false;

    if (ewald::phys::active) {
        for (int blk = 0; blk < NBLK; blk++) {
	    if (blk == WALL && INITED) continue;
	    slist_phys[blk].build(faceList[blk]);
	}
    }

    if (ewald::four::active) {
        for (int blk = 0; blk < NBLK; blk++) {
	    if (blk == WALL && INITED) continue;
            ewald::four::setSources(faceList[blk], slist_four[blk]);
            ewald::four::setTargets(vertList[blk], tlist_four[blk]);
        }
    }

    // Set mesh ownership
    for (int blk = 0; blk < NBLK; blk++) {
        if (blk == WALL && INITED) continue;

	int nmesh = meshList[blk].size();
	if (nmesh == 0) continue;

	// Init
        for (int imesh = 0; imesh < meshList[blk].size(); imesh++) {
	    Mesh *mesh = meshList[blk][imesh];
	    mesh->isActive = 0;
	    mesh->isPrivate = 0;
	}

	// Active flag
	for (int i = 0; i < slist_phys[blk].numFaces(); i++) {
	    Mesh *mesh = slist_phys[blk].faces[i]->mesh;
	    mesh->isActive = 1;
	}

	for (int i = 0; i < slist_four[blk].size(); i++) {
	    Mesh *mesh = slist_four[blk][i]->mesh;
	    mesh->isActive = 1;
	}

	// Private flag
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	MArray<int,1> ownerRank(nmesh), ownerRankMax(nmesh);
	for (int imesh = 0; imesh < nmesh; imesh++) {
	    Mesh *mesh = meshList[blk][imesh];
            ownerRank(imesh) = mesh->isActive ? mpi_rank : -1;
	}

        MPI_Allreduce(ownerRank.data(), ownerRankMax.data(), 
			nmesh, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        for (int imesh = 0; imesh < nmesh; imesh++) {
	    Mesh *mesh = meshList[blk][imesh];
	    mesh->isPrivate = (ownerRankMax(imesh) == mpi_rank);
	}
    } // blk

    INITED = true;
}


/*   Update the geometry of meshes that are active */
void VesWall::updateSourceGeometry()
{
    static bool INITED = false;

    for (int blk = 0; blk < NBLK; blk++) {
        if (blk == WALL && INITED) continue;

        for (int imesh = 0; imesh < meshList[blk].size(); imesh++) {
	    Mesh *mesh = meshList[blk][imesh];
	    mesh->updateGeometry();
        }
    }

    INITED = true;
}


// Update neighbor list
void VesWall::updateNeighborList()
{
    static bool INITED = false;

    if (ewald::phys::active) {
        for (int iblk = 0; iblk < NBLK; iblk++) 
	for (int jblk = 0; jblk < NBLK; jblk++) {
	    if (iblk == WALL && jblk == WALL && INITED) continue;

	    NbrList &nlist = nlist_phys[iblk][jblk];
	    nlist.build(vertList[iblk], slist_phys[jblk]);
        }
    }

    INITED = true;
}


// Calculate BIE matrices
void VesWall::calcBieMat()
{
    static bool INITED = false;

    for (int ib = 0; ib < NBLK; ib++)
    for (int jb = 0; jb < NBLK; jb++) {
	if (ib == WALL && jb == WALL && INITED) continue;

	SrcList &slist = slist_phys[ib];
	NbrList &nlist = nlist_phys[ib][jb];

	int nrow = 3*nlist.verts.size();
	int ncol = 3*vertList[jb].size();

	Mat &SL = matSL[ib][jb];
	Mat &DL = matDL[ib][jb];

	if (SL) MatDestroy(&SL);
	if (DL) MatDestroy(&DL);

	MatCreate(PETSC_COMM_SELF, &SL);
	MatSetType(SL, MATSEQAIJ);
	MatSetSizes(SL, PETSC_DECIDE, PETSC_DECIDE, nrow, ncol);

	MatCreate(PETSC_COMM_SELF, &DL);
	MatSetType(DL, MATSEQAIJ);
	MatSetSizes(DL, PETSC_DECIDE, PETSC_DECIDE, nrow, ncol);

	if (ewald::phys::active) {
	    if (jb == WALL) 
		ewald::phys::calcSurfIntMat(nlist, SL, NULL);
	    else 
		ewald::phys::calcSurfIntMat(nlist, SL, DL);
        }
    }

    INITED = true;
}


/* Add boundary integrals
 * Arguments:
 *   iblk -- the array block where the integrals are added to
 *   vec -- the vector
 *   cf, cg -- single and double layer coefficients 
 * Note:
 *   -- The source terms are stored in vert.f[] and vert.g[] arrays of every mesh, so
 *      they must be set up 
 *   -- For double-layer integral, we use the limiting value from the plus side of the
 *      surface
 *   */
void VesWall::addBoundaryIntegrals(int iblk, Vec vec, const double *cf, const double *cg)
{
    assert(iblk == VESICLE || iblk == WALL);

    const double eps = 1.E-10;

    // Has double layer 
    bool gflag = false;
    for (int jblk = 0; jblk < NBLK; jblk++) {
        if ( fabs(cg[jblk]) > eps ) gflag = true;
    }

    // Add jump terms in the double layer, only for nonzero cg[iblk] 
    bool jump_flag = false;
    if (fabs(cg[iblk]) > eps) jump_flag = true;

    // Physical sum
    if (ewald::phys::active) {
	for (int jblk = 0; jblk < NBLK; jblk++) {
	    if ( fabs(cf[jblk]) < eps && fabs(cg[jblk]) < eps ) continue;

	    int nvert_jblk = vertList[jblk].size();
	    MArray<double,2> f(nvert_jblk,3), g(nvert_jblk,3);
	    f = 0.0;
	    g = 0.0;

	    for (int imesh = 0; imesh < meshList[jblk].size(); imesh++) {
		Mesh &mesh = *meshList[jblk][imesh];
		if (!mesh.isActive ) continue;

		for (int ivert = 0; ivert < mesh.numVerts(); ivert++) {
		    Point &vert = mesh.verts[ivert];

		    if (fabs(cf[jblk]) > eps) m_dcopy(3, vert.f, &f(vert.Gindx,0));
		    if (fabs(cg[jblk]) > eps) m_dcopy(3, vert.g, &g(vert.Gindx,0));
		}
	    }

	    NbrList &nlist = nlist_phys[iblk][jblk];
	    vector<Point*> &tlist = nlist.verts;
	    int n = tlist.size();
	    MArray<double,2> v(n,3);
	    v = 0.0;

	    if (fabs(cf[jblk]) > eps) myMatMultAdd(cf[jblk], matSL[iblk][jblk], f.data(), v.data());
	    if (fabs(cg[jblk]) > eps) myMatMultAdd(cg[jblk], matDL[iblk][jblk], g.data(), v.data());

	    MArray<int,2> idx(n,3);
	    for (int i = 0; i < n; i++) {
		int offset = 3*tlist[i]->Gindx;
		FOR_J3 idx(i,j) = offset++;
	    }
	    VecSetValues(vec, 3*n, idx.data(), v.data(), ADD_VALUES);
	}
    }

    // Fourier sum
    if (ewald::four::active) {
        ewald::four::clear_source();
	for (int jblk = 0; jblk < NBLK; jblk++) {
	    ewald::four::add_source(cf[jblk], cg[jblk], slist_four[jblk]);
        }
        ewald::four::transform();

	vector<Point*> &tlist = tlist_four[iblk];
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

    if (gflag) {
	// Linear term in the double layer potential
	double vlin[3] = {0.0};
	for (int jblk = 0; jblk < NBLK; jblk++) {
	    if ( fabs(cg[jblk]) < eps ) continue;

	    double my_vlin[3] = { 0.0 };
	    double my_vlin_sum[3] = { 0.0 };

	    for (int imesh = 0; imesh < meshList[jblk].size(); imesh++) {
		Mesh &mesh = *meshList[jblk][imesh];
		if (!mesh.isPrivate) continue;

		double tmp[3];
		mesh.calcDoubleLayerLinearTerm(tmp);
		FOR_I3 my_vlin[i] += cg[jblk]*tmp[i];
	    }

	    MPI_Allreduce(my_vlin, my_vlin_sum, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    FOR_I3 vlin[i] += my_vlin_sum[i];
	}

	// Add vlin[] as well as double-layer jump term if it exists
	double beta = 4.0*M_PI*cg[iblk];

	for (int imesh = 0; imesh < meshList[iblk].size(); imesh++) {
	    Mesh &mesh = *meshList[iblk][imesh];
	    if (!mesh.isPrivate) continue;

	    int nvert = mesh.numVerts();
	    MArray<double,2> v(nvert,3);
	    v = 0.0;

	    for (int ivert = 0; ivert < nvert; ivert++) {
		Point &vert = mesh.verts[ivert];
	        if (vert.pp) continue;

	        FOR_I3 v(ivert,i) = vlin[i];

		if (jump_flag) {
		    FOR_I3 FOR_J3 v(ivert,i) += beta*vert.C[i][j]*vert.g[j];
		}
	    }

	    MArray<int,2> idx(nvert,3);
	    for (int ivert = 0; ivert < nvert; ivert++) {
	        int offset = 3*mesh.verts[ivert].Gindx;
		FOR_J3 idx(ivert,j) = offset++;
	    }

	    VecSetValues(vec, 3*nvert, idx.data(), v.data(), ADD_VALUES);
	}
    }

    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);
}



// Calculate the boundary integrals at a set of random points
void VesWall::addBoundaryIntegrals(const MArray<double,2> &x, Vec vec, const double *cf, const double *cg)
{
    const double eps = 1.E-10;

    bool gflag = false;
    for (int jblk = 0; jblk < NBLK; jblk++) {
        if ( fabs(cg[jblk]) > eps ) gflag = true;
    }

    // Set probe points
    int nglb = x.size(0);
    vector<Point> probes(nglb);
    vector<Point*> pprobes(nglb);
    for (int i = 0; i < nglb; i++) {
        FOR_D3 probes[i].x[d] = x(i,d);
	probes[i].Gindx = i;
	pprobes[i] = &probes[i];
    }

    if (ewald::phys::active) {
        for (int blk = 0; blk < NBLK; blk++) {

            NbrList nlist;
            nlist.build(pprobes, slist_phys[blk]);

	    int n = nlist.verts.size();
	    MArray<double,2> v(n,3);
	    v = 0.0;
            ewald::phys::addSurfInt(nlist, cf[blk], cg[blk], v.data());

	    MArray<int,2> idx(n,3);
            for (int i = 0; i < n; i++) {
                int offset = 3*nlist.verts[i]->Gindx;
                FOR_J3 idx(i,j) = offset++;
            }
            VecSetValues(vec, 3*n, idx.data(), v.data(), ADD_VALUES);
        }
    }

    // Fourier sum
    if (ewald::four::active) {
        ewald::four::clear_source();
	for (int jblk = 0; jblk < NBLK; jblk++) {
	    ewald::four::add_source(cf[jblk], cg[jblk], slist_four[jblk]);
        }
        ewald::four::transform();

	vector<Point*> tlist;
	ewald::four::setTargets(pprobes, tlist);
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

    if (gflag) {
	// Linear term in the double layer potential
	double vlin[3] = {0.0};

	for (int jblk = 0; jblk < NBLK; jblk++) {
	    if ( fabs(cg[jblk]) < eps ) continue;

	    double my_vlin[3] = { 0.0 };
	    double my_vlin_sum[3] = { 0.0 };

	    for (int imesh = 0; imesh < meshList[jblk].size(); imesh++) {
		Mesh &mesh = *meshList[jblk][imesh];
		if (!mesh.isPrivate) continue;

		double tmp[3];
		mesh.calcDoubleLayerLinearTerm(tmp);
		FOR_I3 my_vlin[i] += cg[jblk]*tmp[i];
	    }

	    MPI_Allreduce(my_vlin, my_vlin_sum, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    FOR_I3 vlin[i] += my_vlin_sum[i];
	}

	// Add vlin[]
	int idx_low, idx_high;
	VecGetOwnershipRange(vec, &idx_low, &idx_high);
	int n = idx_high - idx_low;
	MArray<double,1> v(n);
	MArray<int,1> idx(n);
	for (int i = 0; i < n; i++) {
	    int d = (idx_low + i)%3;
	    v(i) = vlin[d];
	    idx(i) = idx_low + i;
	}
	VecSetValues(vec, n, idx.data(), v.data(), ADD_VALUES);
    }

    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);
}


/* Calculate the residual in the vesicle BIE
 * Arguments:
 *   rhs -- rhs vector of the vesicle BIE */
void VesWall::calcRhs_vesicleBie(Vec rhs)
{
    // Set up source term
    for (int ives = 0; ives < numVesicles(); ives++) {
	Vesicle &vesicle = vesicles[ives];
	if ( ! vesicle.isActive ) continue;

	int nvert = vesicle.numVerts();
	for (int ivert = 0; ivert < nvert; ivert++) {
	    Point &vert = vesicle.verts[ivert];

	    double ftmp[3];
	    FOR_I3 ftmp[i] = vesicle.fbend(ivert,i) + vesicle.ftens(ivert,i);

//	    // debug
//	    ftmp[0] += sin(2*M_PI*ewald::iL[1]*vesicle.center[1]);
//	    ftmp[1] += sin(2*M_PI*ewald::iL[2]*vesicle.center[2]);
//	    ftmp[2] += sin(2*M_PI*ewald::iL[0]*vesicle.center[0]);
//	    // end debug

	    m_dcopy(3, ftmp, vert.f);
	    m_dclear(3, vert.g);
	}
    }

    for (int iwall = 0; iwall < numWalls(); iwall++) {
        Wall &wall = walls[iwall];
	if (! wall.isActive ) continue;

	int nvert = wall.numVerts();
	for (int ivert = 0; ivert < nvert; ivert++) {
	    Point &vert = wall.verts[ivert];

	    FOR_I3 vert.f[i] = wall.f(ivert,i);
	}
    }

    // Add boundary integrals
    double cf[NBLK] = { -1.0/(8*M_PI), -1.0/(8*M_PI) };
    double cg[NBLK] = { 0.0, 0.0 };
    VecZeroEntries(rhs);
    addBoundaryIntegrals(VESICLE, rhs, cf, cg);

    // Add background velocity
    for (int ives = 0; ives < numVesicles(); ives++) {
        Vesicle &vesicle = vesicles[ives];
	if (!vesicle.isPrivate) continue;

	int nvert = vesicle.numVerts();
	MArray<double,2> v(nvert,3);
	v = 0.0;

	for (int ivert = 0; ivert < nvert; ivert++) {
	    Point &vert = vesicle.verts[ivert];
	    if (vert.pp) continue;

	    FOR_D3 v(ivert,d) = vbkg[d];
	    v(ivert,0) += shRate*(vert.x[2] - 0.5*ewald::L[2]);
	}

	MArray<int,2> idx(nvert,3);
	for (int ivert = 0; ivert < nvert; ivert++) {
	    int offset = 3*vesicle.verts[ivert].Gindx;
	    FOR_J3 idx(ivert,j) = offset++;
	}

	VecSetValues(rhs, 3*nvert, idx.data(), v.data(), ADD_VALUES);
    }
    VecAssemblyBegin(rhs);
    VecAssemblyEnd(rhs);
}


/* Calculate the residual in the wall BIE
 * Argument:
 *   rhs -- */
void VesWall::calcRhs_wallBie(Vec rhs)
{
    // Set up source term
    for (int ives = 0; ives < numVesicles(); ives++) {
	Vesicle &vesicle = vesicles[ives];
	if ( ! vesicle.isActive ) continue;

	int nvert = vesicle.numVerts();
	for (int ivert = 0; ivert < nvert; ivert++) {
	    Point &vert = vesicle.verts[ivert];

	    double ftmp[3];
	    FOR_D3 ftmp[d] = vesicle.fbend(ivert,d) + vesicle.ftens(ivert,d);

//	    // debug
//	    ftmp[0] += sin(2*M_PI*ewald::iL[1]*vesicle.center[1]);
//	    ftmp[1] += sin(2*M_PI*ewald::iL[2]*vesicle.center[2]);
//	    ftmp[2] += sin(2*M_PI*ewald::iL[0]*vesicle.center[0]);
//	    // end debug

	    m_dcopy(3, ftmp, vert.f);
	    m_dcopy(3, &vesicle.v(ivert,0), vert.g);
	}
    }

    // Add boundary integrals
    VecZeroEntries(rhs);
    double cf[NBLK] = { -1.0/(8*M_PI), 0.0 };
    double cg[NBLK] = { (1.0 - viscRat)/(8*M_PI), 0.0 };
    addBoundaryIntegrals(WALL, rhs, cf, cg);

    // Add background velocity
    {
	int idx_low, idx_high;
	VecGetOwnershipRange(rhs, &idx_low, &idx_high);
	int n = idx_high - idx_low;
	MArray<double,1> v(n);
	MArray<int,1> idx(n);
	for (int i = 0; i < n; i++) {
	    int d = (idx_low + i)%3;
	    v(i) = vbkg[d];
	    idx(i) = idx_low + i;
	}
	VecSetValues(rhs, n, idx.data(), v.data(), ADD_VALUES);
    }
    VecAssemblyBegin(rhs);
    VecAssemblyEnd(rhs);
}

// Time integration
void VesWall::timeInt()
{
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    char token[256], tensionfn[256];
    int nout;
    double facelocalpressure;
	double stresstau[3][3];


    //if(mpi_rank == 0) {cout << vesicles[0].verts[0].x[0] << "a\n";}


    // Zero out vesicle.sigma initially
    for (int ives = 0; ives < numVesicles(); ives++) {
	Vesicle &vesicle = vesicles[ives];
	vesicle.sigma = 0.0;
    }

    //if(mpi_rank == 0) {cout << vesicles[0].verts[0].x[0] << "b\n";}


	if(mpi_rank == 0) //Log simulation parameters as header to output summary file
	{
		printf("Output directory = %s\n", outfiledirectory.c_str());
		printf("vbkg = %.5f, %.5f, %.5f\n", vbkg[0], vbkg[1], vbkg[2]);
		printf("EB = %.5f\n", vesicles[0].EB);
		printf("Viscosity ratio = %.5f\n", viscRat);
		printf("Surface Area = %.5f\n", vesicles[0].areaTar);
		printf("Volume = %.5f\n", vesicles[0].volTar);
		printf("Reduced Volume = %.5f\n", 6*sqrt(3.14159265359)*vesicles[0].volTar*pow(vesicles[0].areaTar, -1.5));
	}







    for (; lt <= Nt; lt++) {
        double wtimeBgn = MPI_Wtime();

	if (mpi_rank == 0) printf("lt = %9d time = %.5f\n", lt, time);

    //if(mpi_rank == 0) {cout << vesicles[0].verts[0].x[0] << "c\n";}

	// Domain decomposition
	for (int ives = 0; ives < numVesicles(); ives++) {
	    Vesicle &vesicle = vesicles[ives];
	    vesicle.updateGeometry();
	    vesicle.calcDoubleLayerJump();
	    vesicle.bendForce(vesicle.fbend);
	    vesicle.tensionForce(vesicle.sigma, vesicle.ftens);
	}

    //if(mpi_rank == 0) {cout << vesicles[0].verts[0].x[0] << "d\n";}


	for (int iwall = 0; iwall < numWalls(); iwall++) {
	    Wall &wall = walls[iwall];
	    wall.updateGeometry();
	}
    //if(mpi_rank == 0) {cout << vesicles[0].verts[0].x[0] << "e\n";}

	domainDecomp();
    //if(mpi_rank == 0) {cout << vesicles[0].verts[0].x[0] << "f\n";}

	updateNeighborList();
   // if(mpi_rank == 0) {cout << vesicles[0].verts[0].x[0] << "g\n";}

	calcBieMat();
    //if(mpi_rank == 0) {cout << vesicles[0].verts[0].x[0] << "h\n";}

	// Solve the wall friction force
	Vec wall_rhs = getWallBieRhs();
    //if(mpi_rank == 0) {cout << vesicles[0].verts[0].x[0] << "i\n";}

	calcRhs_wallBie(wall_rhs);
    //if(mpi_rank == 0) {cout << vesicles[0].verts[0].x[0] << "j\n";}

	solveWallBie();
    //if(mpi_rank == 0) {cout << vesicles[0].verts[0].x[0] << "k\n";}

	// Solve the vesicle velocity
	{
	    // Prediction
	    Vec rhs_vel = getVesicleBieRhs();
	        //if(mpi_rank == 0) {cout << vesicles[0].verts[0].x[0] << "l\n";}

	    Vec sol_vel = getVesicleBieSol();
    //if(mpi_rank == 0) {cout << vesicles[0].verts[0].x[0] << "m\n";}

	    VecScatter scatter_vel = getVesicleBieScatter();
    //if(mpi_rank == 0) {cout << vesicles[0].verts[0].x[0] << "n\n";}

	    calcRhs_vesicleBie(rhs_vel);
    //if(mpi_rank == 0) {cout << vesicles[0].verts[0].x[0] << "o\n";}

	    solveVesicleBie(Ts);
   // if(mpi_rank == 0) {cout << vesicles[0].verts[0].x[0] << "p\n";}

	    // Projection
	    int n = vertList[VESICLE].size();
	    MArray<double,2> vold(n,3), vnew(n,3);
	    MArray<double,1> divTar(n), lbd(n);

	    for (int p = 0, ives = 0; ives < numVesicles(); ives++) {
		Vesicle &vesicle = vesicles[ives];
		int nvert = vesicle.numVerts();

		const double tau = param::getDoubleValue("AREA_RELAX_TIME");
		double s = (vesicle.areaTar/vesicle.area - 1.0)/tau;

		for (int ivert = 0; ivert < nvert; ivert++) {
		    divTar(p) = s*vesicle.vertArea(ivert);
		    p++;
		}
	    }
   // if(mpi_rank == 0) {cout << vesicles[0].verts[0].x[0] << "q\n";}

	    myVecScatter(scatter_vel, sol_vel, vold.data(), INSERT_VALUES, SCATTER_FORWARD);
	    projectVel(vold, divTar, vnew, lbd);

	    for (int p = 0, ives = 0; ives < numVesicles(); ives++) {
		Vesicle &vesicle = vesicles[ives];
		int nvert = vesicle.numVerts();

		m_dcopy(3*nvert, &vnew(p,0), vesicle.v.data());
		m_dadd(nvert, &lbd(p), vesicle.sigma.data());

		p += nvert;
	    }
	        //if(mpi_rank == 0) {cout << vesicles[0].verts[0].x[0] << "r\n";}

        }

	// Probe velocity
	calcProbeVel(xprb, vprb);
    //if(mpi_rank == 0) {cout << vesicles[0].verts[0].x[0] << "s\n";}

	// Output the state before system evolves
	writeAll();

    // New feature: write Tension and Pressure file
    if(mpi_rank == 0) 
		{
    		strcpy(token, "VESICLE_OUT");
		    nout = param::exist(token)? param::getIntValue(token) : -1;
		    if (lt%nout == 0) 
		    {
	    	    
                /////////////////////////////
                //Write pressure at each face
				sprintf(tensionfn, "%s%6.6d%s", "D/pressure", lt, ".dat");
				FILE *pressfile = fopen(tensionfn, "w");

		    	for (int ives = 0; ives < numVesicles(); ives++) 
		    	{
		    		Vesicle &vesicle = vesicles[ives];		
					fprintf(pressfile, "Pressure at each face\n");

					for(int ifa = 0; ifa < vesicle.faces.size(); ifa++)
					{   //This loop is a localized version of calcStress()
					    Tri &face = vesicle.faces[ifa];
					    m_dclear(9, *stresstau);

					    double xtri[3][3], ftri[3][3], vtri[3][3]; 

					    for (int l = 0; l < 3; l++) 
					    {
							int ivert = face.ivert[l];
						    m_dcopy(3, vesicle.verts[ivert].x, xtri[l]);

							m_dclear(3, ftri[l]);
							m_dadd(3, &vesicle.fbend(ivert,0), ftri[l]);
							m_dadd(3, &vesicle.ftens(ivert,0), ftri[l]);

							m_dcopy(3, &vesicle.v(ivert,0), vtri[l]);
						}

					    Quad2D &Q = quadrature::select_rule_2d("TRI_3");
					    for (int iq = 0; iq < Q.n(); iq++) 
					    {
						    double s = Q.x(iq), t = Q.y(iq);
							double r = 1.0 - s - t;

						    double xq[3], fq[3], vq[3], dA;
							FOR_I3 {
							    xq[i] = r*xtri[0][i] + s*xtri[1][i] + t*xtri[2][i] - vesicle.center[i];
							    fq[i] = r*ftri[0][i] + s*ftri[1][i] + t*ftri[2][i];
							    vq[i] = r*vtri[0][i] + s*vtri[1][i] + t*vtri[2][i];
						        }
							dA = Q.w(iq)*face.detJ;


                            //Technically, the below statement could be a single loop over the diagonal, but I leave in all of the stress tensor for clarity 
			                FOR_I3
			                FOR_J3 {  
				                    stresstau[i][j] += 0.5*(fq[i]*xq[j] + fq[j]*xq[i])*dA;
								    stresstau[i][j] += (viscRat - 1.0)*
						    		(vq[i]*face.normal[j] + vq[j]*face.normal[i])*dA;
					            }

					            
					    } // iq

                        facelocalpressure = -(stresstau[0][0] + stresstau[1][1] + stresstau[2][2])*0.33333333333333333;
						fprintf(pressfile, "%.6f\n", facelocalpressure);
					} //ifa
					fprintf(pressfile, "\n");
					fclose(pressfile);
			    } //ives
		    }//endif nout
		}//endif mpi_rank


	// Evolve the system
	for (int ives = 0; ives < numVesicles(); ives++) {
	    Vesicle &vesicle = vesicles[ives];
	    int nvert = vesicle.numVerts();

	    MArray<double,2> vrelax(nvert,3);
	    vesicle.calcMeshRelaxVel(vrelax);

	    double UT[3], OMG[3];
	    vesicle.calcTransRotatVel(vesicle.v, UT, OMG);

	    for (int ivert = 0; ivert < nvert; ivert++) {
	        Point &vert = vesicle.verts[ivert];

		double v[3];
		FOR_D3 v[d] = vesicle.v(ivert,d) - UT[d];

		const double *nrml = &vesicle.vertNrml(ivert,0);
		double vn = m_ddot(3, nrml, v);

		FOR_D3 v[d] = UT[d] + vn*nrml[d] + vrelax(ivert,d);

		FOR_D3 vert.x[d] += Ts*v[d];
	    }

	    // Adjust volume
	    vesicle.updateAVC();
	    double Vtar = max(0.99*vesicle.vol, min(1.01*vesicle.vol, vesicle.volTar));
	    double s = pow(Vtar/vesicle.vol, 1.0/3.0);
	    for (int ivert = 0; ivert < nvert; ivert++) {
	        Point &vert = vesicle.verts[ivert];

		FOR_D3 {
		    double xx = vert.x[d] - vesicle.center[d];
		    vert.x[d] = vesicle.center[d] + s*xx;
	        }
	    }
	}



	//======================================================================
	// Post-process
	
	// Prevent cell overlapping
	updateSourceGeometry();
	vector<NbrList*> nlists;
	nlists.push_back(&nlist_phys[VESICLE][VESICLE]);
	collision::forceSeparation(vertList[VESICLE], nlists, 0.03);

	rebox();
	syncCoord();
	time += Ts;

	// Time cost
        double wtimeEnd = MPI_Wtime();
	if (mpi_rank == 0) {
	    printf("    total wtime = %7.2f s\n", wtimeEnd - wtimeBgn);
	    printf("S Area = %.5f\n", vesicles[0].areaTar);
		printf("Volume = %.5f\n", vesicles[0].volTar);
	}

	// Check whether to kill job
	int killjob = 0;
	if (mpi_rank == 0) {
	    const char fn[] = "killjob";
	    struct stat fstat;
	    killjob = (stat(fn, &fstat) == 0);
	    if (killjob) {
	        printf("Found file killjob, terminate job\n");
	        remove(fn);
	    }
	}
	MPI_Bcast(&killjob, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (killjob) break;
    } // for
} 


// Rebox vesicles
void VesWall::rebox()
{
    using ewald::L;
    using ewald::iL;

    vector<Mesh*> meshes;
    meshes.insert(meshes.end(), meshList[VESICLE].begin(), meshList[VESICLE].end());

    for (int imesh = 0; imesh < meshes.size(); imesh++) {
        Mesh *mesh = meshes[imesh];

        // Determine translating distance
        double xmin[3], xmax[3], xx[3];
        mesh->getCoordRange(xmin, xmax);

        FOR_I3  {
            double eps = 0.01*L[i];
            double xc = 0.5*(xmin[i] + xmax[i]);
            xx[i] = (xc > -eps && xc < L[i] + eps) ?  0.0 : floor(xc*iL[i])*L[i];
        }

        // Translate vesicle
        for (int ivert = 0; ivert < mesh->numVerts(); ivert++) {
	    Point &vert = mesh->verts[ivert];
	    FOR_I3 vert.x[i] -= xx[i];
	}
    } // imesh
}


// Synchronize vesicle coordinates
void VesWall::syncCoord()
{
    // Assemble global coordinate arrays
    vector<Mesh*> meshes;
    meshes.insert(meshes.end(), meshList[VESICLE].begin(), meshList[VESICLE].end());

    int npoint = 0;
    for (int imesh = 0; imesh < meshes.size(); imesh++)
        npoint += meshes[imesh]->numVerts();

    double (*x)[3] = new double[npoint][3];

    for (int imesh = 0, p = 0; imesh < meshes.size(); imesh++) {
        Mesh *mesh = meshes[imesh];
	for (int ivert = 0; ivert < mesh->numVerts(); ivert++) {
	    m_dcopy(3, mesh->verts[ivert].x, x[p]);
	    p++;
	}
    }

    // Broadcast from root node
    MPI_Bcast(*x, 3*npoint, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Update vertex coordinates
    for (int imesh = 0, p = 0; imesh < meshes.size(); imesh++) {
        Mesh *mesh = meshes[imesh];
	for (int ivert = 0; ivert < mesh->numVerts(); ivert++) {
	    m_dcopy(3, x[p], mesh->verts[ivert].x);
	    p++;
	}
    }

    // Dealloc temp arrays
    delete [] x;
}


/* Compute the particles extra stress
 * Arguments:
 *   tau -- the particle stress, normalized by vesicle volume fraction
 * Note:
 *   -- Must calculate fbend and ftens for each vesicle before calling
 *   -- Only useful for Couette flow calculation */
void VesWall::calcStress(double (*tau)[3])
{
    // Init
    m_dclear(9, *tau);
    double vol_sum = 0.0;

    for (int ives = 0; ives < vesicles.size(); ives++) 
    {
        Vesicle &vesicle = vesicles[ives];

		for (int ifa = 0; ifa < vesicle.faces.size(); ifa++) 
		{
		    Tri &face = vesicle.faces[ifa];

		    double xtri[3][3], ftri[3][3], vtri[3][3]; 

		    for (int l = 0; l < 3; l++) 
		    {
				int ivert = face.ivert[l];
			        m_dcopy(3, vesicle.verts[ivert].x, xtri[l]);

				m_dclear(3, ftri[l]);
				m_dadd(3, &vesicle.fbend(ivert,0), ftri[l]);
				m_dadd(3, &vesicle.ftens(ivert,0), ftri[l]);

				m_dcopy(3, &vesicle.v(ivert,0), vtri[l]);
			}

		    Quad2D &Q = quadrature::select_rule_2d("TRI_3");
		    for (int iq = 0; iq < Q.n(); iq++) 
		    {
			    double s = Q.x(iq), t = Q.y(iq);
				double r = 1.0 - s - t;

			    double xq[3], fq[3], vq[3], dA;
				FOR_I3 {
				    xq[i] = r*xtri[0][i] + s*xtri[1][i] + t*xtri[2][i] - vesicle.center[i];
				    fq[i] = r*ftri[0][i] + s*ftri[1][i] + t*ftri[2][i];
				    vq[i] = r*vtri[0][i] + s*vtri[1][i] + t*vtri[2][i];
			        }
				dA = Q.w(iq)*face.detJ;

                FOR_I3
                FOR_J3 {
	                    tau[i][j] += 0.5*(fq[i]*xq[j] + fq[j]*xq[i])*dA;
					    tau[i][j] += (viscRat - 1.0)*
			    		(vq[i]*face.normal[j] + vq[j]*face.normal[i])*dA;
		            }
		    } // iq
		} // ifa

	    vol_sum += vesicle.vol;
    } // ives

    double s = 1.0/vol_sum;
    s /= fabs(shRate) + 1.E-10;
    m_dscal(9, s, *tau);
} 


/* Calculate probe velocity 
 * Arguments:
 *  x -- coordinates of probes
 *  v -- velocities at probe points 
 * Note:
 *  - Assume matched viscosity
 *  - Neglect the surface integral on particles 
 *  - Assume that the vesicle and wall surface force density are pre-computed 
 *    and stored in f[][3] member arrays */
void VesWall::calcProbeVel(MArray<double,2> &x, MArray<double,2> &v)
{
    int n = x.size(0);
    if (n <= 0) return;

    // Create a global PETSC vector for velocity calculation
    static Vec vec;
    static VecScatter scatt;
    static bool vec_inited = false;

    // First check if n == (size of v_vec)
    if (vec_inited) {
        int size_vec;
	VecGetSize(vec, &size_vec);
	if (size_vec != 3*n) {
	    VecDestroy(&vec);
	    VecScatterDestroy(&scatt);
	    vec_inited = false;
	}
    }

    if (! vec_inited) {
        VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, 3*n, &vec);
        myVecScatterCreateToAll(vec, scatt);
	vec_inited = true;
    }

    // Set sources
    for (int ives = 0; ives < numVesicles(); ives++) {
	Vesicle &vesicle = vesicles[ives];
	for (int ivert = 0; ivert < vesicle.numVerts(); ivert++) {
	    Point &vert = vesicle.verts[ivert];
	    FOR_D3 {
		vert.f[d] = vesicle.fbend(ivert,d) + vesicle.ftens(ivert,d);
		vert.g[d] = vesicle.v(ivert,d);
	    }
	}
    }

    for (int iwall = 0; iwall < numWalls(); iwall++) {
	Wall &wall = walls[iwall];
	for (int ivert = 0; ivert < wall.numVerts(); ivert++) {
	    Point &vert = wall.verts[ivert];
	    FOR_D3 {
		vert.f[d] = wall.f(ivert,d);
		vert.g[d] = 0.0;
	    }
	}
    }

    // Add boundary integrals
    VecZeroEntries(vec);      
    double cf[NBLK] = { -1.0/(8*M_PI), -1.0/(8*M_PI) };
    double cg[NBLK] = { (1.0 - viscRat)/(8*M_PI), 0.0 };
    addBoundaryIntegrals(x, vec, cf, cg);
    myVecScatter(scatt, vec, v.data(), INSERT_VALUES, SCATTER_FORWARD);

    // Add background velocity
    for (int i = 0; i < n; i++) {
        FOR_D3 v(i,d) += vbkg[d];
	v(i,0) += shRate*(x(i,2) - 0.5*ewald::L[2]);
    }
}


// Write everything
void VesWall::writeAll()
{
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    const char FN_FMT[] = "%s%6.6d%s";
    char token[256], fn[256];
    int nout;
    const int LT_CHUNK = 10000;

    // Vesicles
    strcpy(token, "VESICLE_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	sprintf(fn, FN_FMT, "D/vesicle", lt, ".dat");
	if (mpi_rank == 0) writeVesicles(fn);
    }

    // Walls
    strcpy(token, "WALL_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	sprintf(fn, FN_FMT, "D/wall", lt, ".dat");
	if (mpi_rank == 0) writeWalls(fn);
    }

    // Center of vesicles
    strcpy(token, "VESICLE_CENTER_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
        sprintf(fn, FN_FMT, "D/vesicle_center", (lt/LT_CHUNK)*LT_CHUNK, ".dat");
	if (mpi_rank == 0) writeVesicleCenter(fn);
    }

    // Stress
    strcpy(token, "STRESS_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	sprintf(fn, FN_FMT, "D/stress", (lt/LT_CHUNK)*LT_CHUNK, ".dat");
	if (mpi_rank == 0) writeStress(fn);
    }

    // Restart
    strcpy(token, "RESTART_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	sprintf(fn, FN_FMT, "D/restart", lt, ".dat");
	if (mpi_rank == 0) writeRestart(fn);
    }
}


/* Write all vesicles to file
 * Arguments:
 *   fn -- file name */
void VesWall::writeVesicles(const char *fn)
{
    if (numVesicles() <= 0) return;

    FILE *file = fopen(fn, "w");
    fprintf(file, "# TIME = %f\n", time);

    fprintf(file, "variables = x, y, z, un, vn, wn, tension\n");

    for (int ives = 0; ives < numVesicles(); ives++) {
        Vesicle &vesicle = vesicles[ives];
	int nvert = vesicle.numVerts();
	int nface = vesicle.numFaces();
	fprintf(file, "zone N=%d E=%d F=FEPOINT ET=TRIANGLE\n", nvert, nface);

	// coordinates
	for (int ivert = 0; ivert < nvert; ivert++) 
	{
	    Point &vert = vesicle.verts[ivert];
	    fprintf(file, " %10.7f %10.7f %10.7f", vert.x[0], vert.x[1], vert.x[2]);

		
        //Write normal velocity
        const double *nrml = &vesicle.vertNrml(ivert,0);
		double vn = m_ddot(3, nrml, &vesicle.v(ivert,0));
	    fprintf(file, " %.7f %.7f %.7f", vn*vesicle.vertNrml(ivert,0), vn*vesicle.vertNrml(ivert,1), vn*vesicle.vertNrml(ivert,2));



		//Write tension at each vertex
		fprintf(file, " %.7f\n", vesicle.sigma(ivert));

    }




	// connectivity
	for (int iface = 0; iface < nface; iface++) {
	    Tri &face = vesicle.faces[iface];
	    fprintf(file, " %d %d %d\n", face.ivert[0]+1, face.ivert[1]+1, face.ivert[2]+1);
        }
    }

    fclose(file);
}

/* Write all vesicles to file, just positions, no velocity or tension
 * Arguments:
 *   fn -- file name */
void VesWall::writeVesiclesplain(const char *fn)
{
    if (numVesicles() <= 0) return;

    FILE *file = fopen(fn, "w");
    fprintf(file, "# TIME = %f\n", time);

    fprintf(file, "variables = x, y, z\n");

    for (int ives = 0; ives < numVesicles(); ives++) {
        Vesicle &vesicle = vesicles[ives];
	int nvert = vesicle.numVerts();
	int nface = vesicle.numFaces();
	fprintf(file, "zone N=%d E=%d F=FEPOINT ET=TRIANGLE\n", nvert, nface);

	// coordinates
	for (int ivert = 0; ivert < nvert; ivert++) 
	{
	    Point &vert = vesicle.verts[ivert];
	    fprintf(file, " %10.7f %10.7f %10.7f\n", vert.x[0], vert.x[1], vert.x[2]);
    }




	// connectivity
	for (int iface = 0; iface < nface; iface++) {
	    Tri &face = vesicle.faces[iface];
	    fprintf(file, " %d %d %d\n", face.ivert[0]+1, face.ivert[1]+1, face.ivert[2]+1);
        }
    }

    fclose(file);
}




/* Write all walls to file
 * Arguments:
 *   fn -- file name */
void VesWall::writeWalls(const char *fn)
{
    if (numWalls() <= 0) return;

    FILE *file = fopen(fn, "w");
    fprintf(file, "# TIME = %f\n", time);

    fprintf(file, "VARIABLES = X, Y, Z, FX, FY, FZ\n");

    for (int iwall = 0; iwall < numWalls(); iwall++) {
        Wall &wall = walls[iwall];
	int nvert = wall.numVerts();
	int nface = wall.numFaces();
	bool write_force = (wall.f.size(0) == nvert);

	fprintf(file, "ZONE N=%d E=%d F=FEPOINT ET=TRIANGLE\n", nvert, nface);

	// coordinates
	for (int ivert = 0; ivert < nvert; ivert++) {
	    Point &vert = wall.verts[ivert];

	    fprintf(file, " %10.3E %10.3E %10.3E", vert.x[0], vert.x[1], vert.x[2]);

	    if (write_force) {
	        fprintf(file, " %10.7E %10.7E %10.7E", 
		    wall.f(ivert,0), wall.f(ivert,1), wall.f(ivert,2));
	    } else {
	        fprintf(file, " %10.3E %10.3E %10.3E", 0.0, 0.0, 0.0);
	    }

	    fprintf(file, "\n");
	}

	// connectivity
	for (int iface = 0; iface < nface; iface++) {
	    Tri &face = wall.faces[iface];
	    fprintf(file, "%d  %d  %d\n", 
	    		face.ivert[0]+1, face.ivert[1]+1, face.ivert[2]+1);
	}
    }

    fclose(file);
}


void VesWall::writeRestart(const char *fn)
{
    hid_t fid;
    hsize_t dims[2];
    char token[256];

    fid = H5Fcreate(fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    dims[0] = 3;
    H5LTmake_dataset_double(fid, "EWALD_L", 1, dims, ewald::L);

    dims[0] = 1;
    H5LTmake_dataset_double(fid, "TIME", 1, dims, &time);

    dims[0] = 1;
    H5LTmake_dataset_int(fid, "LT", 1, dims, &lt);

    dims[0] = 3;
    H5LTmake_dataset_double(fid, "VBKG", 1, dims, vbkg);

    dims[0] = 1;
    H5LTmake_dataset_double(fid, "SHRATE", 1, dims, &shRate);

    int nvesicle = numVesicles();
    if (nvesicle > 0) {
	dims[0] = 1;
	H5LTmake_dataset_int(fid, "NVESICLE", 1, dims, &nvesicle);
    }

    int nwall = numWalls();
    if (nwall > 0) {
	dims[0] = 1;
	H5LTmake_dataset_int(fid, "NWALL", 1, dims, &nwall);
    }

    // Vesicles 
    for (int ives = 0; ives < numVesicles(); ives++) {
        Vesicle &vesicle = vesicles[ives];
	int nvert = vesicle.numVerts();
	int nface = vesicle.numFaces();
	double (*x)[3] = new double[nvert][3];
	int (*f2v)[3] = new int[nface][3];
	vesicle.getCoords(x);
	vesicle.getConnectivities(f2v);

        sprintf(token, "VESICLE%d", ives);
	hid_t gid = H5Gcreate(fid, token, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	dims[0] = nvert;
	dims[1] = 3;
        H5LTmake_dataset_double(gid, "X", 2, dims, *x);

	dims[0] = nface;
	dims[1] = 3;
	H5LTmake_dataset_int(gid, "F2V", 2, dims, *f2v);

	H5Gclose(gid);
	delete []x;
	delete []f2v;
    }

    // walls
    for (int iwall = 0; iwall < numWalls(); iwall++) {
        Wall &wall = walls[iwall];
	int nvert = wall.numVerts();
	int nface = wall.numFaces();
	double (*x)[3] = new double[nvert][3];
	int (*f2v)[3] = new int[nface][3];

	wall.getCoords(x);
	wall.getConnectivities(f2v);

        sprintf(token, "WALL%d", iwall);
	hid_t gid = H5Gcreate(fid, token, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	dims[0] = nvert;
	dims[1] = 3;
        H5LTmake_dataset_double(gid, "X", 2, dims, *x);

	dims[0] = nface;
	dims[1] = 3;
	H5LTmake_dataset_int(gid, "F2V", 2, dims, *f2v);

	H5Gclose(gid);
	delete [] x;
	delete [] f2v;
    }

    // Probes
    int nprb = xprb.size(0);
    if (nprb > 0) {
	dims[0] = nprb;
	dims[1] = 3;
	H5LTmake_dataset_double(fid, "PROBE", 2, dims, xprb.data());
    }

    H5Fclose(fid);
}


/* Read restart file in HDF5 format
 * Arguments:
 *  fn -- file name */
void VesWall::readRestart(const char *fn)
{
    const MPI_Comm mpi_comm = MPI_COMM_WORLD;
    int mpi_size = 1;
    int mpi_rank = 0;

    int mpi_inited;
    MPI_Initialized(&mpi_inited);
    if (mpi_inited) {
	MPI_Comm_size(mpi_comm, &mpi_size);
	MPI_Comm_rank(mpi_comm, &mpi_rank);
    }

    hid_t fid, gid;
    hsize_t dims[2];
    H5T_class_t class_id;
    size_t type_size;
    char token[256];

    if (mpi_rank == 0) {
	fid = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (fid < 0)  printf("Can not open %s\n", fn);
    }
    if (mpi_size > 1) {
	MPI_Bcast(&fid, sizeof(hid_t), MPI_BYTE, 0, mpi_comm);
    }
    if (fid < 0) exit(1);

    if (mpi_rank == 0) {
	H5LTread_dataset_double(fid, "EWALD_L", ewald::L);
	H5LTread_dataset_double(fid, "TIME", &time);
	H5LTread_dataset_int(fid, "LT", &lt);
	H5LTread_dataset_double(fid, "VBKG", vbkg);
	H5LTread_dataset_double(fid, "SHRATE", &shRate);

	printf("L = %.3f %.3f %.3f\n", ewald::L[0], ewald::L[1], ewald::L[2]);
	printf("time = %.3f\n", time);
	printf("lt = %d\n", lt);
	printf("vbkg = %.3E %.3E %.3E\n", vbkg[0], vbkg[1], vbkg[2]);
	printf("shRate = %.3E\n", shRate);
    }
    if (mpi_size > 1) {
	MPI_Bcast(ewald::L, 3, MPI_DOUBLE, 0, mpi_comm);
	MPI_Bcast(&time, 1, MPI_DOUBLE, 0, mpi_comm);
	MPI_Bcast(&lt, 1, MPI_INT, 0, mpi_comm);
	MPI_Bcast(vbkg, 3, MPI_DOUBLE, 0, mpi_comm);
	MPI_Bcast(&shRate, 1, MPI_DOUBLE, 0, mpi_comm);
    }


    int nvesicle =0;
    int nwall =0; 
    if (mpi_rank == 0) {
	if (H5LTfind_dataset(fid, "NVESICLE")) H5LTread_dataset_int(fid, "NVESICLE", &nvesicle);
	if (H5LTfind_dataset(fid, "NWALL")) H5LTread_dataset_int(fid, "NWALL", &nwall);

	printf("nvesicle = %d\n", nvesicle);
	printf("nwall = %d\n", nwall);
    }
    if (mpi_size > 1) {
	MPI_Bcast(&nvesicle, 1, MPI_INT, 0, mpi_comm);
	MPI_Bcast(&nwall, 1, MPI_INT, 0, mpi_comm);
    }

    vesicles.resize(nvesicle);
    walls.resize(nwall);

    // Read vesicles
    for (int ives = 0; ives < nvesicle; ives++) {
        Vesicle &vesicle = vesicles[ives];
        int nvert, nface;

	if (mpi_rank == 0) {
	    sprintf(token, "VESICLE%d", ives);
	    gid = H5Gopen(fid, token, H5P_DEFAULT);

	    H5LTget_dataset_info(gid, "X", dims, &class_id, &type_size);
	    nvert = dims[0];

	    H5LTget_dataset_info(gid, "F2V", dims, &class_id, &type_size);
	    nface = dims[0];
        }
	if (mpi_size > 1) {
	    MPI_Bcast(&nvert, 1, MPI_INT, 0, mpi_comm);
	    MPI_Bcast(&nface, 1, MPI_INT, 0, mpi_comm);
	}

        MArray<double,2> x(nvert,3);
	MArray<int,2> f2v(nface,3);

	if (mpi_rank == 0) {
	    H5LTread_dataset_double(gid, "X", x.data());
	    H5LTread_dataset_int(gid, "F2V", f2v.data());
        }
	if (mpi_size > 1) {
	    MPI_Bcast(x.data(), x.size(), MPI_DOUBLE, 0, mpi_comm);
	    MPI_Bcast(f2v.data(), f2v.size(), MPI_INT, 0, mpi_comm);
	}

	vesicle.verts.resize(nvert);
	vesicle.setCoords(x);

	vesicle.faces.resize(nface);
	vesicle.setConnectivities(f2v);

	if (mpi_rank == 0) H5Gclose(gid);
    } 

    // wall
    for (int iwall = 0; iwall < nwall; iwall++) {
        Wall &wall = walls[iwall];
        int nvert, nface;

	if (mpi_rank == 0) {
	    sprintf(token, "WALL%d", iwall);
	    gid = H5Gopen(fid, token, H5P_DEFAULT);

	    H5LTget_dataset_info(gid, "X", dims, &class_id, &type_size);
	    nvert = dims[0];

	    H5LTget_dataset_info(gid, "F2V", dims, &class_id, &type_size);
	    nface = dims[0];
        }
	if (mpi_size > 1) {
	    MPI_Bcast(&nvert, 1, MPI_INT, 0, mpi_comm);
	    MPI_Bcast(&nface, 1, MPI_INT, 0, mpi_comm);
	}

	MArray<double,2> x(nvert,3);
	MArray<int,2> f2v(nface,3);

	if (mpi_rank == 0) {
	    H5LTread_dataset_double(gid, "X", x.data());
	    H5LTread_dataset_int(gid, "F2V", f2v.data());
        }
	if (mpi_size > 1) {
	    MPI_Bcast(x.data(), x.size(), MPI_DOUBLE, 0, mpi_comm);
	    MPI_Bcast(f2v.data(), f2v.size(), MPI_INT, 0, mpi_comm);
	}

	wall.verts.resize(nvert);
	wall.setCoords(x);

	wall.faces.resize(nface);
	wall.setConnectivities(f2v);

	if (mpi_rank == 0) H5Gclose(gid);
    }

    // Probes
    int nprb = 0;

    if (mpi_rank == 0) {
	strcpy(token, "PROBE");
	if (H5LTfind_dataset(fid, token)) {
	    H5LTget_dataset_info(fid, token, dims, &class_id, &type_size);
	    nprb = dims[0];
        }

	printf("nprobe = %d\n", nprb);
    }
    if (mpi_size > 1) {
	MPI_Bcast(&nprb, 1, MPI_INT, 0, mpi_comm);
    }

    xprb.resize(nprb,3);
    vprb.resize(nprb,3);

    if (nprb > 0) {
	if (mpi_rank == 0) {
	    H5LTread_dataset_double(fid, token, xprb.data());
        }
	if (mpi_size > 1) {
	    MPI_Bcast(xprb.data(), 3*nprb, MPI_DOUBLE, 0, mpi_comm);
	}
    }

    if (mpi_rank == 0) H5Fclose(fid);
}


/* Write probes */
void VesWall::writeProbe(const char *fn)
{
    int nprb = xprb.size(0);
    if (nprb <= 0) return;

    hid_t fid, dset;
    int ierr;
    hsize_t dims[2];
    char token[256];

    // Disable error message
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    if (H5Fis_hdf5(fn) > 0) {
        fid = H5Fopen(fn, H5F_ACC_RDWR, H5P_DEFAULT);
    }
    else {
        fid = H5Fcreate(fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }

    sprintf(token, "X");
    if (! H5LTfind_dataset(fid, token)) {
        dims[0] = nprb;
        dims[1] = 3;
        H5LTmake_dataset_double(fid, token, 2, dims, xprb.data());
    } else {
        dset = H5Dopen(fid, token, H5P_DEFAULT);
        H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xprb.data());
        H5Dclose(dset);
    }

    sprintf(token, "V%6.6d", lt);
    if (! H5LTfind_dataset(fid, token)) {
        dims[0] = nprb;
        dims[1] = 3;
        H5LTmake_dataset_double(fid, token, 2, dims, vprb.data());
    }
    else {
        dset = H5Dopen(fid, token, H5P_DEFAULT);
        H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vprb.data());
        H5Dclose(dset);
    }

    H5Fclose(fid);
}


/* Write the center and velocity of vesicles */
void VesWall::writeVesicleCenter(const char *fn)
{
    int nvesicle = numVesicles();
    if (nvesicle <= 0) return;

    MArray<double,2> x(nvesicle,3), v(nvesicle,3);

    for (int ives = 0; ives < nvesicle; ives++) {
        Vesicle &vesicle = vesicles[ives];
	FOR_I3 x(ives,i) = vesicle.center[i];
	vesicle.calcCentroidVel(vesicle.v, &v(ives,0));
    }

    // Write to file
    hid_t fid, dset;
    int ierr;
    hsize_t dims[2];
    char token[256];

    // No error message
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    if (H5Fis_hdf5(fn) > 0) {
	fid = H5Fopen(fn, H5F_ACC_RDWR, H5P_DEFAULT);
    }
    else {
	fid = H5Fcreate(fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }

    if (! H5LTfind_dataset(fid, "NVESICLE")) {
        dims[0] = 1;
	H5LTmake_dataset_int(fid, "NVESICLE", 1, dims, &nvesicle);
    }

    sprintf(token, "X%6.6d", lt);
    if (! H5LTfind_dataset(fid, token)) {
        dims[0] = nvesicle;
	dims[1] = 3;
	H5LTmake_dataset_double(fid, token, 2, dims, x.data());
    }
    else {
        dset = H5Dopen(fid, token, H5P_DEFAULT);
	H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, x.data());
	H5Dclose(dset);
    }

    sprintf(token, "V%6.6d", lt);
    if (! H5LTfind_dataset(fid, token)) {
        dims[0] = nvesicle;
	dims[1] = 3;
	H5LTmake_dataset_double(fid, token, 2, dims, v.data());
    }
    else {
        dset = H5Dopen(fid, token, H5P_DEFAULT);
	H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
	H5Dclose(dset);
    }
    
    H5Fclose(fid);
}


// Write stress
void VesWall::writeStress(const char *fn)
{
    // Need surface force distribution first
    for (int ives = 0; ives < numVesicles(); ives++) {
        Vesicle &vesicle = vesicles[ives];
	vesicle.bendForce(vesicle.fbend);
	vesicle.tensionForce(vesicle.sigma, vesicle.ftens);
    }

    double tau[3][3];
    calcStress(tau);

    // Write to file
    hid_t fid, dset;
    int ierr;
    hsize_t dims[2];
    char token[256];

    // No error message
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    if (H5Fis_hdf5(fn) > 0) {
	fid = H5Fopen(fn, H5F_ACC_RDWR, H5P_DEFAULT);
    }
    else {
	fid = H5Fcreate(fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }

    sprintf(token, "TAU%6.6d", lt);
    if (! H5LTfind_dataset(fid, token)) {
        dims[0] = 3;
        dims[1] = 3;
	H5LTmake_dataset_double(fid, token, 2, dims, *tau);
    }
    else {
        dset = H5Dopen(fid, token, H5P_DEFAULT);
	H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *tau);
	H5Dclose(dset);
    }

    H5Fclose(fid);
}

// Hematocrit, i.e. vesicle volume fraction
double VesWall::volFraction()
{
    double vol = 0.0;
    for (int ives = 0; ives < numVesicles(); ives++) {
        vol += vesicles[ives].vol;
    }

    return vol/ewald::vol;
}
