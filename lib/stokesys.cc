#include "cxxheaders.h"
#include "stokesys.h"
#include "ewald.h"
#include "mathfunc.h"
#include "quadrature.h"
#include "param.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "debugfunc.h"


/* StokeSys::StokeSys 
 * StokeSys::~StokeSys 
 * StokeSys::initialHook 
 * StokeSys::calcFixedSTNList
 * StokeSys::domainDecomp
 * StokeSys::buildNeighborList
 * StokeSys::setMeshOwnership
 * StokeSys::timeInt
 * StokeSys::rebox
 * StokeSys::syncCoord
 * StokeSys::calcPressGrad
 * StokeSys::calcCellVolFrac
 * StokeSys::calcParticleDensity
 * StokeSys::calcProbeVel 
 * StokeSys::writeAll 
 * StokeSys::writeCells 
 * StokeSys::writeRigids 
 * StokeSys::writeWalls 
 * StokeSys::writeRestart 
 * StokeSys::readRestart 
 * StokeSys::writeProbe
 * StokeSys::writeCellCenter 
 * StokeSys::writeRigidCenter 
 * StokeSys::cellVolFraction
 * StokeSys::calcAreaChange */

// Subsystem identifies
namespace {
    const int NBLK = 3;
};


StokeSys::StokeSys()
{
    nprb = 0;
    xprb = NULL;
    vprb = NULL;

    ntrac = 0;
    xtrac = NULL;
    vtrac = NULL;
}


StokeSys::~StokeSys()
{
    delete [] xprb;
    delete [] vprb;

    delete [] xtrac;
    delete [] vtrac;
}

/* -- Index meshes, vertices, and faces
 * -- Build vertex lists and face lists
 * -- Compute geometries of reference cells */
void StokeSys::initialHook()
{
    // Build mesh list and index meshes
    for (int i = 0; i < numCells(); i++) {
        Vesicle &cell = cells[i];

        meshList[CELL].push_back(&cell);

	cell.Gindx = i;
	cell.setInternalPointers();
	cell.connectPeriodicBoundaries();

	/*Cell *cellRef = cell.cellRef;  //RBC reference state
	if (cellRef) {
	    cellRef->setInternalPointers();
	    cellRef->connectPeriodicBoundaries();
	    // Only need to do it once
	    cellRef->buildEdgeList();
	    cellRef->updateGeometry();
	}*/
    }

    for (int i = 0; i < numRigids(); i++) {
        Rigid &rigid = rigids[i];

	meshList[RIGID].push_back(&rigid);

	rigid.Gindx = i;
	rigid.setInternalPointers();
	rigid.connectPeriodicBoundaries();
    }

    for (int i = 0; i < numWalls(); i++) {
	Wall &wall = walls[i];

	meshList[WALL].push_back(&wall);

	wall.Gindx = i;
	wall.setInternalPointers();
	wall.connectPeriodicBoundaries(ewald::L);
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
    }


    //======================================================================
    // Allocate f and v array for each mesh
    for (int icell = 0; icell < numCells(); icell++) {
        Vesicle &cell = cells[icell];

	int nvert = cell.numVerts();
	cell.f.resize(nvert,3);
	cell.v.resize(nvert,3);
    }

    for (int irigid = 0; irigid < numRigids(); irigid++) {
        Rigid &rigid = rigids[irigid];

	int nvert = rigid.numVerts();
	rigid.f.resize(nvert,3);
	rigid.v.resize(nvert,3);
    }

    for (int iwall = 0; iwall < numWalls(); iwall++) {
        Wall &wall = walls[iwall];

	int nvert = wall.numVerts();
	wall.f.resize(nvert,3);
    }
}


// Pre-compute the source, target, and neighbor lists that does not vary with time
void StokeSys::calcFixedSTNList()
{
    static bool inited = false;
    if (inited) return;
    inited = true;

    // Rigid-self
    for (int irigid = 0; irigid < numRigids(); irigid++) {
	Rigid &rigid = rigids[irigid];
	rigid.updateGeometry();
    }

    if (ewald::phys::active) {
	slist_RigidSelf.build(faceList[RIGID]);
	// NO_SELF = false, ONLY_SELF = true
	nlist_RigidSelf.build(vertList[RIGID], slist_RigidSelf, false, true);
    }

    // Wall
    for (int iwall = 0; iwall < numWalls(); iwall++) {
	Wall &wall = walls[iwall];
	wall.updateGeometry();
    }

    if (ewald::phys::active) {
        slist_phys[WALL].build(faceList[WALL]);
	nlist_phys[WALL][WALL].build(vertList[WALL], slist_phys[WALL]);
    }

    if (ewald::four::active) {
        ewald::four::setSources(faceList[WALL], slist_four[WALL]);
	ewald::four::setTargets(vertList[WALL], tlist_four[WALL]);
    }
}



/* Build source and target lists
 * Note:
 *   - The fixed lists are taken care by calcFixedSTNList */
void StokeSys::domainDecomp()
{
    // Physical 
    if (ewald::phys::active) {
        for (int blk = 0; blk < NBLK; blk++) {
	    if (blk == WALL) continue;
	    slist_phys[blk].build(faceList[blk]);
	}
    }

    // Fourier
    if (ewald::four::active) {
        for (int blk = 0; blk < NBLK; blk++) {
	    if (blk == WALL) continue;
            ewald::four::setSources(faceList[blk], slist_four[blk]);
            ewald::four::setTargets(vertList[blk], tlist_four[blk]);
        }
    }
}


// Set mesh active and private flags
void StokeSys::setMeshOwnership()
{
    for (int blk = 0; blk < NBLK; blk++) {
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

	// Need to include rigid self interactions
	if (blk == RIGID) {
	    for (int i = 0; i < slist_RigidSelf.numFaces(); i++) {
	        Mesh *mesh = slist_RigidSelf.faces[i]->mesh;
	        mesh->isActive = 1;
	    }
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
}


/*   Update the geometry of meshes that are active */
void StokeSys::updateSourceGeometry()
{
    for (int icell = 0; icell < numCells(); icell++) {
        Vesicle &cell = cells[icell];
	if (cell.isActive) cell.updateGeometry();
    }

    for (int irigid = 0; irigid < numRigids(); irigid++) {
        Rigid &rigid = rigids[irigid];
	if (rigid.isActive) rigid.updateGeometry();
    }

    // Assume that Wall geometry is calcualted in calcFixedSTNList
}


// Update neighbor list
void StokeSys::buildNeighborList()
{
    if (ewald::phys::active) {
        for (int iblk = 0; iblk < NBLK; iblk++) 
	for (int jblk = 0; jblk < NBLK; jblk++) {
	    NbrList &nlist = nlist_phys[iblk][jblk];

	    if (iblk == WALL && jblk == WALL) {
	        continue;
	    }
	    else if (iblk == RIGID && jblk == RIGID) {
		// Exclude rigid self interactions
	        nlist.build(vertList[iblk], slist_phys[jblk], true);	// NO_SELF = true
		continue;
            } 
	    else {
		nlist.build(vertList[iblk], slist_phys[jblk]);
	    }
        }
    }
}


// Time integration
void StokeSys::timeInt()
{
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    // Init
    for (int ives = 0; ives < numCells(); ives++) {
	Vesicle &vesicle = cells[ives];
	vesicle.v = 0.0;          // initialize velocity and surface tension to zero
	vesicle.sigma = 0.0;
    }

    for (; lt <= Nt; lt++) {
        double wtimeBgn = MPI_Wtime();

	if (mpi_rank == 0) printf("lt = %9d time = %.5f\n", lt, time);

	// Domain decomposition
	for (int icell = 0; icell < numCells(); icell++) {
	    Vesicle &cell = cells[icell];
	    cell.updateGeometry();
	    cell.calcDoubleLayerJump();
	    cell.bendForce(cell.fbend);
	    cell.tensionForce(cell.sigma, cell.ftens);

//	    // debug
//	    // Stir the system
//	    for (int ivert = 0; ivert < cell.numVerts(); ivert++) {
//	        Point &vert = cell.verts[ivert];
//		cell.f(ivert,0) += sin(2*M_PI*ewald::iL[1]*vert.x[1]);
//		cell.f(ivert,1) += sin(2*M_PI*ewald::iL[2]*vert.x[2]);
//		cell.f(ivert,2) += sin(2*M_PI*ewald::iL[0]*vert.x[0]);
//	    }
//	    // end debug
	}

	for (int irigid = 0; irigid < numRigids(); irigid++) {
	    Rigid &rigid = rigids[irigid];
	    rigid.updateGeometry();
	    rigid.calcDoubleLayerJump();
	}

	calcFixedSTNList();
	domainDecomp();
	setMeshOwnership();
	buildNeighborList();

	solveFlow(); // this still needs to be changed (7/25/14)

	if (nprb > 0) calcProbeVel(nprb, xprb, vprb);
	if (ntrac > 0) calcProbeVel(ntrac, xtrac, vtrac);

	// Output the state before system evolves
	writeAll();

	// Diagnose
	double Ht = cellVolFraction();

	//double dAmin, dAmax, dAmean;
	//calcAreaChange(dAmin, dAmax, dAmean);

	if (mpi_rank == 0) {
	    printf("    cell vol frac = %.3f\n", Ht);
	    //printf("    area change (min,max,mean) = %.3f %.3f %.3f\n", 
	    	//	dAmin, dAmax, dAmean);
	}

	// Evolve the system
	/*for (int icell = 0; icell < numCells(); icell++) {
	    Vesicle &cell = cells[icell];
	    for (int ivert = 0; ivert < cell.numVerts(); ivert++) {
	        m_daxpy(3, Ts, &cell.v(ivert,0), cell.verts[ivert].x);
	    }
	}*/

	//Move vesicles
	for (int ives = 0; ives < numCells(); ives++) {
	    Vesicle &vesicle = cells[ives];
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


	}


	for (int irigid = 0; irigid < numRigids(); irigid++) {
	    Rigid &rigid = rigids[irigid];
	    rigid.move(rigid.v, Ts);
	}

	for (int itrac = 0; itrac < ntrac; itrac++) {
	    FOR_I3 xtrac[itrac][i] += Ts*vtrac[itrac][i];
	}

	//======================================================================
	// Post-process
	// Collision detect
	updateSourceGeometry();
	cellNoContact(0.03);

	// Correct cell volume
	for (int ives = 0; ives < numCells(); ives++) {
	    Vesicle &vesicle = cells[ives];
	    	    // Adjust volume
	    vesicle.updateAVC();
	    int nvert = vesicle.numVerts();
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

	rebox();
	syncCoord();
	time += Ts;

	// Time cost
        double wtimeEnd = MPI_Wtime();
	if (mpi_rank == 0) {
	    printf("    total wtime = %7.2f s\n", wtimeEnd - wtimeBgn);
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


// Rebox cells and rigid particles
void StokeSys::rebox()
{
    using ewald::L;
    using ewald::iL;

    vector<Mesh*> meshes;
    meshes.insert(meshes.end(), meshList[CELL].begin(), meshList[CELL].end());
    meshes.insert(meshes.end(), meshList[RIGID].begin(), meshList[RIGID].end());

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

        // Translate cell
        for (int ivert = 0; ivert < mesh->numVerts(); ivert++) {
	    Point &vert = mesh->verts[ivert];
	    FOR_I3 vert.x[i] -= xx[i];
	}
    } // imesh

    for (int itrac = 0; itrac < ntrac; itrac++) {
        FOR_I3  {
            double eps = 0.01*L[i];

	    if (xtrac[itrac][i] < eps  || xtrac[itrac][i] > L[i] + eps) {
	        double xx = floor(xtrac[itrac][i]*iL[i])*L[i];
		xtrac[itrac][i] -= xx;
	    }
        }
    }
}


// Synchronize cell and rigid particle coordinates
void StokeSys::syncCoord()
{
    // Assemble global coordinate arrays
    vector<Mesh*> meshes;
    meshes.insert(meshes.end(), meshList[CELL].begin(), meshList[CELL].end());
    meshes.insert(meshes.end(), meshList[RIGID].begin(), meshList[RIGID].end());

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

    // Synchronize tracer point coordinates
    MPI_Bcast(*xtrac, 3*ntrac, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Dealloc temp arrays
    delete [] x;
}


/* Compute the pressure gradient
 * Arguments:
 *  pgrad -- the pressure gradient */
void StokeSys::calcPressGrad(double *pgrad)
{
    // init
    pgrad[0] = pgrad[1] = pgrad[2] = 0.0;

    for (int iwall = 0; iwall < numWalls(); iwall++) {
        Wall &wall = walls[iwall];

	for (int iface = 0; iface < wall.numFaces(); iface++) {
	    Tri &face = wall.faces[iface];
	    double f[3] = {0.0};

	    for (int l = 0; l < 3; ++l)  {
	        int ivert = face.ivert[l];
	        m_dadd(3, &wall.f(ivert,0), f);
            }

	    m_daxpy(3, THRD*face.area, f, pgrad);
	}
    }

    m_dscal(3, ewald::iVol, pgrad);  //Note to jmb: may need to change use of ewald::iVol for non-flat channel geometries
}


/* Compute the cell volume fraction
 * Arguments:
 *   zmin, zmax -- domain limit
 *   nz -- number of intervals
 *   vf -- volume fraction ratio
 * Note:
 *   1. The force density array must be preset */
void StokeSys::calcCellVolFrac(double zmin, double zmax, int nz, double *vf)
{
    // Generate cutting planes
    double dz = (zmax - zmin)/nz;
    vector<double> zcut(nz);
    for (int iz = 0; iz < nz; iz++) zcut[iz] = zmin + (iz + 0.5)*dz;

    double znormal[3] = {0.0, 0.0, 1.0};

    // Init
    m_dclear(nz, vf);

    for (int icell = 0; icell < cells.size(); icell++) {
        Vesicle &cell = cells[icell];

	// Slice the cell
	// error check, VCell = cell.vol after iteration
	double VCell = 0.0;	

	for (int ifa = 0; ifa < cell.faces.size(); ifa++) {
	    Tri &face = cell.faces[ifa];

	    double xtri[3][3]; 
	    for (int l = 0; l < 3; l++) {
		int ivert = face.ivert[l];
	        m_dcopy(3, cell.verts[ivert].x, xtri[l]);
	    }

	    int nsub;
	    int *npsub = new int[nz];
	    double (*xsub)[5][3] = new double[nz][5][3];
	    slice_triangle(xtri, znormal, zcut, nsub, npsub, xsub);

	    for (int isub = 0; isub < nsub; isub++) {
	        double areaSub, xcSub[3];
		calcPolygonAreaCenter(npsub[isub], xsub[isub], areaSub, xcSub);

		// Contribution to the cell volume
		double dV = areaSub*xcSub[1]*face.normal[1];

		// Error checking
		VCell += dV;

	        // Add volume fraction and osmotic in the slab
		int iz = (int)floor((xcSub[2] - zmin)/dz);
		if (iz < 0 || iz >= nz) continue;
		double Vslab = ewald::vol/nz;

		vf[iz] += dV/Vslab;
	    } // isub

	    // Clean up
	    delete [] npsub;
	    delete [] xsub;
	} // ifa

	// Error checking to see wheather the slicing is done correctly
	if (fabs(VCell - cell.vol) > 1.E-10*cell.vol) {
	    cout << "StokeSys::calcCellVolFrac: ";
	    cout << "cell is not sliced properly" << endl;
	}
    } // icell
} 


/* Compute the particles extra stress
 * Arguments:
 *   zmin, zmax -- domain limit
 *   nz -- number of intervals
 *   Sxz, Sxx, Syy, Szz -- extra particle stress
 *   fx, fy, fz -- extra residual force
 * Note:
 *   1. The force density array must be preset */
void StokeSys::calcStress(double zmin, double zmax, int nz, 
		double *Sxz, double *Sxx, double *Syy, double *Szz,
		double *fx, double *fy, double *fz)
{
    // Generate cutting planes
    double dz = (zmax - zmin)/nz;
    vector<double> zcut(nz);
    for (int iz = 0; iz < nz; iz++) zcut[iz] = zmin + (iz + 0.5)*dz;

    double znormal[3] = {0.0, 0.0, 1.0};

    // Init
    m_dclear(nz, Sxz);
    m_dclear(nz, Sxx);
    m_dclear(nz, Syy);
    m_dclear(nz, Szz);

    m_dclear(nz, fx);
    m_dclear(nz, fy);
    m_dclear(nz, fz);

    for (int icell = 0; icell < cells.size(); icell++) {
        Vesicle &cell = cells[icell];
        cell.vesaddtension();
	double tau[3][3] = {0.0};

	// Compute the total osmotic pressure from this cell
	double opCell = 0.0;

	for (int ifa = 0; ifa < cell.faces.size(); ifa++) {
	    Tri &face = cell.faces[ifa];

	    double xtri[3][3], ftri[3][3]; 
	    for (int l = 0; l < 3; l++) {
		int ivert = face.ivert[l];
	        m_dcopy(3, cell.verts[ivert].x, xtri[l]);
		m_dcopy(3, &cell.f(ivert,0), ftri[l]);
	    }

	    Quad2D &Q = quadrature::select_rule_2d("TRI_3");
	    for (int iq = 0; iq < Q.n(); iq++) {
	        double s = Q.x(iq), t = Q.y(iq);
		double r = 1.0 - s - t;

	        double xq[3], fq[3], dA;
		FOR_I3 {
		    xq[i] = r*xtri[0][i] + s*xtri[1][i] + t*xtri[2][i] - cell.center[i];
		    fq[i] = r*ftri[0][i] + s*ftri[1][i] + t*ftri[2][i];
	        }
		dA = Q.w(iq)*face.detJ;

                FOR_I3
                FOR_J3 
                    tau[i][j] += 0.5*(fq[i]*xq[j] + fq[j]*xq[i])*dA;

	    } // iq
	} // ifa

	// Slice the cell
	double VCell = 0.0;	// error check, VCell = cell.vol after iteration

	for (int ifa = 0; ifa < cell.faces.size(); ifa++) {
	    Tri &face = cell.faces[ifa];

	    double xtri[3][3], ftri[3][3]; 
	    for (int l = 0; l < 3; l++) {
		int ivert = face.ivert[l];
	        m_dcopy(3, cell.verts[ivert].x, xtri[l]);
		m_dcopy(3, &cell.f(ivert,0), ftri[l]);
	    }

	    int nsub;
	    int *npsub = new int[nz];
	    double (*xsub)[5][3] = new double[nz][5][3];
	    slice_triangle(xtri, znormal, zcut, nsub, npsub, xsub);

	    for (int isub = 0; isub < nsub; isub++) {
	        double areaSub, xcSub[3];
		calcPolygonAreaCenter(npsub[isub], xsub[isub], areaSub, xcSub);

		// Contribution to the cell volume
		double dV = areaSub*xcSub[1]*face.normal[1];
		VCell += dV;

		int iz = (xcSub[2] - zmin)/dz;
		if (iz < 0 || iz >= nz) continue;

	        // Accumulate the stresslet
		double Vslab = ewald::vol/nz;

		double wght = (dV/cell.vol)/Vslab;
		Sxz[iz] += tau[0][2]*wght;
		Sxx[iz] += tau[0][0]*wght;
		Syy[iz] += tau[1][1]*wght;
		Szz[iz] += tau[2][2]*wght;

		// Accumulate the force
		double fmean[3];
		FOR_I3 fmean[i] = THRD*(ftri[0][i] + ftri[1][i] + ftri[2][i]);

		wght = areaSub/Vslab;
		fx[iz] += fmean[0]*wght;
		fy[iz] += fmean[1]*wght;
		fz[iz] += fmean[2]*wght;
	    } // isub

	    // Clean up
	    delete []npsub;
	    delete []xsub;
	} // ifa

	// Error checking to see wheather the slicing is done correctly
	if (fabs(VCell - cell.vol) > 1.E-10*cell.vol) {
	    cout << "StokeSys::calcOsmoticPress: ";
	    cout << "cell is not sliced properly" << endl;
	}
    } // icell
} 

/* Assign the while cell stresslet to the z-interval where the cell centroid
 * lies in
 * Arguments:
 *   zmin, zmax, nz
 *   SxzSum -- the sum of total stresslet in each interval
 *   SxxSum, SyySum, SzzSum -- same as SxzSum
 *   VolSum -- total cell volume in each interval */
void StokeSys::calcStressLocal(double zmin, double zmax, int nz, 
		double *SxzSum, double *SxxSum, double *SyySum, double *SzzSum,
		double *VolSum)
{
    // Generate cutting planes
    double dz = (zmax - zmin)/nz;

    // Init
    m_dclear(nz, SxzSum);
    m_dclear(nz, SxxSum);
    m_dclear(nz, SyySum);
    m_dclear(nz, SzzSum);
    m_dclear(nz, VolSum);

    for (int icell = 0; icell < cells.size(); icell++) {
        Vesicle &cell = cells[icell];
        cell.vesaddtension();
	double tau[3][3] = {0.0};

	// Extra particle stress 
	for (int ifa = 0; ifa < cell.faces.size(); ifa++) {
	    Tri &face = cell.faces[ifa];

	    double xtri[3][3], ftri[3][3]; 
	    for (int l = 0; l < 3; l++) {
		int ivert = face.ivert[l];
	        m_dcopy(3, cell.verts[ivert].x, xtri[l]);
		m_dcopy(3, &cell.f(ivert,0), ftri[l]);
	    }

	    Quad2D &Q = quadrature::select_rule_2d("TRI_3");
	    for (int iq = 0; iq < Q.n(); iq++) {
	        double s = Q.x(iq), t = Q.y(iq);
		double r = 1.0 - s - t;

	        double xq[3], fq[3], dA;
		FOR_I3 {
		    xq[i] = r*xtri[0][i] + s*xtri[1][i] + t*xtri[2][i] - cell.center[i];
		    fq[i] = r*ftri[0][i] + s*ftri[1][i] + t*ftri[2][i];
	        }
		dA = Q.w(iq)*face.detJ;

                FOR_I3
                FOR_J3 
                    tau[i][j] += 0.5*(fq[i]*xq[j] + fq[j]*xq[i])*dA;

	    } // iq
	} // ifa

	// Assign the stresslet to the slice where the center of mass lies in
	int iz = (int)floor((cell.center[2] - zmin)/dz);
	if (iz < 0 || iz >= nz) continue;

	VolSum[iz] += cell.vol;

	SxzSum[iz] += tau[0][2];
	SxxSum[iz] += tau[0][0];
	SyySum[iz] += tau[1][1];
	SzzSum[iz] += tau[2][2];
    } // icell
} 


/* Calculate rigid particles' number density profile
 * Arguments:
 *   zmin, zmax -- domain limit
 *   nz -- number of intervals
 *   dens -- rigid particle number densities */
void StokeSys::calcParticleDensity(double zmin, double zmax, int nz, double *dens)
{
    double idz = double(nz)/(zmax - zmin);

    // Initialize
    for (int iz = 0; iz < nz; iz++) dens[iz] = 0.0;

    for (int irigid = 0; irigid < rigids.size(); irigid++) {
        Rigid &rigid = rigids[irigid];
        double zc = rigid.center[2];

        int iz = (int)floor((zc - zmin)*idz);
        if (iz >= 0 && iz < nz) dens[iz]++;
    }

    for (int iz = 0; iz < nz; iz++) dens[iz] *= idz;
}


/* Calculate probe velocity 
 * Arguments:
 *  n -- number of probes
 *  x -- coordinates of probes
 *  v -- velocities at probe points 
 * Note:
 *  - Assume matched viscosity
 *  - Neglect the surface integral on particles 
 *  - Assume that the cell and wall surface force density are pre-computed 
 *    and stored in f[][3] member arrays */
void StokeSys::calcProbeVel(int n, const double (*x)[3], double (*v)[3])
{
    // Create a global PETSC vector for velocity calculation
    static Vec v_vec;
    static VecScatter scatt;
    static bool vec_inited = false;

    // First check if n == (size of v_vec)
    if (vec_inited) {
        int size_vec;
	VecGetSize(v_vec, &size_vec);
	if (size_vec != 3*n) {
	    VecDestroy(&v_vec);
	    VecScatterDestroy(&scatt);
	    vec_inited = false;
	}
    }

    if (! vec_inited) {
        VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, 3*n, &v_vec);
        myVecScatterCreateToAll(v_vec, scatt);
	vec_inited = true;
    }

    // Init
    VecZeroEntries(v_vec);      

    // Set probe points
    vector<Point> probes(n);
    vector<Point*> pprobes(n);
    for (int i = 0; i < n; i++) {
        FOR_D3 probes[i].x[d] = x[i][d];
	probes[i].Gindx = i;
	pprobes[i] = &probes[i];
    }

    // Set sources
    for (int icell = 0; icell < numCells(); icell++) {
	Vesicle &cell = cells[icell];
	for (int ivert = 0; ivert < cell.numVerts(); ivert++) {
	    Point &vert = cell.verts[ivert];
	    FOR_D3 {
		vert.f[d] = cell.f(ivert,d);
		vert.g[d] = 0.0;
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

    // Calculate surface integrals
    if (ewald::phys::active) {
        for (int blk = 0; blk < NBLK; blk++) {
            if (blk != CELL && blk != WALL) continue;

            NbrList nlist;
            nlist.build(pprobes, slist_phys[blk]);

            int nloc = nlist.verts.size();
            double *vloc = new double[3*nloc];
            int *ix = new int[3*nloc];

            for (int i = 0; i < nloc; i++) {
                int iglb = nlist.verts[i]->Gindx;
                FOR_D3 ix[3*i+d] = 3*iglb + d;
            }
            std::fill_n(vloc, 3*nloc, 0.0);
            ewald::phys::addSurfInt(nlist, -1.0/(8.0*M_PI), 0.0, vloc);
            VecSetValues(v_vec, 3*nloc, ix, vloc, ADD_VALUES);

            delete [] vloc;
            delete [] ix;
        }
    }

    if (ewald::four::active) {
        vector<Point*> tgts;
        ewald::four::setTargets(pprobes, tgts);

        int nloc = tgts.size();
        double *vloc = new double[3*nloc];
        int *ix = new int[3*nloc];
        for (int i = 0; i < nloc; i++) {
            int iglb = tgts[i]->Gindx;
            FOR_D3 ix[3*i+d] = 3*iglb + d;
        }

        ewald::four::clear_source();
        ewald::four::add_source(-1.0/(8*M_PI), 0.0, slist_four[CELL]);
        ewald::four::add_source(-1.0/(8*M_PI), 0.0, slist_four[WALL]);
        ewald::four::transform();

        std::fill_n(vloc, 3*nloc, 0.0);
        ewald::four::add_interp_vel(tgts, vloc);
        VecSetValues(v_vec, 3*nloc, ix, vloc, ADD_VALUES);

        delete [] vloc;
        delete [] ix;
    }

    VecAssemblyBegin(v_vec);
    VecAssemblyEnd(v_vec);
    myVecScatter(scatt, v_vec, *v, INSERT_VALUES, SCATTER_FORWARD);

    // Add background velocity
    for (int i = 0; i < n; i++) {
        FOR_D3 v[i][d] += vbkg[d];
    }
}


// Write everything
void StokeSys::writeAll()
{
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    const char FN_FMT[] = "%s%s%6.6d%s";
    char token[256], fn[256];
    int nout;
    const int LT_CHUNK = 10000;

    // Cells
    strcpy(token, "CELL_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	sprintf(fn, FN_FMT, outfiledirectory.c_str(), "/cell", lt, ".dat");
	if (mpi_rank == 0) writeCells(fn);
    }

    // Rigids
    strcpy(token, "RIGID_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	sprintf(fn, FN_FMT, outfiledirectory.c_str(), "/rigid", lt, ".dat");
	if (mpi_rank == 0) writeRigids(fn);
    }

    // Walls
    strcpy(token, "WALL_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	sprintf(fn, FN_FMT, outfiledirectory.c_str(), "/wall", lt, ".dat");
	if (mpi_rank == 0) writeWalls(fn);
    }

    // Restart
    strcpy(token, "RESTART_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	sprintf(fn, FN_FMT, outfiledirectory.c_str(), "/restart", lt, ".dat");
	if (mpi_rank == 0) writeRestart(fn);
    }

    // Center of cells
    strcpy(token, "CELL_CENTER_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
        sprintf(fn, FN_FMT, outfiledirectory.c_str(), "/cell_center", (lt/LT_CHUNK)*LT_CHUNK, ".dat");
        if (mpi_rank == 0) writeCellCenter(fn);
    }

    // Center of rigids
    strcpy(token, "RIGID_CENTER_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
        sprintf(fn, FN_FMT, outfiledirectory.c_str(), "/rigid_center", (lt/LT_CHUNK)*LT_CHUNK, ".dat");
        if (mpi_rank == 0) writeRigidCenter(fn);
    }

    // Probe
    strcpy(token, "PROBE_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	sprintf(fn, FN_FMT, outfiledirectory.c_str(), "/probe", (lt/LT_CHUNK)*LT_CHUNK, ".dat");
	if (mpi_rank == 0) writeProbe(fn);
    }

    // Tracer
    strcpy(token, "TRACER_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	sprintf(fn, FN_FMT, outfiledirectory.c_str(), "/tracer", (lt/LT_CHUNK)*LT_CHUNK, ".dat");
	if (mpi_rank == 0) writeTracer(fn);
    }
}


/* Write all cells to file
 * Arguments:
 *   fn -- file name */
void StokeSys::writeCells(const char *fn)
{
    if (numCells() <= 0) return;

    FILE *file = fopen(fn, "w");
//    fprintf(file, "variables = x, y, z, u, v, w\n");
    fprintf(file, "variables = x, y, z\n");

    for (int icell = 0; icell < numCells(); icell++) {
        Vesicle &cell = cells[icell];
	int nvert = cell.numVerts();
	int nface = cell.numFaces();
	fprintf(file, "zone N=%d E=%d F=FEPOINT ET=TRIANGLE\n", nvert, nface);
        
	// coordinates
	for (int ivert = 0; ivert < nvert; ivert++) {
	    Point &vert = cell.verts[ivert];
	    fprintf(file, "%10.3e %10.3e %10.3e\n", vert.x[0], vert.x[1], vert.x[2]);
 	    //fprintf(file, "%10.3e %10.3e %10.3e\n", 
	    //		cell.v[ivert][0], cell.v[ivert][1], cell.v[ivert][2]);
        }

	// connectivity
	for (int iface = 0; iface < nface; iface++) {
	    Tri &face = cell.faces[iface];
	    fprintf(file, "%d  %d  %d\n", face.ivert[0]+1, face.ivert[1]+1, face.ivert[2]+1);
        }
    }

    fclose(file);
}



/* Write all rigids to file
 * Arguments:
 *  fn -- file name */
void StokeSys::writeRigids(const char *fn)
{
    if (numRigids() <= 0) return;

    FILE *file = fopen(fn, "w");
    fprintf(file, "VARIABLES = X, Y, Z\n");

    for (int irigid = 0; irigid < numRigids(); irigid++) {
        Rigid &rigid = rigids[irigid];
	int nvert = rigid.numVerts();
	int nface = rigid.numFaces();
	fprintf(file, "ZONE N=%d E=%d F=FEPOINT ET=TRIANGLE\n", nvert, nface);
        
	// coordinates
	for (int ivert = 0; ivert < nvert; ivert++) {
	    Point &vert = rigid.verts[ivert];
	    fprintf(file, "%.5f  %.5f  %.5f\n", vert.x[0], vert.x[1], vert.x[2]);
        }

	// connectivity
	for (int iface = 0; iface < nface; iface++) {
	    Tri &face = rigid.faces[iface];
	    fprintf(file, "%d  %d  %d\n", face.ivert[0]+1, face.ivert[1]+1, face.ivert[2]+1);
        }
    }

    fclose(file);
}


/* Write all walls to file
 * Arguments:
 *   fn -- file name */
void StokeSys::writeWalls(const char *fn)
{
    if (numWalls() <= 0) return;

    FILE *file = fopen(fn, "w");
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

	    fprintf(file, "%.5f %.5f %.5f", vert.x[0], vert.x[1], vert.x[2]);

	    if (write_force)
	        fprintf(file, " %.5f %.5f %.5f\n", 
		    wall.f(ivert,0), wall.f(ivert,1), wall.f(ivert,2));
	    else
	        fprintf(file, " %.5f %.5f %.5f\n", 0.0, 0.0, 0.0);
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


void StokeSys::writeRestart(const char *fn)
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

    int ncell = numCells();
    if (ncell > 0) {
	dims[0] = 1;
	H5LTmake_dataset_int(fid, "NCELL", 1, dims, &ncell);
    }

    int nrigid = numRigids();
    if (nrigid > 0) {
	dims[0] = 1;
	H5LTmake_dataset_int(fid, "NRIGID", 1, dims, &nrigid);
    }

    int nwall = numWalls();
    if (nwall > 0) {
	dims[0] = 1;
	H5LTmake_dataset_int(fid, "NWALL", 1, dims, &nwall);
    }

    // Cells 
    for (int icell = 0; icell < numCells(); icell++) {
        Vesicle &cell = cells[icell];
	int nvert = cell.numVerts();
	int nface = cell.numFaces();
	double (*x)[3] = new double[nvert][3];
	int (*f2v)[3] = new int[nface][3];
	cell.getCoords(x);
	cell.getConnectivities(f2v);

        sprintf(token, "CELL%d", icell);
	hid_t gid = H5Gcreate(fid, token, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	dims[0] = nvert;
	dims[1] = 3;
        H5LTmake_dataset_double(gid, "X", 2, dims, *x);

	/*if (cell.cellRef != NULL) {
	    cell.cellRef->getCoords(x);
	    dims[0] = nvert;
	    dims[1] = 3;
	    H5LTmake_dataset_double(gid, "XREF", 2, dims, *x);
	}*/

	dims[0] = nface;
	dims[1] = 3;
	H5LTmake_dataset_int(gid, "F2V", 2, dims, *f2v);

	H5Gclose(gid);
	delete []x;
	delete []f2v;
    }

    // Rigids
    for (int irigid = 0; irigid < numRigids(); irigid++) {
        Rigid &rigid = rigids[irigid];

	int nvert = rigid.numVerts();
	int nface = rigid.numFaces();
	double (*x)[3] = new double[nvert][3];
	int (*f2v)[3] = new int[nface][3];

	rigid.getCoords(x);
	rigid.getConnectivities(f2v);

        sprintf(token, "RIGID%d", irigid);
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
    if (nprb > 0) {
	dims[0] = nprb;
	dims[1] = 3;
	H5LTmake_dataset_double(fid, "PROBE", 2, dims, *xprb);
    }

    // Tracers
    if (ntrac > 0) {
	dims[0] = ntrac;
	dims[1] = 3;
	H5LTmake_dataset_double(fid, "TRACER", 2, dims, *xtrac);
    }

    H5Fclose(fid);
}


/* Read restart file in HDF5 format
 * Arguments:
 *  fn -- file name */
void StokeSys::readRestart(const char *fn)
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

    //    if(mpi_rank ==0) cout << "In Stokesys:readRestart()\n";

    if (mpi_rank == 0) {
	fid = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT);
        //cout << "fid = " << fid << endl;
	if (fid < 0)  printf("Can not open %s\n", fn);
    }
    if (mpi_size > 1) {
	MPI_Bcast(&fid, sizeof(hid_t), MPI_BYTE, 0, mpi_comm);
    }
    if (fid < 0) exit(1);

    if (mpi_rank == 0) {
      //cout << "Reading file\n";
	H5LTread_dataset_double(fid, "EWALD_L", ewald::L);
	H5LTread_dataset_double(fid, "TIME", &time);
	H5LTread_dataset_int(fid, "LT", &lt);
	H5LTread_dataset_double(fid, "VBKG", vbkg);

	printf("L = %.3f %.3f %.3f\n", ewald::L[0], ewald::L[1], ewald::L[2]);
	printf("time = %.3f\n", time);
	printf("lt = %d\n", lt);
	printf("vbkg = %.3E %.3E %.3E\n", vbkg[0], vbkg[1], vbkg[2]);
    }
    if (mpi_size > 1) {
	MPI_Bcast(ewald::L, 3, MPI_DOUBLE, 0, mpi_comm);
	MPI_Bcast(&time, 1, MPI_DOUBLE, 0, mpi_comm);
	MPI_Bcast(&lt, 1, MPI_INT, 0, mpi_comm);
	MPI_Bcast(&vbkg, 3, MPI_DOUBLE, 0, mpi_comm);
    }


    int ncell =0;
    int nrigid =0; 
    int nwall =0; 
    if (mpi_rank == 0) {
	if (H5LTfind_dataset(fid, "NCELL")) H5LTread_dataset_int(fid, "NCELL", &ncell);
	if (H5LTfind_dataset(fid, "NRIGID")) H5LTread_dataset_int(fid, "NRIGID", &nrigid);
	if (H5LTfind_dataset(fid, "NWALL")) H5LTread_dataset_int(fid, "NWALL", &nwall);

	printf("ncell = %d\n", ncell);
	printf("nrigid = %d\n", nrigid);
	printf("nwall = %d\n", nwall);
    }
    if (mpi_size > 1) {
	MPI_Bcast(&ncell, 1, MPI_INT, 0, mpi_comm);
	MPI_Bcast(&nrigid, 1, MPI_INT, 0, mpi_comm);
	MPI_Bcast(&nwall, 1, MPI_INT, 0, mpi_comm);
    }

    cells.resize(ncell);
    rigids.resize(nrigid);
    walls.resize(nwall);

    // Read cells
    for (int icell = 0; icell < ncell; icell++) {
        Vesicle &cell = cells[icell];
        int nvert, nface;//, xref_flag;

	if (mpi_rank == 0) {
	    sprintf(token, "CELL%d", icell);
	    gid = H5Gopen(fid, token, H5P_DEFAULT);

	    H5LTget_dataset_info(gid, "X", dims, &class_id, &type_size);
	    nvert = dims[0];

	    H5LTget_dataset_info(gid, "F2V", dims, &class_id, &type_size);
	    nface = dims[0];

	    //xref_flag = H5LTfind_dataset(gid, "XREF");
        }
	if (mpi_size > 1) {
	    MPI_Bcast(&nvert, 1, MPI_INT, 0, mpi_comm);
	    MPI_Bcast(&nface, 1, MPI_INT, 0, mpi_comm);
	    //MPI_Bcast(&xref_flag, 1, MPI_INT, 0, mpi_comm);
	}

        MArray<double,2> x(nvert,3);//, xref(nvert,3);
	MArray<int,2> f2v(nface,3);

	if (mpi_rank == 0) {
	    H5LTread_dataset_double(gid, "X", x.data());
	    H5LTread_dataset_int(gid, "F2V", f2v.data());
	    //if (xref_flag) H5LTread_dataset_double(gid, "XREF", xref.data());
        }
	if (mpi_size > 1) {
	    MPI_Bcast(x.data(), x.size(), MPI_DOUBLE, 0, mpi_comm);
	    MPI_Bcast(f2v.data(), f2v.size(), MPI_INT, 0, mpi_comm);
	    //if (xref_flag) {
		//MPI_Bcast(xref.data(), xref.size(), MPI_DOUBLE, 0, mpi_comm);
	    //}
	}

	cell.verts.resize(nvert);
	cell.setCoords(x);

	cell.faces.resize(nface);
	cell.setConnectivities(f2v);

      /*  if (xref_flag) {
	    cell.cellRef = new Cell();

	    cell.cellRef->verts.resize(nvert);
	    cell.cellRef->setCoords(xref);

	    cell.cellRef->faces.resize(nface);
	    cell.cellRef->setConnectivities(f2v);
	}*/

	if (mpi_rank == 0) H5Gclose(gid);
    } 

    // Rigid
    for (int irigid = 0; irigid < nrigid; irigid++) {
        Rigid &rigid = rigids[irigid];
        int nvert, nface;

	if (mpi_rank == 0) {
	    sprintf(token, "RIGID%d", irigid);
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

	rigid.verts.resize(nvert);
	rigid.setCoords(x);

	rigid.faces.resize(nface);
	rigid.setConnectivities(f2v);

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
    nprb = 0;

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

    xprb = new double[nprb][3];
    vprb = new double[nprb][3];

    if (nprb > 0) {
	if (mpi_rank == 0) {
	    H5LTread_dataset_double(fid, token, *xprb);
        }
	if (mpi_size > 1) {
	    MPI_Bcast(xprb, 3*nprb, MPI_DOUBLE, 0, mpi_comm);
	}
    }

    // Tracers
    ntrac = 0;

    if (mpi_rank == 0) {
	strcpy(token, "TRACER");
	if (H5LTfind_dataset(fid, token)) {
	    H5LTget_dataset_info(fid, token, dims, &class_id, &type_size);
	    ntrac = dims[0];
        }
	printf("ntrace = %d\n", ntrac);
    }
    if (mpi_size > 1) {
	MPI_Bcast(&ntrac, 1, MPI_INT, 0, mpi_comm);
    }

    xtrac = new double[ntrac][3];
    vtrac = new double[ntrac][3];

    if (ntrac > 0) {
	if (mpi_rank == 0) {
	    H5LTread_dataset_double(fid, token, *xtrac);
	}
	if (mpi_size > 1) {
	    MPI_Bcast(xtrac, 3*ntrac, MPI_DOUBLE, 0, mpi_comm);
	}
    }

    if (mpi_rank == 0) H5Fclose(fid);
}


/* Read restart file in HDF5 format without outputting anything
 * Arguments:
 *  fn -- file name */
void StokeSys::readRestart_quiet(const char *fn)
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

//	printf("L = %.3f %.3f %.3f\n", ewald::L[0], ewald::L[1], ewald::L[2]);
//	printf("time = %.3f\n", time);
//	printf("lt = %d\n", lt);
//	printf("vbkg = %.3E %.3E %.3E\n", vbkg[0], vbkg[1], vbkg[2]);
    }
    if (mpi_size > 1) {
	MPI_Bcast(ewald::L, 3, MPI_DOUBLE, 0, mpi_comm);
	MPI_Bcast(&time, 1, MPI_DOUBLE, 0, mpi_comm);
	MPI_Bcast(&lt, 1, MPI_INT, 0, mpi_comm);
	MPI_Bcast(&vbkg, 3, MPI_DOUBLE, 0, mpi_comm);
    }


    int ncell =0;
    int nrigid =0; 
    int nwall =0; 
    if (mpi_rank == 0) {
	if (H5LTfind_dataset(fid, "NCELL")) H5LTread_dataset_int(fid, "NCELL", &ncell);
	if (H5LTfind_dataset(fid, "NRIGID")) H5LTread_dataset_int(fid, "NRIGID", &nrigid);
	if (H5LTfind_dataset(fid, "NWALL")) H5LTread_dataset_int(fid, "NWALL", &nwall);

//	printf("ncell = %d\n", ncell);
//	printf("nrigid = %d\n", nrigid);
//	printf("nwall = %d\n", nwall);
    }
    if (mpi_size > 1) {
	MPI_Bcast(&ncell, 1, MPI_INT, 0, mpi_comm);
	MPI_Bcast(&nrigid, 1, MPI_INT, 0, mpi_comm);
	MPI_Bcast(&nwall, 1, MPI_INT, 0, mpi_comm);
    }

    cells.resize(ncell);
    rigids.resize(nrigid);
    walls.resize(nwall);

    // Read cells
    for (int icell = 0; icell < ncell; icell++) {
        Vesicle &cell = cells[icell];
        int nvert, nface, xref_flag;

	if (mpi_rank == 0) {
	    sprintf(token, "CELL%d", icell);
	    gid = H5Gopen(fid, token, H5P_DEFAULT);

	    H5LTget_dataset_info(gid, "X", dims, &class_id, &type_size);
	    nvert = dims[0];

	    H5LTget_dataset_info(gid, "F2V", dims, &class_id, &type_size);
	    nface = dims[0];

	    //xref_flag = H5LTfind_dataset(gid, "XREF");
        }
	if (mpi_size > 1) {
	    MPI_Bcast(&nvert, 1, MPI_INT, 0, mpi_comm);
	    MPI_Bcast(&nface, 1, MPI_INT, 0, mpi_comm);
	    //MPI_Bcast(&xref_flag, 1, MPI_INT, 0, mpi_comm);
	}

        MArray<double,2> x(nvert,3);//, xref(nvert,3);
	MArray<int,2> f2v(nface,3);

	if (mpi_rank == 0) {
	    H5LTread_dataset_double(gid, "X", x.data());
	    H5LTread_dataset_int(gid, "F2V", f2v.data());
	    //if (xref_flag) H5LTread_dataset_double(gid, "XREF", xref.data());
        }
	if (mpi_size > 1) {
	    MPI_Bcast(x.data(), x.size(), MPI_DOUBLE, 0, mpi_comm);
	    MPI_Bcast(f2v.data(), f2v.size(), MPI_INT, 0, mpi_comm);
	    //if (xref_flag) {
		//MPI_Bcast(xref.data(), xref.size(), MPI_DOUBLE, 0, mpi_comm);
	    //}
	}

	cell.verts.resize(nvert);
	cell.setCoords(x);

	cell.faces.resize(nface);
	cell.setConnectivities(f2v);

        /*if (xref_flag) {
	    cell.cellRef = new Cell();

	    cell.cellRef->verts.resize(nvert);
	    cell.cellRef->setCoords(xref);

	    cell.cellRef->faces.resize(nface);
	    cell.cellRef->setConnectivities(f2v);
	}*/

	if (mpi_rank == 0) H5Gclose(gid);
    } 

    // Rigid
    for (int irigid = 0; irigid < nrigid; irigid++) {
        Rigid &rigid = rigids[irigid];
        int nvert, nface;

	if (mpi_rank == 0) {
	    sprintf(token, "RIGID%d", irigid);
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

	rigid.verts.resize(nvert);
	rigid.setCoords(x);

	rigid.faces.resize(nface);
	rigid.setConnectivities(f2v);

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
    nprb = 0;

    if (mpi_rank == 0) {
	strcpy(token, "PROBE");
	if (H5LTfind_dataset(fid, token)) {
	    H5LTget_dataset_info(fid, token, dims, &class_id, &type_size);
	    nprb = dims[0];
        }

//	printf("nprobe = %d\n", nprb);
    }
    if (mpi_size > 1) {
	MPI_Bcast(&nprb, 1, MPI_INT, 0, mpi_comm);
    }

    xprb = new double[nprb][3];
    vprb = new double[nprb][3];

    if (nprb > 0) {
	if (mpi_rank == 0) {
	    H5LTread_dataset_double(fid, token, *xprb);
        }
	if (mpi_size > 1) {
	    MPI_Bcast(xprb, 3*nprb, MPI_DOUBLE, 0, mpi_comm);
	}
    }

    // Tracers
    ntrac = 0;

    if (mpi_rank == 0) {
	strcpy(token, "TRACER");
	if (H5LTfind_dataset(fid, token)) {
	    H5LTget_dataset_info(fid, token, dims, &class_id, &type_size);
	    ntrac = dims[0];
        }
//	printf("ntrace = %d\n", ntrac);
    }
    if (mpi_size > 1) {
	MPI_Bcast(&ntrac, 1, MPI_INT, 0, mpi_comm);
    }

    xtrac = new double[ntrac][3];
    vtrac = new double[ntrac][3];

    if (ntrac > 0) {
	if (mpi_rank == 0) {
	    H5LTread_dataset_double(fid, token, *xtrac);
	}
	if (mpi_size > 1) {
	    MPI_Bcast(xtrac, 3*ntrac, MPI_DOUBLE, 0, mpi_comm);
	}
    }

    if (mpi_rank == 0) H5Fclose(fid);
}

/* Write probes */
void StokeSys::writeProbe(const char *fn)
{
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
        H5LTmake_dataset_double(fid, token, 2, dims, *xprb);
    } else {
        dset = H5Dopen(fid, token, H5P_DEFAULT);
        H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *xprb);
        H5Dclose(dset);
    }

    sprintf(token, "V%6.6d", lt);
    if (! H5LTfind_dataset(fid, token)) {
        dims[0] = nprb;
        dims[1] = 3;
        H5LTmake_dataset_double(fid, token, 2, dims, *vprb);
    }
    else {
        dset = H5Dopen(fid, token, H5P_DEFAULT);
        H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *vprb);
        H5Dclose(dset);
    }

    H5Fclose(fid);
}


void StokeSys::writeTracer(const char *fn)
{
    if (ntrac <= 0) return;

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

    sprintf(token, "X%6.6d", lt);
    if (! H5LTfind_dataset(fid, token)) {
        dims[0] = ntrac;
        dims[1] = 3;
        H5LTmake_dataset_double(fid, token, 2, dims, *xtrac);
    }
    else {
        dset = H5Dopen(fid, token, H5P_DEFAULT);
        H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *xtrac);
        H5Dclose(dset);
    }

    sprintf(token, "V%6.6d", lt);
    if (! H5LTfind_dataset(fid, token)) {
        dims[0] = ntrac;
        dims[1] = 3;
        H5LTmake_dataset_double(fid, token, 2, dims, *vtrac);
    }
    else {
        dset = H5Dopen(fid, token, H5P_DEFAULT);
        H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *vtrac);
        H5Dclose(dset);
    }

    H5Fclose(fid);
}


/* Write the center and velocity of cells */
void StokeSys::writeCellCenter(const char *fn)
{
    int ncell = numCells();
    if (ncell <= 0) return;

    MArray<double,2> x(ncell,3), v(ncell,3);

    for (int icell = 0; icell < ncell; icell++) {
        Vesicle &cell = cells[icell];
	FOR_I3 x(icell,i) = cell.center[i];
	cell.calcCentroidVel(cell.v, &v(icell,0));
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

    if (! H5LTfind_dataset(fid, "NCELL")) {
        dims[0] = 1;
	H5LTmake_dataset_int(fid, "NCELL", 1, dims, &ncell);
    }

    sprintf(token, "X%6.6d", lt);
    if (! H5LTfind_dataset(fid, token)) {
        dims[0] = ncell;
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
        dims[0] = ncell;
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


/* Write the center and velocity of rigid particles */
void StokeSys::writeRigidCenter(const char *fn)
{
    int nrigid = numRigids();
    if (nrigid <= 0) return;

    double (*x)[3], (*v)[3], (*omg)[3];
    x = new double[nrigid][3];
    v = new double[nrigid][3];
    omg = new double[nrigid][3];

    for (int irigid = 0; irigid < nrigid; irigid++) {
        Rigid &rigid = rigids[irigid];

	FOR_I3 x[irigid][i] = rigid.center[i];
	rigid.calcTransRotatVel(rigid.v, v[irigid], omg[irigid]);
    }

    // Write to file
    hid_t fid, dset;
    int ierr;
    hsize_t dims[2];
    char token[256];

    // Turn off error message
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    if (H5Fis_hdf5(fn) > 0) {
	fid = H5Fopen(fn, H5F_ACC_RDWR, H5P_DEFAULT);
    }
    else {
	fid = H5Fcreate(fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }

    if (! H5LTfind_dataset(fid, "NRIGID")) {
        dims[0] = 1;
	H5LTmake_dataset_int(fid, "NRIGID", 1, dims, &nrigid);
    }

    sprintf(token, "X%6.6d", lt);
    if (! H5LTfind_dataset(fid, token)) {
        dims[0] = nrigid;
	dims[1] = 3;
	H5LTmake_dataset_double(fid, token, 2, dims, *x);
    }
    else {
        dset = H5Dopen(fid, token, H5P_DEFAULT);
	H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *x);
	H5Dclose(dset);
    }

    sprintf(token, "V%6.6d", lt);
    if (! H5LTfind_dataset(fid, token)) {
        dims[0] = nrigid;
	dims[1] = 3;
	H5LTmake_dataset_double(fid, token, 2, dims, *v);
    }
    else {
        dset = H5Dopen(fid, token, H5P_DEFAULT);
	H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *v);
	H5Dclose(dset);
    }
    
    H5Fclose(fid);
    delete [] x;
    delete [] v;
    delete [] omg;
}


// Hematocrit, i.e. cell volume fraction
double StokeSys::cellVolFraction()
{
    double vol = 0.0;
    for (int icell = 0; icell < numCells(); icell++) {
        vol += cells[icell].vol;
    }

    return vol/ewald::vol;
}


// Max area change
/*
void StokeSys::calcAreaChange(double &dAmin, double &dAmax, double &dAmean)
{
    // Init
    dAmin = 1.E10;
    dAmax = -1.E10;

    double dAsum = 0.0, Asum = 0.0;

    for (int icell = 0; icell < numCells(); icell++) {
	Vesicle &cell = cells[icell];
	Vesicle *cellRef = cell.cellRef;
	if (cellRef == NULL) continue;

	for (int iface = 0; iface < cell.numFaces(); iface++) {
	    Tri &face = cell.faces[iface];
	    Tri &faceRef = cellRef->faces[iface];

	    double dA = face.area/faceRef.area - 1.0;

	    dAmin = min(dAmin, dA);
	    dAmax = max(dAmax, dA);

	    dAsum += fabs(face.area - faceRef.area);
	    Asum += faceRef.area;
        } // iface
    } // icell

    dAmean = dAsum/Asum;

    // Special case
    if (dAmin > dAmax) dAmin = dAmax = 0.0;
}*/
