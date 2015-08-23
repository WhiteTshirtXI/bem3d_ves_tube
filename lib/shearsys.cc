#include "cxxheaders.h"
#include "shearsys.h"
#include "ewald.h"
#include "mathfunc.h"
#include "quadrature.h"
#include "param.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "miscUtils.h"
#include "debugfunc.h"


/* ShearSys::ShearSys 
 * ShearSys::~ShearSys 
 * ShearSys::initialHook 
 * ShearSys::domainDecomp
 * ShearSys::updateSourceGeometry
 * ShearSys::buildNeighborList
 * ShearSys::timeInt
 * ShearSys::rebox
 * ShearSys::syncCoord
 * ShearSys::calcStress
 * ShearSys::calcProbeVel
 * ShearSys::writeAll 
 * ShearSys::writeCells 
 * ShearSys::writeRigids 
 * ShearSys::writeRestart 
 * ShearSys::readRestart 
 * ShearSys::writeProbe
 * ShearSys::writeCellCenter
 * ShearSys::writeRigidCenter 
 * ShearSys::writeStress 
 * ShearSys::cellVolFraction
 * ShearSys::calcMaxPerturbVel
 * ShearSys::calcAreaChange */

// Subsystem identifies
namespace {
    const int NBLK = 2;
};


ShearSys::ShearSys()
{
    nprb = 0;
    xprb = NULL;
    vprb = NULL;

    ntrac = 0;
    xtrac = NULL;
    vtrac = NULL;
}


ShearSys::~ShearSys()
{
    delete [] xprb;
    delete [] vprb;

    delete [] xtrac;
    delete [] vtrac;
}

/* -- Index meshes, vertices, and faces
 * -- Build vertex lists and face lists
 * -- Compute geometries of reference cells */
void ShearSys::initialHook()
{
    // Build mesh list and index meshes
    for (int i = 0; i < numCells(); i++) {
        Cell &cell = cells[i];
        meshList[CELL].push_back(&cell);

	cell.Gindx = i;
	cell.setInternalPointers();
	cell.noPeriodicBoundaries();

	Cell *cellRef = cell.cellRef;
	if (cellRef) {
	    cellRef->setInternalPointers();
	    cellRef->noPeriodicBoundaries();
	    cellRef->updateGeometry();
	}
    }

    for (int i = 0; i < numRigids(); i++) {
        Rigid &rigid = rigids[i];
	meshList[RIGID].push_back(&rigid);

	rigid.Gindx = i;
	rigid.setInternalPointers();
	rigid.noPeriodicBoundaries();
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
    // Allocate f and v arrays
    for (int icell = 0; icell < numCells(); icell++) {
        Cell &cell = cells[icell];
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
}


void ShearSys::domainDecomp()
{
    // Build source list
    if (ewald::phys::active) {
        for (int blk = 0; blk < NBLK; blk++) {
            slist_phys[blk].build(faceList[blk]);
        }
    }

    if (ewald::four::active) {
        for (int blk = 0; blk < NBLK; blk++) {
            ewald::four::setSources(faceList[blk], slist_four[blk]);
        }
    }

    // Build target list
    if (ewald::four::active) {
        for (int blk = 0; blk < NBLK; blk++) {
            ewald::four::setTargets(vertList[blk], tlist_four[blk]);
        }
    }

    // Set mesh active and private flags
    for (int blk = 0; blk < NBLK; blk++) {
	int nmesh = meshList[blk].size();
	if (nmesh == 0) return;

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


/*   Update the geometry of meshes that are active (i.e. containing source elements) */
void ShearSys::updateSourceGeometry()
{
    for (int icell = 0; icell < numCells(); icell++) {
        Cell &cell = cells[icell];
	if (cell.isActive) cell.updateGeometry();
    }

    for (int irigid = 0; irigid < numRigids(); irigid++) {
        Rigid &rigid = rigids[irigid];
	if (rigid.isActive) rigid.updateGeometry();
    }
}


/* Build neighbor list */
void ShearSys::buildNeighborList()
{
    if (ewald::phys::active) {
        for (int iblk = 0; iblk < NBLK; iblk++) 
	for (int jblk = 0; jblk < NBLK; jblk++) {
	    NbrList &nlist = nlist_phys[iblk][jblk];

	    if (iblk == RIGID && jblk == RIGID) {
	        // Exclude rigid self interactions
		// NO_SELF = true
		nlist.build(vertList[iblk], slist_phys[jblk], true);	
	    } else {
		nlist.build(vertList[iblk], slist_phys[jblk]);
	    }
	}
    }
}


// Time integration
void ShearSys::timeInt()
{
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    Nt0 = lt;

    for (; lt <= Nt; lt++) {
        double wtimeBgn = MPI_Wtime();
	if (mpi_rank == 0) printf("lt = %9d time = %.5f\n", lt, time);

        // Compute cell surface force
	for (int icell = 0; icell < numCells(); icell++) {
	    Cell &cell = cells[icell];
	    cell.updateGeometry();
	    cell.calcDoubleLayerJump();
	    cell.calcSurfForce(cell.f);

//	    // debug
//	    // Stir the system
//	    for (int ivert = 0; ivert < cell.numVerts(); ivert++) {
//		Point &vert = cell.verts[ivert];
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

	domainDecomp();
	buildNeighborList();

	solveFlow();
	
	if (nprb > 0) calcProbeVel(nprb, xprb, vprb);
	if (ntrac > 0) calcProbeVel(ntrac, xtrac, vtrac);

	// Output the state before system evolves
	writeAll();

	// Diagnose
	double Ht = cellVolFraction();

	double dAmin, dAmax, dAmean;
	calcAreaChange(dAmin, dAmax, dAmean);

	double dvmax[3], UTmax[3], OMGmax[3];
	calcMaxPerturbVel(dvmax, UTmax, OMGmax);

	if ( m_dmax(3, dvmax) > 100*max(1.0, fabs(shRate)*ewald::L[2])) {
	    if (mpi_rank == 0) {
	        printf("Numerical instability detected, stop\n");
	    }

	    break;
	}


	if (mpi_rank == 0) {
	    if(numCells() > 0) {
		printf("    cell vol frac   = %10.3f\n", Ht);
		printf("    area change (min,max,mean) = %10.3f %10.3f %10.3f\n", 
			    dAmin, dAmax, dAmean);
		printf("    max cell perturb vel = %10.3f %10.3f %10.3f\n", 
			    dvmax[0], dvmax[1], dvmax[2]);
            }

	    if (numRigids() > 0) {
		printf("    max rigid UT =  %10.3f %10.3f %10.3f\n", UTmax[0], UTmax[1], UTmax[2]);
		printf("    max rigid omg = %10.3f %10.3f %10.3f\n", OMGmax[0], OMGmax[1], OMGmax[2]);
	    }
	}


	// Evolve the system
	for (int icell = 0; icell < numCells(); icell++) {
	    Cell &cell = cells[icell];
	    for (int ivert = 0; ivert < cell.numVerts(); ivert++) {
	        m_daxpy(3, Ts, &cell.v(ivert,0), cell.verts[ivert].x);
	    }
	}

	for (int irigid = 0; irigid < numRigids(); irigid++) {
	    Rigid &rigid = rigids[irigid];
	    rigid.move(rigid.v, Ts);
	}


	for (int itrac = 0; itrac < ntrac; itrac++) {
	    FOR_I3 xtrac[itrac][i] += Ts*vtrac[itrac][i];
	}

	ewald::strain += Ts*shRate;
	ewald::updateThreeAxes();

	//======================================================================
	// Post-process
	// Prevent cell-cell collision
	updateSourceGeometry();
	cellNoContact(0.03);

	// Correct cell volume
	for (int icell = 0; icell < cells.size(); ++icell) {
	    Cell &cell = cells[icell];
	    cell.updateAVC();

	    double s = 1 + Ts*(cell.volTar/cell.vol - 1);
	    s = pow(s, 1.0/3.0);

	    for (int ivert = 0; ivert < cell.numVerts(); ++ivert) {
		Point &vert = cell.verts[ivert];
		FOR_I3 vert.x[i] = cell.center[i] + s*(vert.x[i] - cell.center[i]);
	    } // ivert
	}

	rebox();
	syncCoord();
	time += Ts;

	// Time cost
        double wtimeEnd = MPI_Wtime();
	if (mpi_rank == 0) {
	    printf("    wall time = %7.2f s\n", wtimeEnd - wtimeBgn);
	}

	// Check whether to kill job
	int killjob = 0;
	if (mpi_rank == 0) {
	    const char fn[] = "killjob";
	    killjob = miscUtils::fileExists(fn);
	    if (killjob) {
	        printf("Found file %s, terminate job\n", fn);
	        remove(fn);
	    }
	}
	MPI_Bcast(&killjob, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (killjob) break;
    } // lt
}


// Rebox cells and rigid particles
void ShearSys::rebox()
{
    using ewald::L;
    using ewald::iL;

    vector<Mesh*> meshes;
    meshes.insert(meshes.end(), meshList[CELL].begin(), meshList[CELL].end());
    meshes.insert(meshes.end(), meshList[RIGID].begin(), meshList[RIGID].end());

    for (int imesh = 0; imesh < meshes.size(); imesh++) {
        Mesh *mesh = meshes[imesh];

        // Determine translating distance
        double xmin[3], xmax[3], xc[3];

        mesh->getCoordRange(xmin, xmax);
	FOR_D3 xc[d] = 0.5*(xmin[d] + xmax[d]);

	bool move_mesh = false;
	FOR_D3 {
            double eps = 0.01*L[d];
	    if (xc[d] < -eps || xc[d] > L[d] + eps) move_mesh = true;
	}

	if (move_mesh) {
	    double origin[3];
	    FOR_D3 origin[d] = 0.5*ewald::L[d];

	    double xcnew[3];
	    FOR_D3 xcnew[d] = xc[d];
	    ewald::to_CloseBy(origin, xcnew);

	    double xx[3];
	    FOR_D3 xx[d] = xcnew[d] - xc[d];

	    for (int ivert = 0; ivert < mesh->numVerts(); ivert++) {
		Point &vert = mesh->verts[ivert];
		FOR_I3 vert.x[i] += xx[i];
	    }
	}
    } // imesh

    for (int itrac = 0; itrac < ntrac; itrac++) {
	bool move_mesh = false;
        FOR_D3  {
            double eps = 0.1*L[d];
	    if (xtrac[itrac][d] < -eps  || xtrac[itrac][d] > L[d] + eps) move_mesh = true;
        }

	if (move_mesh) {
	    double origin[3];
	    FOR_D3 origin[d] = 0.5*L[d];
	    ewald::to_CloseBy(origin, xtrac[itrac]); 
	}
    }
}


// Synchronize coordinates
void ShearSys::syncCoord()
{
    // Assemble global coordinate arrays
    vector<Mesh*> meshes;
    meshes.insert(meshes.end(), meshList[CELL].begin(), meshList[CELL].end());
    meshes.insert(meshes.end(), meshList[RIGID].begin(), meshList[RIGID].end());

    int n = 0;
    for (int imesh = 0; imesh < meshes.size(); imesh++)
        n += meshes[imesh]->numVerts();

    MArray<double,2> x(n,3);

    int p = 0;
    for (int imesh = 0; imesh < meshes.size(); imesh++) {
        Mesh *mesh = meshes[imesh];
	int nvert = mesh->numVerts();

	for (int ivert = 0; ivert < mesh->numVerts(); ivert++) {
	    Point &vert = mesh->verts[ivert];
	    FOR_J3 x(p,j) = vert.x[j];
	    p++;
	}
    }

    MPI_Bcast(x.data(), x.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    p = 0;
    for (int imesh = 0, p = 0; imesh < meshes.size(); imesh++) {
        Mesh *mesh = meshes[imesh];
	for (int ivert = 0; ivert < mesh->numVerts(); ivert++) {
	    Point &vert = mesh->verts[ivert];
	    FOR_J3 vert.x[j] = x(p,j);
	    p++;
	}
    }

    // Synchronize tracer point coordinates
    MPI_Bcast(*xtrac, 3*ntrac, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

/* Compute the particle stress
 * Arguments:
 *   tau -- added stress, normalized by (cell vol) * (shear rate)  */
void ShearSys::calcStress(double (*tau)[3])
{
    // Init
    double vol = 0.0;
    m_dclear(9, *tau);


    for (int icell = 0; icell < numCells(); icell++) {
        Cell &cell = cells[icell];

	for (int ifa = 0; ifa < cell.faces.size(); ifa++) {
	    Tri &face = cell.faces[ifa];

	    double xtri[3][3], ftri[3][3], vtri[3][3]; 
	    for (int l = 0; l < 3; l++) {
		int ivert = face.ivert[l];
	        m_dcopy(3, cell.verts[ivert].x, xtri[l]);
		m_dcopy(3, &cell.f(ivert,0), ftri[l]);
		m_dcopy(3, &cell.v(ivert,0), vtri[l]);
	    }

	    Quad2D &Q = quadrature::select_rule_2d("TRI_3");
	    for (int iq = 0; iq < Q.n(); iq++) {
	        double s = Q.x(iq), t = Q.y(iq);
		double r = 1.0 - s - t;

	        double xq[3], fq[3], vq[3], dA;
		FOR_I3 {
		    xq[i] = r*xtri[0][i] + s*xtri[1][i] + t*xtri[2][i] - cell.center[i];
		    fq[i] = r*ftri[0][i] + s*ftri[1][i] + t*ftri[2][i];
		    vq[i] = r*vtri[0][i] + s*vtri[1][i] + t*vtri[2][i];
	        }

		dA = Q.w(iq)*face.detJ;

		FOR_I3
		FOR_J3 {
		    tau[i][j] += 0.5*(fq[i]*xq[j] + fq[j]*xq[i])*dA;
		               + (viscRat - 1)*(vq[i]*face.normal[j] 
			                      + vq[j]*face.normal[i])*dA;
		}
	    } // iq
	} // ifa

	vol += cell.vol;
    } // icell

    m_dscal(9, 1.0/(vol*shRate), *tau);
} 


/* Calculate probe velocity 
 * Arguments:
 *  n -- number of probes
 *  x -- coordinates of probes
 *  v -- velocities at probe points 
 * Note:
 *  - Assume matched viscosity
 *  - Neglect the surface integral on particles 
 *  - Assume that the cell surface force density are pre-computed 
 *    and stored in f(:,3) */
void ShearSys::calcProbeVel(int n, const double (*x)[3], double (*v)[3])
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
	Cell &cell = cells[icell];
	for (int ivert = 0; ivert < cell.numVerts(); ivert++) {
	    Point &vert = cell.verts[ivert];
	    m_dcopy(3, &cell.f(ivert,0), vert.f);
	    m_dclear(3, vert.g);
	}
    }

    // Calculate surface integrals
    if (ewald::phys::active) {
        for (int blk = 0; blk < NBLK; blk++) {
            if (blk != CELL) continue;

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
        double vbkg[3];
	calcBkgVel(x[i], vbkg);
	m_dadd(3, vbkg, v[i]);
    }
}


// Write everything
void ShearSys::writeAll()
{
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    const char FN_FMT[] = "%s%6.6d%s";
    char token[256], fn[256];
    int nout;
    const int LT_CHUNK = 10000;

    // Cells
    strcpy(token, "CELL_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	sprintf(fn, FN_FMT, "D/cell", lt, ".dat");
	if (mpi_rank == 0) writeCells(fn);
    }

    // Rigids
    strcpy(token, "RIGID_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	sprintf(fn, FN_FMT, "D/rigid", lt, ".dat");
	if (mpi_rank == 0) writeRigids(fn);
    }

    // Restart
    strcpy(token, "RESTART_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	sprintf(fn, FN_FMT, "D/restart", lt, ".dat");
	if (mpi_rank == 0) writeRestart(fn);
    }

    // Center of cells
    strcpy(token, "CELL_CENTER_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
        sprintf(fn, FN_FMT, "D/cell_center", (lt/LT_CHUNK)*LT_CHUNK, ".dat");
        if (mpi_rank == 0) writeCellCenter(fn);
    }

    // Center of rigids
    strcpy(token, "RIGID_CENTER_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
        sprintf(fn, FN_FMT, "D/rigid_center", (lt/LT_CHUNK)*LT_CHUNK, ".dat");
        if (mpi_rank == 0) writeRigidCenter(fn);
    }

    // Probe
    strcpy(token, "PROBE_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	sprintf(fn, FN_FMT, "D/probe", (lt/LT_CHUNK)*LT_CHUNK, ".dat");
	if (mpi_rank == 0) writeProbe(fn);
    }

    // Tracer
    strcpy(token, "TRACER_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	sprintf(fn, FN_FMT, "D/tracer", (lt/LT_CHUNK)*LT_CHUNK, ".dat");
	if (mpi_rank == 0) writeTracer(fn);
    }

    // Stress
    strcpy(token, "STRESS_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	double tau[3][3];
	calcStress(tau);

	sprintf(fn, FN_FMT, "D/stress", Nt0, ".dat");
	if (mpi_rank == 0) writeStress(fn, tau);
    }
}


/* Write all cells to file
 * Arguments:
 *   fn -- file name */
void ShearSys::writeCells(const char *fn)
{
    if (numCells() <= 0) return;

    FILE *file = fopen(fn, "w");
//    fprintf(file, "variables = x, y, z, u, v, w\n");
    fprintf(file, "variables = x, y, z\n");

    for (int icell = 0; icell < numCells(); icell++) {
        Cell &cell = cells[icell];
	int nvert = cell.numVerts();
	int nface = cell.numFaces();
	fprintf(file, "zone N=%d E=%d F=FEPOINT ET=TRIANGLE\n", nvert, nface);
        
	// coordinates
	for (int ivert = 0; ivert < nvert; ivert++) {
	    Point &vert = cell.verts[ivert];
	    fprintf(file, "%10.3e %10.3e %10.3e\n", vert.x[0], vert.x[1], vert.x[2]);
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
void ShearSys::writeRigids(const char *fn)
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


void ShearSys::writeRestart(const char *fn)
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

    dims[0] = 1;
    H5LTmake_dataset_double(fid, "STRAIN", 1, dims, &ewald::strain);

    dims[0] = 1;
    H5LTmake_dataset_double(fid, "SHRATE", 1, dims, &shRate);

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

    // Cells 
    for (int icell = 0; icell < numCells(); icell++) {
        Cell &cell = cells[icell];
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

	if (cell.cellRef != NULL) {
	    cell.cellRef->getCoords(x);
	    dims[0] = nvert;
	    dims[1] = 3;
	    H5LTmake_dataset_double(gid, "XREF", 2, dims, *x);
	}

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
void ShearSys::readRestart(const char *fn)
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

	// For compatibility
	if (H5LTfind_dataset(fid, "STRAIN")) {
	    H5LTread_dataset_double(fid, "STRAIN", &ewald::strain);
	}
	else if (H5LTfind_dataset(fid, "XSHEAR")) {
	    double xshear;
	    H5LTread_dataset_double(fid, "XSHEAR", &xshear);
	    ewald::strain = xshear/ewald::L[2];
	}

	H5LTread_dataset_double(fid, "SHRATE", &shRate);

	printf("time = %.3f\n", time);
	printf("lt = %d\n", lt);
	printf("L = %.3f %.3f %.3f\n", ewald::L[0], ewald::L[1], ewald::L[2]);
	printf("strain = %.3f\n", ewald::strain);
	printf("shear rate = %.3E\n", shRate);
    }
    if (mpi_size > 1) {
	MPI_Bcast(ewald::L, 3, MPI_DOUBLE, 0, mpi_comm);
	MPI_Bcast(&time, 1, MPI_DOUBLE, 0, mpi_comm);
	MPI_Bcast(&lt, 1, MPI_INT, 0, mpi_comm);
	MPI_Bcast(&ewald::strain, 1, MPI_DOUBLE, 0, mpi_comm);
	MPI_Bcast(&shRate, 1, MPI_DOUBLE, 0, mpi_comm);
    }

    int ncell =0;
    int nrigid =0; 
    if (mpi_rank == 0) {
	if (H5LTfind_dataset(fid, "NCELL")) H5LTread_dataset_int(fid, "NCELL", &ncell);
	if (H5LTfind_dataset(fid, "NRIGID")) H5LTread_dataset_int(fid, "NRIGID", &nrigid);

	printf("ncell = %d\n", ncell);
	printf("nrigid = %d\n", nrigid);
    }
    if (mpi_size > 1) {
	MPI_Bcast(&ncell, 1, MPI_INT, 0, mpi_comm);
	MPI_Bcast(&nrigid, 1, MPI_INT, 0, mpi_comm);
    }

    cells.resize(ncell);
    rigids.resize(nrigid);

    // Read cells
    for (int icell = 0; icell < ncell; icell++) {
        Cell &cell = cells[icell];
        int nvert, nface, xref_flag;

	if (mpi_rank == 0) {
	    sprintf(token, "CELL%d", icell);
	    gid = H5Gopen(fid, token, H5P_DEFAULT);

	    H5LTget_dataset_info(gid, "X", dims, &class_id, &type_size);
	    nvert = dims[0];

	    H5LTget_dataset_info(gid, "F2V", dims, &class_id, &type_size);
	    nface = dims[0];

	    xref_flag = H5LTfind_dataset(gid, "XREF");
        }
	if (mpi_size > 1) {
	    MPI_Bcast(&nvert, 1, MPI_INT, 0, mpi_comm);
	    MPI_Bcast(&nface, 1, MPI_INT, 0, mpi_comm);
	    MPI_Bcast(&xref_flag, 1, MPI_INT, 0, mpi_comm);
	}

        MArray<double,2> x(nvert,3), xref(nvert,3);
	MArray<int,2> f2v(nface,3);

	if (mpi_rank == 0) {
	    H5LTread_dataset_double(gid, "X", x.data());
	    H5LTread_dataset_int(gid, "F2V", f2v.data());
	    if (xref_flag) H5LTread_dataset_double(gid, "XREF", xref.data());
        }
	if (mpi_size > 1) {
	    MPI_Bcast(x.data(), x.size(), MPI_DOUBLE, 0, mpi_comm);
	    MPI_Bcast(f2v.data(), f2v.size(), MPI_INT, 0, mpi_comm);
	    if (xref_flag) {
		MPI_Bcast(xref.data(), xref.size(), MPI_DOUBLE, 0, mpi_comm);
	    }
	}

	cell.verts.resize(nvert);
	cell.setCoords(x);

	cell.faces.resize(nface);
	cell.setConnectivities(f2v);

        if (xref_flag) {
	    cell.cellRef = new Cell();

	    cell.cellRef->verts.resize(nvert);
	    cell.cellRef->setCoords(xref);

	    cell.cellRef->faces.resize(nface);
	    cell.cellRef->setConnectivities(f2v);
	}

	if (mpi_rank == 0) H5Gclose(gid);
    } // icell 


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

    // Probes
    strcpy(token, "PROBE");
    int probe_flag;

    if (mpi_rank == 0) {
        probe_flag = H5LTfind_dataset(fid, token);
    }
    if (mpi_size > 1) {
	MPI_Bcast(&probe_flag, 1, MPI_INT, 0, mpi_comm);
    }

    if (probe_flag) {
        if (mpi_rank == 0) {
	    H5LTget_dataset_info(fid, token, dims, &class_id, &type_size);
	    nprb = dims[0];
	    printf("nprobe = %d\n", nprb);
	}
	if (mpi_size > 1) {
	    MPI_Bcast(&nprb, 1, MPI_INT, 0, mpi_comm);
	}

	xprb = new double[nprb][3];
	vprb = new double[nprb][3];

	if (mpi_rank == 0) {
	    H5LTread_dataset_double(fid, token, *xprb);
	}
	if (mpi_size > 1) {
	    MPI_Bcast(*xprb, 3*nprb, MPI_DOUBLE, 0, mpi_comm);
        }
    } else {
        nprb = 0;
	xprb = vprb = NULL;
    }

    // Tracers
    strcpy(token, "TRACER");
    int tracer_flag;

    if (mpi_rank == 0) {
        tracer_flag = H5LTfind_dataset(fid, token);
    }
    if (mpi_size > 1) {
	MPI_Bcast(&tracer_flag, 1, MPI_INT, 0, mpi_comm);
    }

    if (tracer_flag) {
        if (mpi_rank == 0) {
	    H5LTget_dataset_info(fid, token, dims, &class_id, &type_size);
	    ntrac = dims[0];
	    printf("ntracer = %d\n", ntrac);
        }
	if (mpi_size > 1) {
	    MPI_Bcast(&ntrac, 1, MPI_INT, 0, mpi_comm);
	}

	xtrac = new double[ntrac][3];
	vtrac = new double[ntrac][3];

	if (mpi_rank == 0) {
	    H5LTread_dataset_double(fid, token, *xtrac);
	}
	if (mpi_size > 1) {
	    MPI_Bcast(*xtrac, 3*ntrac, MPI_DOUBLE, 0, mpi_comm);
        }
    }
    else {
        ntrac = 0;
	xtrac = vtrac = NULL;
    }

    if (mpi_rank == 0) H5Fclose(fid);
}


/* Write probes */
void ShearSys::writeProbe(const char *fn)
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


void ShearSys::writeTracer(const char *fn)
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


void ShearSys::writeCellCenter(const char *fn)
{
    int ncell = numCells();
    if (ncell <= 0) return;

    MArray<double,2> x(ncell,3), v(ncell,3);

    for (int icell = 0; icell < ncell; icell++) {
        Cell &cell = cells[icell];
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
void ShearSys::writeRigidCenter(const char *fn)
{
    int nrigid = numRigids();
    if (nrigid <= 0) return;

    MArray<double,2> x(nrigid,3), v(nrigid,3);

    for (int irigid = 0; irigid < nrigid; irigid++) {
        Rigid &rigid = rigids[irigid];
	FOR_I3 x(irigid,i) = rigid.center[i];
	rigid.calcCentroidVel(rigid.v, &v(irigid,0));
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
	H5LTmake_dataset_double(fid, token, 2, dims, x.data());
    }
    else {
        dset = H5Dopen(fid, token, H5P_DEFAULT);
	H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, x.data());
	H5Dclose(dset);
    }

    sprintf(token, "V%6.6d", lt);
    if (! H5LTfind_dataset(fid, token)) {
        dims[0] = nrigid;
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


/* Write the add stress */
void ShearSys::writeStress(const char *fn, const double (*tau)[3])
{
    static bool first_write = true;
    FILE *file;

    if (first_write) {
        file = fopen(fn, "w");
        first_write = false;
    }
    else {
        file = fopen(fn, "a");
    }

    fprintf(file, "%12.5E %12.5E %12.5E %12.5E %12.5E\n",
                time, tau[0][2], tau[0][0], tau[1][1], tau[2][2]);

    fclose(file);
}


/* Cell volume fraction 
 * Note:
 *   -- Collective */
double ShearSys::cellVolFraction()
{
    double vol = 0.0;

    for (int icell = 0; icell < numCells(); icell++) {
        Cell &cell = cells[icell];
        vol += cell.vol;
    }

    return vol/ewald::vol;
}


/* Calculate max perturbation velocity 
 * Argument:
 *   dvmax -- max cell surface perturbation velcoity
 *   UTmax -- max rigid center translation velocity
 *   OMGmax -- max rigid angular velocity
 * Note:
 *   -- Collective */
void ShearSys::calcMaxPerturbVel(double *dvmax, double *UTmax, double *OMGmax)
{
    // Init
    FOR_D3 dvmax[d] = 0.0;
    FOR_D3 UTmax[d] = 0.0;
    FOR_D3 OMGmax[d] = 0.0;

    // Cells
    for (int icell = 0; icell < numCells(); icell++) {
	Cell &cell = cells[icell];

	for (int ivert = 0; ivert < cell.numVerts(); ivert++) {
	    Point &vert = cell.verts[ivert];

	    double vbkg[3];
	    calcBkgVel(vert.x, vbkg);

	    FOR_D3 dvmax[d] = max(dvmax[d], fabs(cell.v(ivert,d) - vbkg[d]));
	}
    } 

    // Rigid
    for (int irigid = 0; irigid < numRigids(); irigid++) {
	Rigid &rigid = rigids[irigid];

	double UT[3], OMG[3];
	rigid.calcTransRotatVel(rigid.v, UT, OMG);

	double UT0[3];
	calcBkgVel(rigid.center, UT0);
	FOR_I3 UT[i] -= UT0[i];

	FOR_I3 UTmax[i] = max( UTmax[i], fabs(UT[i]));
	FOR_I3 OMGmax[i] = max( OMGmax[i], fabs(OMG[i]));
    }
}


/* Area change per element
 * Arguments:
 *   dAmin, dAmax, dAmean -- min, max and mean relative area change
 * Note:
 *   -- Collective */
void ShearSys::calcAreaChange(double &dAmin, double &dAmax, double &dAmean)
{
    // Init
    dAmin = FLT_MAX;
    dAmax = -FLT_MAX;
    double dAsum = 0.0, Asum = 0.0;

    for (int icell = 0; icell < numCells(); icell++) {
	Cell &cell = cells[icell];
	Cell *cellRef = cell.cellRef;
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

    // Special case
    if (dAmin > dAmax) dAmin = dAmax = 0.0;
    dAmean = dAsum/Asum;
}
